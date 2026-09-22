"""Read ordinary HDF5 arrays through h5py or an existing HDF5 C library.

The C-library fallback is read-only and uses the documented H5F/H5D/H5S/H5T
API. It supports numeric and string datasets needed for deposited 10X matrices.
No custom HDF5 binary decoder is used. https://support.hdfgroup.org/documentation/
"""
import ctypes as C
import ctypes.util
import os
from pathlib import Path
import numpy as np


class H5Reader:
    def __init__(self, path):
        self.path = Path(path)
        try:
            import h5py
            self.file = h5py.File(self.path, 'r')
            self.backend = 'h5py ' + h5py.__version__
            self.lib = None
            return
        except ImportError:
            pass
        paths = [os.environ.get('HDF5_SHARED_LIBRARY')]
        for n in ['hdf5_serial', 'hdf5']:
            paths.append(C.util.find_library(n))
        root = os.environ.get('CODEX_PRIMARY_RUNTIME_ROOT')
        if root:
            paths.extend(str(p) for p in Path(root).glob('dependencies/python/lib/python*/site-packages/vtkmodules/libvtkhdf5-*.so'))
        err = None
        for p in filter(None, paths):
            try:
                lib = C.CDLL(p)
                prefix = 'vtkhdf5_' if hasattr(lib, 'vtkhdf5_H5Fopen') else ''
                if hasattr(lib, prefix + 'H5Fopen'):
                    self.lib, self.prefix, self.backend = lib, prefix, str(p)
                    break
            except OSError as e:
                err = e
        else:
            raise ImportError('h5py or an existing HDF5 shared library is required') from err
        hid = C.c_longlong
        self._bind('H5open', [], C.c_int)()
        self.fopen = self._bind('H5Fopen', [C.c_char_p, C.c_uint, hid], hid)
        self.fclose = self._bind('H5Fclose', [hid], C.c_int)
        self.dopen = self._bind('H5Dopen2', [hid, C.c_char_p, hid], hid)
        self.dclose = self._bind('H5Dclose', [hid], C.c_int)
        self.dspace = self._bind('H5Dget_space', [hid], hid)
        self.dtype = self._bind('H5Dget_type', [hid], hid)
        self.ndim = self._bind('H5Sget_simple_extent_ndims', [hid], C.c_int)
        self.dims = self._bind('H5Sget_simple_extent_dims', [hid, C.POINTER(C.c_ulonglong), C.POINTER(C.c_ulonglong)], C.c_int)
        self.tclass = self._bind('H5Tget_class', [hid], C.c_int)
        self.tsize = self._bind('H5Tget_size', [hid], C.c_size_t)
        self.tsign = self._bind('H5Tget_sign', [hid], C.c_int)
        self.torder = self._bind('H5Tget_order', [hid], C.c_int)
        self.tvariable = self._bind('H5Tis_variable_str', [hid], C.c_int)
        self.read = self._bind('H5Dread', [hid, hid, hid, hid, hid, C.c_void_p], C.c_int)
        self.sclose = self._bind('H5Sclose', [hid], C.c_int)
        self.tclose = self._bind('H5Tclose', [hid], C.c_int)
        self.reclaim = self._bind('H5Dvlen_reclaim', [hid, hid, hid, C.c_void_p], C.c_int)
        self.fid = self.fopen(os.fsencode(self.path), 0, 0)
        if self.fid < 0:
            raise OSError(f'Cannot open {self.path}')

    def _bind(self, name, args, result):
        f = getattr(self.lib, self.prefix + name)
        f.argtypes, f.restype = args, result
        return f

    def array(self, name):
        if self.lib is None:
            a = self.file[name][:]
            if a.dtype.kind in ['S', 'O']:
                return np.array([x.decode() if isinstance(x, bytes) else x for x in a])
            return a
        did = self.dopen(self.fid, name.encode(), 0)
        if did < 0:
            raise KeyError(name)
        sid, tid = self.dspace(did), self.dtype(did)
        try:
            ndim = self.ndim(sid)
            dims = (C.c_ulonglong * ndim)()
            self.dims(sid, dims, None)
            shape = tuple(dims)
            n = int(np.prod(shape))
            cls, size = self.tclass(tid), self.tsize(tid)
            if cls == 3 and self.tvariable(tid) > 0:
                buf = (C.c_char_p * n)()
                if self.read(did, tid, 0, 0, 0, buf) < 0:
                    raise OSError(name)
                a = np.array([x.decode('utf-8') if x is not None else '' for x in buf]).reshape(shape)
                self.reclaim(tid, sid, 0, buf)
                return a
            if cls == 3:
                dt = np.dtype(f'S{size}')
            elif cls in [0, 1]:
                endian = '<' if self.torder(tid) == 0 else '>'
                kind = 'f' if cls == 1 else ('i' if self.tsign(tid) else 'u')
                dt = np.dtype(f'{endian}{kind}{size}')
            else:
                raise TypeError(f'Unsupported HDF5 class {cls} at {name}')
            a = np.empty(shape, dtype=dt)
            if self.read(did, tid, 0, 0, 0, a.ctypes.data_as(C.c_void_p)) < 0:
                raise OSError(name)
            return np.char.decode(a, 'utf-8') if cls == 3 else a
        finally:
            self.tclose(tid)
            self.sclose(sid)
            self.dclose(did)

    def close(self):
        self.file.close() if self.lib is None else self.fclose(self.fid)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
