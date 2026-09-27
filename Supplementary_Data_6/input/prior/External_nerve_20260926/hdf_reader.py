"""Read the seven CellBender datasets with h5py or an installed HDF5 library.

The fallback calls the official HDF5 C API. It does not implement HDF5 parsing.
Only integer, floating-point and string arrays are supported.
"""
import ctypes as C
import ctypes.util
import importlib.util
import io
from pathlib import Path
import tempfile
import numpy as np

PATHS=['features/name','barcodes','data','indices','indptr','shape']


def read_cellbender(raw):
    try:
        import h5py
    except ImportError:
        return _read_c_api(raw)
    with h5py.File(io.BytesIO(raw),'r') as f:
        return {p:f['matrix/'+p][:] for p in PATHS}


def _read_c_api(raw):
    library=ctypes.util.find_library('hdf5_serial') or ctypes.util.find_library('hdf5')
    prefix=''
    if library is None:
        spec=importlib.util.find_spec('vtkmodules')
        if spec is None: raise ImportError('Install h5py in the analysis environment to read CellBender files.')
        files=list(Path(next(iter(spec.submodule_search_locations))).glob('libvtkhdf5-*.so'))
        if len(files)!=1: raise ImportError('No usable installed HDF5 library found.')
        library=str(files[0]); prefix='vtkhdf5_'
    lib=C.CDLL(library)
    hid=C.c_int64
    def fn(name,restype,argtypes):
        f=getattr(lib,prefix+name); f.restype=restype; f.argtypes=argtypes
        return f
    init=fn('H5open',C.c_int,[])
    fopen=fn('H5Fopen',hid,[C.c_char_p,C.c_uint,hid])
    fclose=fn('H5Fclose',C.c_int,[hid])
    dopen=fn('H5Dopen2',hid,[hid,C.c_char_p,hid])
    dclose=fn('H5Dclose',C.c_int,[hid])
    space=fn('H5Dget_space',hid,[hid])
    sclose=fn('H5Sclose',C.c_int,[hid])
    npoints=fn('H5Sget_simple_extent_npoints',C.c_int64,[hid])
    dtype=fn('H5Dget_type',hid,[hid])
    tclose=fn('H5Tclose',C.c_int,[hid])
    tclass=fn('H5Tget_class',C.c_int,[hid])
    tsize=fn('H5Tget_size',C.c_size_t,[hid])
    tsign=fn('H5Tget_sign',C.c_int,[hid])
    torder=fn('H5Tget_order',C.c_int,[hid])
    tvlen=fn('H5Tis_variable_str',C.c_int,[hid])
    dread=fn('H5Dread',C.c_int,[hid,hid,hid,hid,hid,C.c_void_p])
    reclaim=fn('H5Dvlen_reclaim',C.c_int,[hid,hid,hid,C.c_void_p])
    assert init()>=0
    result={}
    with tempfile.NamedTemporaryFile(suffix='.h5') as tmp:
        tmp.write(raw); tmp.flush()
        file=fopen(tmp.name.encode(),0,0)
        assert file>=0
        try:
            for path in PATHS:
                d=dopen(file,('matrix/'+path).encode(),0)
                assert d>=0,path
                s=space(d); t=dtype(d); n=npoints(s); cls=tclass(t); size=tsize(t)
                try:
                    assert n>=0
                    if cls==3 and tvlen(t)>0:
                        buf=(C.c_char_p*n)()
                        assert dread(d,t,0,0,0,buf)>=0,path
                        result[path]=np.array([buf[i] for i in range(n)],dtype=object)
                        assert reclaim(t,s,0,buf)>=0
                    else:
                        if cls==3: dt=np.dtype('S'+str(size))
                        elif cls in [0,1]:
                            order=torder(t)
                            assert order in [0,1],(path,order)
                            end='<' if order==0 else '>'
                            kind='f' if cls==1 else ('i' if tsign(t)==1 else 'u')
                            dt=np.dtype(end+kind+str(size))
                        else: raise TypeError((path,cls))
                        a=np.empty(n,dtype=dt)
                        assert dread(d,t,0,0,0,a.ctypes.data)>=0,path
                        result[path]=a
                finally:
                    tclose(t); sclose(s); dclose(d)
        finally: fclose(file)
    return result
