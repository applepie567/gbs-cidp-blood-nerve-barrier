"""Retrieve source metadata that define vector slots, filters and sample links."""
from pathlib import Path
import base64
import json
from public_sqlite import CachedReader

ROOT=Path(__file__).resolve().parent
ACCESS=ROOT/'access'
SOURCE='https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD056286/QExactiveHF02_18988_N1530.'

def enc(x):return {'base64':base64.b64encode(x).decode()}

def main():
    jobs={
      'pdResult':['DataDistributionMaps','DataDistributionLevels','DataDistributionBoxes','DataDistributionBoxExtendedData','DataTypesColumns','StudyInformation','AnalysisDefinition','WorkflowInputFiles','StudyInformationWorkflowInputFiles','EnumDataTypes','EnumDataTypeValues','ProcessingNodeCustomData','AbundanceCorrectionItems','SettingsValues','Workflows','ResultStatistics'],
      'pdResultView':['DataTypeSubsets','DataTypeSubsetDetails','SettingsValues','SettingsOwner','DataTypeAdditionalColumns','TargetProteins_Subset_4','TargetProteins_Subset_5'],
      'msfView':['DataTypeSubsets','DataTypeSubsetDetails','SettingsValues','SettingsOwner','DataTypeAdditionalColumns']}
    for suffix,names in jobs.items():
        cache=ACCESS/('pdResult_blocks64' if suffix=='pdResult' else suffix+'_blocks')
        cache.mkdir(exist_ok=True)
        r=CachedReader(SOURCE+suffix,cache)
        (ACCESS/(suffix+'_schema.json')).write_text(json.dumps(r.schema,indent=2))
        for name in names:
            out=ACCESS/(suffix+'_'+name+'_with_ids.json')
            if out.exists():continue
            if not any(x[0]=='table' and x[1]==name for x in r.schema):continue
            rows=list(r.table(name,with_rowid=True))
            out.write_text(json.dumps(rows,indent=2,default=enc))
            print(suffix,name,len(rows),flush=True)
        (ACCESS/(suffix+'_metadata_ranges.json')).write_text(json.dumps(r.manifest(),indent=2))

if __name__=='__main__':main()
