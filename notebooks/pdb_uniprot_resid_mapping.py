import pandas as pd
import numpy as np
from Bio.PDB import *
import os 

os.chdir(r"path")
data = pd.read_excel("data.xlsx")

def pdb_uniprot_mapping(pdb_id):
    
    try:
        pdb_id = pdb_id+".cif"
        mmcif_dict = MMCIF2Dict.MMCIF2Dict(pdb_id)
    except:
        pass
    
    try:
        pdb = mmcif_dict["_struct_ref_seq.pdbx_PDB_id_code"]
        uniprot = mmcif_dict["_struct_ref_seq.pdbx_db_accession"]
        chain_id = mmcif_dict["_struct_ref_seq.pdbx_strand_id"] 
        sp_beg = mmcif_dict["_struct_ref_seq.db_align_beg"] #SP_BEG
        sp_end = mmcif_dict["_struct_ref_seq.db_align_end"] #SP_END
        res_beg = mmcif_dict["_struct_ref_seq.seq_align_beg"] #RES_BEG
        res_end = mmcif_dict["_struct_ref_seq.seq_align_end"] #RES_END
        
    except:
        pdb = None
        uniprot = None
        chain_id = None
        sp_beg = None
        sp_end = None
        res_beg = None
        res_end = None
        
    return pdb, chain_id, uniprot, res_beg, res_end, sp_beg, sp_end


#Call Function

func = list(map(pdb_uniprot_mapping,data["PDB"]))

# Write new DF 

columns= ["PDB", "CHAIN", "UNIPROT_ID", "RES_BEG", "RES_END", "SP_BEG", "SP_END"]
data = pd.DataFrame(func,columns=columns)
data= data.apply(pd.Series.explode)


