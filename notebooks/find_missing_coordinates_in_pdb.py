from Bio.PDB import *
import pandas as pd
import numpy as np
import os
from collections import defaultdict

path=r"enter path"
os.chdir(path)

def find_missing(pdb_id, chain_id):
    
    try:
        pdb_id = pdb_id+".cif"
        chain_id = chain_id
        mmcif_dict = MMCIF2Dict.MMCIF2Dict(pdb_id)
    except:
        pass
    
    try:
        #translate PDB id(id in right column) to uniprot id (id in left column)
        missing_coord_chain_pdb = mmcif_dict["_pdbx_unobs_or_zero_occ_residues.auth_asym_id"]
        missing_coord_chain_uniprot_ori = mmcif_dict["_pdbx_unobs_or_zero_occ_residues.label_asym_id"]
    
        #no duplicates in chain id
        missing_coord_chain_pdb = list(dict.fromkeys(missing_coord_chain_pdb))
        missing_coord_chain_uniprot =list(dict.fromkeys(missing_coord_chain_uniprot_ori))
    
        zip_chains = zip(missing_coord_chain_pdb, missing_coord_chain_uniprot)
    
        d1 = defaultdict(list)
        for k, v in zip_chains:
            d1[k].append(v)
        chain_conversion = d1.get(chain_id)[0]
        
    except:
        missing_coord_chain_pdb = None
        missing_coord_chain_uniprot_ori = None
        missing_coord_chain_uniprot = None
        chain_conversion = None

    try: 
        missing_coord_uniprot = mmcif_dict["_pdbx_unobs_or_zero_occ_residues.label_seq_id"] 

        missing_info_list = zip(missing_coord_chain_uniprot_ori, missing_coord_uniprot)
        d2 = defaultdict(list)

        for k, v in missing_info_list:
            d2[k].append(v)
        missing_info = (dict((k, list(v)) for k, v in d2.items())).get(chain_conversion)
        
    except:
        
        missing_coord_uniprot = None
        missing_info = None

    return [missing_info]
	

# Call function
	
human_proteome_pdb = pd.read_excel(file)
pdb_id = human_proteome_pdb["PDB"]
chain = human_proteome_pdb["CHAIN"]
get_pdb_res_interval_output = list(map(find_missing,pdb_id,chain))
human_proteome_pdb["Missing"] = get_pdb_res_interval_output

	
