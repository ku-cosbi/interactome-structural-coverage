import pandas as pd
import numpy as np


def remove_duplicates(int_dtb, int_dtb_ids=['protein1_uniprot', 'protein2_uniprot']):
    int_dtb_interface = int_dtb[int_dtb_ids]
    dups_removed = int_dtb_interface.apply(lambda x: tuple(sorted([x[int_dtb_ids[0]], x[int_dtb_ids[1]]])), axis=1)
    dups_removed_df = pd.DataFrame(list(dict.fromkeys(dups_removed)), columns=int_dtb_ids)
    return dups_removed_df

def get_coverage(interactome, model_dtb, int_dtb_ids=['protein1_uniprot', 'protein2_uniprot']):
    "Get statistics of an interactome database wrt 3D structure source"
    interactome_proteins = np.unique(interactome[int_dtb_ids].to_numpy().flatten()) # all unique proteins in the interactome dtb
    model_dtb = model_dtb.drop_duplicates(['Entry']) # protein IDs covered by the structure/model dtb
    num_nodes = model_dtb[model_dtb['Entry'].isin(interactome_proteins)] # how many of them are in the pdb/homology dtb
    cov_ints_both_sides = interactome[(interactome[int_dtb_ids[0]].isin(model_dtb['Entry'])) & (interactome[int_dtb_ids[1]].isin(model_dtb['Entry']))]
    return len(num_nodes), len(cov_ints_both_sides)

def get_coverage_both_sides_df(interactome, model_dtb, int_dtb_ids=['protein1_uniprot', 'protein2_uniprot']):
    interactome_proteins = np.unique(interactome[int_dtb_ids].to_numpy().flatten())
    model_dtb = model_dtb.drop_duplicates(['Entry'])
    num_nodes = model_dtb[model_dtb['Entry'].isin(interactome_proteins)]
    cov_ints_both_sides = interactome[(interactome[int_dtb_ids[0]].isin(model_dtb['Entry'])) & (interactome[int_dtb_ids[1]].isin(model_dtb['Entry']))]
    return cov_ints_both_sides

def get_ready_for_upset(df, df_name, int_dtb_ids=['protein1_uniprot', 'protein2_uniprot']):
    df2 = df.apply(lambda x: tuple(sorted([x[int_dtb_ids[0]], x[int_dtb_ids[1]]])), axis=1)
    df2_dict = {key: 1 for key in dict.fromkeys(df2)}
    df3 = pd.DataFrame(list(df2_dict.items()))
    df3 = df3.rename(columns={0: 'interaction', 1: df_name})
    return df3


def get_ready_for_upset_str(interactome_mapped, structure_dtb, df_name, int_dtb_ids=['protein1_uniprot', 'protein2_uniprot']):
    #dups_removed = remove_duplicates(interactome_mapped, int_dtb_ids=int_dtb_ids)
    df_both_sides = get_coverage_both_sides_df(interactome=interactome_mapped, model_dtb=structure_dtb, int_dtb_ids=int_dtb_ids)
    df_both_sides = df_both_sides.reset_index(drop=True)
    df_upset_ready = get_ready_for_upset(df_both_sides, df_name=df_name)
    return df_upset_ready


def get_ready_for_upset_proteome(df, df_name, col='Entry'):
    df2 = df[col].drop_duplicates()
    df2_dict = {key: 1 for key in dict.fromkeys(df2)}
    df3 = pd.DataFrame(list(df2_dict.items()))
    df3 = df3.rename(columns={0: 'Entry', 1: df_name})
    return df3
