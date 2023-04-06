import numpy as np
import pandas as pd
from collections import Counter

def calculate_coverage(begin_end, prot_length): 
    covered_resid1 = set() 
    for pair in begin_end:
        p1, p2 = pair[0],pair[1]
        p_min, p_max = sorted([int(p1), int(p2)])
        covered_resid1.update(range(p_min, p_max + 1)) 
    coverage_percent = len(covered_resid1) * 100 / prot_length
    return coverage_percent


def get_res_coverage_by_threshold(df, threshold, gb="Entry", upper=True):
    dfs = [x for _, x in df.groupby(by=gb)]
    coverage_percents = []
    for df in dfs:
        covered_resid = set()
        prot_length = list(set(df["Length"].to_list()))[0]
        begin_end = df[['TargetBeg', 'TargetEnd']].to_numpy()
        for pair in begin_end:
            p1, p2 = pair[0], pair[1]
            p_min, p_max = sorted([int(p1), int(p2)])
            covered_resid.update(range(p_min, p_max + 1))
        coverage_percents.append(len(covered_resid) * 100 / prot_length)
    num_proteins = 0
    
    for i in coverage_percents:
        if upper:
            if i >= threshold:
                num_proteins += 1
        else:
            if i < threshold:
                num_proteins += 1
    # pct_threshold = num_threshold/len(coverage_percents) * 100 
    return num_proteins


def get_res_coverage_by_threshold_extra(df, threshold, gb="Entry"):
    dfs = [x for _, x in df.groupby(by=gb)]
    coverage_percents = []
    for df in dfs:
        covered_resid = set()
        prot_length = list(set(df["Length"].to_list()))[0]
        for model in range(len(df)):
            covered_resid.update(df.iloc[model].difference)
        coverage_percents.append(len(covered_resid) * 100 / prot_length)
    
    num_proteins = 0
    for i in coverage_percents:
        if i >= threshold:
            num_proteins += 1
    return num_proteins
def get_protein_based_coverage(num_proteins, total_prot_num):
    prot_cov = num_proteins/total_prot_num * 100
    return np.round(prot_cov, 2)


def get_coverage_percent(df, proteome_length, gb="Entry"):
    dfs = [x for _, x in df.groupby(by=gb)]
    coverage_percent = 0
    for df in dfs:
        covered_resid = set()
        begin_end = df[['TargetBeg','TargetEnd']].to_numpy()
        for pair in begin_end:
            p1, p2 = pair[0], pair[1]
            p_min, p_max = sorted([int(p1), int(p2)])
            covered_resid.update(range(p_min, p_max + 1))
        coverage_percent += len(covered_resid) * 100 / proteome_length
    return coverage_percent


def get_coverage_percent_extra(df, proteome_length, gb="Entry"):
    dfs = [x for _, x in df.groupby(by=gb)]
    coverage_percent = 0
    for df in dfs:
        covered_resid = set()
        for model in range(len(df)):
            covered_resid.update(df.iloc[model].difference)
        coverage_percent += len(covered_resid) * 100 / proteome_length
    return coverage_percent

def get_residue_coverage(df, gb="Entry", length="Length", res_beg_end=['TargetBeg','TargetEnd']):
    dfs = [x for _, x in df.groupby(by=gb)]
    coverage_percents = []
    for df in dfs:
        covered_resid = set()
        prot_length = list(set(df[length].to_list()))[0]
        begin_end = df[res_beg_end].to_numpy()
        for pair in begin_end:
            p1, p2 = pair[0],pair[1]
            p_min, p_max = sorted([int(p1), int(p2)])
            covered_resid.update(range(p_min, p_max + 1))
        coverage_percents.append(len(covered_resid) * 100 / prot_length)
    return coverage_percents


def get_uncommon_ind(data, ind):
    uncommon_entry_ind = []
    for index in list(data.index):
        if index not in ind:
            uncommon_entry_ind.append(index)
    return uncommon_entry_ind


def get_common_uncommon_elements(data1, data2):
    list1 = data1.Entry.tolist()
    list2 = data2.Entry.tolist()
    
    common_entries = []
    list1_ind = []
    list2_ind = []
    for i, element in enumerate(list1):
        if element in list2:
            common_entries.append(element)
            list1_ind.append(i)
    for j, element in enumerate(list2):
        if element in list1:
            list2_ind.append(j)
    
    data1_common = data1.iloc[list1_ind]
    data2_common = data2.iloc[list2_ind]
    
    commons = pd.concat([data1_common, data2_common], ignore_index=True)
    
    un_data1_ind = get_uncommon_ind(data1, list1_ind)
    un_data2_ind = get_uncommon_ind(data2, list2_ind)
    
    data1_uncommon = data1.iloc[un_data1_ind]
    data2_uncommon = data2.iloc[un_data2_ind]
    
    
    return commons, data1_common, data2_common, data1_uncommon, data2_uncommon