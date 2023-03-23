import pandas as pd
import numpy as np
from collections import Counter


data = pd.read_excel(r"data.xlsx")

data["Missing"] = data["Missing"].fillna("[-100000]")
data["Missing"] = list(map(lambda s: list(map(int,s[1:-1].split(","))), data["Missing"]))

data["dif"] = data["SP_BEG"] - data["RES_BEG"]
data["MissingInterval"] = data.apply(lambda x: list(map(lambda y: y+x["dif"], x["Missing"])), axis=1)
data["MissingInterval"] = data.apply(lambda x: None if x["MissingInterval"][0] < -9999 else x["MissingInterval"], axis=1)

data['combined'] = data.apply(lambda x: list([x["SP_BEG"],x["SP_END"]]),axis=1)


def exclude_out_missing(missing_all, interval):
    try:
        filtered = filter(lambda num: num >= interval[0] and num <= interval[1] , missing_all)
    except:
        filtered = []
    return list(filtered)

data["MissingTrueInterval"] = data.apply(lambda x: exclude_out_missing(missing_all = x["MissingInterval"], interval = x["combined"]), axis=1)


def calc_cov(begin_end, missing, prot_length):

    covered_resid1 = []
    covered_resid2 = set()
 
    for pair in begin_end:
        
        p1, p2 = pair[0],pair[1]
        p_min, p_max = sorted([int(p1), int(p2)])
        
        covered_resid1.extend(range(p_min, p_max + 1))
    
    covered_resid1_exclude_missing=list((Counter(covered_resid1) - Counter(missing)).elements())
    

    covered_resid2 = set(covered_resid1_exclude_missing)
    coverage_percent = len(covered_resid2) * 100 / prot_length
    
    return coverage_percent
	

data["Coverage"] = data.apply(lambda x: calc_cov(begin_end = x["combined"], 
                                                              missing = x["MissingTrueInterval"],
                                                              prot_length = x["Length"]), axis=1)
