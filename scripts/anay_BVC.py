#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the bilevel vertex cover results 

@author: Xueyu shi
"""


from collections import defaultdict 
import numpy as np
import re,os
data = defaultdict(dict)


BVC_TYPE = "Symm" # "ASymm" for asymmetric case. "Symm" for symmetric case


for n in range(10, 55, 5):
# for n in [25, 35, 40, 45, 50]:
    deg = b = int(np.ceil(n/2))
    
    mibs_log = os.path.join("Log/Log_BVC_{}".format(BVC_TYPE), "BVC_n{}d{}b{}".format(str(n), str(deg), str(deg)))

    for j in range(10):
        log_name = mibs_log + 'j{}_mibs.log'.format(str(j))
        if not os.path.exists(log_name): continue
        data[(n,deg, b)][j] = {"mibs": np.array([0, 0, 0], dtype='float')}
        
        with open(log_name, 'r') as f:
            for line in f.readlines():
                if "Cost" in line:
                    data[(n,deg, b)][j]["mibs"][0] = -float(line.strip().split()[-1])
                elif "optimality gap" in line:
                    data[(n,deg, b)][j]["mibs"][1] = 1 / ( 1 - float(line.strip().split()[-1][:-1]) / 100)
                    


with open("BVC_{}_results_mibs.log".format(BVC_TYPE), 'r') as f:
    for line in f.readlines():
        
        line = line.strip().split()
        ins = line[0] 
        
        n = int(re.search('n[0-9]+', ins).group(0)[1:])
        deg = int(re.search('d[0-9]+', ins).group(0)[1:])
        b = int(re.search('b[0-9]+', ins).group(0)[1:])
        j = int(re.search('j[0-9]+', ins).group(0)[1:])
        data[(n,deg, b)][j]["mibs"][-1] = min(1800, float(line[-1]))
        
        


with open("BVC_{}_results.txt".format(BVC_TYPE), 'r') as f:
    line_num = 0
    for line in f.readlines():
        line_num += 1
        if line_num <= 2: continue
        
        line = line.strip().split()
        ins = line[0]
        
        n = int(re.search('n[0-9]+', ins).group(0)[1:])
        deg = int(re.search('d[0-9]+', ins).group(0)[1:])
        b = int(re.search('b[0-9]+', ins).group(0)[1:])
        j = int(re.search('j[0-9]+', ins).group(0)[1:])
        
        if (n, deg, b) not in data: continue

        data[(n,deg, b)][j][0] = np.array([float(x) for x in line[1:4]])
        for k in range(1,4):
            data[(n,deg, b)][j][k] = np.array([float(x) for x in line[4*k:4*(k+1)]])
        
        
        base_val = data[(n,deg, b)][j]["mibs"][0]
        for k in range(4):
            tmp = (data[(n,deg, b)][j][k][0] + 0.1) / (base_val + 0.1)
            data[(n,deg, b)][j][k][0] = (data[(n,deg, b)][j][k][1] + 0.1) / (base_val + 0.1)
            data[(n,deg, b)][j][k][1] = tmp

avg_data = {}

for key, val in data.items():
    avg_val = [0]*17
    num_ins = len(val)
    for _, ins in val.items():
        
        avg_val[:2] += ins['mibs'][1:] / num_ins
        avg_val[2:5] +=  ins[0] / num_ins
        for k in range(1, 4):
            avg_val[1 + 4*k:1+4*(k+1)] += ins[k] / num_ins
    
    avg_data[key] = [round(x, 2) for x in avg_val]


import pandas as pd
dataframe = pd.DataFrame.from_dict(avg_data, orient='index')
dataframe.to_excel('BVC_{}_avg_results_new.xlsx'.format(BVC_TYPE))
        
        
        
        
