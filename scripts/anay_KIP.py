#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the knapsack interdiction results 

@author: Xueyu shi
"""

from collections import defaultdict 
import numpy as np
import re,os
data = defaultdict(dict)


for n in range(10, 55, 10):
    for r in range(10):
        mibs_log = os.path.join("Log/Log_KIP", "DNEG_n{}k{}".format(str(n), str(r)))
    
        for j in range(10):
            log_name = mibs_log + 'j{}_mibs.log'.format(str(j))
            data[(n, r)][j] = {"mibs": np.array([0, 0, 0], dtype='float')}
            with open(log_name, 'r') as f:
                for line in f.readlines():
                    if "Cost" in line:
                        data[(n, r)][j]["mibs"][0] = float(line.strip().split()[-1])
                    elif "optimality gap" in line:
                        data[(n, r)][j]["mibs"][1] = 1.0 + float(line.strip().split()[-1][:-1]) / 100.0
                        
            data[(n, r)][j]["mibs"][1] = data[(n, r)][j]["mibs"][0] / data[(n, r)][j]["mibs"][1]
    


with open("DNeg_results_mibs.log", 'r') as f:
    for line in f.readlines():
        
        line = line.strip().split()
        ins = line[0] 
        
        n = int(re.search('n[0-9]+', ins).group(0)[1:])
        r = int(re.search('k[0-9]+', ins).group(0)[1:])
        j = int(re.search('j[0-9]+', ins).group(0)[1:])
        data[(n,r)][j]["mibs"][-1] = min(600, float(line[-1]))




with open("DNeg_results.txt", 'r') as f:
    line_num = 0
    for line in f.readlines():
        line_num += 1
        if line_num <= 2: continue
        
        line = line.strip().split()
        ins = line[0]
        
        n = int(re.search('n[0-9]+', ins).group(0)[1:])
        r = int(re.search('k[0-9]+', ins).group(0)[1:])
        j = int(re.search('j[0-9]+', ins).group(0)[1:])
        
        # data[(n,r)][j]['CCLW'] = np.array([float(x) for x in line[1:3]])
        data[(n,r)][j][0] = np.array([float(x) for x in line[5:8]])
        for k in range(1,4):
            data[(n,r)][j][k] = np.array([float(x) for x in line[4+k*6:8+6*k]])
        
        base_val = data[(n,r)][j]['mibs'][0]
        
        for k in ['mibs'] + list(range(4)):
            data[(n,r)][j][k][0] = (data[(n,r)][j][k][0]+1) / (base_val + 1)
            data[(n,r)][j][k][1] = (data[(n,r)][j][k][1]+1) / (base_val + 1)
            
        

avg_data = {'hard': {}, 'easy': {}}

for key, val in data.items():
    n, r = key
    avg_val = [0]*17
    num_ins = len(val)
    for _, ins in val.items():
        avg_val[0:2] += ins['mibs'][1:] / num_ins
        avg_val[2:5] += ins[0] / num_ins
        
        for k in range(1, 4):
            avg_val[1+4*k:5+4*k] += ins[k] / num_ins

    if int(r) < 5:
        
        avg_data['hard'][(n, r+1)] = [round(x, 2) for x in avg_val]
    else:
        avg_data['easy'][(n, r+1)] = [round(x, 2) for x in avg_val]


import pandas as pd
with pd.ExcelWriter('DNEG_avg_results_new2.xlsx', engine='xlsxwriter') as writer:
    for key, val in avg_data.items():
        df = pd.DataFrame.from_dict(val, orient='index')
        df.to_excel(writer, sheet_name=key)





        
        
        

