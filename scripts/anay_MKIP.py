#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the multi-knapsack interdiction results 

@author: Xueyu shi
"""

from collections import defaultdict 
import numpy as np
import re,os
data = defaultdict(dict)


mibs_result = "MKIP_results_mibs_m3.log"
bpk_result = "MKIP_results_m3.txt"
write_file = 'MKIP_avg_results_m3.xlsx'


m1 = 1
m2 = 3
for n in range(10, 55, 10):
    for r in range(2, 7):
        mibs_log = os.path.join("Log_MKIP", "MKIP_n{}k{}".format(str(n), str(r)))
    
        for j in range(10):
            log_name = mibs_log + 'j{}_m1_{}_m2_{}_mibs.log'.format(str(j), str(m1), str(m2))
            data[(n, r)][j] = {"mibs": np.array([0, 0, 0], dtype='float')}
            with open(log_name, 'r') as f:
                for line in f.readlines():
                    if "Cost" in line:
                        data[(n, r)][j]["mibs"][0] = float(line.strip().split()[-1])
                    elif "optimality gap" in line:
                        data[(n, r)][j]["mibs"][1] = 1. + float(line.strip().split()[-1][:-1]) / 100

            data[(n, r)][j]["mibs"][1] = data[(n, r)][j]["mibs"][0] / data[(n, r)][j]["mibs"][1]
    

with open(mibs_result, 'r') as f:
    for line in f.readlines():
        
        line = line.strip().split()
        ins = line[0] 
        
        n = int(re.search('n[0-9]+', ins).group(0)[1:])
        r = int(re.search('k[0-9]+', ins).group(0)[1:])
        j = int(re.search('j[0-9]+', ins).group(0)[1:])
        if (n, r) not in data: continue
        
        data[(n,r)][j]["mibs"][-1] = min(600, float(line[-1]))


with open(bpk_result, 'r') as f:
    line_num = 0
    for line in f.readlines():
        line_num += 1
        
        line = line.strip().split()
        ins = line[0]
        
        n = int(re.search('n[0-9]+', ins).group(0)[1:])
        r = int(re.search('k[0-9]+', ins).group(0)[1:])
        j = int(re.search('j[0-9]+', ins).group(0)[1:])
        
        # data[(n,r)][j]['CCLW'] = np.array([float(x) for x in line[1:3]])
        data[(n,r)][j][0] = np.array([float(x) for x in line[1:5]])
        for k in range(1,4):
            data[(n,r)][j][k] = np.array([float(x) for x in line[k*8-3:8*k+1]])
            data[(n,r)][j][str(k) + 'ext'] = np.array([float(x) for x in line[k*8+1:8*k+5]])
        
        base_val = data[(n,r)][j]['mibs'][0]
        data[(n,r)][j]['mibs'][1] = (data[(n,r)][j]['mibs'][1]+1) / (base_val + 1)
        
        for k in list(range(4)) + [ str(k) + 'ext' for k in range(1,4) ]:
            if  len(data[(n,r)][j][k]) == 0:
                data[(n,r)][j][k] = data[(n,r)][j][3]
                continue
            data[(n,r)][j][k][2] = (data[(n,r)][j][k][2]+1) / (base_val + 1)
            data[(n,r)][j][k][1] = (data[(n,r)][j][k][1]+1) / (base_val + 1)
            
            data[(n,r)][j][k][-1] = min(600, data[(n,r)][j][k][-1])
            
        

avg_data = {'hard': {}, 'easy': {}}

for key, val in data.items():
    n, r = key
    if r > 6 or r<2: continue
    avg_val = [0]*17
    num_ins = len(val)
    for _, ins in val.items():
        avg_val[0:2] += ins['mibs'][1:] / num_ins
        avg_val[2:5] += ins[0][1:] / num_ins
        
        for k in range(1, 4):
            if k == 3 and len(ins[k]) == 0: continue
            avg_val[4*k+1] += max(ins[k][1], ins[str(k) + 'ext'][1]) / num_ins
            avg_val[4*k+2] += min(ins[k][2], ins[str(k) + 'ext'][2]) / num_ins
            avg_val[4*k+3] += ins[k][3] / num_ins
            avg_val[4*k+4] += ins[str(k) + 'ext'][-1] / num_ins
            
   
    avg_data['hard'][(n, r+1)] = [round(x, 2) for x in avg_val]
    
import pandas as pd
with pd.ExcelWriter(write_file, engine='xlsxwriter') as writer:
    for key, val in avg_data.items():
        df = pd.DataFrame.from_dict(val, orient='index')
        df.to_excel(writer, sheet_name=key)





        
        
        

