#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the bilevel minimum spanning tree results 

@author: Xueyu shi
"""


filename = "BMST_results.txt"

from collections import defaultdict 
import numpy as np
data = defaultdict(dict)
with open(filename, 'r') as f:
    line_num = 0
    for line in f.readlines():
        line_num += 1
        if line_num <= 2: continue
        line = line.strip().split()
        ins = line[0].split('_')
        
        n = ins[0]
        r = int(ins[2]) / 100
        j = int(ins[-1])
        
        if r not in data[n]:
            data[n][r] = defaultdict(dict)
        
        data[n][r][j] = np.array([float(line[5]), float(line[6]), float(line[9]), float(line[10])])
        
        if 'meta' not in data[n]:
            data[n]['meta'] = [int(line[1]), int(line[2])]
        
        
        
        

avg_data = {}

for key, val in data.items():
    avg_val = [0]*12
    
    for i, r in enumerate([0.05, 0.1, 0.15]):
        num_ins = len(val[r])
        for _, ins in val[r].items():
            avg_val[4*i:4*(i+1)] += ins / num_ins
    
    avg_data[key] = val['meta'] + [round(x, 2) for x in avg_val]


import pandas as pd
dataframe = pd.DataFrame.from_dict(avg_data, orient='index')
dataframe.to_excel('BMST_avg_results.xlsx')
        
        
        
