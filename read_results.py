import os
import sys
import math
import mmap
import argparse
import collections
import pandas as pd

# ''' Merge csv files '''
# list = []
# dirs = [d for d in os.listdir('.') if os.path.isdir(os.path.join('.',d))]
# for d in dirs:
# 	list += [os.path.join(d,file) for file in os.listdir(d) if file.endswith(".csv")]
# # print(list)

# fout=open("results.csv","a")

# # first file:
# for line in open(list[0]):
#     fout.write(line)
# # now the rest:
# for num in range(1,len(list)):
#     f = open(list[num])
#     f.__next__() # skip the header
#     for line in f:
#          fout.write(line)
#     f.close() # not really needed


''' Process results '''
fields = ['structure', 'method', 'opt_time']

df = pd.read_csv('test0/results.csv', skipinitialspace=True, usecols=fields)


# for method in list(df['method'].unique()):
# 	rows = df.loc[df['method'] == method]['opt_time']
# 	times = [t for t in rows if not math.isnan(t)]
# 	avg = sum(times) / len(times)
# 	print(avg)

# fout.close()
