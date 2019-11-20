import os
import sys
import math
import mmap
import argparse
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


''' Merge csv files '''
# list = []
# dirs = [d for d in os.listdir('.') if os.path.isdir(os.path.join('.',d))]
# for d in dirs:
# 	list += [os.path.join(d,file) for file in os.listdir(d) if file.endswith(".csv")]

# temp_list = [list[2], list[3]]

# fout=open("results.csv","a")

# # first file:
# for line in open(temp_list[0]):
#     fout.write(line)
# # now the rest:
# for num in range(1,len(temp_list)):
#     f = open(temp_list[num])
#     f.__next__() # skip the header
#     for line in f:
#          fout.write(line)
#     f.close() # not really needed


''' Process results '''
df   = pd.read_csv('results.csv', skipinitialspace=True)
# print(df[['structure', 'method', 'opt_succ', 'failure']])

labels = list(df['method'].unique())

dicts = {}

for method in labels:
	random  = df[df['structure'].str.contains("rat") == False].loc[df['method'] == method]
	rattled = df[df['structure'].str.contains("rat")].loc[df['method'] == method]

	ran_succ = len([b for b in random['opt_succ'] if b])
	rat_succ = len([b for b in rattled['opt_succ'] if b])

	dicts[method] = {'ran_total' : len(random), 'rat_total' : len(rattled), 
						'ran_succ' : ran_succ, 'rat_succ' : rat_succ}

# print(dicts)

# c_succ = [20, 34]
# women_means = [25, 32]

# x = np.arange(len(labels))  # the label locations
# width = 0.35  # the width of the bars

# fig, ax = plt.subplots()
# rects1 = ax.bar(x - width/2, men_means, width, label='Men')
# rects2 = ax.bar(x + width/2, women_means, width, label='Women')

# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
# ax.set_xticks(x)
# ax.set_xticklabels(labels)
# ax.legend()


# def autolabel(rects):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for rect in rects:
#         height = rect.get_height()
#         ax.annotate('{}'.format(height),
#                     xy=(rect.get_x() + rect.get_width() / 2, height),
#                     xytext=(0, 3),  # 3 points vertical offset
#                     textcoords="offset points",
#                     ha='center', va='bottom')


# autolabel(rects1)
# autolabel(rects2)

# fig.tight_layout()

# plt.show()

# for method in list(df['method'].unique()):
# 	rows = df.loc[df['method'] == method]['opt_time']
# 	times = [t for t in rows if not math.isnan(t)]
# 	avg = sum(times) / len(times)
# 	print(avg)

# fout.close()

# fig = plt.figure()  # an empty figure with no axes
# fig.suptitle('No axes on this figure')  # Add a title so we know which it is
