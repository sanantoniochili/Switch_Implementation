import os
import sys
import math
import mmap
import argparse
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

''' Attach a text label above each bar in *rects*,
	displaying its height. '''


def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


''' Merge csv files '''

# list = []
# dirs = [d for d in os.listdir('.') if os.path.isdir(os.path.join('.',d))]
# for d in dirs:
# 	list += [os.path.join(d,file) for file in os.listdir(d) if file.endswith(".csv")]

# temp_list = [list[2], list[3]]

# fout=open("results.csv","a")
# fout=open("results_test.csv","a")

# # first file:
# for line in open(temp_list[0]):
#     fout.write(line)
# # now the rest:
# for num in range(1,len(list)):
#     f = open(list[num])
#     f.__next__() # skip the header
#     for line in f:
#          fout.write(line)
    # f.close() # not really needed


''' Process results '''

# ''' No. of optimised structures per method (same no. of total structures) '''

# ran_succ_list = []
# rat_succ_list = []

# df = pd.read_csv('results_test.csv', skipinitialspace=True)
# labels = list(df['method'].unique())
# dicts = {}

# for method in labels:
#     random = df[df['structure'].str.contains(
#         "rat") == False].loc[df['method'] == method]
#     rattled = df[df['structure'].str.contains(
#         "rat")].loc[df['method'] == method]

#     ran_succ = len([b for b in random['opt_succ'] if b])
#     rat_succ = len([b for b in rattled['opt_succ'] if b])

#     ran_succ_list.append(ran_succ)
#     rat_succ_list.append(rat_succ)

# total = len(df.loc[df['method'] == method]) / 2

# x = np.arange(len(labels))  # the label locations
# width = 0.35  # the width of the bars

# fig, ax = plt.subplots()
# ax.set_ylim(0,total+10) # set total samples as y limit (add space for view)
# rects1 = ax.bar(x - width/2, ran_succ_list, width, label='Random')
# rects2 = ax.bar(x + width/2, rat_succ_list, width, label='Rattled')

# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Structures')
# ax.set_title('Successfully optimised structures')
# ax.set_xticks(x)
# ax.set_xticklabels(labels)
# ax.legend()

# autolabel(rects1)
# autolabel(rects2)

# fig.tight_layout()

# plt.show()


''' No. of optimised structures per method (diff no. of total structures) '''

ran_succ_list = []
rat_succ_list = []

df = pd.read_csv('results_test.csv', skipinitialspace=True)
labels = list(df['method'].unique())
dicts = {}

for method in labels:
    random = df[df['structure'].str.contains(
        "rat") == False].loc[df['method'] == method]
    rattled = df[df['structure'].str.contains(
        "rat")].loc[df['method'] == method]

    ran_succ = len([b for b in random['opt_succ'] if b])
    rat_succ = len([b for b in rattled['opt_succ'] if b])

    ran_succ_list.append(ran_succ)
    rat_succ_list.append(rat_succ)

total = len(df.loc[df['method'] == method]) / 2

y = np.arange(1)  # the label locations
width = 0.35  # the width of the bars

fig, axs = plt.subplots(len(labels))
fig.suptitle('Successfully optimised structures')

count = 0
for ax in axs:
	# ax.set_ylim(0,total+10) # set total samples as y limit (add space for view)
	rects1 = ax.barh(y - width/2, ran_succ_list[count], width, label='Random')
	rects2 = ax.barh(y + width/2, rat_succ_list[count], width, label='Rattled')

	# Add some text for labels, title and custom x-axis tick labels, etc.
	# ax.set_xlabel(labels[count])
	ax.set_yticks(y)
	ax.set_yticklabels([labels[count]])
	ax.invert_yaxis()  # labels read top-to-bottom
	ax.legend()

	# autolabel(rects1)
	# autolabel(rects2)

	count += 1

# fig.tight_layout()

plt.show()