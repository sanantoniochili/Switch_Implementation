import os
import sys
import math
import mmap
import argparse
import collections
import pandas as pd

class Info:
	def __init__(self, file, catg):  # create dict
		self.file = file
		self.catg = catg

if __name__ == "__main__":
	''' Get input file and method to use from user '''
	parser = argparse.ArgumentParser(
		description='Define input')
	# parser.add_argument(
	# 	'ofilename', metavar='--output', type=str,
	# 	help='.csv file to produce')
	parser.add_argument(
		'ifilename', metavar='--input', type=str,
		help='.got file to read')
	parser.add_argument(
		'method', metavar='--method', type=str,
		help='Method used')
	parser.add_argument(
		'test_dir', metavar='--test_folder', type=str,
		help='Define folder with samples')
	args = parser.parse_args()

	test_dir = args.test_dir+'/'+args.method+'/output'
	flist = os.listdir(test_dir)	

	count=0
	ofilename = args.test_dir+args.ofilename
	''' Get energies and gnorms for all files 
		in method directory in one dataframe '''
	for filename in flist:
		with open(args.ifilename, 'r') as file:
			info = Info(file, {})
			info.catg['structure'] = [
				args.ifilename.split('/')[-1].split('.')[0]]
			info.catg['method'] = args.method

			df = pd.DataFrame.from_dict(info.catg, orient='columns')

			Es = [] # lists for each file
			Gns = []

			for line in file:
				if "Cycle" in line:
					energy = line.split(':')[2].split(' ')[-3]
					if "**" not in energy:
						Es.append(energy)	

			dfe = pd.DataFrame(Es).T
			df = df.join(dfe)
			df = df.set_index(['structure', 'method'])

		''' Merge dataframes '''
		if not count:
			dfes = 
		
		count += 1

		# 	try:
		# 		with open(args.ofilename, 'w') as f:
		# 			df.to_csv(f, header=True)
		# 	finally:
		# 		f.close()
		# else:
		# 	try:
		# 		with open(args.ofilename, 'a') as f:
		# 			df.to_csv(f, header=False)
		# 	finally:
		# 		f.close()
