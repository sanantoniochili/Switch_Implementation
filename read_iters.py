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
	parser.add_argument(
		'ofilename', metavar='--output', type=str,
		help='.csv file to produce')
	parser.add_argument(
		'test_dir', metavar='--test_folder', type=str,
		help='Define folder with samples')
	args = parser.parse_args()

	flist = []
	dirs = [d for d in os.listdir(args.test_dir) # find all objects
			if os.path.isdir(os.path.join(args.test_dir,d))] # check if is directory
	for d in dirs:
		d += '/output'
		if not os.path.exists(os.path.join(args.test_dir,d)):
			continue
		d = os.path.abspath(args.test_dir+'/'+d)
		flist += [os.path.join(d,file) for file in os.listdir(d) if file.endswith(".got")]

	count=0
	ofilename = args.test_dir+'/'+args.ofilename
	# Get energies and gnorms for all files 
	#	in method directory in one dataframe
	for filename in flist:
		with open(filename, 'r') as file:
			info = Info(file, {})
			info.catg['structure'] = [
				filename.split('/')[-1].split('.')[0]]
			info.catg['method'] = filename.split('/')[-3]
			df = pd.DataFrame.from_dict(info.catg, orient='columns')

			Es = [] # lists for each file
			Gns = []

			for line in file:
				if "Cycle" in line:
					energy = line.split(':')[2].split(' ')[-3]
					gnorm = line.split(':')[3].split(' ')[-3]
					if "**" not in energy:
						Es.append(energy)
					if "**" not in gnorm:
						Gns.append(gnorm)	

			dfe = pd.DataFrame(Es).T
			dfe_ = df.join(dfe)
			dfe_ = dfe_.set_index(['structure', 'method'])

			dfg = pd.DataFrame(Gns).T
			dfg_ = df.join(dfg)
			dfg_ = dfg_.set_index(['structure', 'method'])

		''' Merge dataframes '''
		if count:
			dfes = pd.concat([dfes,dfe_], axis=0, sort=False)
			dfgs = pd.concat([dfgs,dfg_], axis=0, sort=False)			
		else:
			dfes = dfe_
			dfgs = dfg_
		count += 1
	print(count)
	try:
		with open(args.test_dir+'/'+args.ofilename+'_energy.csv', 'w') as f:
			dfes.to_csv(f, header=True)
	finally:
		f.close()

	try:
		with open(args.test_dir+'/'+args.ofilename+'_gnorm.csv', 'w') as f:
			dfgs.to_csv(f, header=True)
	finally:
		f.close()
