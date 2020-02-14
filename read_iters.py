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
		help='Define test environment folder')
	args = parser.parse_args()

	flist = []
	dirs = [d for d in os.listdir(args.test_dir) # find all directory objects
			if os.path.isdir(os.path.join(args.test_dir,d))] # check if is directory
	for d in dirs: # for every method folder
		d += '/output'
		if not os.path.exists(os.path.join(args.test_dir,d)):
			continue
		path = args.test_dir+'/'+d
		sampledirs = [d_ for d_ in os.listdir(path) # find all directory objects
			if os.path.isdir(os.path.join(path,d_))] # check if is directory
		for r in sampledirs: # rattled or random
			rpath = os.path.join(path,r)
			flist += [os.path.join(rpath,file) # list all files
						for file in os.listdir(rpath) if file.endswith(".got")]

	count=0
	ofilename = args.test_dir+'/'+args.ofilename
	# Get energies and gnorms for all files 
	#	in method directory in one dataframe
	for filename in flist:
		count_stepc=0 # count step size iters
		with open(filename, 'r') as file:
			info = Info(file, {})
			info.catg['structure'] = [
				filename.split('/')[-1].split('.')[0]]
			info.catg['method'] = filename.split('/')[-4]
			info.catg['folder'] = filename.split('/')[-2]
			df = pd.DataFrame.from_dict(info.catg, orient='columns')

			Es = [] # lists for each file
			Gns = []
			Steps = []

			for line in file:
				if "new    ..." in line:
					line_ = line.split('...')[1].rstrip().lstrip(' ').split(' ') # remove first part of line
					nline_ = list(filter(None, line_)) # remove blanks
					step = nline_[0]
					energy = nline_[-1]
					if ("*" in energy[1:]):
						step = energy.split('*')[0]
						energy = -math.inf
					elif ("-" in energy[1:]):
						both = energy.lstrip('-')
						step = both.split('-')[0]
						if energy[0] == "-":
							step = "-"+step
						energy = "-"+both.split('-')[-1]

					Es.append(energy)
					Steps.append(step)
					count_stepc += 1

				if "Cycle" in line:
					energy = line.split(':')[2].split(' ')[-3]
					gnorm = line.split(':')[3].split(' ')[-3]
					# if "**" not in energy:
						# Es.append(energy)
					if "**" not in gnorm:
						Gns.append(gnorm)
					count_stepc=0 # step is stabilized

			dfe = pd.DataFrame(Es).T
			dfe_ = df.join(dfe)
			dfe_ = dfe_.set_index(['structure', 'method'])

			dfg = pd.DataFrame(Gns).T
			dfg_ = df.join(dfg)
			dfg_ = dfg_.set_index(['structure', 'method'])

			dfs = pd.DataFrame(Steps).T
			dfs_ = df.join(dfs)
			dfs_ = dfs_.set_index(['structure', 'method'])

			# dfsi = pd.DataFrame(Stepi).T
			# dfsi_ = df.join(dfsi)
			# dfsi_ = dfsi_.set_index(['structure', 'method'])

		''' Merge dataframes '''
		if count:
			dfes = pd.concat([dfes,dfe_], axis=0, sort=False)
			dfgs = pd.concat([dfgs,dfg_], axis=0, sort=False)
			dfss = pd.concat([dfss,dfs_], axis=0, sort=False)
			# dfsis = pd.concat([dfsis,dfsi_], axis=0, sort=False)												
		else: # initialise
			dfes = dfe_
			dfgs = dfg_
			dfss = dfs_
			# dfsis = dfsi_
		count += 1
	for d in dirs:
		with open(args.test_dir+'/'+args.ofilename+'_energy.csv', 'w') as f:
			dfes.to_csv(f, header=True)

		with open(args.test_dir+'/'+args.ofilename+'_gnorm.csv', 'w') as f:
			dfgs.to_csv(f, header=True)

		with open(args.test_dir+'/'+args.ofilename+'_step.csv', 'w') as f:
			dfss.to_csv(f, header=True)

		# with open(args.test_dir+'/'+args.ofilename+'_stepi.csv', 'w') as f:
		# 	dfsis.to_csv(f, header=True)
