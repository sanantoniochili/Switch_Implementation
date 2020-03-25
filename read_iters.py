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
	parser.add_argument(
		'-e', '--energy', action='store_true',
		help='Records of final energy')
	parser.add_argument(
		'-g', '--gnorm', action='store_true',
		help='Records of final gradient norm')
	parser.add_argument(
		'-s', '--step', action='store_true',
		help='Records of step sizes')
	parser.add_argument(
		'-i', '--interatomic', action='store_true',
		help='Records of interatomic energy values')
	parser.add_argument(
		'--all', action='store_true',
		help='Records of all values')
	args = parser.parse_args()

	if args.all:
		args.energy = True
		args.gnorm  = True
		args.step   = True
		args.interatomic  = True

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
			Inters = []
			iflag = True # keep only first interatomic potential

			for line in file:
				if args.step: # keeping records of step sizes and respective energy levels
					if "new    ..." in line:
						line_ = line.split('...')[1].rstrip().lstrip(' ').split(' ') # remove first part of line
						nline_ = list(filter(None, line_)) # remove blanks
						step = nline_[0]
						energy = nline_[-1]
						if ("*" in energy[1:]):
							step = energy.split('*')[0]
							energy = -math.inf
						elif ("-" in energy[1:]): # take care of overlapping values
							both = energy.lstrip('-')
							step = both.split('-')[0]
							if energy[0] == "-":
								step = "-"+step
							energy = "-"+both.split('-')[-1]

						Es.append(energy)
						Steps.append(step)
						count_stepc += 1

				if "Cycle" in line:
					if args.energy:
						energy = line.split(':')[2].split(' ')[-3]
					if args.gnorm:
						gnorm = line.split(':')[3].split(' ')[-3]
					# if args.energy:
						# if "**" not in energy:
							# Es.append(energy)
					if args.gnorm:
						if "**" not in gnorm:
							Gns.append(gnorm)
					count_stepc=0 # step is stabilized

				if "Interatomic potentials     =" in line:
					if args.interatomic & iflag:
						Inters.append(line.split()[-2])
						iflag = False

			if args.energy:
				dfe = pd.DataFrame(Es).T
				dfe_ = df.join(dfe)
				dfe_ = dfe_.set_index(['structure', 'method'])

			if args.gnorm:
				dfg = pd.DataFrame(Gns).T
				dfg_ = df.join(dfg)
				dfg_ = dfg_.set_index(['structure', 'method'])

			if args.step:
				dfs = pd.DataFrame(Steps).T
				dfs_ = df.join(dfs)
				dfs_ = dfs_.set_index(['structure', 'method'])

			if args.interatomic:
				dfsei = pd.DataFrame(Inters).T
				dfsei_ = df.join(dfsei)
				dfsei_ = dfsei_.set_index(['structure', 'method'])

			# dfsi = pd.DataFrame(Stepi).T
			# dfsi_ = df.join(dfsi)
			# dfsi_ = dfsi_.set_index(['structure', 'method'])

		''' Merge dataframes '''
		if count:
			if args.energy:
				dfes = pd.concat([dfes,dfe_], axis=0, sort=False)
			if args.gnorm:
				dfgs = pd.concat([dfgs,dfg_], axis=0, sort=False)
			if args.step:
				dfss = pd.concat([dfss,dfs_], axis=0, sort=False)
			if args.interatomic:
				dfseis = pd.concat([dfseis,dfsei_], axis=0, sort=False)
			# dfsis = pd.concat([dfsis,dfsi_], axis=0, sort=False)												
		else: # initialise
			if args.energy:
				dfes = dfe_
			if args.gnorm:
				dfgs = dfg_
			if args.step:
				dfss = dfs_
			if args.interatomic:
				dfseis = dfsei_
			# dfsis = dfsi_
		count += 1

	# for d in dirs:
	if args.energy:
		with open(args.test_dir+'/'+args.ofilename+'_energy.csv', 'w') as f:
			dfes.to_csv(f, header=True)

	if args.gnorm:
		with open(args.test_dir+'/'+args.ofilename+'_gnorm.csv', 'w') as f:
			dfgs.to_csv(f, header=True)

	if args.step:
		with open(args.test_dir+'/'+args.ofilename+'_step.csv', 'w') as f:
			dfss.to_csv(f, header=True)

	if args.interatomic:
		with open(args.test_dir+'/'+args.ofilename+'_interatomic.csv', 'w') as f:
			dfseis.to_csv(f, header=True)

		# with open(args.test_dir+'/'+args.ofilename+'_stepi.csv', 'w') as f:
		# 	dfsis.to_csv(f, header=True)
