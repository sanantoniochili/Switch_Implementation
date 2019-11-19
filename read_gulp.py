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
		'ifilename', metavar='--input', type=str,
		help='.got file to read')
	parser.add_argument(
		'ofilename', metavar='--output', type=str,
		help='.csv file to produce')
	parser.add_argument(
		'method', metavar='--method', type=str,
		help='Method used')
	parser.add_argument(
		'-c', '--categories',
		action='store_true',
		help='Add columns to file')
	args = parser.parse_args()

	# Initialisation
	switch_flag = False  	# check if method switches
	gnorm_flag  = False  	# check if gnorm value is the criterion
	switch      = ""  		# criterion in case there is method switch
	fail 		= ""  		# reason of failure
	error 		= ""  		# look for errors
	c_cnt 		= 0  		# count iterations
	H_cnt 		= 0  		# count Hessian calculations
	switch 		= -1  		# count iterations until switching method
	options 	= ""  		# add options used

	with open('temp.txt', 'r') as temp:
		for line in temp:
			options += line
			if "switch_minimiser" in line:
				switch_flag = True  # will switch method
				mod, value = line.split(' ')[-2], float(line.split(' ')[-1])
				if mod == "gnorm":  # switch based on gnorm
					gnorm_flag = True

	with open(args.ifilename, 'r') as file:
		info = Info(file, {})
		info.catg['structure'] = [
			args.ifilename.split('/')[-1].split('.')[0]]
		info.catg['method'] = args.method

		# Initialisation
		info.catg['energy'] 	= [math.inf]
		info.catg['gnorm'] 		= [0]
		info.catg['opt_time'] 	= [0]
		info.catg['peak_mem'] 	= [0]
		info.catg['cpu_time'] 	= [0]

		for line in file:
			if "Cycle" in line:
				# no. of iterations without cycle 0
				c_cnt += 1
				if switch_flag and gnorm_flag:  # check if method switched
					str = line.split(' ')[-7].split(":")[-1]  # with gnorm value
					if "**" not in str and float(str) < value:
						# check if has been error and if method switched
						switch = c_cnt-2  # remove cycle 0 and current
						gnorm_flag = False # count only first time
			elif "Hessian calculated" in line:
				# no. of Hessian calculations
				H_cnt += 1
			elif " **** " in line:
				if "Optimisation achieved" in line:
					# Check if optimisation succeeded
					info.catg['opt_succ'] = [True]
				else:
					# Check reason of failure
					info.catg['opt_succ'] = [False]
					fail = line
					line = file.readline()
					while " **** " in line:
						fail += line
						line = file.readline()
			elif "ERROR" in line:
				# Look for errors
				error += line
			elif "Final energy" in line:  # Optimised energy
				if "**" not in line:
					info.catg['energy'] = [float(line.split(" ")[-2].split(":")[-1])]
			elif "Final Gnorm" in line:  # Final gradient norm
				if "**" not in line:
					info.catg['gnorm'] = [float(
						line.split(" ")[-1].split(":")[-1].rstrip('\n'))]
			elif "Time to end of optimisation" in line:
				# Optimisation duration
				if "**" not in line:
					info.catg['opt_time'] = [float(line.split(" ")[-2])]
			elif "Peak dynamic memory used" in line:  # Most memory used
				if "**" not in line:
					info.catg['peak_mem'] = [float(line.split(" ")[-3])]
			elif "Total CPU time" in line:  # Total CPU time
				if "**" not in line:
					info.catg['cpu_time'] = [float(
						line.split(" ")[-1].rstrip('\n'))]

		info.catg['cycles'] 	= [c_cnt-1]  # Remove Cycle 0
		info.catg['hessian'] 	= [H_cnt]
		info.catg['failure'] 	= [fail]
		info.catg['options'] 	= [options]
		info.catg['switch'] 	= [switch]

		df = pd.DataFrame.from_dict(info.catg, orient='columns')
		df = df.set_index(['structure', 'method'])

	if args.categories:
		try:
			with open(args.ofilename, 'w') as f:
				df.to_csv(f, header=True)
		finally:
			f.close()
	else:
		try:
			with open(args.ofilename, 'a') as f:
				df.to_csv(f, header=False)
		finally:
			f.close()
