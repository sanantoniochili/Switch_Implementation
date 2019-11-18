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
	switch = ""  # criterion in case there is method switch
	fail = ""  # reason of failure
	c_cnt = 0  # count iterations
	H_cnt = 0  # count Hessian calculations
	options = ""  # add options used
	with open('temp.txt', 'r') as temp:
		for line in temp:
			options += line

	with open(args.ifilename, 'r') as file:
		info = Info(file, {})
		info.catg['structure'] = [
			args.ifilename.split('/')[-1].split('.')[0]]
		info.catg['method'] = args.method

		# Initialisation
		info.catg['energy'] = [math.inf]
		info.catg['gnorm'] = [0]
		info.catg['opt_time'] = [0]
		info.catg['peak_mem'] = [0]
		info.catg['cpu_time'] = [0]
	
		for line in file:
			if "Cycle" in line:
				# no. of iterations without cycle 0
				c_cnt += 1
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
			elif "Final energy" in line:  # Optimised energy
				info.catg['energy'] = [float(line.split(" ")[-2])]
			elif "Final Gnorm" in line:  # Final gradient norm
				info.catg['gnorm'] = [float(
					line.split(" ")[-1].rstrip('\n'))]
			elif "Time to end of optimisation" in line:
									# Optimisation duration
				info.catg['opt_time'] = [float(line.split(" ")[-2])]
			elif "Peak dynamic memory used" in line:  # Most memory used
				info.catg['peak_mem'] = [float(line.split(" ")[-3])]
			elif "Total CPU time" in line:  # Total CPU time
				info.catg['cpu_time'] = [float(
					line.split(" ")[-1].rstrip('\n'))]
			# elif "Minimiser to switch" in line:
			# 						# Switch of minimisers criterion
			# 	switch += line.rstrip('\n') + file.readline().rstrip('\n')

		info.catg['cycles'] = [c_cnt-1]  # Remove Cycle 0
		info.catg['hessian'] = [H_cnt]
		info.catg['failure'] = [fail]
		info.catg['options'] = [options]

		df = pd.DataFrame.from_dict(info.catg, orient='columns')
		df = df.set_index(['structure','method'])

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
