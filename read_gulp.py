import sys
import argparse


class Info:
    def __init__(self, file, catg):
        self.file = file
        self.catg = catg

    # def get_jfinish(self, str):
    #     self.jfinish = str

    # def get_sw_criterion(self, str):
    #     self.sw_criterion = str

    # def get_fenergy(self, energy):
    #     self.energy = energy

    # def get_fgnorm(self, gnorm):
    #     self.gnorm = gnorm

    # def get_opti_time(self, time):
    # 	self.opti_time = time

    # def get_mem(self, mem):
    # 	self.mem = mem

    # def get_cpu_time(self, time):
    # 	self.cpu_time = time


if __name__ == "__main__":
    ''' Get input file and method to use from user '''
    parser = argparse.ArgumentParser(
        description='Define input')
    parser.add_argument(
        'ifilename', metavar='--input', type=str,
        help='.got file to read')
    args = parser.parse_args()
    try:
        with open(args.ifilename, 'r') as file:
            info = Info(file, {})

            c_cnt = 0
            H_cnt = 0
            for line in file:
                if "Cycle" in line:
                    c_cnt += 1  # Counts also Cycle 0
                elif "Hessian calculated" in line:
                    H_cnt += 1
                elif "Final energy" in line:
                    info.catg['fenergy'] = float(line.split(" ")[-2])
                elif "Final Gnorm" in line:
                    info.catg['gnorm'] = float(
                        line.split(" ")[-1].rstrip('\n'))
                elif "Time to end of optimisation" in line:
                    info.catg['opt_time'] = float(line.split(" ")[-2])
                elif "Peak dynamic memory used" in line:
                    info.catg['peak_mem'] = float(line.split(" ")[-3])
                elif "Total CPU time" in line:
                    info.catg['cpu_time'] = float(line.split(" ")[-1].rstrip('\n'))
                elif "Job Finished at" in line:
                      info.catg['jfinish'] = (line.rstrip('\n'))
                elif "Minimiser to switch" in line:
                    line = line + file.readline()
                    info.catg['sw_criterion'] = line

            print(info.catg)

    finally:
        file.close()
