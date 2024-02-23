# usage: python <this_script.py> <syc_file_for_filtering> <Tlower> <Tupper> > output_filtered.sync
import fileinput
import sys

def SumCell(s):
    r = s.split(":")[:4]
    results = sum(map(int, r))
    return results

# Tlower: lower coverage filter, provided as a third argument in command line; similarly for Tupper (forth argument)
Tlower = int(sys.argv[2])
Tupper = int(sys.argv[3])

for line in fileinput.input(sys.argv[1]):
    a = line.strip().split()
    a1 = a[:3] 
    a2 = a[3:]
    r = list(map(SumCell, a2))
    AvCov = float(sum(r))/len(r)
    if Tlower < AvCov:
    	if AvCov > Tupper:
        print(line.strip())
