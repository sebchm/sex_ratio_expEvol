# filter positions based on mean coverage
# adjust Tlower and Tupper
# usage: python <this_script.py> <syc_file_for_filtering> <Tlower> <Tupper> > output_filtered.sync
# author: Mateusz Konczal, Evolutionary Biology group, mateusz.konczal(_a_)amu.edu.pl
import fileinput
import sys

def SumCell(s):
    r = s.split(":")[:4]
    results = sum(map(int, r))
    return results

# Tlower: lower coverage filter, provided as a third argument in command line; similarly for Tupper (forth argument)
Tlower = int(sys.argv[2])
Tupper = int(sys.argv[3])

for line in fileinput.input():
        a = line.strip().split()
        a1 = a[:3]
        a2 = a[3:]
        r = list(map(SumCell, a2))
        AvCov = float(sum(r))/len(r)
        if AvCov < Tupper:
                if AvCov > Tlower:
                        print(line.strip())
