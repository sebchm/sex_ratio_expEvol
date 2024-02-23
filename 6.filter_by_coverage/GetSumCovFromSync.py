# using python scripts the sync file was sampled every 10,000 position, based on this the sync file was filtered to contained positions that were within a range from the expected coverage, the pileup file was then filtered to the same positions. The first script was then repeated using the newly filtered sync file.

# usage:
        # awk 'NR % 10000  == 0' sync_file.sync | python this_script.py > SumShuf.OneEverytenthous.txt


import fileinput

def SumCell(s):
        r = s.split(":")[:4]
        results = sum(map(int, r))
        return str(results)

for line in fileinput.input():
        a = line.strip().split()
        a1 = a[:3]
        a2 = a[3:]
        r = list(map(SumCell, a2))
        print("\t".join(a1 + r))
