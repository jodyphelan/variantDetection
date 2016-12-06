import sys
from collections import defaultdict
from tqdm import tqdm
import argparse


######### functions #########
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def loadBed(bedFile):
	bed = defaultdict(set)
	for i in [x.rstrip().split() for x in open(bedFile).readlines()]:
		for j in range(int(i[1]),int(i[2])+1):
			bed[i[0]].add(str(j))
	return bedFile

####### main functions ########
def main_filter(args):
	bed = loadBed(args.bed_file)
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			o.write(f.readline())
			for i in tqdm(range(file_len(sys.argv[2])-1)):
				l = f.readline()
				arr = l.rstrip().split()
				if arr[1] in bed[arr[0]]:
					o.write(l)



parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('filter', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('bed_file',help='BED File')
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.set_defaults(func=main_filter)

args = parser.parse_args()
args.func(args)

