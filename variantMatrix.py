#! /home/jody/software/anaconda2/bin/python
from __future__ import division
import sys
from collections import defaultdict
from tqdm import tqdm
import argparse
from scipy.stats import chisquare
from subprocess import Popen,PIPE

######### global variables #########
allele_start_col = 5


######### functions #########
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def bed2set(bedFile):
	bed = defaultdict(set)
	for i in [x.rstrip().split() for x in open(bedFile).readlines()]:
		for j in range(int(i[1]),int(i[2])+1):
			bed[i[0]].add(str(j))
	return bed

def loadSampleFile(sampleFile):
	with open(sampleFile,"r") as f:
		return map(lambda x:x.rstrip(),f.readlines())

def loadMeta(metaFile):
	temp = {}
	for i in [x.rstrip().split() for x in open(metaFile).readlines()]:
		temp[i[0]] = i[1]
	return temp

def chisq(tab):
	
	num_samp_mut = tab["1-0"]+tab["1-1"]
	num_samp_wt = tab["0-0"]+tab["0-1"]
	num_samp_m1 = tab["0-1"]+tab["1-1"]
	num_samp_m0 = tab["0-0"]+tab["1-0"]
	pct_samp_m1 = num_samp_m1/(num_samp_m1+num_samp_m0)
	pct_samp_m0 = 1-pct_samp_m1
	exp = {}
	exp["0-0"] = num_samp_m0-(num_samp_mut*pct_samp_m0)
	exp["0-1"] = num_samp_m1-(num_samp_mut*pct_samp_m1)
	exp["1-0"] = num_samp_wt-(num_samp_wt*pct_samp_m0)
	exp["1-1"] = num_samp_mut-(num_samp_wt*pct_samp_m1)
	if num_samp_mut==0:
		return 1.0
	elif num_samp_wt==0:
		return 1.0
        elif num_samp_wt/(num_samp_wt+num_samp_mut)<0.05:
		return 1.0
#	print "%s\t%s" % (tab["0-0"],tab["0-1"])
#        print "%s\t%s" % (tab["1-0"],tab["1-1"])
	
	csq_pval = chisquare([tab["0-0"],tab["1-0"],tab["0-1"],tab["1-1"]],f_exp=[exp["0-0"],exp["1-0"],exp["0-1"],exp["1-1"]]).pvalue
	return csq_pval

####### main functions ########
def main_filter(args):
	bed = bed2set(args.bed_file)
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			o.write(f.readline())
			for i in tqdm(range(file_len(args.matrix)-1)):
				l = f.readline()
				arr = l.rstrip().split()
				if arr[1] in bed[arr[0]]:
					o.write(l)

def main_subset(args):
	samples = set(loadSampleFile(args.sample_file))
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			l = f.readline()
			header = l.rstrip().split()
			o.write(l)
			for _ in tqdm(range(file_len(args.matrix)-1)):
				l = f.readline()
				arr = l.rstrip().split()
				new_arr = []
				for i in range(allele_start_col,len(arr)):
					if header[i] in samples:
						new_arr.append(arr[i])
				o.write("\t".join(arr[:allele_start_col])+"\t"+"\t".join(new_arr)+"\n")

def main_chisq(args):
	meta_dict = loadMeta(args.meta_file)
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			header =  f.readline().rstrip().split()
			for i in tqdm(range(file_len(args.matrix)-1)):
				arr = f.readline().rstrip().split()
				table_dict = defaultdict(int)
				table_dict["0-0"] = 0
				table_dict["0-1"] = 0
				table_dict["1-0"] = 0
				table_dict["1-1"] = 0

				for j in range(allele_start_col,len(arr)):
					table_dict[arr[j]+"-"+meta_dict[header[j]]] += 1
				pval = chisq(table_dict)
				o.write("%s\t%s\t%s\n" % (arr[0],arr[1],pval))
			
def main_var2tab(args):
	chrom,pos = args.pos.split(":")
	meta_dict = loadMeta(args.meta_file)

	with open(args.matrix) as f:
		header =  f.readline().rstrip().split()
		for i in (range(file_len(args.matrix)-1)):
			arr = f.readline().rstrip().split()
			if arr[0]==chrom and arr[1]==pos:
				table_dict = defaultdict(int)
				table_dict["0-0"] = 0
				table_dict["0-1"] = 0
				table_dict["1-0"] = 0
				table_dict["1-1"] = 0
				for j in range(allele_start_col,len(arr)):
					table_dict[arr[j]+"-"+meta_dict[header[j]]] += 1
				if args.verbose:
					print table_dict
				if args.pretty:
					print "\n"
					print "|"+args.pos+(20-len(args.pos))*" "+"|"+"Sensitive"+(10-len("Sensitive"))*" "+"|"+"Resistant"+(10-len("Resistant"))*" "+"|"
					print "-"*44
					print "|"+"Wild-Type"+(20-len("Wild-Type"))*" "+"|"+str(table_dict["0-0"])+(10-len(str(table_dict["0-0"])))*" "+"|"+str(table_dict["0-1"])+(10-len(str(table_dict["0-1"])))*" "+"|"
					print "-"*44
					print "|"+"Mutant"+(20-len("Mutant"))*" "+"|"+str(table_dict["1-0"])+(10-len(str(table_dict["1-0"])))*" "+"|"+str(table_dict["1-1"])+(10-len(str(table_dict["1-1"])))*" "+"|"
					print "\n"
				else:
					print "%s\t%s" % (table_dict["0-0"],table_dict["0-1"])
				 	print "%s\t%s" % (table_dict["1-0"],table_dict["1-1"])
				quit()


parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('filter', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('bed_file',help='BED File')
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.set_defaults(func=main_filter)

parser_sub = subparsers.add_parser('subset', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('sample_file',help='Sample File')
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.set_defaults(func=main_filter)

parser_sub = subparsers.add_parser('chisq', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('meta_file',help='Sample File')
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.set_defaults(func=main_chisq)

parser_sub = subparsers.add_parser('pos2tab', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('meta_file',help='Sample File')
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('pos',help='Output File')
parser_sub.add_argument('--pretty',action='store_true')
parser_sub.add_argument('--verbose',action='store_true')
parser_sub.set_defaults(func=main_var2tab)


args = parser.parse_args()
args.func(args)

