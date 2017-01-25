#! /home/jody/software/anaconda2/bin/python
from __future__ import division
import sys
from collections import defaultdict
from tqdm import tqdm
import argparse
from scipy.stats import chisquare,fisher_exact
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

def index_matrix(infile):
	import os
	import subprocess
	tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"
	gzipped = infile+".gz"
	tbi = infile+".gz.tbi"
	if not os.path.isfile(gzipped):
		subprocess.call("cat %s | bgzip -c >%s" % (infile,gzipped),shell=True)
	if not os.path.isfile(tbi):
		subprocess.call("%s -b 2 -e 2 -s 1 %s" % (tabix,gzipped),shell=True)

def extract_from_index(infile,chrom_pos):
	tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"
	cmd = Popen("%s %s %s -h" % (tabix,infile,chrom_pos),shell=True,stdout=PIPE)
	results = []
	for l in cmd.stdout:
		results.append(l)
	return results[0].rstrip().split(),results[1].rstrip().split()

def load_bed(bedfile):
	temp = defaultdict(lambda :defaultdict(str))
	loci = []
	starts = {}
	for l in open(bedfile):
		arr = l.rstrip().split()
		loci.append(arr[3])
		starts[arr[3]] = (arr[0],arr[1])
		for i in range(int(arr[1]),int(arr[2])):
			temp[arr[0]][str(i)] = arr[3]
	return temp,loci,starts

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
			samples_idx  = []
			new_arr = []
			for i in range(allele_start_col,len(header)):
				if header[i] in samples:
					samples_idx.append(i)
					new_arr.append(header[i])
			o.write("\t".join(header[:allele_start_col])+"\t"+"\t".join(new_arr)+"\n")
			for _ in tqdm(range(file_len(args.matrix)-1)):
				l = f.readline()
				arr = l.rstrip().split()
				new_arr = []
				for i in samples_idx:
					if header[i] in samples:
						new_arr.append(arr[i])
				o.write("\t".join(arr[:allele_start_col])+"\t"+"\t".join(new_arr)+"\n")


def main_binarise(args):
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			l = f.readline()
			header = l.rstrip().split()
			o.write(l)
			for _ in tqdm(range(file_len(args.matrix)-1)):
				l = f.readline()
				arr = l.rstrip().split()
				for i in range(allele_start_col,len(arr)):
					if arr[i]=="NA":
						continue
					elif arr[i]!=arr[2]:
						arr[i] = "1"
					else:	
						arr[i] = "0"
				o.write("\t".join(arr)+"\n")	

def main_snps2fasta(args):
	print "Loading Matrix"
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			l = f.readline()
			header = l.rstrip().split()
			fa_dict = defaultdict(list)
			for _ in tqdm(range(file_len(args.matrix)-1)):
				l = f.readline()
				arr = l.rstrip().split()
				if arr[4]!="SNP":
					continue
				for i in range(allele_start_col,len(arr)):
					if arr[i]=="NA":
						arr[i]="N"
					fa_dict[header[i]].append(arr[i][0]) # [0] is because sometimes you can have an indel in a "snp" position
				if args.ref:
					fa_dict["ref"].append(arr[2])
			print "Writing Fasta"
			for s in tqdm(header[allele_start_col:]):
				o.write(">%s\n"%s)
				o.write("".join(fa_dict[s])+"\n")
			if args.ref:
				o.write(">ref\n")
				o.write("".join(fa_dict["ref"])+"\n")

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
	index_matrix(args.matrix)
	chrom,pos = args.pos.split(":")
	chrom_pos = chrom+":"+pos+"-"+pos
	header,arr = extract_from_index(args.matrix+".gz",chrom_pos)
	meta_dict = loadMeta(args.meta_file)
	table_dict = defaultdict(int)
	table_dict["0-1"] = 0
	table_dict["1-0"] = 0
	table_dict["1-1"] = 0
	table_dict["1-NA"] = 0

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
	elif args.mutants:
		print "%s\t%s\t%s\t%s\t%s" % (arr[0],arr[1],table_dict["1-0"],table_dict["1-1"],table_dict["1-NA"])
	else:
		print "%s\t%s" % (table_dict["0-0"],table_dict["0-1"])
	 	print "%s\t%s" % (table_dict["1-0"],table_dict["1-1"])
	if args.samples:
		temp_arr = []
		for j in range(allele_start_col,len(arr)):
			if arr[j]=="1" and meta_dict[header[j]]=="1":
				temp_arr.append(header[j])
		print temp_arr
	if args.af:
		print "Allele frequency: %s" % (table_dict["1-1"]/(table_dict["0-1"]+table_dict["1-1"]))

def main_2by2(args):
	index_matrix(args.matrix)
	chrom1,pos1 = args.pos1.split(":")
	chrom_pos1 = chrom1+":"+pos1+"-"+pos1
	chrom2,pos2 = args.pos2.split(":")
	chrom_pos2 = chrom2+":"+pos2+"-"+pos2
	header,arr1 = extract_from_index(args.matrix+".gz",chrom_pos1)
	header,arr2 = extract_from_index(args.matrix+".gz",chrom_pos2)

	table_dict = defaultdict(int)
	table_dict["0-1"] = 0
	table_dict["1-0"] = 0
	table_dict["1-1"] = 0
	table_dict["1-NA"] = 0
	for j in range(allele_start_col,len(arr1)):
		table_dict[arr1[j]+"-"+arr2[j]] += 1
	if args.verbose:
		print table_dict
	if args.pretty:
		label = args.pos1+"-"+args.pos2
		print "\n"
		print "|"+label+(50-len(label))*" "+"|"+"0"+(10-len("0"))*" "+"|"+"1"+(10-len("1"))*" "+"|"
		print "-"*74
		print "|"+"0"+(50-len("0"))*" "+"|"+str(table_dict["0-0"])+(10-len(str(table_dict["0-0"])))*" "+"|"+str(table_dict["0-1"])+(10-len(str(table_dict["0-1"])))*" "+"|"
		print "-"*74
		print "|"+"1"+(50-len("1"))*" "+"|"+str(table_dict["1-0"])+(10-len(str(table_dict["1-0"])))*" "+"|"+str(table_dict["1-1"])+(10-len(str(table_dict["1-1"])))*" "+"|"
		print "\n"
	else:
		print "%s\t%s" % (table_dict["0-0"],table_dict["0-1"])
	 	print "%s\t%s" % (table_dict["1-0"],table_dict["1-1"])
	if args.fisher:
		print fisher_exact([[table_dict["0-0"],table_dict["0-1"]],[table_dict["1-0"],table_dict["1-1"]]])

def main_af(args):
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			header =  f.readline().rstrip().split()
			for _ in tqdm(range(file_len(args.matrix)-1)):
				arr = f.readline().rstrip().split()
				num_samps = 0
				for x in arr[5:]:
					if x=="1":
						num_samps+=1
				af = num_samps/(len(arr)-5)
				o.write("%s\t%s\t%s\t%s\t%s\n" % (arr[0],arr[1],arr[4],num_samps,af))
					

def main_locus_agg(args):
	bed,loci,starts = load_bed(args.bed_file)
	results = defaultdict(dict)
	with open(args.out_file,"w") as o:
		with open(args.matrix) as f:
			line = f.readline()
			header =  line.rstrip().split()
			o.write(line)
			for l in loci:
				for s in header[5:]:
					results[l][s] = "0"
			for _ in tqdm(range(file_len(args.matrix)-1)):
				line = f.readline()
				if args.concat:
					out.write(line)	
				arr = line.rstrip().split()
				chrom,pos = arr[:2]
				for i in range(5,len(arr)):
					if arr[i]=="1":
						results[bed[chrom][pos]][header[i]] = "1"
		for l in loci:
			o.write(starts[l][0]+"\t"+starts[l][1]+"\t.\t"+l+"\t"+"locus_sum"+"\t"+"\t".join([results[l][x] for x in header[5:]])+"\n")

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
parser_sub.set_defaults(func=main_subset)

parser_sub = subparsers.add_parser('binary', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.set_defaults(func=main_binarise)

parser_sub = subparsers.add_parser('locus_sum', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('bed_file',help='BED File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.add_argument('--concat',action='store_true')
parser_sub.set_defaults(func=main_locus_agg)


parser_sub = subparsers.add_parser('af', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.set_defaults(func=main_af)



parser_sub = subparsers.add_parser('snps2fa', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('out_file',help='Output File')
parser_sub.add_argument('--ref',action='store_true')
parser_sub.set_defaults(func=main_snps2fasta)


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
parser_sub.add_argument('--mutants',action='store_true')
parser_sub.add_argument('--verbose',action='store_true')
parser_sub.add_argument('--samples',action='store_true')
parser_sub.add_argument('--af',action='store_true')
parser_sub.set_defaults(func=main_var2tab)

parser_sub = subparsers.add_parser('2by2', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('matrix',help='Variant Matrix File')
parser_sub.add_argument('pos1',help='Output File')
parser_sub.add_argument('pos2',help='Output File')
parser_sub.add_argument('--pretty',action='store_true')
parser_sub.add_argument('--fisher',action='store_true')
parser_sub.add_argument('--verbose',action='store_true')
parser_sub.set_defaults(func=main_2by2)

args = parser.parse_args()
args.func(args)

