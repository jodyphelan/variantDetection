#! /home/jody/software/anaconda2/bin/python
import sys
from tqdm import tqdm
from psGeneral import file_len


def matrix_subset(infile,outfile,samples):
	with open(outfile,"w") as o:
		with open(infile) as f:
			header = f.readline().split()
			sample_idx_arr = []
			sample_list = []
			for i,s in enumerate(header):
				if s in samples:
					sample_idx_arr.append(i) 
					sample_list.append(s)
			print header[:4]
			print sample_list
			o.write("\t".join(header[:4]+sample_list)+"\n")
			print "Filtering variants"
			for i in tqdm(range(file_len(infile)-1)):
				line = f.readline()
				arr = line.rstrip().split("\t")
				new_arr = arr[:4]
				for i in sample_idx_arr:
					new_arr.append(arr[i])
				o.write("\t".join(new_arr)+"\n")


infile = sys.argv[1]
samples = set([x.rstrip() for x in open(sys.argv[2]).readlines()])
outfile = sys.argv[3]
matrix_subset(infile,outfile,samples)
