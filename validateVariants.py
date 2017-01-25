#! /home/jody/software/anaconda2/bin/python
import json
import sys
from subprocess import Popen,PIPE
from collections import defaultdict
import os
if len(sys.argv)<6:
	print "validateVariants.py <sample> <base_dir> <min_cov> <pct_reads> <outfile>"
	quit()

script,sample,base_dir,min_cov,pct,out = sys.argv
scriptDir = os.path.dirname(os.path.realpath(__file__))
tabix = scriptDir+"htslib/tabix"
min_cov = int(min_cov)
pct = float(pct)
covFile = base_dir+"/pileup/"+sample+".pileup.gz"

def index_file(infile):
	tbi = infile+".tbi"
	if not os.path.isfile(tbi):
		import subprocess
		subprocess.call("%s -b 2 -e 2 -s 1 %s" % (tabix,infile),shell=True)

index_file(covFile)




variants = json.loads(open("variants.json").readline())
covCMD = Popen([tabix,covFile,"-T","variants.bed"],stdout=PIPE)


final_calls = defaultdict(dict)
chromosomes = []
positions = defaultdict(list)
for chrom in variants:
	chromosomes.append(chrom)
	for var in variants[chrom]:
		positions[chrom].append(var["pos"])
		final_calls[chrom][var["pos"]] = "NA"


for l in covCMD.stdout:
	chrom,pos,ref,temp_alts,temp_meta =  l.rstrip().split()
	cov = temp_meta.split(":")[1]
	alt_list = temp_alts.split(",")
	cov_list = cov.split(",")
	call_dict = {}
	for i in range(len(alt_list)):
		call_dict[alt_list[i]] = int(cov_list[i])

	total_cov = sum(call_dict.values())
	if total_cov<min_cov:
		final_calls[chrom][pos] = "NA"
		continue
	cutoff = 0
	if total_cov==min_cov:
		cutoff = min_cov
	elif total_cov>min_cov:
		cutoff = total_cov*pct
	allele = ""
	for call in call_dict:
		if call_dict[call]>=cutoff:
			allele = allele+call
	if allele=="":
		allele = "NA"
	final_calls[chrom][pos] = allele

with open(out,"w") as o:
	o.write(sample+"\n")
	for chrom in chromosomes:
		for pos in positions[chrom]:
			o.write(final_calls[chrom][pos]+"\n")
