#! /home/jody/software/anaconda2/bin/python
from sys import argv
import gzip
import json
from subprocess import PIPE,Popen
import os
import re 
from collections import Counter

script,sample,refFile,baseDir,cutoff,pct,outFile = argv

scriptDir = os.path.dirname(os.path.realpath(__file__))


pct = float(pct)
FNULL = open(os.devnull, 'w') # for the stderr
cutoff = float(cutoff)
tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"

indels = json.loads(open("indels.json","r").readline().rstrip())
covFile = baseDir+"/coverage/"+sample+".coverage.gz"

ofCalls = open(outFile,"w")
ofCalls.write("%s\n" % sample)

def callDel(res):
	arr = res.rstrip().split("\t")

	if len(arr)<2:
		return "NA"
	pileup = arr[4]
	if len(pileup)<cutoff:
		return "NA"
	insArr = re.findall("\+[0-9]+[ACGTNacgtn]+",pileup)
	insReads = len(insArr)
	refReads = len(re.findall("[.,]",pileup)) - insReads

	
	if (insReads==0) and (refReads>=cutoff):
		return arr[2]
	elif (insReads==0) and (refReads<cutoff):
		return "NA"
	
	cnt = Counter()
	for s in [x.upper() for  x in insArr]:
		cnt[s]+=1
	indel_seq = cnt.most_common(1)[0][0]

	test = float(insReads)/(float(insReads)+float(refReads))

	if test > pct:
		return indel_seq
	elif (1-test) > pct:
		return arr[2]
	else:
		return "NA"

		
def callInd(res):
	arr = res.rstrip().split("\t")
	

	if len(arr)<2:
		return "NA"
	pileup = arr[4]
	if len(pileup)<cutoff:
		return "NA"
	insArr = re.findall("\+[0-9]+[ACGTNacgtn]+",pileup)
	insReads = len(insArr)
	refReads = len(re.findall("[.,]",pileup)) - insReads

	if (insReads==0) and (refReads>=cutoff):
		return arr[2]
	elif (insReads==0) and (refReads<cutoff):
		return "NA"

	cnt = Counter()
	for s in [x.upper() for  x in insArr]:
		cnt[s]+=1
	indel_seq = cnt.most_common(1)[0][0]

	test = float(insReads)/(float(insReads)+float(refReads))
	
	if test > pct:
		return indel_seq
	elif (1-test) > pct:
		return arr[2]
	else:
		return "NA"

dict_pileup = {}
chr_set = set()
bamFile = baseDir+"/bam/"+sample+".bam"
p1 = Popen(["samtools", "view", "-b",bamFile,"-L","indels.bed"],stderr=FNULL,stdout=PIPE)
p2 = Popen(["samtools","mpileup","-f",refFile,"-"],stdin=p1.stdout,stderr=FNULL,stdout=PIPE)
p1.stdout.close()
p3 = Popen(["python",scriptDir+"/filterPileup.py","indels.json"],stdin=p2.stdout,stderr=FNULL,stdout=PIPE)
p2.stdout.close()
for line in p3.stdout:
	chr,pos = line.split()[:2]
	if chr not in chr_set:
		chr_set.add(chr)
		dict_pileup[chr] = {}
	dict_pileup[chr][pos] = line.rstrip()

for chrom in indels:
	for var in indels[chrom]:
		if var["pos"] not in dict_pileup[chrom]:
			ofCalls.write("%s\n" % "NA")
			continue
		res = dict_pileup[chrom][var["pos"]]
		if var["type"]=="deletion":
			ofCalls.write("%s\n" % callDel(res))
		elif var["type"]=="insertion":
			ofCalls.write("%s\n" % callInd(res))


