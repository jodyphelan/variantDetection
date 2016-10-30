#! /usr/bin/python
from sys import argv
import gzip
import json
from subprocess import PIPE,Popen
import os
import re 

script,sample,refFile,baseDir,cutoff,pct,outFile = argv

pct = float(pct)
FNULL = open(os.devnull, 'w') # for the stderr
cutoff = float(cutoff)
tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"

indels = json.loads(open("indels.json","r").readline().rstrip())
covFile = baseDir+"/coverage/"+sample+".coverage.gz"

ofCalls = open(outFile,"w")
ofCalls.write("%s\n" % sample)

def callDel(res):
	arr = []
	for l in res.split("\n"):
		if l=="":
			break
		tot = l.split("\t")[2]
		indel = l.split("\t")[2]
		arr.append(int(tot))
	meanCov = reduce(lambda x, y: x + y, arr) / float(len(arr))
	if meanCov>cutoff:
		return "0"
	else:
		return "1"
		
def callInd(res):
	arr = res.rstrip().split("\t")
	if len(arr)<2:
		return "N"
	pileup = arr[4]
	print pileup
	if len(pileup)<10:
		print "NAN"
		return "N"
	insReads = len(re.findall("\+[0-9]+[ACGTNacgtn]+",pileup))
	refReads = len(re.findall("[.,]",pileup)) - insReads

	if (insReads==0) and (refReads>10):
		return "1"
	elif (insReads==0) and (refReads<10):
		return "N"

	print insReads
	print refReads
	test = float(insReads)/(float(insReads)+float(refReads))
	if test > pct:
		print "True"
		return "1"
	else:
		print "False"
		return "0"

for chr in indels:
	for var in indels[chr]:
		print "-"*80
		if var["type"]=="deletion":
			indelLen = len(var["ref"])
			posEnd = int(var["pos"])+indelLen
			locus = "%s:%s-%s" % (chr,var["pos"],str(posEnd))
			print locus
			covCMD = Popen([tabix,covFile,locus],stdout=PIPE)
			res = covCMD.communicate()[0]
			ofCalls.write("%s\n" % callDel(res))
		elif var["type"]=="insertion":
			bamFile = baseDir+"/bam/"+sample+".bam"
			locus = chr+":"+var["pos"]+"-"+var["pos"]
			print locus
			p1 = Popen(["samtools", "view", "-b",bamFile,locus],stderr=FNULL,stdout=PIPE)
			p2 = Popen(["samtools","mpileup","-f",refFile,"-"],stdin=p1.stdout,stderr=FNULL,stdout=PIPE)
			p1.stdout.close()
			p3 = Popen(["grep",var["pos"]],stdin=p2.stdout,stderr=FNULL,stdout=PIPE)
			p2.stdout.close()
			res = p3.communicate()[0]
			ofCalls.write("%s\n" % callInd(res))

#for l in covCMD.communicate()[0].split("\n"):
#	if l=="":
#		continue
#	chr,pos  = l.split("\t")[:2]	
#	allele = callAllele(minCov,pct,l)
#	ofCalls.write("%s\n" % (allele))
