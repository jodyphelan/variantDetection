#! /usr/bin/python
from sys import argv
import gzip
import json
import os
from subprocess import PIPE,Popen
import subprocess

script,sample,baseDir,minCov,pct,outFile = argv

minCov = int(minCov)
pct = float(pct)
tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"
covFile = baseDir+"/coverage/"+sample+".coverage.gz"
if not os.path.isfile(covFile+".tbi"):
	tabix_cmd = tabix + " -0 -b 2 -e 2 -s 1 -S 1 "+covFile
	subprocess.call(tabix_cmd,shell=True)
covCMD = Popen([tabix,covFile,"-T","snps.bed"],stdout=PIPE)

ofCalls = open(outFile,"w")
ofCalls.write("%s\n" % sample)

def callAllele(minCov,pct,l):
	chrom,pos,tot,a,c,g,t,deletion,refskip,sample = l.split("\t")
	pos,tot,a,c,g,t = map(lambda x:int(x),[pos,tot,a,c,g,t])
	if tot<minCov:
		return "NA"
	cutoff = 0
	if tot==minCov:
		cutoff = minCov
	elif tot>minCov:
		cutoff = tot*pct
	nucObj = {"A":a,"C":c,"G":g,"T":t}
	allele = ""
	for nuc in nucObj:
		if nucObj[nuc]>=cutoff:
			allele = allele+nuc
	if allele=="":
		allele = "NA"
	return allele

for l in covCMD.communicate()[0].split("\n"):
	if l=="":
		continue
	chrom,pos  = l.split("\t")[:2]	
	allele = callAllele(minCov,pct,l)
	ofCalls.write("%s\n" % (allele))
