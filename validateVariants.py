#! /usr/bin/python
from sys import argv
import gzip
import json
from subprocess import PIPE,Popen


script,sample,baseDir,minCov,pct,outFile = argv

minCov = int(minCov)
pct = float(pct)
tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"
snps = json.loads(open("snps.json","r").readline().rstrip())
covFile = baseDir+"/coverage/"+sample+".coverage.gz"
covCMD = Popen([tabix,covFile,"-T","snps.bed"],stdout=PIPE)

ofCalls = open(outFile,"w")
ofCalls.write("%s\n" % sample)

def callAllele(minCov,pct,l):
	chr,pos,tot,a,c,g,t,deletion,refskip,sample = l.split("\t")
	pos,tot,a,c,g,t = map(lambda x:int(x),[pos,tot,a,c,g,t])
	if tot<minCov:
		return "-"
	cutoff = 0
	if tot==minCov:
		cutoff = minCov
	elif tot>minCov:
		cutoff = tot*pct
	nucObj = {"A":a,"C":c,"G":g,"T":t}
	allele = ""
	for nuc in nucObj:
		if nucObj[nuc]>=cutoff:
			print("%s %s %s %s" % (chr,pos,nuc,+nucObj[nuc]))
			allele = allele+nuc
	if allele=="":
		allele = "N"
	return allele

for l in covCMD.communicate()[0].split("\n"):
	if l=="":
		continue
	chr,pos  = l.split("\t")[:2]	
	allele = callAllele(minCov,pct,l)
	ofCalls.write("%s\n" % (allele))
