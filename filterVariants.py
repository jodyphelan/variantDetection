#! /usr/bin/python
import os
import argparse
import sys
import gzip
import json

scriptDir = os.path.dirname(os.path.realpath(__file__))

script,sampleFile,outFile = sys.argv


def loadSampleFile(sampleFile):
        f = open(sampleFile,"r")
        return map(lambda x:x.rstrip(),f.readlines())
	f.close()

def parseVCF(variants,baseDir,sample):
        f = gzip.open(baseDir+"/vcf/"+sample+".filt.vcf.gz","rb")
	for line in f.readlines():
		if line[0] == "#":
			continue
		chr,pos,id,ref,alt,qual,filter,info,format,sample = line.rstrip().split("\t")
		if chr not in variants:
			variants[chr] = {pos:ref}
		elif pos not in variants[chr]:
			variants[chr][pos] = ref
	f.close()

def splitVariants(obj):
	snps = {}
	indels = {}
	for chr in obj:
		for pos in obj[chr]:
			if len(obj[chr][pos])==1:
				if chr not in snps:
					snps[chr] = {pos:obj[chr][pos]}
				elif pos not in snps[chr]:
					snps[chr][pos] = obj[chr][pos]
			elif len(obj[chr][pos])>1:
				if chr not in indels:
					indels[chr] = {pos:obj[chr][pos]}
				elif pos not in indels[chr]:
					indels[chr][pos] = obj[chr][pos]
	return [snps,indels]

def snps2bed(snps):
	ofBED = open("snps.bed","w")
	for chr in snps:
		snpPos = [int(x) for x in snps[chr]]
		snpPos.sort()
		for pos in snpPos:
			ofBED.write("%s\t%s\t%s\n" % (chr,pos,pos))

def snps2txt(snps):
	ofTxt = open("snps.txt","w")
	ofTxt.write("chr\tpos\tref\n")
	for chr in snps:
		snpPos = [int(x) for x in snps[chr]]
		snpPos.sort()
		for pos in snpPos:
			pos = str(pos)
			ofTxt.write("%s\t%s\t%s\n" % (chr,pos,snps[chr][pos]))

def getCalls():
	os.system("cat samples | xargs -i -P20 sh -c \"./validateVariants.py {} . 10 0.7 {}.calls\"")

def calls2mat(samples):
	tempSamples = samples
	j = 0
	while len(tempSamples)>0:
		j = j+1
		arr = []
		for i in (range(2)):
			if len(tempSamples)>0:
				arr.append(tempSamples.pop(0))
		files = ".calls ".join(arr)+".calls" 
		os.system("paste %s > %s.tempMat" % (files,j))
	arr = [str(x) for x in list(range(1,j+1))]
	
	files = ".tempMat ".join(arr)+".tempMat"	
	files = "snps.txt "+files 
	os.system("paste %s > %s" % (files,"unfiltered.mat"))
	

def main():
	samples = loadSampleFile(sampleFile)
	print samples
	variants = {}
	for s in samples:
		parseVCF(variants,".",s)
	snps,indels = splitVariants(variants)
	ofINDELS = open("indels.json","w")
	ofINDELS.write(json.dumps(indels))
	ofINDELS.close()
	ofSNPs = open("snps.json","w")
	ofSNPs.write(json.dumps(snps))
	ofSNPs.close()
	snps2bed(snps)
	snps2txt(snps)
	getCalls()
	calls2mat(samples)

main()
