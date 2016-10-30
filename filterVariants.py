#! /usr/bin/python
import os
import argparse
import sys
import gzip
import json
import re 

scriptDir = os.path.dirname(os.path.realpath(__file__))



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
		type = ""
		if re.search("INDEL",line)==None:
			type = "SNP"
		elif len(ref)<len(alt):
			type = "insertion"
		elif len(ref)>len(alt):
			type = "deletion"
		else:	
			type = "unknown"
		
		if chr not in variants:
			variants[chr] = [{"pos":pos,"ref":ref,"alt":alt,"type":type}]
		elif pos not in map(lambda x: x["pos"], variants[chr]):
			variants[chr].append({"pos":pos,"ref":ref,"alt":alt,"type":type})
	f.close()

def sortVariantHash(obj):
	newObj = {}
	for chr in obj:
		newObj[chr] = []
		sortedPos = sorted(map(lambda x:int(x["pos"]),obj[chr]))
		for pos in sortedPos:
			newObj[chr].append(obj[chr][map(lambda x: x["pos"],obj[chr]).index(str(pos))])
	return newObj

def splitVariants(obj):
	snps = {}
	indels = {}
	for chr in obj:
		for var in obj[chr]:
			if var["type"]=="SNP":
				if chr not in snps:
					snps[chr] = [var]
				else:
					snps[chr].append(var)
			elif var["type"]=="insertion" or var["type"]=="deletion":
				if chr not in indels:
					indels[chr] = [var]
				else:
					indels[chr].append(var)
	snps = sortVariantHash(snps)
	indels = sortVariantHash(indels)
	return [snps,indels]

def snps2bed(snps):
	ofBED = open("snps.bed","w")
	for chr in snps:
		for var in snps[chr]:
			ofBED.write("%s\t%s\t%s\n" % (chr,var["pos"],var["pos"]))

def variants2txt(outname,vars):
	ofTxt = open(outname,"w")
	ofTxt.write("chr\tpos\tref\n")
	for chr in vars:
		for var in vars[chr]:
			ofTxt.write("%s\t%s\t%s\t%s\n" % (chr,var["pos"],var["ref"],var["type"]))

def getCalls(sampleFile,baseDir,refFile):
	os.system("cat %s | xargs -i -P20 sh -c \"%s/validateSNPs.py {} %s 10 0.7 {}.snpCalls\"" % (sampleFile,scriptDir,baseDir))
	os.system("cat %s | xargs -i -P20 sh -c \"%s/validateIndels.py {} %s %s 10 0.7 {}.indelCalls\"" % (sampleFile,scriptDir,refFile,baseDir))

def calls2mat(samples):
	tempSamples = samples
	j = 0
	while len(tempSamples)>0:
		j = j+1
		arr = []
		for i in (range(500)):
			if len(tempSamples)>0:
				arr.append(tempSamples.pop(0))
		snpFiles = ".snpCalls ".join(arr)+".snpCalls" 
		indelFiles = ".indelCalls ".join(arr)+".indelCalls"
		os.system("paste %s > %s.snpTempMat" % (snpFiles,j))
		os.system("paste %s > %s.indelTempMat" % (indelFiles,j))
	arr = [str(x) for x in list(range(1,j+1))]
	
	snpFiles = ".snpTempMat ".join(arr)+".snpTempMat"
	indelFiles = ".indelTempMat ".join(arr)+".indelTempMat"
		
	snpFiles = "snps.txt "+snpFiles 
	indelFiles = "indels.txt "+indelFiles
	os.system("paste %s > %s" % (snpFiles,"raw.snp.mat"))
	os.system("paste %s > %s" % (indelFiles,"raw.indel.mat"))
	matRemoveMono("raw.snp.mat","snps.unfiltered.mat")
	matRemoveMono("raw.indel.mat","indels.unfiltered.mat")

def testMulti(line):
	print line
	tempArr = line.split("\t")
	print tempArr
	for i in range(4):
		tempArr.pop(0)
	if len(set(tempArr) - set(["-","N"]))>1:
		return True
	else:
		return False

def matRemoveMono(inFile,outFile):
	o = open(outFile,"w")
	with open(inFile) as f:
		o.write(f.readline())
		for i in range(file_len(inFile)-1): 
		        l = f.readline()
			print l
			if testMulti(l.rstrip()):
				o.write(l)
			
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def main(args):
	samples = loadSampleFile(args.sampleFile)
	print samples
	variants = {}
	for s in samples:
		parseVCF(variants,args.baseDir,s)
	snps,indels = splitVariants(variants)
	ofINDELS = open("indels.json","w")
	ofINDELS.write(json.dumps(indels))
	ofINDELS.close()
	ofSNPs = open("snps.json","w")
	ofSNPs.write(json.dumps(snps))
	ofSNPs.close()
	snps2bed(snps)
	variants2txt("snps.txt",snps)
	variants2txt("indels.txt",indels)
	getCalls(args.sampleFile,args.baseDir,args.ref)
	calls2mat(samples)
	

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_raw = subparsers.add_parser('raw', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sampleFile',help='Sample File')
parser_raw.add_argument('baseDir',help='Base directory')
parser_raw.add_argument('ref',help='RefFile')
parser_raw.set_defaults(func=main)


args = parser.parse_args()
args.func(args)
