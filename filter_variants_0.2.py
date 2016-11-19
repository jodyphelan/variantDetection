#! /Users/jody/anaconda2/bin/python
from __future__ import division
import os
import argparse
import sys
import gzip
import json
import re
from tqdm import tqdm
from psGeneral import file_len
import matplotlib.pyplot as plt
import numpy as np
import subprocess


####### Global Variables #########
scriptDir = os.path.dirname(os.path.realpath(__file__))



####### Functions ###########


def loadSampleFile(sampleFile):
        f = open(sampleFile,"r")
        return map(lambda x:x.rstrip(),f.readlines())
	f.close()

def parseVCF(samples,baseDir):
	print "Loading VCFs"
	setChr = set()
	setPos = {}
	variants = {}
	for i in tqdm(range(len(samples))):
		sample = samples[i]
	        f = gzip.open(baseDir+"/vcf/"+sample+".vcf.gz","rb")
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

			if chr not in setChr:
				setChr.add(chr)
				setPos[chr] = set()
				variants[chr] = [{"pos":pos,"ref":ref,"alt":alt,"type":type}]
			elif pos not in setPos[chr]:
				setPos[chr].add(pos)
				variants[chr].append({"pos":pos,"ref":ref,"alt":alt,"type":type})
		f.close()
	return variants


def variants2txt(outname,vars):
	ofTxt = open(outname,"w")
	ofTxt.write("chr\tpos\tref\tid\ttype\n")
	for chr in vars:
		for var in vars[chr]:
			var_id = "."
			ofTxt.write("%s\t%s\t%s\t%s\t%s\n" % (chr,var["pos"],var["ref"][0],var_id,var["type"]))

def variants2bed(fname,variants):
	ofBED = open(fname,"w")
	for chr in variants:
		for var in variants[chr]:
			ofBED.write("%s\t%s\t%s\n" % (chr,var["pos"],var["pos"]))

def getCalls(sampleFile,baseDir,depth_cut,pct_cut):
	subprocess.call("cat %s | xargs -i -P20 sh -c \"%s/validateVariants.py {} %s %s %s {}.variantCalls\"" % (sampleFile,scriptDir,baseDir,depth_cut,pct_cut),shell=True)

####################### Main functions ############################

def main_raw(args):
    samples = loadSampleFile(args.sampleFile)
#    variants = parseVCF(samples,args.baseDir)
#    open("variants.json","w").write(json.dumps(variants))
    variants = json.loads(open("variants.json").readline())
    variants2txt("variants.txt",variants)
    variants2bed("variants.bed",variants)
    getCalls(args.sampleFile,args.baseDir,args.depth,args.percent)

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('raw', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('sampleFile',help='Sample File')
parser_sub.add_argument('baseDir',help='Base directory')
parser_sub.add_argument('depth',help='RefFile')
parser_sub.add_argument('percent',help='RefFile')
parser_sub.set_defaults(func=main_raw)

args = parser.parse_args()
args.func(args)
