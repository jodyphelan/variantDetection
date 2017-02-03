#! /home/jody/software/anaconda2/bin/python
from __future__ import division
import os
import argparse
import sys
import gzip
import json
import re
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from collections import defaultdict

####### Global Variables #########
scriptDir = os.path.dirname(os.path.realpath(__file__))
pltvals = []
meta_col_num = 5



####### Functions ###########
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def loadSampleFile(sampleFile):
	with open(sampleFile,"r") as f:
		return map(lambda x:x.rstrip(),f.readlines())
	
def check_for_files(samples,base_dir):
	for s in samples:
		vcf_file = base_dir+"/vcf/"+s+".vcf.gz"
		pileup_file = base_dir+"/pileup/"+s+".pileup.gz"
		if not os.path.isfile(vcf_file):
			print "Could not find %s" % vcf_file
			quit()
		elif not os.path.isfile(pileup_file):
			print "Could not find %s" % pileup_file
			quit()

def parseVCF(samples,baseDir):
	print "Loading VCFs"
	temp_variants = defaultdict(dict)
	for i in tqdm(range(len(samples))):
		sample = samples[i]
		f = gzip.open(baseDir+"/vcf/"+sample+".vcf.gz","rb")
		for line in f.readlines():
			if line[0] == "#":
				continue
			chr,pos,id,ref,alt,qual,filter,info,format,sample = line.rstrip().split("\t")
			temp_variants[chr][int(pos)] = {"pos":pos,"ref":ref,"alt":alt,"type":"."}
		f.close()
	variants = defaultdict(list)
	for chrom in temp_variants:
		for pos in sorted(temp_variants[chrom].keys()):
			variants[chrom].append(temp_variants[chrom][pos]) 
	return variants


def variants2txt(outname,vars):
	ofTxt = open(outname,"w")
	ofTxt.write("#chr\tpos\tref\tid\ttype\n")
	for chr in vars:
		for var in vars[chr]:
			var_id = "."
			ofTxt.write("%s\t%s\t%s\t%s\t%s\n" % (chr,var["pos"],var["ref"][0],var_id,var["type"]))

def variants2bed(fname,variants):
	ofBED = open(fname,"w")
	for chr in variants:
		for var in variants[chr]:
			ofBED.write("%s\t%s\t%s\n" % (chr,var["pos"],var["pos"]))

def matRemoveMono(inFile,outFile):
	def testMulti(line):
		tempArr = line.split("\t")
		for i in range(meta_col_num):
			tempArr.pop(0)
		if len(set(tempArr) - set(["-","N","NA"]))>1:
			return True
		else:
			return False

	print "Removing Monomorphic variants"
	o = open(outFile,"w")
	with open(inFile) as f:
		o.write(f.readline())
		for i in tqdm(range(file_len(inFile)-1)): 
			l = f.readline()
			if testMulti(l.rstrip()):
				o.write(l)


def getCalls(sampleFile,baseDir,depth_cut,pct_cut):
	print "Calling variants"
	subprocess.call("cat %s | xargs -i -P20 sh -c \"%s/validateVariants.py {} %s %s %s {}.variantCalls\" " % (sampleFile,scriptDir,baseDir,depth_cut,pct_cut),shell=True)


def calls2mat(samples,prefix):
	print "Running calls2mat"
	tempSamples = samples
	j = 0
	while len(tempSamples)>0:
		j = j+1
		arr = []
		for i in range(500):
			if len(tempSamples)>0:
				arr.append(tempSamples.pop(0))
		callFiles = ".variantCalls ".join(arr)+".variantCalls" 
		os.system("paste %s > %s.varTempMat" % (callFiles,j))
	arr = [str(x) for x in list(range(1,j+1))]
	
	varFiles = ".varTempMat ".join(arr)+".varTempMat"
	
	varFiles = "variants.txt "+varFiles 
	os.system("paste %s > %s" % (varFiles,"raw.mat"))
	matRemoveMono("raw.mat","%s.mat" % prefix)
	subprocess.call("rm *TempMat",shell=True)	

def plotStats(dat,num):
	pltvals = []
	def onclick(event):
		global pltvals
		pltvals = []
		pltvals.append(event.ydata)
		if len(pltvals)>=num:
			plt.close()
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(dat)
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	plt.show()

def sampleStats(infile):
	if os.path.isfile("sampleStats.json"):
		temp_stats = json.load(open("sampleStats.json"))
		dict_pct_NA = temp_stats["na"]
		dict_pct_MX = temp_stats["mx"]
		arr_pct_NA = dict_pct_NA.values()
		arr_pct_MX = dict_pct_MX.values()
	else:
		print "Computing Sample stats"
		dict_num_NA = {}
		dict_num_MX ={}
		header = []
		indel_re = re.compile("[\+-][0-9]+[ACGTNacgtn]+")
		with open(infile) as f:
			header = f.readline().rstrip().split("\t")
			for i in range(meta_col_num,len(header)):
				dict_num_NA[header[i]] = 0
				dict_num_MX[header[i]] = 0
			for i in tqdm(range(file_len(infile)-1)):
				arr = f.readline().rstrip().split("\t")
				for j in range(meta_col_num,len(header)):
					if arr[j] == "NA":
						dict_num_NA[header[j]] += 1
					elif len(indel_re.findall(arr[j]))>0: ## For indels
						pass
					elif arr[j] != "NA" and len(arr[j])>1:
						dict_num_MX[header[j]] += 1
		dict_pct_NA = {}
		arr_pct_NA = []
		dict_pct_MX = {}
		arr_pct_MX = []
		i += 1 ## because the i enumerate starts at 0 we need to add one
		for sample in dict_num_NA:
			dict_pct_NA[sample] = dict_num_NA[sample]/i
			arr_pct_NA.append(dict_num_NA[sample]/i)
			dict_pct_MX[sample] = dict_num_MX[sample]/i
			arr_pct_MX.append(dict_num_MX[sample]/i)
	
		open("sampleStats.json","w").write(json.dumps({"na":dict_pct_NA,"mx":dict_pct_MX}))

	arr_pct_NA.sort()
	arr_pct_MX.sort()
	plotStats(arr_pct_NA,1)
	miss_cut = pltvals[0]
	plotStats(arr_pct_MX,1)
	mix_cut = pltvals[0]
	return (miss_cut,mix_cut)

def sampleFilter(na_cut,mx_cut):
	print "Filtering Samples"
	header = []
	samples_fail = set()
	sample_stats = json.loads(open("sampleStats.json").readline())
	samples = set()
	filtering_results = {"pass":[],"na":[],"mx":[]}
	for s in sample_stats["na"]:
		samples.add(s)
		if float(sample_stats["na"][s])>float(na_cut):
			samples_fail.add(s)
			filtering_results["na"].append(s)
		elif float(sample_stats["mx"][s])>float(mx_cut):
			samples_fail.add(s)
			filtering_results["mx"].append(s)
		samples_pass = samples-samples_fail
	for s in samples_pass:
		filtering_results["pass"].append(s)
	open("sampleFiltering.json","w").write(json.dumps(filtering_results))
	calls2mat(sorted(list(samples_pass)),"sample.filt")

def loadMappability(infile):
	dict_mappability = {}
	with open(infile) as f:
		print "Loading Mappability"
		for i in tqdm(range(file_len(infile)-1)):
			chrom,start,end,id,map_val = f.readline().rstrip().split()
			if float(map_val)==1:
				continue
			if chrom not in dict_mappability:
				dict_mappability[chrom] = set()
			for i in range(int(start),int(end)):
				dict_mappability[chrom].add(i)
	return dict_mappability


def varStats(infile):
	print "Computing variant stats"
	dict_pct_NA = defaultdict(dict)
	dict_pct_MX = defaultdict(dict)
	arr_pct_NA = []
	arr_pct_MX = []
	indel_re = re.compile("[\+-][0-9]+[ACGTNacgtn]+")
	with open(infile) as f:
		f.readline()
		for i in tqdm(range(file_len(infile)-1)):
			line = f.readline().rstrip()
			arr = line.split("\t")
			pct_NA = len(filter(lambda x: x=="NA", arr[meta_col_num:]))/len(arr[meta_col_num:])
			num_indels = len(indel_re.findall(line))
			pct_MX = (len(filter(lambda x: x!="NA" and  len(x)>1, arr[meta_col_num:]))-num_indels)/len(arr[meta_col_num:])
			if pct_MX>0:
				print arr
			dict_pct_NA[arr[0]][arr[1]] = pct_NA
			dict_pct_MX[arr[0]][arr[1]] = pct_MX
			arr_pct_NA.append(pct_NA)
			arr_pct_MX.append(pct_MX)
	arr_pct_NA.sort()
	arr_pct_MX.sort()
	plotStats(arr_pct_NA,1)
	na_cut = pltvals[0]
	plotStats(arr_pct_MX,1)
	mx_cut = pltvals[0]
	open("varStats.json","w").write(json.dumps({"na":dict_pct_NA,"mx":dict_pct_MX}))
	print "NA cutoff: "+str(na_cut)
	return (na_cut,mx_cut)

def filterVars(infile,outfile,na_cut,mx_cut):
	print "Filtering variants"
	filtering_results = {"pass":[],"na":[],"mx":[],"map":[]}
	dict_map = loadMappability("mappability.bed")
	dict_var_stats = json.loads(open("varStats.json").readline())
	o = open(outfile,"w")
	with open(infile) as f:
		o.write(f.readline())
		print "Filtering variants"
		for i in tqdm(range(file_len(infile)-1)):
			line = f.readline()
			arr = line.rstrip().split("\t")
			pct_NA = dict_var_stats["na"][arr[0]][arr[1]]
			pct_MX = dict_var_stats["mx"][arr[0]][arr[1]]
			if int(arr[1]) in dict_map[arr[0]]:
				filtering_results["map"].append((arr[0],arr[1]))
			elif float(pct_NA)>float(na_cut):
				filtering_results["na"].append((arr[0],arr[1]))
			elif float(pct_MX)>float(mx_cut):
				filtering_results["mx"].append((arr[0],arr[1]))
			else:
				filtering_results["pass"].append((arr[0],arr[1]))
				o.write(line)
	open("variantFiltering.json","w").write(json.dumps(filtering_results))

####################### Main functions ############################

def main_raw(args):
	samples = loadSampleFile(args.sampleFile)
	check_for_files(samples,args.baseDir)
	variants = parseVCF(samples,args.baseDir)
	open("variants.json","w").write(json.dumps(variants))
	variants2txt("variants.txt",variants)
	variants2bed("variants.bed",variants)
	getCalls(args.sampleFile,args.baseDir,args.depth,args.read_proportion)
	calls2mat(samples,"unfiltered")

def main_filter(args):
	sample_na_cut,sample_mx_cut = sampleStats("unfiltered.mat")
	sampleFilter(sample_na_cut,sample_mx_cut)
	var_na_cut,var_mx_cut = varStats("sample.filt.mat")
	filterVars("sample.filt.mat","variant.sample.filt.mat",var_na_cut,var_mx_cut)	
####################### Argument Parser ###########################

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('raw', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('sampleFile',help='Sample File')
parser_sub.add_argument('baseDir',help='Base directory')
parser_sub.add_argument('depth',help='minimum depth required to call alleles')
parser_sub.add_argument('read_proportion',help='minimum proportion of all reads to call alleles')
parser_sub.set_defaults(func=main_raw)

parser_sub = subparsers.add_parser('filter', help='Generate filtered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.set_defaults(func=main_filter)


args = parser.parse_args()
args.func(args)
