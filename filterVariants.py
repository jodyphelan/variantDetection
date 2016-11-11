#! /home/jody/software/anaconda2/bin/python
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


### Global variables
scriptDir = os.path.dirname(os.path.realpath(__file__))
pltvals = []
meta_col_num = 5


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
			
			if chr not in setChr:
				setChr.add(chr)
				setPos[chr] = set()
				variants[chr] = [{"pos":pos,"ref":ref,"alt":alt,"type":type}]
			elif pos not in setPos[chr]:
				setPos[chr].add(pos)
				variants[chr].append({"pos":pos,"ref":ref,"alt":alt,"type":type})
		f.close()
	return variants

def sortVariantHash(obj):
	newObj = {}
	for chr in obj:
		newObj[chr] = []
		sortedPos = sorted(map(lambda x:int(x["pos"]),obj[chr]))
		for pos in sortedPos:
			newObj[chr].append(obj[chr][map(lambda x: x["pos"],obj[chr]).index(str(pos))])
	return newObj

def splitVariants(obj):
	snp_set_chr = set()
	snp_set_pos = {}
	indel_set_chr = set()
	indel_set_pos = {}
	snps = {}
	indels = {}
	for chr in obj:
		for var in obj[chr]:
			if var["type"]=="SNP":
				if chr not in snp_set_chr:
					snp_set_chr.add(chr)
					snps[chr] = [var]
				else:
					snps[chr].append(var)
			elif var["type"]=="insertion" or var["type"]=="deletion":
				if chr not in indel_set_chr:
					indel_set_chr.add(chr)
					indels[chr] = [var]
				else:
					indels[chr].append(var)

	for chr in snps:
		snps[chr].sort(key=lambda x: int(x["pos"]))
	for chr in indels:
		indels[chr].sort(key=lambda x: int(x["pos"]))
	return [snps,indels]

def snps2bed(snps):
	ofBED = open("snps.bed","w")
	for chr in snps:
		for var in snps[chr]:
			ofBED.write("%s\t%s\t%s\n" % (chr,var["pos"],var["pos"]))
def indels2bed(indels):
	ofBED = open("indels.bed","w")
	for chr in indels:
		for var in indels[chr]:
			end = int(var["pos"])+1
			ofBED.write("%s\t%s\t%s\n" % (chr,var["pos"],end))

def variants2txt(outname,vars):
	ofTxt = open(outname,"w")
	ofTxt.write("chr\tpos\tref\tid\ttype\n")
	for chr in vars:
		for var in vars[chr]:
			var_id = "."
			ofTxt.write("%s\t%s\t%s\t%s\t%s\n" % (chr,var["pos"],var["ref"][0],var_id,var["type"]))

def getCalls(sampleFile,baseDir,refFile,depth_cut,pct_cut):
	os.system("cat %s | xargs -i -P20 sh -c \"%s/validateSNPs.py {} %s %s %s {}.snpCalls\"" % (sampleFile,scriptDir,baseDir,depth_cut,pct_cut))
	os.system("cat %s | xargs -i -P20 sh -c \"%s/validateIndels.py {} %s %s %s %s {}.indelCalls\"" % (sampleFile,scriptDir,refFile,baseDir,depth_cut,pct_cut))

def calls2mat(samples,prefix):
	tempSamples = samples
	j = 0
	while len(tempSamples)>0:
		j = j+1
		arr = []
		for i in range(500):
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
#	mat2bin("%s.snps.mat" % prefix,"%s.snps.mat.bin" % prefix)
#	mat2bin("raw.snp.mat","raw.snp.mat.bin")
	matRemoveMono("raw.snp.mat","%s.snps.mat" % prefix)
	matRemoveMono("raw.indel.mat","%s.indels.mat" % prefix)

def testMulti(line):
	tempArr = line.split("\t")
	for i in range(meta_col_num):
		tempArr.pop(0)
	if len(set(tempArr) - set(["-","N","NA"]))>1:
		return True
	else:
		return False

def matRemoveMono(inFile,outFile):
	print "Removing Monomorphic variants"
	o = open(outFile,"w")
	with open(inFile) as f:
		o.write(f.readline())
		for i in tqdm(range(file_len(inFile)-1)): 
		        l = f.readline()
			if testMulti(l.rstrip()):
				o.write(l)
			


def raw(args):
	samples = loadSampleFile(args.sampleFile)
	variants = parseVCF(samples,args.baseDir)
	snps,indels = splitVariants(variants)
	open("indels.json","w").write(json.dumps(indels))
	open("snps.json","w").write(json.dumps(snps))
	
	snps2bed(snps)
	indels2bed(indels)
	variants2txt("snps.txt",snps)
	variants2txt("indels.txt",indels)
	print("Filtering variants")
	getCalls(args.sampleFile,args.baseDir,args.ref,args.depth,args.percent)
	print("Generating variant matrix")
	calls2mat(samples,"unfiltered")

def mat2bin(inFile,outFile):
	print "Creating binary matrix"
	o = open(outFile,"w")
	with open(inFile) as f:
		o.write(f.readline())
		for i in tqdm(range(file_len(inFile)-1)): 
		        l = f.readline().rstrip()
			arr = l.split("\t")
			for i in range(meta_col_num,len(arr)):
				if arr[i]=="-" or arr[i]=="N" or arr[i]=="NA":
					arr[i] = "NA"
				elif len(arr[i])>1:
					arr[i] = "0.5"
				elif arr[i]==arr[2]:
					arr[i] = "0"
				elif arr[i] in (set(["A","C","G","T"])-set([arr[2]])):
					arr[i] = "1"
				else:
					print arr[i]
					raise ValueError("Invalid nucleotide found: " % str(arr[i]))
			o.write("\t".join(arr))
			o.write("\n")

def wrapMat2Bin(args):
	mat2bin(args.infile,args.outfile)

def varStats(infile):
	print "Computing variant stats"
	arr_pct_NA = []
	arr_pct_MX = []

	with open(infile) as f:
		f.readline()
		for i in tqdm(range(file_len(infile)-1)):
			arr = f.readline().rstrip().split("\t")
			pct_NA = len(filter(lambda x: x=="NA", arr[meta_col_num:]))/len(arr[meta_col_num:])
			if arr[4]=="SNP":
				pct_MX = len(filter(lambda x: x!="NA" and  len(x)>1, arr[meta_col_num:]))/len(arr[meta_col_num:])
			else:
				pct_MX = 0 ######## Dont filter mixed values based in indels #############
			arr_pct_NA.append(pct_NA)
			arr_pct_MX.append(pct_MX)
	arr_pct_NA.sort()
	arr_pct_MX.sort()
	plotStats(arr_pct_NA,1)
	na_cut = pltvals[0]
	plotStats(arr_pct_MX,1)
	mx_cut = pltvals[0]
	open("varStats.json","w").write(json.dumps({"na":arr_pct_NA,"mx":arr_pct_MX}))
	return (na_cut,mx_cut)
		
def sampleStats(infile):
	print "Computing Sample stats"
	dict_num_NA = {}
	dict_num_MX ={}
	header = []
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
				elif arr[j]!="NA" and len(arr[j])>1:
					dict_num_MX[header[j]] += 1
	dict_pct_NA = {}
	arr_pct_NA = [] 
	dict_pct_MX = {}
	arr_pct_MX = []
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


def plotStats(dat,num):
	pltvals = []
	def onclick(event):
		global pltvals
		pltvals.append(event.ydata)	
		if len(pltvals)>=num:
			plt.close()
				
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(dat)
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	plt.show()
	
def filterSNPs(infile,outfile,na_cut,mx_cut):
	filtering_results = {"pass":[],"na":[],"mx":[],"map":[]}
	dict_map = loadMappability("mappability.bed")
	o = open(outfile,"w")
	with open(infile) as f:
		o.write(f.readline())
		print "Filtering variants"
		for i in tqdm(range(file_len(infile)-1)):
			line = f.readline()
			arr = line.rstrip().split("\t")
			pct_NA = len(filter(lambda x: x=="NA", arr))/len(arr[meta_col_num:])
			pct_MX = len(filter(lambda x: x=="0.5", arr))/len(arr[meta_col_num:])
			if float(pct_NA)>float(na_cut):
				filtering_results["na"].append((arr[0],arr[1]))
			elif float(pct_MX)>float(mx_cut):
				filtering_results["mx"].append((arr[0],arr[1]))
			elif arr[1] in dict_map[arr[0]]:
				filtering_results["map"].append((arr[0],arr[1]))
			else:
				filtering_results["pass"].append((arr[0],arr[1]))
				o.write(line)
	open("variantFiltering.json","w").write(json.dumps(filtering_results))

def sampleFilter(na_cut,mx_cut):
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
	print sorted(list(samples_pass))
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

def wrapMap(args):
	dict_map = loadMappability(args.infile)
	print dict_map

def wrapStats(args):
	print "Computing sample stats"
#	sample_na_cut,sample_mx_cut = 0.05,0.05
	sample_na_cut,sample_mx_cut = sampleStats("unfiltered.snps.mat")
	print "Filtering with sample Missing-cutoff:%s and Mixed-cutoff:%s " % (sample_na_cut,sample_mx_cut)
	sampleFilter(sample_na_cut,sample_mx_cut)
	mergeMatrices("sample.filt.indels.mat","sample.filt.snps.mat","sample.filt.mat.bin")
#	var_na_cut,var_mx_cut =0.05,0.05
	var_na_cut,var_mx_cut = varStats("sample.filt.mat.bin")

	print "Filtering with variant Missing-cutoff:%s and Mixed-cutoff:%s " % (var_na_cut,var_mx_cut)
	filterSNPs("sample.filt.mat.bin","variant.sample.filt.mat",var_na_cut,var_mx_cut)

def mergeMatrices(file1,file2,file3):
	subprocess.call("cat %s %s | sort -k1,1 -k2,2n  | uniq > %s" % (file1,file2,file3),shell=True)

def wrapMerge(args):
	mergeMatrices(args.f1,args.f2)

def plotData(args):
	results = {}
	sample_stats = json.loads(open("sampleStats.json").readline())
	results["sampleStats"] = {"na":[],"mx":[]}
	for s in sample_stats["na"]:
		results["sampleStats"]["na"].append({"id":s,"val":sample_stats["na"][s]})
		results["sampleStats"]["mx"].append({"id":s,"val":sample_stats["mx"][s]})
	results["varStats"] = {"na":[],"mx":[]}
	results["varStats"]["na"] = json.loads(open("varStats.json").readline())["na"]
	results["varStats"]["mx"] = json.loads(open("varStats.json").readline())["mx"]
	var_filter = json.loads(open("variantFiltering.json").readline())

	results["varFilter"] = {}
	results["varFilter"]["pass"] = len(var_filter["pass"])
	results["varFilter"]["na"] = len(var_filter["na"])
	results["varFilter"]["mx"] = len(var_filter["mx"])
	results["varFilter"]["map"] = len(var_filter["map"])
	
	sample_filter = json.loads(open("sampleFiltering.json").readline())
	results["sampleFilter"] = {}
	results["sampleFilter"]["pass"] = len(sample_filter["pass"])
	results["sampleFilter"]["na"] = len(sample_filter["na"])
	results["sampleFilter"]["mx"] = len(sample_filter["mx"])

	open("plot_data.json","w").write(json.dumps(results))

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_raw = subparsers.add_parser('raw', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sampleFile',help='Sample File')
parser_raw.add_argument('baseDir',help='Base directory')
parser_raw.add_argument('ref',help='RefFile')
parser_raw.add_argument('depth',help='RefFile')
parser_raw.add_argument('percent',help='RefFile')
parser_raw.set_defaults(func=raw)

parser_raw = subparsers.add_parser('binary', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('infile',help='RefFile')
parser_raw.add_argument('outfile',help='RefFile')
parser_raw.set_defaults(func=wrapMat2Bin)

parser_raw = subparsers.add_parser('stats', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.set_defaults(func=wrapStats)

parser_raw = subparsers.add_parser('plots', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.set_defaults(func=plotData)

parser_raw = subparsers.add_parser('merge', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('f1',help='RefFile')
parser_raw.add_argument('f2',help='RefFile')
parser_raw.set_defaults(func=wrapMerge)

parser_raw = subparsers.add_parser('map', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('infile',help='RefFile')
parser_raw.set_defaults(func=wrapMap)


args = parser.parse_args()
args.func(args)
