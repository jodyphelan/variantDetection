import sys
from collections import defaultdict
from tqdm import tqdm
import json
import os
import subprocess

if len(sys.argv)!=6:
	print "script.py <samples> <base_dir> <window> <min_cov> <pct>"
	quit()

####### Global Variables #########
sample_file = sys.argv[1]
base_dir = sys.argv[2]
window = int(sys.argv[3])
min_cov = sys.argv[4]
pct = sys.argv[5]
meta_col_num = 5
script_dir = os.path.dirname(os.path.realpath(__file__))

####### Functions ###########
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def loadSampleFile(sampleFile):
	with open(sampleFile,"r") as f:
		return map(lambda x:x.rstrip(),f.readlines())


def getDellyCalls(sampleFile):
	cmd = "cat %s | xargs -i -P20 sh -c \"python %s/delly_cov2calls.py {} %s %s %s %s\"" % (sampleFile,script_dir,base_dir,window,min_cov,pct)
	subprocess.call(cmd,shell=True)



def fa2dict(filename):
        fa_dict = {}
        seq_name = ""
        with open(filename) as f:
                for i in range(file_len(filename)):
                        line = f.readline().rstrip()
                        if line[0] == ">":
                                seq_name = line[1:].split()[0]
                                fa_dict[seq_name] = []
                        else:
                                fa_dict[seq_name].append(line)
        result = {}
        for seq in fa_dict:
                result[seq] = "".join(fa_dict[seq])
        return result


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
		callFiles = ".delly_calls ".join(arr)+".delly_calls" 
		subprocess.call("paste %s > %s.varTempMat" % (callFiles,j),shell=True)
	arr = [str(x) for x in list(range(1,j+1))]
	
	varFiles = ".varTempMat ".join(arr)+".varTempMat"
	
	varFiles = "variants.txt "+varFiles 
	subprocess.call("paste %s > %s" % (varFiles,"raw.mat"),shell=True)
	matRemoveMono("raw.mat","%s.mat" % prefix)
	subprocess.call("rm *TempMat",shell=True)	


ref_dict = fa2dict("/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa")


with open("variants.txt","w") as o:
	o.write("chr\tpos\tref\tinfo\ttype\n")
	for chrom in ref_dict:
		for i in range(window,len(ref_dict[chrom]),window):
					
			o.write("%s\t%s\t%s\t.\t.\n" % (chrom,i-window,i))


samples = loadSampleFile(sample_file)
getDellyCalls(sample_file)
calls2mat(samples,"unfiltered")
