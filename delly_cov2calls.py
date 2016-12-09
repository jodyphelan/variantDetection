from __future__ import division
import sys
from collections import defaultdict
from math import ceil
import json

if len(sys.argv)!=6:
	print "delly_cov2calls.py <sample> <base_dir> <window> <min_cov> <pct>"
	quit() 

sample = sys.argv[1]
base_dir = sys.argv[2]
cov_file = base_dir + "/delly_cov/" + sample + ".delly_cov"
window = int(sys.argv[3])
min_cov = int(sys.argv[4])
pct = float(sys.argv[5])
cov = defaultdict(lambda : defaultdict(list))
out = sample+".delly_calls"



def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


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

ref_dict = fa2dict("/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa")




for l in open(cov_file):
	arr = l.rstrip().split()
	w = int(ceil(int(arr[1])/window)*window)
	dr,dv,rr,rv = [int(x) for x in arr[3:7]]
	paired_cov = dr+dv
	split_cov = rr+rv
	if paired_cov<min_cov and split_cov<min_cov:
		continue
	paired_test=0
	split_test=0
	if paired_cov>min_cov:
		paired_test = dv/(dr+dv)
	if split_cov>min_cov:
		split_test = rv/(rr+rv)
	if paired_test>pct or split_test>pct:
		cov[arr[0]][w].append({"chr":arr[0],"start":arr[1],"end":arr[2],"dr":dr,"dv":dv,"rr":rr,"rv":rv})

calls = []
for chrom in ref_dict:
	for w in cov[chrom]:
		for var in cov[chrom][w]:
			calls.append(var)
json.dump(calls,open(sample+".indels","w"))


with open(out,"w") as o:
	o.write(sample+"\n")
	for chrom in ref_dict:
		for i in range(window,len(ref_dict[chrom]),window):
			if i in cov[chrom]:
				call="1"
			else:
				call="0"
			o.write("%s\n" % (call))
