from __future__ import division
import sys
import json
from collections import defaultdict,Counter

def fa2dict(filename):
        print "Loading sequences"
        fa_dict = {}
        seq_name = ""
        for l in open(filename):
                line = l.rstrip()
                if line[0] == ">":
                        seq_name = line[1:].split()[0]
                        fa_dict[seq_name] = []
                else:
                        fa_dict[seq_name].append(line)
        result = {}
        for seq in fa_dict:
                result[seq] = "".join(fa_dict[seq])
        return result



sample_file = sys.argv[1]
json_dir = sys.argv[2]
ref_file = sys.argv[3]
outfile = sys.argv[4]


samp_miss_cut = 0.1
var_miss_cut = 0.1

fa_dict = fa2dict(ref_file)
samples = [x.rstrip() for x in open(sample_file).readlines()]
all_vars = defaultdict(lambda:defaultdict(dict))
loci = set()
for s in samples:
	var_file = "%s/%s.var.json" % (json_dir,s)
	var = json.load(open(var_file))
	for chrom in var:
		for pos in var[chrom]:
			loci.add((chrom,int(pos)))
			all_vars[chrom][int(pos)][s] = var[chrom][pos]

const_miss1 = set(["-","N"])
const_miss2 = set(["-"])
const_miss3 = set(["N"])

for chrom in all_vars:
	for pos in sorted(all_vars[chrom]):
		tv = set(all_vars[chrom][pos].values())
		if tv==const_miss2 or tv==const_miss1 or tv==const_miss3: 
			del all_vars[chrom][pos]
			loci.remove((chrom,pos))

samp_missing = defaultdict(int)
samp_miss_pass = []
samp_miss_fail = []  
for s in samples:
	calls = Counter([all_vars[c][p][s] for c,p in loci if s in all_vars[c][p]])
	miss = (calls["N"]+calls["-"])/len(loci)
	samp_missing[s] = miss
	if miss<samp_miss_cut:
		samp_miss_pass.append(s)
	else:
		samp_miss_fail.append(s)

print "%s/%s samples pass missingness cutoff" % (len(samp_miss_pass),len(samples))

for s in samp_miss_fail:
	for chrom,pos in loci:
		if s in all_vars[chrom][pos]:
			del all_vars[chrom][pos][s]

samples = samp_miss_pass

var_missing = defaultdict(dict)
var_miss_pass = []
var_miss_fail = []
for chrom,pos in loci:
	calls = Counter(all_vars[chrom][pos].values())
	miss = (calls["N"]+calls["-"])/len(samples)
	var_missing[(chrom,pos)] = miss
	if miss<var_miss_cut:
		var_miss_pass.append((chrom,pos))
	else:
		var_miss_fail.append((chrom,pos))
		del all_vars[chrom][pos]

print "%s/%s variants pass missingess cutoff" % (len(var_miss_pass),len(loci))

loci = var_miss_pass

json.dump(all_vars,open(outfile,"w"))
