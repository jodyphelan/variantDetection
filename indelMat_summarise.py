import sys
from tqdm import tqdm
from collections import defaultdict
import re
import numpy as np

indel_re = re.compile("[ACGTNacgtn]([\+-])([0-9]+)([ACGTNacgtn]+)")

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1


indels = []
with open(sys.argv[1]) as f:
	f.readline()
	for x in tqdm(range(file_len(sys.argv[1])-1)):
		l = f.readline()
		arr = l.rstrip().split()
		calls = defaultdict(int)
		for i in range(5,len(arr)):
			calls[arr[i]] += 1
		indels.append({"chr":arr[0],"pos":arr[1],"type":arr[4],"calls":calls})

o = open(sys.argv[2]+".indelReport.txt","w")
all_indel_size = []
for var in indels:
	indel_len = []
	indel_freq = []
	for call in var["calls"]:
		if indel_re.match(call)!=None:
			size = indel_re.match(call).group(2)
			indel_len.append(size)
			indel_freq.append(str(var["calls"][call]))
			for i in range(var["calls"][call]):
				all_indel_size.append(int(size))

	o.write("%s\t%s\t%s\t%s\t%s" % (var["chr"],var["pos"],var["type"],";".join(indel_len),";".join(indel_freq))+"\n")

print "Min\tMedian\tMax"
print "%s\t%s\t%s" % (min(all_indel_size),np.median(all_indel_size),max(all_indel_size))
