import sys
from tqdm import tqdm
from collections import defaultdict

genotypes_start_col = 3

if len(sys.argv)<3:
	print "indels2genesum <indels.mat.bin> <out.genesum>"
	quit()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

print "Loading Annotation"
gene_dict = {}
gene_set = set()
with open("/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.ann") as f:
	for i in tqdm(range(file_len("/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.ann"))):
		l = f.readline()
		arr = l.split()
		gene_dict[arr[1]] = arr[16]


mut_dict = defaultdict(lambda: defaultdict(int))
print "Loading indels"
with open(sys.argv[1]) as f:
	header = f.readline().rstrip().split()
	for x in tqdm(range(file_len(sys.argv[1])-1)):
		arr = f.readline().rstrip().split()
		for i in range(genotypes_start_col,len(arr)):
			if arr[i]=="1":
				mut_dict[gene_dict[arr[1]]][header[i]] += 1

print "Computing genesum" 
with open(sys.argv[2],"w") as o:
	o.write(" \"")
	o.write("\" \"".join(header[genotypes_start_col:]))
	o.write("\"\n")
	for gene in tqdm(sorted(mut_dict.keys())):
		o.write("\""+gene+"\"")
		for s in header[genotypes_start_col:]:
			o.write(" "+str(mut_dict[gene][s]))
		o.write("\n")
