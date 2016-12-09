import sys

gene_dict = {}
for l in open("/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.ann"):
	arr = l.split()
	gene_dict[arr[1]] = arr[16]

for l in open(sys.argv[1]):
	line = l.rstrip()
	arr = line.split()
	if arr[0]=="chr":
		print line+"\tgene"
		continue
	if arr[1] in gene_dict:
		line = line + "\t"+gene_dict[arr[1]]
	print line
