import sys
from collections import defaultdict
from subprocess import Popen,PIPE
from tqdm import tqdm

if len(sys.argv)<5:
	print "validateVariants.py <bed> <samples> <base_dir> <outfile>"
	quit()

tabix = "/home/jody/software/samtools-1.3.1/htslib-1.3.1/tabix"

def loadBed(bedFile):
	bed = defaultdict(set)
	for i in [x.rstrip().split() for x in open(bedFile).readlines()]:
	        for j in range(int(i[1])+1,int(i[2])+1):
	                bed[i[0]].add(str(j))
	return bed

def loadSampleFile(sampleFile):
        with open(sampleFile) as f:
                return map(lambda x:x.rstrip(),f.readlines())


bed_file = sys.argv[1]
bed = loadBed(bed_file)




base_dir = sys.argv[3]
samples = loadSampleFile(sys.argv[2])
cov_dict = defaultdict(lambda :defaultdict(dict))
for sample in tqdm(samples):
	covFile = base_dir+"/pileup/"+sample+".pileup.gz"
	covCMD = Popen([tabix,covFile,"-R",bed_file],stdout=PIPE)
	for l in covCMD.stdout:
		chrom,pos,ref,temp_alts,temp_meta =  l.rstrip().split()
		cov = temp_meta.split(":")[1]
		cov_list = cov.split(",")
		total_cov = sum([int(x) for x in cov_list])
		cov_dict[sample][chrom][pos] = total_cov
	
with open(sys.argv[4],"w") as o:
	o.write("chr\tpos"+"\t"+"\t".join(samples)+"\n")
	for chrom in bed:
		for pos in sorted(list(bed[chrom])):
			arr = []
			for s in samples:
				if pos in cov_dict[s][chrom]:
					arr.append(cov_dict[s][chrom][pos])
				else:
					arr.append(0)
			o.write(chrom+"\t"+pos+"\t"+"\t".join([str(x) for x in arr])+"\n")

