import sys
import json

if len(sys.argv)!=6:
	print "script.py <sample> <base_dir> <chrom> <start> <end>"
	quit()

script,sample,base_dir,chrom,start,end = sys.argv

start = int(start)
end = int(end)

indel_file = base_dir+"/"+sample+".indels"


indels = json.load(open(indel_file))
for indel in indels:
	if indel["chr"]!=chrom:
		continue
	indel["start"] = int(indel["start"])
	indel["end"] = int(indel["end"])
	if indel["start"]>start and indel["end"]<end:
		size = indel["end"]-indel["start"]
		sys.stdout.write("%s\t%s\t%s\t%s\t%s\n" % (sample,indel["chr"],indel["start"],indel["end"],size))
