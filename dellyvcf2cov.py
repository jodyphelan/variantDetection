import sys
import re

if len(sys.argv) != 3:
	print "dellyvcf2cov.py <in.vcf> <out.cov>"
	quit()

with open(sys.argv[2],"w") as o:
	for l in open(sys.argv[1]):
		if l[0]=="#":
			continue
		arr = l.rstrip().split()
		stats_arr = arr[9].split(":")
		dr,dv,rr,rv = stats_arr[8:12]
		re_obj = re.search("END=(\d+)",arr[7])
		o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (arr[0],arr[1],re_obj.group(1),dr,dv,rr,rv))
		o.write("\n")
		
