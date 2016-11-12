import json
import sys


variants = json.loads(open(sys.argv[1]).readline())
set_pos = {}
for chr in variants:
	set_pos[chr] = set()
	for var in variants[chr]:
		set_pos[chr].add(var["pos"])

for line in sys.stdin:
	chr,pos = line.rstrip().split()[:2]
	if pos in set_pos[chr]:
		sys.stdout.write(line)
