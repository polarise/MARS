import sys

def current_line(f,g):
	while f:
		yield f.readline(),g.readline()
	
f = open(sys.argv[1])
g = open(sys.argv[2])
for f_line,g_line in current_line(f,g):
	print f_line.strip()
	print g_line.strip()
	print
