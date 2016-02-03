import time
f = open("paulsYRI.txt")

t = time.time()
K =  list()
for row in f:
	K.append(row)
print time.time() - t
print

f.close()

f = open("paulsYRI.txt")
t = time.time()
K = f.readlines()
print time.time() - t
print

f.close()

f = open("paulsYRI.txt")
t = time.time()
K = [row for row in f if sldfjsdkl else ljsldjflsjk]
print time.time() - t
print

f.close()

