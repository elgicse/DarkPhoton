import math

num = []
ct = []
gamma = []
mass = []
eps = []

res = []


with open("fastDataTable_20140502.dat","r") as f:
    for line in f:
        line=line.split("\t")
        ct.append(pow(10.,-6.)*float(line[0]))
        gamma.append(float(line[1]))
        num.append(float(line[2]))
        eps.append(float(line[3]))
        mass.append(float(line[4]))

for i in xrange(len(num)):
    res.append(  math.exp((-60.)/(gamma[i]*ct[i])) - math.exp((-100.)/(gamma[i]*ct[i])) )

print "mass\t eps\t ctau\t esp\t expected"
for i in xrange(len(res)):
     print "%s\t%s\t%s\t%s\t%s"%(mass[i], eps[i], ct[i], res[i], num[i]*res[i])

