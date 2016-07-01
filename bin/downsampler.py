#!/usr/bin/python3

# https://github.com/iskold
# Created Tue Jun 28 21:21:01 CEST 2016

import argparse,sys,os,re,math,time #standard libs I frequently use
from datetime import datetime as dt

parser = argparse.ArgumentParser(prog="downsample.py",usage="python3 downsample.py -f -t <target-read-count> -i <fastq-file> -o <fastq-file.output>",epilog="Written by https://github.com/iskold, on Tue Jun 28 21:21:01 CEST 2016",description="Description: Randomly downsamples a fastq-file to a specific count. NOTE! Does not normalize. If you want to normalize, use BBmap, and the tool called BBNorm.sh. You should downsample if you want to simulate lower sequencing depth. For example if you want to do gene counts matrices based on reads across samples with different depths.")
parser.add_argument("-i",metavar="fastq", help="Input file",required="True",nargs="*")
parser.add_argument("-o",metavar="fasta", help="Output file. Default: STDOUT. If paired end data, STDOUT will be disabled and output will automatically have 'DS' added to the input file name.")
parser.add_argument("-t",metavar="number",help="Target read count",required="True",type=int,nargs=1)
parser.add_argument("-a", help="Low-memory mode (uses basically no memory, but results in the target read counts +/- a few reads. Default: Off.",action="store_true")
parser.add_argument("-v",help="Verbose. Prints helpful progress message to STDERR. Default is off.",action="store_true")
args = parser.parse_args()

if args.v:
	def eprint(*args,**kwargs):
		for arg in args:
			print(arg, file=sys.stderr, **kwargs)
			sys.stdout.flush()
else:
	def eprint(*args,**kwarfs):
		pass
	#eprint = lambda *a: None


if args.o and len(args.i) == 1:
	def oprint(*args,**kwargs):
		if args[0] == "tmp2":
			pass
		else:
			writefastqs(args[0],args[1])
#			writefastqs(args,kwargs)
elif not args.o and len(args.i) == 1:
	def oprint(*args,**kwargs):
		if args[0] == "tmp":
			pass
		else:
			print("{0}".format(args[1].rstrip()))
else:
	def oprint(*args,**kwargs):
		writefastqs(args[0],args[1],args[2],args[3])

if len(args.i) == 1: #syntax is writefastqs(fid1,entry1,fid2,entry2)
	fastq1 = args.i[0]
	if not args.o:
		fastq1o = "tmp2"
	else:
		fastq1o = args.o
	fastq2 = "tmp1"
	fastq2o = "tmp2"
	def writefastqs(*args,**kwargs):
		sup = args[0].write(args[1])

elif len(args.i) == 2:
	fastq1 = args.i[0]
	fastq1o = ".".join(args.i[0].strip(".gz").strip(".gzip").split(".")[0:-1])+".DS.fq"
	fastq2 = args.i[1]
	fastq2o = ".".join(args.i[1].strip(".gz").strip(".gzip").split(".")[0:-1])+".DS.fq"
	def writefastqs(*args,**kwargs):
		sup = args[0].write(args[1])
		sup = args[2].write(args[3])

else:
	sys.exit("Number of input arguments ({0}) not understood".format(len(args.i)))


eprint("#Provided arguments: {0}".format(args))		
start = time.time()
d_ = dt.today()
timestarted = d_.strftime("%Y-%m-%d %H:%M")
eprint("#Started at: {0}".format(timestarted))
eprint("\n#Loading modules...",end="")
import numpy as np
import matplotlib.pyplot as plt
import random
import gzip
eprint("Done!")
sys.stdout.flush()
n = args.t[0]

eprint("#Counting fastq entries...",end="")
def fastqItr(myfile):
	if myfile == "tmp1":
		while 1:
			yield ""
	ext = myfile.split(".")[-1]
	if ext == "gz" or ext == "gzip":
		with gzip.open(myfile,"rt") as fid:
			entry = ""
			line = fid.readline()
			while line:
				entry += line
				for i in range(3):
					entry += fid.readline()
					yield entry
				entry = ""
				line = fid.readline()	
	else:
		with open(myfile) as fid:
			entry = ""
			line = fid.readline()
			while line:
				entry += line
				for i in range(3):
					entry += fid.readline()
					yield entry
				entry = ""
				line = fid.readline()



def fastqItr_generator(fid): #Had to split this up because of gzip, and "with"-statement messing up
	entry = ""
	line = fid.readline()
	while line:
		entry += line
		for i in range(3):
			entry += fid.readline()
		yield entry
		entry = ""
		line = fid.readline()

c = 0
for entry in fastqItr(fastq1):
	c+=1
eprint("Done! (Entries: {0})".format(c))

if n > c:
	sys.exit("#ERROR!!!! Target number ({0}) is higher than the number of fastq-entries ({1})!".format(n,c))

#ratio = int(round(100*n/(c+0.0)))
ratio = n/(c+0.0)

#clist = [1 for x in range(ratio)]+[0 for x in range(100-ratio)]
#random.shuffle(clist)

eprint("#Building random list...",end="")
clist = []
if not args.a: # make a list of numbers
	clist  = set()
	while len(clist) != n:
		clist.add(random.randint(0,c-1))
	clist = sorted(clist)
	eprint("Done!")
else:
	eprint("Skipped!")

eprint("#Choosing fastq entries...",end="")
#with open("random.list") as fid:

def up2date(f,x,cc,c):
#	print(x,cc,c,file=sys.stderr)
	while cc <= x and cc != c+1:
		cc += 1
		entry = next(f)
#	print(x,cc,c,file=sys.stderr)
	return entry,cc



if not clist:
	with open(fastq1o,"w") as fo1, open(fastq2o,"w") as fo2:
		i = 0
		c2 = 0
		k = 0
		new = list()
		f1 = fastqItr(fastq1)
		f2 = fastqItr(fastq2)
		for line in f1:
			if random.random() <= ratio:
				i += 1
				entry,c2 = up2date(f2,k,c2,c)
				oprint(fo1,line,fo2,entry)
			k += 1
else:
	with open(fastq1o,"w") as fo1, open(fastq2o,"w") as fo2:
#		eprint(clist,len(clist))
		i = len(clist)
		c1 = 0
		c2 = 0
		k = 0
		new = list()
		f1 = fastqItr(fastq1)
		f2 = fastqItr(fastq2)
		for x in clist:
			entry,c1 = up2date(f1,x,c1,c)
			entry2,c2 = up2date(f2,x,c2,c)
			oprint(fo1,entry,fo2,entry2)

eprint("Done!\n")
if os.path.isfile("tmp1"):
	os.remove("tmp1")
if os.path.isfile("tmp2"):
	os.remove("tmp2")
timeused = (time.time() - start) / 60.0
#eprint("Time used: "+str(round(timeused)) + " min ("+str(round(timeused/60,1))+" hours)\n",args.quiet)
timefinish = d_.strftime("%Y-%m-%d %H:%M")
eprint("#Finished at: {0}".format(timefinish))
eprint("#Time used: {0} min ({1} hours)".format(round(timeused,2),round(timeused/60.0,1)))


sys.exit() #Legacy code below this line. Mostly used to visualize the randomness of the algorithm
if 1:
 i = 0
 new = list()
 for line in fid:
  #if random.choice(clist) == 1:
  if random.random() <= ratio:
   i += 1
   new.append(line.rstrip())




y = [int(i) for i in new]
x = range(n)
x = [i+1 for i in x]

yy = list()
for i in range(c):
 if i in y:
  yy.append(1)
 else:
  yy.append(0)

zz1 = list()
for i in range(c):
 if i in zz:
  zz1.append(1)
 else:
  zz1.append(0)



plt.plot(x[0:1000], yy[0:1000])
plt.plot(x[0:1000], zz1[0:1000])
plt.show()


timeused = (time.time() - start) / 60.0
#eprint("Time used: "+str(round(timeused)) + " min ("+str(round(timeused/60,1))+" hours)\n",args.quiet)
timefinish = d_.strftime("%Y-%m-%d_%H:%M")
eprint("#Finished at: {0}".format(timefinish))
eprint("#Time used: {0} min ({1} hours)".format(timeused,round(timeused/60.0,1)))
