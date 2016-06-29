#!/bin/python

from __future__ import print_function
import sys
import os

if "-h" in sys.argv:
 print("usage: python idbastat.py [scaffold.fa] [10]")
 sys.exit()
if len(sys.argv) > 2:
 try:
  n = int(sys.argv[2])
 except:
  print("Error in number argument.")
  sys.exit()
 if os.path.isfile(sys.argv[1]):
  myfile = sys.argv[1]
 else:
  print("Error in input file")
  sys.exit()
else:
 n = 10

if len(sys.argv) == 1:
 myfile = "scaffold.fa"
 n = 10

if len(sys.argv) == 2:
 if os.path.isfile(sys.argv[1]):
  myfile = sys.argv[1]
 else:
  if os.path.isfile("scaffold.fa"):
   myfile = "scaffold.fa"
  else:
   print("Can't read input file: {0}".format(sys.argv[1]))
   sys.exit()
  try:
   n = int(sys.argv[1])
  except:
   print("Something wrong with arguments")
   n = 10
   sys.exit()

with open(myfile) as fid:
#with open("scaffold.fa") as fid:
 d = {}
 seq = 0
 entry = fid.readline().rstrip()[1:]
 line = fid.readline().rstrip()
 
 while line:
  if line[0] == ">":
   d[entry] = seq
   seq = 0
   entry = line[1:]
  else:
   seq += len(line)
  line = fid.readline().rstrip()
 else:
  d[entry] = seq
import operator
sorted_d = sorted(d.items(), key=operator.itemgetter(1),reverse=True)

print("Number of scaffolds: {0}\n".format(len(d)))
if n > len(d):
 n = len(d)
for i in range(0,n):
 print("{0}:\t\t{1}".format(sorted_d[i][0],sorted_d[i][1]))
