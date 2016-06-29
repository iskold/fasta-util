#!/bin/python

#Rotate fa
import sys
if len(sys.argv) < 2:
 print("Usage: shiftfasta.py contig.fa 65")

with open(sys.argv[1]) as fid:
 rot = int(sys.argv[2])-1
 seq = ""
 entry = ""
 for line in fid:
  if line[0] == ">":
   entry = line
  else:
   seq += line.rstrip()
 if rot > 0:
  print(entry+seq[rot:]+seq[0:rot])
 else:
  print(entry+seq[rot+1:]+seq[0:rot+1])
