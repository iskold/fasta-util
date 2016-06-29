#!/bin/python

#write output
def writeit(entry,seq):
 print(entry+  "".join(    [A[x] for x in seq.upper()[::-1] ]  ) )

#RevComplt fa
import sys
with open(sys.argv[1]) as fid:
 seq = ""
 A = {"A":"T","T":"A","G":"C","C":"G","N":"N","Y":"Y","W":"W"}
 for line in fid:
  if line[0] == ">":
   if seq:
    writeit(entry,seq)
   entry = line.rstrip()+"_RevCom\n"
  elif line != "\n":
   seq += line.rstrip()
  else:
   continue
 else:
  writeit(entry,seq)
