#!/bin/python
import sys
with open(sys.argv[1]) as fid:
 seq = ""
 head = fid.readline().rstrip()
 line = fid.readline().rstrip()
 while line and line[0] != ">":
  seq += line
  line = fid.readline().rstrip()
 print("{0}|{1}".format(head,len(seq)))
 print(seq)
