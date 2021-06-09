# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:57:11 2021
@author: Sasha

USAGE: python Script2_local8_retrieve_cys_pattern_2.7.py sys.argv[1], sys.argv[2]
sys.argv[1] is a space delimited text file with sequence ID in the first column, and AA sequence in the second one.

"""
import os, sys
matpep = open(sys.argv[1], 'r')
cysfile = open(sys.argv[2], 'w')
for line in matpep:
    data=line.rstrip().split()
    if data[1].count('C') == 0:
        print >>cysfile, data[0]+'\tNoCys'
    else:
        if data[1][0] == 'C':
            CP = '-C'
        else:
            CP = '-'
        for i in data[1][1:]:
            if i == 'C':
                CP = CP + 'C'
            elif i != 'C' and CP[-1] == 'C':
                CP = CP + '-'
            else:
                continue
        if CP[-1] == '-':
            CP = CP[:-1]
        print >>cysfile, data[0]+'\t'+CP[1:]
matpep.close()
cysfile.close()                