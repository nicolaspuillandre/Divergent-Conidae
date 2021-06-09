# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:59:27 2021
@author: Sasha
READS:
    cluster file output by CD-Hit_ss
    a tab-delimmited file with all ORFs containing the expression data
DOES:
    retrieves content and identity data for each cluster
    picks up clusters that contain no less than a pre-defined number of constituient ORFs from no less than 2 taxa
    among these ORFs identifies those with high expression and writes down all clusters containing at least one such ORF in a separate file
    (potentially) writes list of clusters for each taxon in a separate file to be used for Venn diagram
"""

import sys, os
import numpy as np

Clength = 3

CDHITdir = 'E:\\Dropbox\\CONIDAE_venom_evolution\\CD-Hit_ss'
Newwritedir = CDHITdir + '\\SuperNewClusts'
filelist = os.listdir(CDHITdir)
ORFdata = {}
ORFile = open(CDHITdir + '\\allORFdata.txt', 'r')
ORFile.readline()
for line in ORFile:
    line=line.rstrip()
    data = line.split('\t')
    ORFdata[data[0]] = data[1:]
ORFile.close()
print ORFdata.keys()[0], ORFdata[ORFdata.keys()[0]]

keyword = 'Conidae'
for filename in filelist:
    #if filename.split('.')[1] == 'clstr' and keyword in filename:
    if filename == 'Conidae_ss065.clstr':
        clusters = {}
        read1 = open(CDHITdir+'\\'+filename, 'r')
        for line in read1:
            line = line.rstrip()
            if line.startswith('>Cluster'):
                clID = line[1:].replace(' ', '_')
                clusters[clID] = []
            else:
                ORFID = line.split('>')[1].split('...')[0]
                identity = line.split()[-1]
                clusters[clID].append([ORFID, identity])
        read1.close()
        print '\n****************************************************************************'
        print filename, len(clusters), len([key for key in clusters if len(clusters[key]) >1])
        
        write1 = open(CDHITdir+'\\'+filename+'1', 'w')
        for key in clusters:
            if len(clusters[key]) == 1:
                clusters[key][0][1] = '*1'
            for thing in clusters[key]:
                print >>write1, thing[0]+'\t'+key+'\t'+thing[1]
        write1.close()
        
        qclusters = {}
        nclusters = {}
        for key in clusters:
            gens = list(set([thing[0][:3] for thing in clusters[key]]))
            minsim = np.mean([float(thing[1][:-1]) for thing in clusters[key] if '.' in thing[1]])
            if len(clusters[key]) >= Clength and len(gens) >= 2:
                qclusters[key] = [gens, minsim]
                nclusters[key] = len(clusters[key])
        
        write3 = open(CDHITdir+'\\'+filename.split('.')[0]+'_clusters'+str(Clength)+'plus.txt', 'w')
        write4 = open(CDHITdir+'\\'+filename.split('.')[0]+'_clusters'+str(Clength)+'plus_HExpressed.txt', 'w')
        for key in sorted(nclusters, key = nclusters.get, reverse=True):
            print >>write3, key, len(clusters[key]), qclusters[key][0], qclusters[key][1]
            exp = []
            for item in clusters[key]:
                ORF = item[0]
                if ORFdata[ORF][2] == 'q':
                    count = 0
                elif 'tpm:' in ORFdata[ORF][2]:
                    count = int(ORFdata[ORF][2].split(':')[-1].split(',')[0])
                else:
                    count = int(ORFdata[ORF][2])
                exp.append(count)
            if max(exp) >=100:
                alnfile = open(Newwritedir+'\\'+key+'_'+str(max(exp))+'.fas', 'w')
                for thing in clusters[key]:
                    ORF = thing[0]
                    print >>alnfile, '>'+ORF+'|'+ORFdata[ORF][0]+'|'+ORFdata[ORF][2]+'\n'+ORFdata[ORF][3]
                alnfile.close()
                print >>write4, key, len(clusters[key]), qclusters[key][0], qclusters[key][1]
                
        write3.close()
        write4.close()
        """
        write2 = open(CDHITdir+'\\'+filename+'2', 'w')
        clust = {}
        for thing in ['Pfc', 'Cll', 'Clc', 'Pyg']:
            clust[thing] = []

        for cl in clusters.keys():#[:3]:
            if len(clusters[cl]) >= 5:
                thinglist = []
                for item in clusters[cl]:
                    if item[0][:3] not in thinglist:
                        thinglist.append(item[0][:3])
                #print thinglist
                for thing in ['Pfc', 'Cll', 'Clc', 'Pyg']:
                    if thing in thinglist:
                        clust[thing].append(cl)
        for thing in ['Pfc', 'Cll', 'Clc', 'Pyg']:
            print >>write2, '\n*******************' + thing
            for item in clust[thing]:
                print >>write2, item
        write2.close()
        """

                