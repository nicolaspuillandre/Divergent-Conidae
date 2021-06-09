# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 15:31:06 2021
@author: Sasha
READS:
    BLAST database
    for each taxon of the taxalist:
        Tab delimited BLASTp output
        'genes.results' file with TPM expression values from RSEM
        assembly (not interleaved format)
DOES:
    identifies BLAST hits that fulfil the set criteria
    retrieves expression data for these hits and sums them up by gene superfamily
    extracts the contigs with hits and writes down to output files (sorted by gene superfamily) their strings that correspond to the BLAST ccordinates
        
USAGE: runs directly from the Spyder environment 
"""
import matplotlib.pyplot as plt
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

taxalist = ['arenatus', 'coronatus', 'ebraeus', 'ermineus', 'imperialis', 'lividus', 'magus3', 'marmoreus', 'quercinus', 'rattus', 'sponsalis', 'textile', 'varius', 'virgo', 'praecellens', 'tribblei']
print len(taxalist)

lookups = {}
lookups['arenatus'] = '|A|'
lookups['coronatus'] = '|A|'
lookups['ebraeus'] = '|SF-mi2|'
lookups['ermineus'] = '|A|'
lookups['imperialis'] = '|A|'
lookups['lividus'] = '|B2|'
lookups['magus3'] = '|A|'
lookups['marmoreus'] = '|B2|'
lookups['quercinus'] = '|B2|'
lookups['rattus'] = '|B4|'
lookups['textile'] = '|A|'
lookups['sponsalis'] = '|A|'
lookups['varius'] = '|B2|'
lookups['virgo'] = '|A|'
lookups['praecellens'] = '|S|'
lookups['tribblei'] = '|B2|'
ngsfs = ['Divergent_MWSGK', 'New-Geo-1', 'DivConN1', 'DivConN2', 'DivConN3', 'DivConN4', 'DivConN5', 'DivConN6', 'DivConN7', 'DivConN8']
sppexp = []

print ngsfs
### READ IN BLAST DB
tox = open('E:\\Dropbox\\CONIDAE_venom_evolution\\DC_new_clusters_mp20plus3S.fas', 'r')

toxins = {}
for line in tox:
    line=line.rstrip()
    if line.startswith('>'):
        ID = line[1:]
        toxins[ID] = ''
    else:
        toxins[ID] += line
tox.close()
print len(toxins)

ngsfsexp = open('E:\\Dropbox\\CONIDAE_venom_evolution\\NEW_GSFS_summed_tpms_ID03_new.txt', 'w')
print >>ngsfsexp, 'taxon\t'+str(ngsfs).replace(',', '\t')
### READ IN BLAST RESULTS AND FIND BEST HITS FOR EACH DB ENTRY (but each time different contig)
for taxon in taxalist:
    checked_contigs = []
    besthits = {}
    coords = {}
    blast = open('E:\\BLASTx_NewClusters\\'+taxon+'.blastx.out', 'r')
    for line in blast:
        line = line.rstrip()
        data = line.split('\t')
        contig = data[0]
        hit = data[1] 
        PID = float(data[2])
        length = int(data[3])
        score = float(data[-1])
        goodhit = False
        if lookups[taxon] in hit and PID >=85 and length >= 0.6*len(toxins[hit]):
            goodhit = True
        if taxon == 'tribblei' and lookups[taxon] in hit and PID >=60 and length >= 0.6*len(toxins[hit]):
            goodhit = True
        if 'DivCon' in hit and PID >=30 and length >= 0.5*len(toxins[hit]):
            goodhit = True
        if 'New-Geo-1' in hit and PID >=30 and length >= 0.5*len(toxins[hit]):
            goodhit = True
        if 'Divergent_MWSGK' in hit and PID >=30 and length >= 0.5*len(toxins[hit]):
            goodhit = True
        if contig in checked_contigs:
            goodhit = False
        if goodhit == True:
            if not hit in besthits.keys():
                besthits[hit] = [[contig, hit, PID, length, int(data[4]), int(data[5]), int(data[6]), int(data[7]), int(data[8]), int(data[9]), float(data[10]), float(data[11]), len(toxins[hit])]]
            else:
                besthits[hit].append([contig, hit, PID, length, int(data[4]), int(data[5]), int(data[6]), int(data[7]), int(data[8]), int(data[9]), float(data[10]), float(data[11]), len(toxins[hit])])
            coords[contig] = [int(data[6]), int(data[7]), int(data[8]), hit]
            checked_contigs.append(contig)
            
    blast.close()
    
    ### READ IN EXPRESSION DATA
    tpms = {}
    exp = open('E:\\BLASTx_NewClusters\\'+taxon+'.genes.results', 'r')
    for line in exp:
        line = line.rstrip()
        data = line.split('\t')
        if data[0] in checked_contigs:
            tpms[data[0]] = float(data[-2])
    exp.close()
    
    ### SUMM UP EXPRESSION PER GSF
    scaletpm = 0
    for key in besthits:
        if lookups[taxon] in key:
            for thing in besthits[key]:
                scaletpm += tpms[thing[0]]
    #print scaletpm
    explist = [0]*len(ngsfs)
    for key in besthits:
        if not lookups[taxon] in key:
            for item in besthits[key]:
                tpm = tpms[item[0]]
                for i in range(len(ngsfs)):
                    if ngsfs[i] in key:
                        explist[i] += tpm
    explist.append(scaletpm)
    print taxon, explist
    print >>ngsfsexp, taxon+'\t'+str(explist).replace(',', '\t')
    
    ### READ IN ASSEMBLY, EXTRACT BLASTED CONTIGS, WRITE THE OUTPUTS
    assembly = open('E:\\BLASTx_NewClusters\\'+taxon+'.trinity.clustered.NI.fasta', 'r')
    qseqs = {}
    for line in assembly:
        line = line.strip()
        if line.startswith('>'):
            contID = line.split()[0][1:]
            toprint = False
            if contID in coords.keys() and 'DivCon' in coords[contID][-1] and tpms[contID] > 5:
                toprint = True
            if contID in coords.keys() and 'New-Geo-1' in coords[contID][-1] and tpms[contID] > 5:
                toprint = True
            if contID in coords.keys() and 'Divergent_MWSGK' in coords[contID][-1] and tpms[contID] > 5:
                toprint = True
            if toprint == True:
                if 'New-Geo-1' in coords[contID][-1]:
                    sfID = 'New-Geo-1'
                elif 'Divergent_MWSGK' in coords[contID][-1]:
                    sfID = 'Divergent_MWSGK'
                else:
                    sfID = coords[contID][-1].split('_')[-1].split('-')[0]
                seq = assembly.next().strip()
                if sfID in ngsfs:
                    #seq = Seq(seq, IUPAC.unambiguous_dna)
                    #print len(seq)
                    if coords[contID][0] < coords[contID][1]:
                        #NTseq = seq[coords[contID][0]-1:coords[contID][1]]
                        NTseq = seq[coords[contID][0]-1:]
                        AAseq = Seq(NTseq, IUPAC.unambiguous_dna).translate()
                    else:
                        #NTseq = Seq(seq[coords[contID][1]-1:coords[contID][0]], IUPAC.unambiguous_dna).reverse_complement()
                        NTseq = Seq(seq[:coords[contID][0]], IUPAC.unambiguous_dna).reverse_complement()
                        AAseq = NTseq.translate()
                    newseqsAA = open('E:\\BLASTx_NewClusters\\New_'+sfID+'_AA_30percent_cutoff.fasta', 'a')
                    newseqsNT = open('E:\\BLASTx_NewClusters\\New_'+sfID+'_NT_30percent_cutoff.fasta', 'a')
                    print >>newseqsAA, '>'+taxon+'|'+sfID+'|'+contID+'|tpm-'+str(tpms[contID])+'\n'+AAseq
                    print >>newseqsNT, '>'+taxon+'|'+sfID+'|'+contID+'|tpm-'+str(tpms[contID])+'\n'+NTseq
                    newseqsAA.close()
                    newseqsNT.close()
    assembly.close()
ngsfsexp.close()








