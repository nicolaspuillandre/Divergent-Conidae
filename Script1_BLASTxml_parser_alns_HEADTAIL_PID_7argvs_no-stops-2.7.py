"""
modified version of Neda's script to parse BLAST output xml file in a readable and informative txt file.
#USAGE: python Script1_BLASTxml_parser_alns_HEADTAIL_PID_7argvs_no-stops-2.7.py argv[1], argv[1], ..., argv[7]
argv[1] - BLAST.xml,
argv[2] - BLAST database,
argv[3] - query assembly,
argv[4] - parsed output,
argv[5] - threshold PID,
argv[6] - contigs with BLAST hit,
argv[7] - blast-table like output.
only the alignments with PID >= user-defined threshold (argv[5]), and without stop-codons (line 72) are printed to a human readable file
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from sys import argv

def writeIDline(query, ref):#both are strings of equal length
    IDline = ''
    for i in range(len(query)):
        if query[i] == ref[i]:
            IDline = IDline + '|'
        else:
            IDline = IDline + ' '
    return IDline

PID_threshold = float(argv[5])#0.65

blast_records = NCBIXML.parse(open(argv[1]))
queries = []
refs = []
PIDS = {}
for blast_record in blast_records:
    if len(blast_record.alignments) > 0:
        PIDS[blast_record.query] = 0
        queries.append(blast_record.query)
        for alignment in blast_record.alignments:
            ref = alignment.title.split()[1]
            for hsp in alignment.hsps:
                if float(hsp.identities)/hsp.align_length >= PID_threshold:
                    PIDS[blast_record.query] += 1
            if not ref in refs:
                refs.append(ref)
print 'Blast records:', len(queries)
dbase = open(argv[2], 'r')
refseqs = {}
for line in dbase:
    line = line.strip()
    if line.startswith('>') and line[1:] in refs:
        refseqs[line[1:]] = dbase.next().strip()
    else:
        continue
dbase.close()
print 'Database hits:', len(refseqs)
assembly = open(argv[3], 'r')
qseqs = {}
for line in assembly:
    line = line.strip()
    if line.startswith('>') and line[1:] in queries:
        qseqs[line[1:]] = Seq(assembly.next().strip(), IUPAC.unambiguous_dna)
    else:
        continue
print 'contigs with BLAST hits:', len(qseqs)
assembly.close()
blast_records = NCBIXML.parse(open(argv[1]))
outfile = open(argv[4], 'w')
qcontigs = open(argv[6], 'w')
coordinates = open(argv[7], 'w')
print >> coordinates, 'blast_query'+'\t'+'blast_query_length'+'\t'+'blast_hit'+'\t'+'alignment_length'+'\t'+'PID'+'\t'+'query_start'+'\t'+'query_end'+'\t'+'frame'
for blast_record in blast_records:
    if len(blast_record.alignments) > 0 and PIDS[blast_record.query] > 0:
        #print blast_record.query
        print >> outfile, '*** Transcript name:  %s  Length: %d' % (blast_record.query, blast_record.query_length)
        print >> outfile, '\n'
        print >> qcontigs, '>'+ blast_record.query + '\n' + qseqs[blast_record.query]
        for alignment in blast_record.alignments:
            printing = False
            for hsp in alignment.hsps:
                if float(hsp.identities)/hsp.align_length >= PID_threshold and hsp.query.count('*') == 0:
                    printing = True
            if printing == True:
                print >> coordinates, blast_record.query+'\t'+str(blast_record.query_length)+'\t'+alignment.title+'\t'+str(alignment.length)+'\t'+str(float(hsp.identities)/hsp.align_length)+'\t'+str(hsp.query_start)+'\t'+str(hsp.query_end)+'\t'+str(hsp.frame)
                print >> outfile, ' BLAST Hit:  %s  Length: %d' % (alignment.title, alignment.length)
                for hsp in alignment.hsps:
                    print >> outfile, ' Frame: %s       Query start: %d    Query end: %d    Hit start: %d    Hit end: %d    Aligned Length: %d    Pident: %f    Score: %f    E-value: %f' % \
                                      (hsp.frame, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, float(hsp.identities)/hsp.align_length, hsp.score, hsp.expect)
                    missingQhead = ''
                    emptyhead = ''
                    missingRhead = ''
                    missingQtail = ''
                    emptytail = ''
                    missingRtail = ''
                    if hsp.sbjct_start > 1:
                        missinghead = hsp.sbjct_start -1# -1 corrects for python indexing
                        if hsp.frame[0] in [1,2,3]:
                            queryAAhead = qseqs[blast_record.query][hsp.frame[0]-1: hsp.query_start-1].translate()#all AA add '-1'
                        else:
                            queryAAhead = qseqs[blast_record.query][hsp.query_end-1:].reverse_complement().translate()
                            #print blast_record.query, queryAAhead
                            #queryAAhead = qseqs[blast_record.query][abs(hsp.frame[0])-1: hsp.query_end-1].reverse_complement().translate()
                        if len(queryAAhead) >= missinghead:
                            missingQhead = queryAAhead[-missinghead:]
                        else:
                            missingQhead = '-'*(missinghead-len(queryAAhead)) + queryAAhead
                        missingRhead = refseqs[alignment.title.split()[1]][:hsp.sbjct_start-1]
                        emptyhead = writeIDline(missingQhead, missingRhead)#' '*missinghead
                    
                    if alignment.length - hsp.sbjct_end > 0:
                        missingtail = alignment.length - hsp.sbjct_end
                        if hsp.frame[0] in [1,2,3]:
                            queryAAtail = qseqs[blast_record.query][hsp.query_end:].translate()
                        else:
                            queryAAtail = qseqs[blast_record.query][abs(hsp.frame[0]+1):hsp.query_start-1].reverse_complement().translate()
                        if len(queryAAtail) >= missingtail:
                            missingQtail = queryAAtail[:missingtail-1]
                        else:
                            missingQtail = queryAAtail + '-'*(missingtail-len(queryAAtail)-1)
                        missingRtail = refseqs[alignment.title.split()[1]][hsp.sbjct_end+1:]
                        emptytail = writeIDline(missingQtail, missingRtail)#' '*(missingtail-1)
                        
                    print >> outfile, '\n'
                    print >> outfile, '    Q:  ' + missingQhead+' %s ' % hsp.query + missingQtail 
                    print >> outfile, '        ' + emptyhead +' %s ' % hsp.match + emptytail 
                    print >> outfile, '    R:  ' + missingRhead+' %s ' % hsp.sbjct + missingRtail 
                    print >> outfile, '\n'
outfile.close()
qcontigs.close()
coordinates.close()
