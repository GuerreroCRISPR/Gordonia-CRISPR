#!/usr/bin/env python2.7


from Bio import pairwise2, SeqIO
from Bio.Seq import reverse_complement
from Bio.pairwise2 import format_alignment
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('spacers', help='spacers file [.fa]')
parser.add_argument('reads', help='comma separeted reads files [1.fastq,2.fastq]') #
parser.add_argument('names', help='Names file of final result')
parser.add_argument('seq', help='.seq file of 1st step')
parser.add_argument('out', help='out name file')
parser.add_argument('-r', type=int, default=3, help="mismatch [Default: 3]")
args = parser.parse_args()

spacersDic = dict() # dict of spacers and seq (from fasta)
readsFilesList = []
readsDict = dict()  # dict of reads 1 or 2 (fastq)
#readsDict1 = dict()  # dict of reads temp (fastq)
readQualDic = dict() # dict of spacersName:readsQual list
readNameDic = dict() # dict of spacerName:readsName list
spacerLen = []
namesDict = dict() # Dict of names from names file {name1: [name1,name2,name3,...]}
seqDict = dict()  # Dict {spacerOcu:readName}
spacerReadsDic = dict() # merge of namesDict and seqDict {spacer:[R1,R2,R3,..]}
empList = []

## Dict of results
aveQBaseDic = dict()
stdQBaseDic = dict()
minQBaseDic = dict()
maxQBaseDic = dict()
aveQPositionDic = dict()

count = 0
nCount = -1
nSumCount = -1

## Make list of reads files:
readsFilesList = args.reads.split(",")
print readsFilesList

## Open and Parse spacers_file
with open(args.spacers,'rt') as spacers_file:
    for n in SeqIO.parse(spacers_file, "fasta"):
        print(n.id),(n.seq),len(n.seq)
        spacersDic[n.id] = str(n.seq)
        readQualDic[n.id] = list()
        readNameDic[n.id] = list()
        spacerLen.append(len(n.seq))

while len(empList) <= max(spacerLen):
    empList.append(0)

## Open and parse names file:
with open(args.names,'rt') as names_file:
    for name in names_file.readlines():
        nameSpl = name[:-1].split("\t")
        namesDict[nameSpl[0]] = nameSpl[1].split(",")

## Open and parse seq file:
with open(args.seq,'rt') as seq_file:
    for seq in seq_file.readlines()[1:]:
        seqSpl = seq[:-1].split("\t")
        seqDict[seqSpl[0]] = seqSpl[2]

## Merge of namesDict and seqDict:
for sp in namesDict:
    spacerReadsDic[sp] = list()
    for spR in namesDict[sp]:
        spacerReadsDic[sp].append(seqDict[spR])  # Dict of spacer:R1,R2,R3,..

namesDict = dict()
seqDict = dict()

for readsin in readsFilesList:
## Open and parse fastq files:
    with open(readsin,'rt') as reads_file1:
        readsDict = SeqIO.to_dict(SeqIO.parse(reads_file1,"fastq"))

    print(len(readsDict))
    for s in spacersDic:
        for rd in spacerReadsDic[s]:
            seqRead = str(readsDict.get(rd,"None"))
            if seqRead != "None":
                seqSpacer = spacersDic[s]
                seqRead = str(readsDict[rd].seq)
                alignment = pairwise2.align.localms(seqRead,seqSpacer,2,-.1,-3,-2, one_alignment_only=True)
                if alignment[0][2] >= (len(seqSpacer)*2)-(args.r*2.1):
                    print(alignment[0]), len(seqSpacer)*2
                    print(format_alignment(*alignment[0]))
                    readNameDic[s].append(readsDict[rd].id)
                    readQualDic[s].append(readsDict[rd].letter_annotations["phred_quality"][alignment[0][3]:alignment[0][4]])
                ## Reverse aligment read 1:
                alignmentR = pairwise2.align.localms(reverse_complement(seqRead),seqSpacer,2,-.1,-3,-2, one_alignment_only=True)
                if alignmentR[0][2] >= (len(seqSpacer)*2)-(args.r*2.1):
                    print(alignmentR[0]), len(seqSpacer)*2
                    print(format_alignment(*alignmentR[0]))
                    readNameDic[s].append(readsDict[rd].id)
                    rR = readsDict[rd].letter_annotations["phred_quality"]
                    rR.reverse()
                    readQualDic[s].append(rR[alignmentR[0][3]:alignmentR[0][4]])

readsDict = dict()  ## clean memory

# Generate out file:
with open(args.out+".resultQ.file.test.csv",'wt') as spacers_out_file:
    spacers_out_file.write("id\tNR\tNR100%\tAvQ100%\tstdQ100%\tAveQ100%List\n")

with open(args.out+".resultQ.file.test.full_report.txt",'wt') as spacers_out_file2:
    spacers_out_file2.write(args.out+"_full_report\n")

with open(args.out+".resultQ.file.test.short_report.txt",'wt') as spacers_out_file3:
    spacers_out_file3.write(args.out+"_short_report\n")

for qval in readQualDic:
    if readQualDic[qval] != []:
        aveQBaseDic[qval] = list()
        stdQBaseDic[qval] = list()
        minQBaseDic[qval] = list()
        maxQBaseDic[qval] = list()
        aveQPositionDic[qval] = list()[0:len(spacersDic[qval])]=empList[0:len(spacersDic[qval])]
        for qqval in readQualDic[qval]:
            if len(qqval) == len(spacersDic[qval]):
                count += 1
                print "\n",qqval, "count:",count
                print "max Q base:",max(qqval), "| min Q base:", min(qqval), "| Average Q base:", round(float(sum(qqval))/len(qqval),2)
                with open(args.out+".resultQ.file.test.full_report.txt",'at') as spacers_out_file2:
                    spacers_out_file2.write(str(qval)+"| ")
                    spacers_out_file2.write("max Q base:" + str(max(qqval)) + "| min Q base:" + str(min(qqval)) + "| Average Q base: " + str(round(float(sum(qqval))/len(qqval),2)) + "\n")
                    spacers_out_file2.write(str(qqval)+" count: "+str(count)+"\n")
                aveQBaseDic[qval].append(round(float(sum(qqval))/len(qqval),2))
                stdQBaseDic[qval].append(round(np.std(qqval),2))
                minQBaseDic[qval].append(min(qqval))
                maxQBaseDic[qval].append(max(qqval))
                for n in qqval:
                    nCount +=1
                    aveQPositionDic[qval][nCount] += n
            nCount = -1
        if count != 0:
            for nSum in aveQPositionDic[qval]:
                nSumCount += 1
                aveQPositionDic[qval][nSumCount] = round(float(nSum)/count,1)
            nSumCount = -1
            print(qval), "| Nro de reads100%len:", count,\
            "| Nro de reads:",len(readQualDic[qval]),\
            "| Q media spacer100%len:",round(float(sum(aveQBaseDic[qval]))/len(aveQBaseDic[qval]),2),\
            "| Q std spacer100%len:", round(np.std(aveQPositionDic[qval]),2)
            print qval,"Q average per base:",aveQPositionDic[qval]

            with open(args.out+".resultQ.file.test.short_report.txt",'at') as spacers_out_file3:
                spacers_out_file3.write(str(qval)+" | Q average per base:"+str(aveQPositionDic[qval])+"\n")

## Output files to visualize on gephi
            with open(args.out+".resultQ.file.test.csv",'at') as spacers_out_file:
                spacers_out_file.write(qval+"\t"+str(len(readQualDic[qval]))+"\t"+str(count)+"\t"+str(round(float(sum(aveQBaseDic[qval]))/len(aveQBaseDic[qval]),2))+"\t"+str(round(np.std(aveQPositionDic[qval]),2))+"\t"+str(aveQPositionDic[qval])+"\n")
            count = 0
