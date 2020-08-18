#############################################
###         Configuration options         ###
#############################################

# Specify directory containing fastq files:
baseDir = "/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/fastq/"

# Specify wildtype (Wt) and mutant (Mut) sequences to look for.
sequenceWt  = "TATTGCCCAATAAACATTG"
sequenceMut = "TATTGCCCAATATACATTG"

# Specify part of the filename that indicates first or second read of paired end reads.
# Script will use these to try and match the paired files together.
# (Usually '_R1_' & '_R2_', include underscores to improve matching).
firstReadOfPair = '_R1_'
secondReadOfPair = '_R2_'

#############################################
###         Load important packages       ###
#############################################

# Import RegEx functionalities.
import re

# Import file system functionalities.
import os
import fileinput
import fnmatch
import gzip

# Import debugging tools
# TODO: Remove for release.
from pprint import pprint


#############################################
###     Define functions for analysis     ###
#############################################

def totalReads(linesFromFastq):
    """
    Return total number of reads in input lines from fastq-file
    :param linesFromFastq: Lines from sam-file
    :return: String
    """
    reads = len(range(1, len(linesFromFastq), 4))
    return reads


def getSeq(fastqInput, wtSequence, mutSequence):
    """
    Return all reads from fastqInput matching any of search input strings from the other 2 arguments
    :param fastqInput:
    :param wtSequence: Wild type DNA sequence of interest
    :param mutSequence: Mutant DNA sequence of interest
    :return: Dict containing a list of tuples for wild type hits, list of tuples for mutant hits and the total number
                of reads in this file. Each tuple is a hit and consists of both reads of the pair.
    """
    # Get reverse complimentary sequence of provided sequence to look for reverse reads too.
    # TODO: Make doing this configurable, it might not be desirable under all circumstances.
    wtSequenceRv = reverseSeq(wtSequence)
    mutSequenceRv = reverseSeq(mutSequence)

    # Create lists to save hits to.
    wtHits = []
    mutHits = []

    for i in range(1, len(fastqInput[1]), 4):
        # Iterate over every sequencing read in the first read file.
        if len(re.findall(wtSequence+"|"+wtSequenceRv, fastqInput[1][i])) > 0:
            # Save read and paired read to list of wild type hits.
            wtHits.append((fastqInput[1][i],fastqInput[0][i]))
        elif len(re.findall(mutSequence+"|"+mutSequenceRv, fastqInput[1][i])) > 0:
            # Save read and paired read to list of mutant hits
            mutHits.append((fastqInput[1][i], fastqInput[0][i]))

    return {'wtHits': wtHits, 'mutHits': mutHits, 'totalReads': i/4}


def reverseSeq(sequence):
    """Returns reverse complementary sequence of provided DNA sequence"""
    sequenceRev = ""
    for i in sequence:
        if i == "A" or i == "a":
            sequenceRev = "T" + sequenceRev
        elif i == "T" or i == "t":
            sequenceRev = "A" + sequenceRev
        elif i == "C" or i == "c":
            sequenceRev = "G" + sequenceRev
        elif i == "G" or i == "g":
            sequenceRev = "C" + sequenceRev
    return (sequenceRev)


def findUniqueReads(pairedReads):
    """
    Filters pairedReads to exclude duplicate reads by looking at the combination of read sequence + UMI
    :param pairedReads: List of tuples of format (read1, read2)
    :return: Filtered copy of pairedReads containing unique reads.

    TODO: Add logic dealing with the uncertainty regarding read accuracy (1 mismatch might mean a read error instead of an actual difference between transcripts).
    """
    return pairedReads


def assignToCells(pairedReads):
    """
    Iterates over pairedReads, splits it into new lists (of tuples) based on the cell barcode
    :param pairedReads: List of tuples of format (read1, read2)
    :return: Dict containing lists of tuples (read1, read2), every list indexed under cell barcode.

    TODO: Add logic dealing with the uncertainty regarding read accuracy (1 mismatch might mean a read error instead of an actually difference in the barcode).
    """
    return {'ATGC': pairedReads[0:99], 'ATGG': pairedReads[100:199]}

#############################################
###                Load files             ###
#############################################

# Obtain directory list of baseDir
baseDirContent = os.listdir(baseDir)

# Discover .fastq(.gz) files that belong to a pair, open them and add their lines as a tuple to dict fileList.
fileList = {}
includedFiles = []
for fileName in baseDirContent:
    # Check whether current file is a fastq-file (.fastq.gz also allowed)
    if not (fileName.endswith('.fastq.gz') | fileName.endswith('.fastq')):
        # File is not a fastq-file, so skip to next iteration of the loop
        continue

    # Check if this file isn't included in the list yet
    if fileName in includedFiles:
        continue

    # Check if this is a part of two paired files
    if re.search(firstReadOfPair, fileName):
        # This is the first read. Generate name of second read.
        firstReadFile = fileName
        secondReadFile = re.sub(firstReadOfPair,secondReadOfPair,fileName)
    elif re.search(secondReadOfPair, fileName):
        # This is the second read. Check whether first read file exists
        secondReadFile = fileName
        firstReadFile = re.sub(secondReadOfPair,firstReadOfPair,fileName)
    else:
        # Files don't match the pattern. Skip to next iteration.
        continue

    # Check whether both filenames exist
    if firstReadFile in baseDirContent and secondReadFile in baseDirContent:
        # Add files to included files
        includedFiles.append(firstReadFile)
        includedFiles.append(secondReadFile)
    else:
        # No matching file is available to make this a pair. Don't include this file. Skip to next iteration.
        continue

    # Open both files of pair and add them to fileList using the filename firstReadFile as index
    if fileName.endswith('.fastq.gz'):
        fileList[firstReadFile] = (gzip.open(baseDir+firstReadFile, 'rn').readlines(), gzip.open(baseDir+secondReadFile, 'rn').readlines())
    else:
        fileList[firstReadFile] = (open(baseDir+firstReadFile, 'rn').readlines(), open(baseDir+secondReadFile, 'rn').readlines())


#############################################
###              Process data             ###
#############################################
readNums = {}
for fileName in fileList:
    readNums[fileName] = totalReads(fileList[fileName])

hits = {}
for fileName in fileList:
    hits[fileName] = getSeq(fileList[fileName], sequenceWt, sequenceMut)


pprint(hits)
