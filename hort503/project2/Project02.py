from sys import argv

import numpy as np
import pandas as pd
import matplotlib as matplotlib
import matplotlib.pyplot as plt

def pass_basic_tests(cols):

    # Check rule number 1.
    # Skip lines with a deletion or insertion. These lines contain one
    # of these characters: *, + or -.  We will not try to deal with insertions
    # or deletions.
    if (cols[4].find('+') >= 0):
        return False
    if (cols[4].find('*') >= 0):
        return False
    if (cols[4].find('-') >= 0):
        return False

    # Check rule #2.
    # Skip lines with a reference skip. These lines contain one of these
    # characters: >, or <.
    if (cols[4].find('>') >= 0):
        return False
    if (cols[4].find('<') >= 0):
        return False

    # Check rule #7.
    # There must be at least 10 reads aligned at a position to call a SNP
    # in Winter Dawn.
    if (int(cols[3]) < 10):
        return False

    return True

def is_SNP(cols):
    depth = int(cols[3])
    quality = get_quality(cols)
    numA = numC = numT = numG = numR = 0
    for i in range(0, depth):

        # We can't count bases with a score less than 30.
        if (quality[i] < 30):
            continue

        # Count the number of variants
        if (cols[4][i] == 'A' or cols[4][i] == 'a'):
            numA = numA + 1
        if (cols[4][i] == 'T' or cols[4][i] == 't'):
            numT = numT + 1
        if (cols[4][i] == 'C' or cols[4][i] == 'c'):
            numC = numC + 1
        if (cols[4][i] == 'G' or cols[4][i] == 'g'):
            numG = numG + 1
        if (cols[4][i] == ',' or cols[4][i] == '.'):
            numR = numR + 1

    retval = { 'A': numA, 'T': numT, 'C': numC, 'G': numG, 'Ref': numR }

    if (numA > 3 or numC > 3 or numT > 3 or numG > 3):
        return retval

    return False

def get_quality(cols):
    quality = cols[5]
    qlen = len(quality)
    q = pd.Series(np.empty(qlen, dtype=int))
    for i in range(0, qlen):
        q[i] = ord(quality[i]) - 33
    return q

def main(argv):
    script, mpileup_in, mpileup_out = argv

    fhi = open(mpileup_in, 'r')
    fho = open(mpileup_out, 'w')
    line = fhi.readline().strip()
    while(line):
        cols = line.split("\t")
        cols[3] = int(cols[3])
        # Make sure the line can pass the basic tests:
        #  1. Not an insertion or deletion
        #  2. Not a reference Skip
        #  3. There must be at least 10 reads aligned.
        if (pass_basic_tests(cols)):

            # Is it a SNP? (i.e. has three of the same variants)
            # If not a SNP then skip this line.
            counts = is_SNP(cols)
            if (counts != False):
                # A position could have more than one SNP so we need to write
                # separate entry to the output file for each one.
                if (counts['A'] > 0):
                    fho.write("{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[1], cols[2], 'A', counts['A'] / cols[3]))
                if (counts['G'] > 0):
                    fho.write("{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[1], cols[2], 'G', counts['G'] / cols[3]))
                if (counts['C'] > 0):
                    fho.write("{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[1], cols[2], 'C', counts['C'] / cols[3]))
                if (counts['T'] > 0):
                    fho.write("{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[1], cols[2], 'T', counts['T'] / cols[3]))
                print(cols)

        line = fhi.readline().strip()

    # Now read in the out file we just created as a data frame and pyplot
    #results = pd.read_csv(mpileup_out, sep="\t")
    results = pd.read_csv('out.txt', header=None, sep="\t", names=('Chr', 'Position', 'Reference', 'SNP', 'Frequency'))
    # Add colors
    results['Colors'] = results.SNP.map({'A': 'red', 'G': 'blue', 'T': 'green', 'C': 'purple'})
    ax = results.plot(x='Position', y='Frequency', kind='bar', legend = False, figsize=(16,4), color=results.Colors)
    xticks = list(range(1, len(results), 100))
    ax.set_xticks(xticks)
    ax.set_xticklabels(results.loc[xticks, 'Position'])

    fig = ax.get_figure()
    fig.savefig("project02.figure.png", format='png', dpi=300)


main(argv)
