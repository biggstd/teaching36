"""
#####################
Hort 503 -- Project 2
#####################

By: Tyler Biggs
"""

# Import standard packages.
from sys import argv
import collections
import re

# Import data science packages.
from tqdm import tqdm
import matplotlib.pyplot as plt



def pass_simple_tests(row):
    """Checks a row for a series of characters, returns False if they are found,
    and True otherwise. Also ensures that there are at least 10 reads.

    :param row:
        A single line from the given pileup file.

    :returns:
        True or False
    """
    # Scrape the carrot and the character immediately following it.
    nucleotide_string = row[4]
    nucleotide_string = re.sub()

    bad_symbols = ('+', '*', '-', '>', '<')
    minimum_reads = 10

    if any(c in bad_symbols for c in nucleotide_string):
        return False
    elif row[3] < minimum_reads:
        return False

    return True


def is_snp(row):

    quality = get_quality(row)
    overlap_reads = row[4]

    counter = collections.Counter()

    for q, r in zip(quality, overlap_reads):

        if q < 30:
            continue
        counter += collections.Counter(r)

    counter['A'] += counter['a']
    counter['T'] += counter['t']
    counter['C'] += counter['c']
    counter['G'] += counter['g']
    counter[','] += counter['.']

    considered_counts = [
        counter['A'],
        counter['T'],
        counter['C'],
        counter['G'],
        # counter[','],
    ]

    if any(c > 3 for c in considered_counts):
        return counter


def get_quality(row):
    chars = row[5]
    return [ord(c) - 33 for c in chars]


def main(argv):
    """Main processing function.
    """

    # Unpack the supplied arguments.
    script, pileup_in, pileup_out = argv

    # Check how many lines are in the file.
    line_count = 0
    for line in open(pileup_in).readlines():
        line_count += 1

    # Create the cool progress bar.
    pbar = tqdm(total=line_count)

    # Open the input and output files.
    with open(pileup_in, 'r') as pin, open(pileup_out, 'w') as pout:

        # Read the input file line by line.
        line = pin.readline().split()

        while line:
            # Update the progress bar.
            pbar.update(1)

            # Ensure the number of reads is an integer.
            line[3] = int(line[3])

            # Run the basic test function.
            if pass_simple_tests(line):

                # Run the snip check function.
                counts = is_snp(line)

                if counts:
                    for x in ['A', 'G', 'C', 'T']:
                        if counts[x] > 0:
                            pout.write(f'{line[0]}\t{line[1]}\t{line[2]}\t{x}\t{counts[x] / line[3]}\n')

            # Read the next line and start again.
            line = pin.readline().split()

"""
Protect the functionsn within the script, this allows us to import them while
still being able to call this script from the terminal.
"""
if __name__ == "__main__":
    main(argv)
