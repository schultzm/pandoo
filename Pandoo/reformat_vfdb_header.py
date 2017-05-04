#!/usr/bin/env python3

'''
    Uses python3.

    This file is a crude python script for modifying the headers
    in a VFDB file.

    Copyright (C) 2017 Mark B Schultz
    https://github.com/schultzm/
    email: dr.mark.schultz@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


# Modules import.
import argparse
import os
import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna

# Set up the arguments PARSER to deal with the command line input.
PARSER = argparse.ArgumentParser(description="replaces fasta header inside \
                                 the file with that file's filename")
PARSER.add_argument('-i', '--infile', help="fasta formatted infile",
                    required=True)
ARGS = PARSER.parse_args()


def readconvert(infile):
    '''
    Convert the fasta headers in the VFDB formatted fasta infile.
    '''
    for record in SeqIO.parse(infile, 'fasta', generic_dna):
        parts = record.description.split()
        record.id = parts[1].replace('(', '').replace(')', '') +\
                    '-'+parts[0].replace('_', '-')
        record.description = ' '.join(parts[2:])
        SeqIO.write(record, sys.stdout, 'fasta')


def main():
    '''
    Main function.
    '''
    if os.path.exists(os.path.abspath(ARGS.infile)):
        readconvert(os.path.abspath(ARGS.infile))
    else:
        sys.stderr.write('\nFile '+os.path.abspath(ARGS.infile)+\
                         ' not found.\n')

if __name__ == '__main__':
    main()
