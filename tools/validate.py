#!/usr/bin/python

#####
#
# Parsing command line arguments
#
#####
import argparse
parser = argparse.ArgumentParser(description='Run helper for rigid fibers')
parser.add_argument('current_matrix',metavar='CURRENT_MATRIX',type=open,help='the current matrix')
parser.add_argument('reference_matrix',metavar='REFERENCE_MATRIX',type=open,help='the reference matrix')

args = parser.parse_args()

#####
#
# Parsing matrices
#
#####
current_matrix = [[float(col) for col in row.strip().split()] for row in args.current_matrix]
reference_matrix = [[float(col) for col in row.strip().split()] for row in args.reference_matrix]
args.current_matrix.close()
args.reference_matrix.close()

#####
#
# Calculating delta
#
#####
import math
import sys

count = 0
total = 0.0
maximum = -1.0
max_location = None

if len(current_matrix) != len(reference_matrix):
    sys.exit(1)

for idx_row in xrange(len(current_matrix)):
    current_row = current_matrix[idx_row]
    reference_row = reference_matrix[idx_row]

    if len(current_row) != len(reference_row):
        sys.exit(2)

    for idx_col in xrange(len(current_row)):
        count += 1

        current_element = current_row[idx_col]
        reference_element = reference_row[idx_col]

        delta = abs(current_element - reference_element)
        total += delta
        if delta > maximum:
            maximum = delta
            max_location = (idx_row / 3, idx_row % 3)

print '**************************************************'
print 'Validation:'
print '  Comparison Matrix :',args.current_matrix.name
print '  Reference Matrix  :',args.reference_matrix.name
print '  Maximum Delta     :',maximum
if maximum > 1e-5:
    print '  Maximum Location  :',max_location
print '  Average Delta     :',total/count
