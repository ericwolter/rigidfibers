#!/usr/bin/python

#####
#
# Parsing command line arguments
#
#####
import argparse
parser = argparse.ArgumentParser(description='Run helper for rigid fibers')
parser.add_argument('current_map',metavar='CURRENT_MAPPING',type=open,help='the current mapping')
parser.add_argument('current_matrix',metavar='CURRENT_MATRIX',type=open,help='the current matrix')
parser.add_argument('reference_map',metavar='REFERENCE_MAPPING',type=open,help='the reference mapping')
parser.add_argument('reference_matrix',metavar='REFERENCE_MATRIX',type=open,help='the reference matrix')

args = parser.parse_args()

#####
#
# Parsing mapping
#
#####

# Indexed by row and column and returns simulation variables
current_mapping = {}
for idx,line in enumerate(args.current_map):
    c = [int(component) for component in line.strip().split("|")[1:]]
    current_mapping[(c[6],c[7])] = (c[0],c[1],c[2],c[3],c[4],c[5])

# Indexed by simulation variables and returns row and column
reference_mapping = {}
for idx,line in enumerate(args.reference_map):
    c = [int(component) for component in line.strip().split("|")[1:]]
    reference_mapping[(c[0],c[1],c[2],c[3],c[4],c[5])] = (c[6],c[7])

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

count = 0
total = 0.0
maximum = -1.0
max_location = None

is_vector = len(current_matrix[0]) == 1

if not is_vector:
    for idx_row in xrange(len(current_matrix)):
        row = current_matrix[idx_row]
        for idx_col in xrange(len(row)):
            count += 1

            current_element = row[idx_col]

            # get simulation variables
            current_location = current_mapping[(idx_row,idx_col)]

            # if element was not explicitly set it should always be zero.
            # this should only occur on the diagonal
            if current_location[0] == -1:
                reference_element = 0.0
            else:
                # translate to row and column in reference matrix
                reference_location = reference_mapping[current_location]
                reference_element = reference_matrix[reference_location[0]][reference_location[1]]

            delta = abs(current_element - reference_element)
            total += delta
            if delta > maximum:
                maximum = delta
                max_location = current_location
else:
    for idx_row in xrange(len(current_matrix)):
        count +=1

        current_element = current_matrix[idx_row][0]
        current_location = current_mapping[(idx_row,idx_row)]

        reference_location = reference_mapping[current_location]
        reference_element = reference_matrix[reference_location[0]][0]

        delta = abs(current_element - reference_element)
        total += delta
        if delta > maximum:
            maximum = delta
            max_location = current_location


print '**************************************************'
print 'Validation:'
print '  Comparison Matrix :',args.current_matrix.name
print '  Reference Matrix  :',args.reference_matrix.name
print '  Maximum Delta     :',maximum
if maximum > 1e-5:
    print '  Maximum Location  :',max_location
print '  Average Delta     :',total/count
