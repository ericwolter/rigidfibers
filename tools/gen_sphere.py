import math
import random
import sys
from time import sleep

N = 100
sys.argv = sys.argv[1:]

if len(sys.argv) >= 1 and sys.argv[0]:
    N = int(sys.argv[0])

FIBER_LENGTH = 2
if len(sys.argv) >= 2 and sys.argv[1]:
    FIBER_LENGTH = float(sys.argv[1])

RADIUS = 1.5 * FIBER_LENGTH
RADIUS_2 = RADIUS * RADIUS

DIM_X = RADIUS
DIM_Y = RADIUS
DIM_Z = RADIUS

MIN_T = -1
MAX_T = 1

# print "domain: ["+str(-DIM_X)+","+str(-DIM_Y)+","+str(-DIM_Z)+"]-["+str(DIM_X)+","+str(DIM_Y)+","+str(DIM_Z)+"]"
#
f = []
while len(f) < N:
    [x,y,z] = [random.uniform(-DIM_X,DIM_X),random.uniform(-DIM_Y,DIM_Y),random.uniform(-DIM_Z,DIM_Z)]
    if x * x + y * y + z * z <= RADIUS_2:
        f.append([x,y,z])

export = open('XcT_sphere'+str(N)+'.in', 'w')
export.write(str(N) + '\n')

for i in xrange(N):
    export.write('\t'.join([str(x) for x in f[i]]) + '\n')
    t = [random.uniform(MIN_T, MAX_T),random.uniform(MIN_T, MAX_T),random.uniform(MIN_T, MAX_T)]
    # normalize t
    t_len = math.sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2])
    t[0] = str(t[0]/t_len);
    t[1] = str(t[1]/t_len);
    t[2] = str(t[2]/t_len);
    export.write('\t'.join(t) + '\n')

export.close()
