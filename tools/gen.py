import random
import math
import sys

NUM_OF_FIBERS = 100
if sys.argv[1]:
    NUM_OF_FIBERS = int(sys.argv[1])

MIN_X = -1000
MAX_X = 1000

MIN_Y = -1000
MAX_Y = 1000

MIN_Z = -1000
MAX_Z = 1000

MIN_T = -1
MAX_T = 1

f = open('XcT_gen'+str(NUM_OF_FIBERS)+'.in', 'w')
f.write(str(NUM_OF_FIBERS) + '\n')

for i in xrange(0,NUM_OF_FIBERS):
    x = [str(random.uniform(MIN_X, MAX_X)),str(random.uniform(MIN_Y, MAX_Y)),str(random.uniform(MIN_Z, MAX_Z))]
    f.write('\t'.join(x) + '\n')
    t = [random.uniform(MIN_T, MAX_T),random.uniform(MIN_T, MAX_T),random.uniform(MIN_Z, MAX_T)]
    # normalize t
    t_len = math.sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2])
    t[0] = str(t[0]/t_len);
    t[1] = str(t[1]/t_len);
    t[2] = str(t[2]/t_len);
    f.write('\t'.join(t) + '\n')

f.close()
