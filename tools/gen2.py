import math
import random
import sys
from time import sleep

def average_distance(f):
    N = len(f)
    count = 0
    total = 0.0
    min_dist = sys.float_info.max

    for i in xrange(N):
        p_i = f[i]

        nearest_dist = sys.float_info.max
        for j in xrange(N):

            if i == j:
                continue

            p_j = f[j]

            dist = (p_i[0] - p_j[0])**2 + (p_i[1] - p_j[1])**2 + (p_i[2] - p_j[2])**2
            dist = math.sqrt(dist)

            if dist < nearest_dist:
                nearest_dist = dist

        total += nearest_dist
        count += 1

        if nearest_dist < min_dist:
            min_dist = nearest_dist

    print "min", min_dist
    print "avg", total/count

N = 1000
if sys.argv[1]:
    N = int(sys.argv[1])

MIN_DISTANCE = 0.2
MIN_DISTANCE2 = MIN_DISTANCE * MIN_DISTANCE
AVG_DISTANCE = MIN_DISTANCE * 1.7 #0.4

DIM_X = (N-1)**(1.0/3) * AVG_DISTANCE
DIM_Y = (N-1)**(1.0/3) * AVG_DISTANCE
DIM_Z = (N-1)**(1.0/3) * AVG_DISTANCE

MIN_T = -1
MAX_T = 1

STEP = 0.1 * 0.2

print "domain: ["+str(-DIM_X)+","+str(-DIM_Y)+","+str(-DIM_Z)+"]-["+str(DIM_X)+","+str(DIM_Y)+","+str(DIM_Z)+"]"

f = []
for i in xrange(N):
    f.append([random.uniform(-DIM_X,DIM_X),random.uniform(-DIM_Y,DIM_Y),random.uniform(-DIM_Z,DIM_Z)])

print "before:"
average_distance(f)

optimal = False
while not optimal:
    optimal = True
    for i in xrange(N):
        p_i = f[i]

        if p_i[0] > DIM_X:
            p_i[0] -= random.uniform(0, STEP)
            optimal=False
        elif p_i[0] < -DIM_X:
            p_i[0] += random.uniform(0, STEP)
            optimal=False

        if p_i[1] > DIM_Y:
            p_i[1] -= random.uniform(0, STEP)
            optimal=False
        elif p_i[1] < -DIM_Y:
            p_i[1] += random.uniform(0, STEP)
            optimal=False

        if p_i[2] > DIM_Z:
            p_i[2] -= random.uniform(0, STEP)
            optimal=False
        elif p_i[2] < -DIM_Z:
            p_i[2] += random.uniform(0, STEP)
            optimal=False

        nearest_dist = sys.float_info.max
        for j in xrange(N):
            if i==j:
                continue

            p_j = f[j]
            dist = (p_i[0] - p_j[0])**2 + (p_i[1] - p_j[1])**2 + (p_i[2] - p_j[2])**2

            if dist < nearest_dist:
                nearest_dist = dist

        if nearest_dist < MIN_DISTANCE2:
            p_i[0] += random.uniform(-STEP,STEP)
            p_i[1] += random.uniform(-STEP,STEP)
            p_i[2] += random.uniform(-STEP,STEP)
            optimal = False

        f[i] = p_i
        #print i, p_i, nearest_dist < MIN_DISTANCE, optimal
        # sleep(0.1)
    # sleep(1)

print "after:"
average_distance(f)

export = open('XcT_gen'+str(N)+'.in', 'w')
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
