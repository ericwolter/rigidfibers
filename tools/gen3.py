import math
import random
import sys
from time import sleep

class Point(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(a,b):
        return Point(a.x+b.x,a.y+b.y,a.z+b.z)

    def __sub__(a,b):
        return Point(a.x-b.x,a.y-b.y,a.z-b.z)

    def __mul__(a,b):
        return Point(a.x*b,a.y*b,a.z*b)

    def __rmul__(a,b):
        return Point(a.x*b,a.y*b,a.z*b)

def dot(p,q):
    return p.x * q.x + p.y * q.y + p.z * q.z

def norm(p):
    return math.sqrt(dot(p,p))

def normalize(p):
    l = norm(p)
    return Point(p.x / l, p.y / l, p.z / l)

# http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
def distanceSegments(p_i, o_i, p_j, o_j):
    S1P0 = p_i - o_i
    S1P1 = p_i + o_i
    S2P0 = p_j - o_j
    S2P1 = p_j + o_j

    # print "segme:", S1P0.x, S1P0.y, S1P0.z
    # print "segme:", S1P1.x, S1P1.y, S1P1.z
    # print "segme:", S2P0.x, S2P0.y, S2P0.z
    # print "segme:", S2P1.x, S2P1.y, S2P1.z

    u = S1P1 - S1P0
    v = S2P1 - S2P0
    w = S1P0 - S2P0

    a = dot(u,u)
    b = dot(u,v)
    c = dot(v,v)
    d = dot(u,w)
    e = dot(v,w)

    X = a * c - b * b
    sc = X
    sN = X
    sD = X
    tc = X
    tN = X
    tD = X

    if X < 1e-5:
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        sN = (b * e - c * d)
        tN = (a * e - b * d)
    if sN < 0.0:
        sN = 0.0
        tN = e
        tD = c
    elif sN > sD:
        sN = sD
        tN = e + b
        tD = c

    if tN < 0.0:
        tN = 0.0
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
          sN = -d
          sD = a
    elif tN > tD:
        tN = tD
        if (-d + b) < 0.0:
            sN = 0.0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = (-d + b)
            sD = a

    if abs(sN) < 1e-5:
        sc = 0.0
    else:
        sc = sN / sD
    if abs(tN) < 1e-5:
        tc = 0.0
    else:
        tc = tN / tD

    dP = w + (sc * u) - (tc * v)
    return (dot(dP,dP), dP)

def stats(N, positions, orientations):
    minimal_distance_segment = sys.float_info.max
    minimal_distance_center = sys.float_info.max
    distance_segment = 0.0
    distance_center = 0.0
    total_segment = 0.0
    total_center = 0.0

    min_domain = Point(sys.float_info.max,sys.float_info.max,sys.float_info.max)
    max_domain = Point(-sys.float_info.max,-sys.float_info.max,-sys.float_info.max)

    for i in xrange(N):
        p_i = positions[i]
        o_i = orientations[i]

        min_domain.x = min(min_domain.x, p_i.x)
        min_domain.y = min(min_domain.y, p_i.y)
        min_domain.z = min(min_domain.z, p_i.z)

        max_domain.x = max(max_domain.x, p_i.x)
        max_domain.y = max(max_domain.y, p_i.y)
        max_domain.z = max(max_domain.z, p_i.z)

        nearest_distance_segment = sys.float_info.max
        nearest_distance_center = sys.float_info.max
        for j in xrange(N):
            if i != j:
                p_j = positions[j]
                o_j = orientations[j]

                distance_segment_2, dP = distanceSegments(p_i, o_i, p_j, o_j)
                distance_segment = math.sqrt(distance_segment_2)
                distance_center = norm(p_i - p_j)

                nearest_distance_segment = min(nearest_distance_segment, distance_segment)
                nearest_distance_center = min(nearest_distance_center, distance_center)

        total_segment += nearest_distance_segment
        total_center += nearest_distance_center

        minimal_distance_segment = min(minimal_distance_segment, nearest_distance_segment)
        minimal_distance_center = min(minimal_distance_center, nearest_distance_center)

    print "Distribution Statistics:"
    print "  - domain: ["+str(min_domain.x)+","+str(min_domain.y)+","+str(min_domain.z)+"]-["+str(max_domain.x)+","+str(max_domain.y)+","+str(max_domain.z)+"]"
    print "  - minimal distance segments: ", minimal_distance_segment
    print "  - minimal distance centers:  ", minimal_distance_center
    print "  - average distance segments: ", total_segment/N
    print "  - average distance centers: ", total_center/N

def sampleBox(side):
    p = Point(random.uniform(-side/2,side/2),random.uniform(-side/2,side/2),random.uniform(-side/2,side/2))
    o = Point(random.uniform(-1,1),random.uniform(-1,1),random.uniform(-1,1))
    o = normalize(o)
    return (p, o)

def sampleSphere(min_radius, max_radius):
    min_radius_2 = min_radius * min_radius
    max_radius_2 = max_radius * max_radius
    hit = False
    while not hit:
        p = Point(random.uniform(-max_radius,max_radius),random.uniform(-max_radius,max_radius),random.uniform(-max_radius,max_radius))

        norm_2 = dot(p,p)
        if norm_2 > min_radius_2 and norm_2 < max_radius_2:
            hit = True

    o = Point(random.uniform(-1,1),random.uniform(-1,1),random.uniform(-1,1))
    o = normalize(o)
    return (p, o)

N = 1000
if len(sys.argv) > 1 and sys.argv[1]:
    N = int(sys.argv[1])

min_distance = 0.2
radius = min_distance
side = min_distance**(1.0/3.0)
min_distance_2 = min_distance * min_distance

# p_i = Point(0.140158015896,-0.0339196372195,0.0525623317914)
# o_i = Point(0.49463284706,-0.275602349505,-0.8242461353)
# p_j = Point(0.00593866115171,-0.044934855849,0.0118355945455)
# o_j = Point(0.195279856614,0.756273544675,-0.624432625048)
#
# print distanceSegments(p_i,o_i,p_j,o_j)

positions = []
orientations = []
# for i in xrange(N):
#     p, o = sampleSphere(RADIUS)
#     positions.append(p)
#     orientations.append(o)

for i in xrange(N):
    # print ".. inserting", i
    inserted = False

    while not inserted:
        retries = 0
        rejected = True
        while rejected and retries < 1000:
            # print ".. try", retries
            retries = retries + 1
            rejected = False
            # ensure sufficient overlap previously generated fibers
            p_i, o_i = sampleSphere(max(0, radius - 4 * min_distance), radius)
            # p_i, o_i = sampleBox(side)

            for j in xrange(len(positions)):
                p_j = positions[j]
                o_j = orientations[j]

                distance_segment_2, dP = distanceSegments(p_i, o_i, p_j, o_j)

                if distance_segment_2 < min_distance_2:
                    rejected = True
                    break

            if not rejected:
                positions.append(p_i)
                orientations.append(o_i)
                inserted = True

        if not inserted:
            radius = radius + min_distance
            print ".. increasing radius", radius
            # side = side + min_distance**(1.0/3.0)

print "after:"
stats(N, positions, orientations)

export = open('XcT_gen'+str(N)+'.in', 'w')
export.write(str(N) + '\n')

for i in xrange(N):
    export.write('\t'.join([str(x) for x in [positions[i].x,positions[i].y,positions[i].z]]) + '\n')
    export.write('\t'.join([str(x) for x in [orientations[i].x,orientations[i].y,orientations[i].z]]) + '\n')

export.close()
