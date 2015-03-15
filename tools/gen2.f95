#define PI (3.1415926535897932D0)

PROGRAM gen
  USE omp_lib
  IMPLICIT NONE

  CHARACTER(LEN=256)::program_name, str_number_of_fibers, str_min_distance, str_avg_distance
  CHARACTER(LEN=6)::domain_type

  INTEGER::N
  REAL*4::min_distance, min_distance_2, distance_segments

  REAL*4,SAVE,DIMENSION(:,:),ALLOCATABLE::positions
  REAL*4,SAVE,DIMENSION(:,:),ALLOCATABLE::orientations

  INTEGER::i,j
  INTEGER::count_pairs
  REAL*4::radius

  REAL*4,DIMENSION(3)::p_i,p_j,o_i,o_j
  LOGICAL::inserted,rejected
  INTEGER::retries

  REAL*4::distanceSegments
  REAL*4,DIMENSION(3)::dP

  CALL init_random_seed()

  CALL GETARG(0, program_name)
  CALL GETARG(1, domain_type)
  CALL GETARG(2, str_number_of_fibers)
  CALL GETARG(3, str_min_distance)
  program_name = trim(adjustl(program_name))
  domain_type = trim(adjustl(domain_type))
  READ(str_number_of_fibers,'(I10)') N
  READ(str_min_distance,'(F16.8)') min_distance

  ALLOCATE(positions(N,3))
  ALLOCATE(orientations(N,3))

  radius = min_distance
  min_distance_2 = min_distance * min_distance
  ! IF (domain_type .EQ. "box") THEN
  !   CALL boxDomain(N, avg_distance, positions, orientations, domain)
  ! ELSE IF (domain_type .EQ. "sphere") THEN
  !   CALL sphereDomain(N, avg_distance, positions, orientations, domain)
  ! ELSE
  !   PRINT *, "[ERROR] Unknown domain type: ", domain_type
  !   CALL EXIT(1)
  ! END IF

  DO i = 1, N
    inserted = .FALSE.

    DO WHILE (inserted .NEQV. .TRUE.)
      retries = 0
      rejected = .TRUE.

      DO WHILE (rejected .EQV. .TRUE. .AND. retries < 1000)
        retries = retries + 1
        rejected = .FALSE.

        CALL sampleSphere(max(0.0, radius - 4 * min_distance), radius, p_i, o_i)

        DO j = 1, i-1
          p_j = positions(j,:)
          o_j = orientations(j,:)

          distance_segments = distanceSegments(p_i, o_i, p_j, o_j, dP)

          IF (distance_segments < min_distance_2) THEN
            rejected = .TRUE.
            EXIT
          END IF

        END DO

        IF (rejected .NEQV. .TRUE.) THEN
          positions(i,:) = p_i
          orientations(i,:) = o_i
          inserted = .TRUE.
        END IF

      END DO

      IF (inserted .NEQV. .TRUE.) THEN
        radius = radius + min_distance
      END IF
    END DO
  END DO

  PRINT *, "---  after   ---"
  CALL stats(N, positions, orientations,min_distance,count_pairs)
  PRINT *, "----------------"

  OPEN(10,file="XcT_gen"//TRIM(str_number_of_fibers)//".in")
  WRITE(10,'(*(I10))') (N)
  DO i=1,N
    WRITE(10,'(*(F32.16))') (positions(i,:))
    WRITE(10,'(*(F32.16))') (orientations(i,:) / SQRT(DOT_PRODUCT(orientations(i,:),orientations(i,:))))
  END DO
  CLOSE(10)

END PROGRAM gen

SUBROUTINE sampleSphere(min_radius, max_radius, p, o)
  IMPLICIT NONE

  REAL*4,DIMENSION(3),INTENT(OUT)::p,o
  REAL*4,INTENT(IN)::min_radius,max_radius
  REAL*4::min_radius_2,max_radius_2
  REAL*4::p_norm_2
  LOGICAL::hit

  min_radius_2 = min_radius * min_radius
  max_radius_2 = max_radius * max_radius

  hit = .FALSE.

  DO WHILE (hit .NEQV. .TRUE.)
    CALL RANDOM_NUMBER(p)
    p = p * (max_radius - (-max_radius)) + (-max_radius)

    p_norm_2 = DOT_PRODUCT(p,p)
    IF ((p_norm_2 > min_radius_2) .AND. (p_norm_2 < max_radius_2)) THEN
      hit = .TRUE.
    END IF
  END DO

  CALL RANDOM_NUMBER(o)
  o = o * (1 - (-1)) + (-1)
  o = o / SQRT(DOT_PRODUCT(o,o))

END SUBROUTINE

!! http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
FUNCTION distanceSegments(p_i, o_i, p_j, o_j, dP) RESULT(distance)
  IMPLICIT NONE

  REAL*4,DIMENSION(3),INTENT(IN)::p_i,p_j,o_i,o_j
  REAL*4::distance

  REAL*4,DIMENSION(3)::S1P1,S1P0,S2P1,S2P0
  REAL*4,DIMENSION(3)::u,v,w
  REAL*4,DIMENSION(3),INTENT(OUT)::dP
  REAL*4::a,b,c,d,e
  REAL*4::X
  REAL*4::sc,sN,sD
  REAL*4::tc,tN,tD

  S1P0 = p_i - o_i
  S1P1 = p_i + o_i
  S2P0 = p_j - o_j
  S2P1 = p_j + o_j

  u = S1P1 - S1P0
  v = S2P1 - S2P0
  w = S1P0 - S2P0

  a = DOT_PRODUCT(u,u)
  b = DOT_PRODUCT(u,v)
  c = DOT_PRODUCT(v,v)
  d = DOT_PRODUCT(u,w)
  e = DOT_PRODUCT(v,w)

  X = a * c - b * b
  sc = X
  sN = X
  sD = X
  tc = X
  tN = X
  tD = X

  IF ( X < 1e-5 ) THEN
    sN = 0.0
    sD = 1.0
    tN = e
    tD = c
  ELSE
    sN = (b * e - c * d)
    tN = (a * e - b * d)
    IF (sN < 0.0) THEN
      sN = 0.0
      tN = e
      tD = c
    ELSE IF (sN > sD) THEN
      sN = sD
      tN = e + b
      tD = c
    END IF
  END IF

  IF ( tN < 0.0 ) THEN
    tN = 0.0
    IF (-d < 0.0) THEN
      sN = 0.0
    ELSE IF (-d > a) THEN
      sN = sD
    ELSE
      sN = -d
      sD = a
    END IF
  ELSE IF ( tN > tD) THEN
    tN = tD
    IF ((-d + b) < 0.0) THEN
      sN = 0.0
    ELSE IF ((-d + b) > a) THEN
      sN = sD
    ELSE
      sN = (-d + b)
      sD = a
    END IF
  END IF

  IF (ABS(sN) < 1e-5) THEN
    sc = 0.0
  ELSE
    sc = sN / sD
  END IF
  IF (ABS(tN) < 1e-5) THEN
    tc = 0.0
  ELSE
    tc = tN / tD
  END IF

  dP = w + (sc * u) - (tc * v)

  distance = DOT_PRODUCT(dp,dp)

END FUNCTION distanceSegments

SUBROUTINE stats(N, positions, orientations, min_distance, count_pairs)
  IMPLICIT NONE

  INTEGER,INTENT(IN)::N
  REAL*4,INTENT(IN),DIMENSION(N,3)::positions
  REAL*4,INTENT(IN),DIMENSION(N,3)::orientations
  REAL*4,INTENT(IN)::min_distance

  INTEGER::i,j,count

  INTEGER,INTENT(OUT)::count_pairs
  REAL*4::distance_segment, distance_center
  REAL*4::nearest_distance_segment, nearest_distance_center
  REAL*4::total_segment, total_center
  REAL*4::minimal_distance_segment, minimal_distance_center
  REAL*4,DIMENSION(3)::p_i,p_j,o_i,o_j
  REAL*4::distanceSegments
  REAL*4,DIMENSION(3)::dP
  REAL*4,DIMENSION(3)::min_domain, max_domain

  min_domain = 99999999.0
  max_domain = -99999999.0

  count = 0
  total_segment = 0.0
  total_center = 0.0
  count_pairs = 0
  minimal_distance_segment = 99999999.0
  minimal_distance_center = 99999999.0

  DO i = 1, N

    p_i = positions(i,:)
    o_i = orientations(i,:)

    IF (p_i(1) < min_domain(1)) THEN
      min_domain(1) = p_i(1)
    END IF
    IF (p_i(2) < min_domain(2)) THEN
      min_domain(2) = p_i(2)
    END IF
    IF (p_i(3) < min_domain(3)) THEN
      min_domain(3) = p_i(3)
    END IF
    IF (p_i(1) > max_domain(1)) THEN
      max_domain(1) = p_i(1)
    END IF
    IF (p_i(2) > max_domain(2)) THEN
      max_domain(2) = p_i(2)
    END IF
    IF (p_i(3) > max_domain(3)) THEN
      max_domain(3) = p_i(3)
    END IF

    nearest_distance_segment = 99999999.0
    nearest_distance_center = 99999999.0
    DO j= 1, N

      IF (i /= j) THEN

        p_j = positions(j,:)
        o_j = orientations(j,:)

        distance_segment = SQRT(distanceSegments(p_i, o_i, p_j, o_j,dP))
        distance_center = SQRT(SUM((p_i - p_j)**2))

        IF (distance_segment < nearest_distance_segment) THEN
          nearest_distance_segment = distance_segment
        END IF
        IF (distance_center < nearest_distance_center) THEN
          nearest_distance_center = distance_center
        END IF
        IF (distance_segment < min_distance) THEN
          count_pairs = count_pairs + 1
        END IF
      END IF

    END DO

    total_segment = total_segment + nearest_distance_segment
    total_center = total_center + nearest_distance_center

    IF (nearest_distance_segment < minimal_distance_segment) THEN
      minimal_distance_segment = nearest_distance_segment
    END IF
    IF (nearest_distance_center < minimal_distance_center) THEN
      minimal_distance_center = nearest_distance_center
    END IF

  END DO

  PRINT *, "Distribution Statistics:"
  PRINT *, "  - number of too close pairs: ", count_pairs
  PRINT *, "  - domain: [", min_domain, ",", max_domain, "]"
  PRINT *, "  - minimal distance segments: ", minimal_distance_segment
  PRINT *, "  - minimal distance centers:  ", minimal_distance_center
  PRINT *, "  - average distance segments: ", total_segment/N
  PRINT *, "  - average distance centers: ", total_center/N

END SUBROUTINE stats

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
