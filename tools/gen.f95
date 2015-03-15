#define PI (3.1415926535897932D0)

PROGRAM gen
  USE omp_lib
  IMPLICIT NONE

  CHARACTER(LEN=256)::program_name, str_number_of_fibers, str_min_distance, str_avg_distance
  CHARACTER(LEN=6)::domain_type

  INTEGER::N
  REAL*4::min_distance, avg_distance, step, rotation_step

  REAL*4,SAVE,DIMENSION(:,:),ALLOCATABLE::positions
  REAL*4,SAVE,DIMENSION(:,:),ALLOCATABLE::orientations

  LOGICAL::optimal,boxBoundary,sphereBoundary
  INTEGER::i,j,optimization_iteration
  REAL*4,DIMENSION(3)::p_i,p_j,o_i,o_j
  REAL*4::r
  REAL*4::distance, nearest_distance
  REAL*4,DIMENSION(3)::nearest_dP
  INTEGER::nearest_index
  INTEGER::count_pairs,old_count_pairs

  REAL*4::distanceSegments
  REAL*4,DIMENSION(3)::dP

  REAL*4::domain

  CALL init_random_seed()

  CALL GETARG(0, program_name)
  CALL GETARG(1, domain_type)
  CALL GETARG(2, str_number_of_fibers)
  CALL GETARG(3, str_min_distance)
  CALL GETARG(4, str_avg_distance)
  program_name = trim(adjustl(program_name))
  domain_type = trim(adjustl(domain_type))
  READ(str_number_of_fibers,'(I10)') N
  READ(str_min_distance,'(F16.8)') min_distance
  READ(str_avg_distance,'(F16.8)') avg_distance

  ! PRINT *, domain_type, N, min_distance, avg_distance

  ! Values for performance charts in paper
  ! min_distance = 0.2
  ! avg_distance = min_distance * 5

  ALLOCATE(positions(N,3))
  ALLOCATE(orientations(N,3))
  step = 1.0 * min_distance
  rotation_step = 1.0

  !! to ensure enough space to account for fiber orientations
  !! increase avg_distance for domain calculation by 2
  !! this factor empirically result in the expected average distance
  !avg_distance = avg_distance * 2

  IF (domain_type .EQ. "box") THEN
    CALL boxDomain(N, avg_distance, positions, orientations, domain)
  ELSE IF (domain_type .EQ. "sphere") THEN
    CALL sphereDomain(N, avg_distance, positions, orientations, domain)
  ELSE
    PRINT *, "[ERROR] Unknown domain type: ", domain_type
    CALL EXIT(1)
  END IF

  PRINT *, "---  before  ---"
  CALL stats(N, positions, orientations, min_distance,count_pairs)
  PRINT *, "----------------"

  !distance = distanceSegments(positions(1,:), orientations(1,:), positions(2,:), orientations(2,:))
  ! PRINT *, distance, SQRT(SUM((positions(1,:) - positions(2,:))**2))

  optimal = .FALSE.

  PRINT *, "... optimizing distribution"
  optimization_iteration = 0
  DO WHILE (optimal .NEQV. .TRUE.)
    optimal = .TRUE.
    optimization_iteration = optimization_iteration + 1

    IF (MODULO(optimization_iteration, 100) == 0) THEN
      PRINT *, "--- optimize ---"
      CALL stats(N, positions, orientations, min_distance, count_pairs)
      IF (old_count_pairs > count_pairs) THEN
        step = 0.98 * step
        rotation_step = 0.98 * rotation_step
      ELSE
        step = 1.01 * step
        rotation_step = 1.01 * rotation_step
      END IF
      old_count_pairs = count_pairs
      PRINT *, "steps: ", step
      PRINT *, "----------------"
    END IF

    DO i = 1, N

      p_i = positions(i,:)
      o_i = orientations(i,:)

      ! force fibers inside domain boundary
      IF (domain_type .EQ. "box") THEN
        optimal = optimal .AND. boxBoundary(p_i, domain, step)
      ELSE IF (domain_type .EQ. "sphere") THEN
        optimal = optimal .AND. sphereBoundary(p_i, domain, step)
      END IF

      nearest_distance = 99999999.0
      DO j = 1, N

        IF (i /= j) THEN

          p_j = positions(j,:)
          o_j = orientations(j,:)

          distance = distanceSegments(p_i, o_i, p_j, o_j, dP)
          distance = SQRT(distance)

          IF (distance < nearest_distance) THEN
            nearest_distance = distance
            nearest_dP = dP
            nearest_index = j
          END IF

        END IF

      END DO

      IF (nearest_distance < min_distance) THEN

        p_j = positions(nearest_index,:)
        o_j = orientations(nearest_index,:)

        ! randomly move fiber
        CALL RANDOM_NUMBER(r)
        p_i(1) = p_i(1) + r * (step - (-step)) + (-step)
        CALL RANDOM_NUMBER(r)
        p_i(2) = p_i(2) + r * (step - (-step)) + (-step)
        CALL RANDOM_NUMBER(r)
        p_i(3) = p_i(3) + r * (step - (-step)) + (-step)

        CALL RANDOM_NUMBER(r)
        p_j(1) = p_j(1) + r * (step - (-step)) + (-step)
        CALL RANDOM_NUMBER(r)
        p_j(2) = p_j(2) + r * (step - (-step)) + (-step)
        CALL RANDOM_NUMBER(r)
        p_j(3) = p_j(3) + r * (step - (-step)) + (-step)

        ! ! exactly seperate fibers
        ! ! normalize difference vector
        ! nearest_dP = nearest_dP / SQRT(DOT_PRODUCT(nearest_dp,nearest_dp))
        ! ! scale difference vector by error
        ! ! add small additional distance for numeric precision
        ! nearest_dP = nearest_dP * (min_distance - nearest_distance) * (1+1e-2)
        ! ! move fibers
        ! p_i = p_i + nearest_dP / 2.0
        ! p_j = p_j - nearest_dP / 2.0

        CALL RANDOM_NUMBER(r)
        o_i(1) = o_i(1) + r * (rotation_step - (-rotation_step)) + (-rotation_step)
        CALL RANDOM_NUMBER(r)
        o_i(2) = o_i(2) + r * (rotation_step - (-rotation_step)) + (-rotation_step)
        CALL RANDOM_NUMBER(r)
        o_i(3) = o_i(3) + r * (rotation_step - (-rotation_step)) + (-rotation_step)

        o_i = o_i / SQRT(DOT_PRODUCT(o_i, o_i))

        CALL RANDOM_NUMBER(r)
        o_j(1) = o_j(1) + r * (rotation_step - (-rotation_step)) + (-rotation_step)
        CALL RANDOM_NUMBER(r)
        o_j(2) = o_j(2) + r * (rotation_step - (-rotation_step)) + (-rotation_step)
        CALL RANDOM_NUMBER(r)
        o_j(3) = o_j(3) + r * (rotation_step - (-rotation_step)) + (-rotation_step)

        o_j = o_j / SQRT(DOT_PRODUCT(o_j, o_j))

        positions(j,:) = p_j
        orientations(j,:) = o_j

        optimal = optimal .AND. .FALSE.
      END IF

      positions(i,:) = p_i
      orientations(i,:) = o_i

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

FUNCTION boxBoundary(p, domain, step) RESULT(inside)
  IMPLICIT NONE

  REAL*4,DIMENSION(3),INTENT(INOUT)::p
  REAL*4,INTENT(IN)::domain
  REAL*4,INTENT(IN)::step
  REAL*4::r
  LOGICAL::inside

  inside = .TRUE.
  IF (p(1) > domain) THEN
    CALL RANDOM_NUMBER(r)
    p(1) = p(1) - r * step
    inside = .FALSE.
  ELSE IF (p(1) < -domain) THEN
    CALL RANDOM_NUMBER(r)
    p(1) = p(1) + r * step
    inside = .FALSE.
  END IF

  IF (p(2) > domain) THEN
    CALL RANDOM_NUMBER(r)
    p(2) = p(2) - r * step
    inside = .FALSE.
  ELSE IF (p(2) < -domain) THEN
    CALL RANDOM_NUMBER(r)
    p(2) = p(2) + r * step
    inside = .FALSE.
  END IF

  IF (p(3) > domain) THEN
    CALL RANDOM_NUMBER(r)
    p(3) = p(3) - r * step
    inside = .FALSE.
  ELSE IF (p(3) < -domain) THEN
    CALL RANDOM_NUMBER(r)
    p(3) = p(3) + r * step
    inside = .FALSE.
  END IF

END FUNCTION

FUNCTION sphereBoundary(p, radius, step) RESULT(inside)
  IMPLICIT NONE

  REAL*4,DIMENSION(3),INTENT(INOUT)::p
  REAL*4::p_length
  REAL*4,INTENT(IN)::radius
  REAL*4,INTENT(IN)::step
  REAL*4::r
  LOGICAL::inside

  inside = .TRUE.
  p_length = SQRT(DOT_PRODUCT(p,p))
  IF (p_length > radius) THEN
    CALL RANDOM_NUMBER(r)
    !normalize distance to 0,0,0
    p = p / p_length
    !shrink length randomly
    p_length = p_length - r * step
    !rescale vector to reduced length
    p = p * p_length
    inside = .FALSE.
  END IF

END FUNCTION

SUBROUTINE boxDomain(N, average_distance, positions, orientations, domain)
  IMPLICIT NONE

  INTEGER,INTENT(IN)::N
  REAL*4,INTENT(IN)::average_distance
  REAL*4,INTENT(INOUT),DIMENSION(N,3)::positions
  REAL*4,INTENT(INOUT),DIMENSION(N,3)::orientations
  REAL*4,INTENT(OUT)::domain
  REAL*4::side_length
  INTEGER::i

  side_length = (N+1)**(1.0/3.0) * average_distance
  domain = side_length;

  PRINT *, "Generating box configuration with ", N, " fibers..."
  PRINT *, "... using bounding box with domain: [",-domain,", ",domain,"]"

  CALL RANDOM_NUMBER(positions)
  CALL RANDOM_NUMBER(orientations)

  positions = positions * (domain - (-domain)) + (-domain)
  orientations = orientations * (1 - (-1)) + (-1)

  ! normalize orientations
  DO i=1,N
    orientations(i,:) = orientations(i,:) / SQRT(DOT_PRODUCT(orientations(i,:),orientations(i,:)))
  END DO

END SUBROUTINE

SUBROUTINE sphereDomain(N, average_distance, positions, orientations, domain)
  IMPLICIT NONE

  INTEGER,INTENT(IN)::N
  REAL*4,INTENT(IN)::average_distance
  REAL*4,INTENT(INOUT),DIMENSION(N,3)::positions
  REAL*4,INTENT(INOUT),DIMENSION(N,3)::orientations
  REAL*4,INTENT(OUT)::domain
  REAL*4::radius, radius_2
  REAL*4::side_length, r
  REAL*4,DIMENSION(3)::p,o
  INTEGER::i

  side_length = (N+1)**(1.0/3.0) * average_distance
  ! to make room for oriented fibers increase box by factor 2
  side_length = 2 * side_length
  ! ensure same volume for both box and sphere
  radius = ((3 * side_length**3.0) / (4 * PI))**(1.0/3.0)
  radius_2 = radius * radius
  domain = radius

  PRINT *, "Generating sphere configuration with ", N, " fibers..."
  PRINT *, "... using bounding sphere with radius: ", radius

  CALL RANDOM_NUMBER(orientations)
  orientations = orientations * (1 - (-1)) + (-1)

  DO i=1,N
    DO WHILE (.TRUE.)
      CALL RANDOM_NUMBER(r)
      p(1) = r * (domain - (-domain)) + (-domain)
      CALL RANDOM_NUMBER(r)
      p(2) = r * (radius - (-domain)) + (-domain)
      CALL RANDOM_NUMBER(r)
      p(3) = r * (radius - (-domain)) + (-domain)

      IF (DOT_PRODUCT(p,p) < radius_2) EXIT
    END DO

    positions(i,:) = p

    ! normalize orientations
    orientations(i,:) = orientations(i,:) / SQRT(DOT_PRODUCT(orientations(i,:),orientations(i,:)))
  END DO

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

  ! S1P0(1) = 1.55296648;
  ! S1P0(2) = -0.639337718;
  ! S1P0(3) = -1.15629041;
  !
  ! S1P1(1) = 0.815716386;
  ! S1P1(2) = -1.94480538;
  ! S1P1(3) = 0.167422503;
  !
  ! S2P0(1) = -1.82547736;
  ! S2P0(2) = 3.73734117;
  ! S2P0(3) = 1.79775524;
  !
  ! S2P1(1) = -1.06064153;
  ! S2P1(2) = 3.02075219;
  ! S2P1(3) = 3.50114202;

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
    count = count + 1

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
  PRINT *, "  - average distance segments: ", total_segment/count
  PRINT *, "  - average distance centers: ", total_center/count

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
