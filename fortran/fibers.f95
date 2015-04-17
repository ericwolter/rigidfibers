#include "constants.incl"

#define PLUS_EQUALS(A, B) A = A + B
#define MINUS_EQUALS(A, B) A = A - B
#define MULTIPLY_EQUALS(A, B) A = A * B
#define DIVIDE_EQUALS(A, B) A = A / B

!#define CONDITION_NUMBER

PROGRAM fibers
  USE omp_lib
  IMPLICIT NONE

  CHARACTER(LEN=256)::program_name, data_name

#if defined(BENCHMARK)
  INTEGER total_count_rate,total_count_max,total_count1,total_count2
  REAL*8 total_CPU_p
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p
#endif

  !--------------------------------------------------
  ! Initalize Memory
  !--------------------------------------------------
  REAL*4,ALLOCATABLE,TARGET,DIMENSION(:)::t_previous_positions, t_current_positions, t_next_positions
  REAL*4,ALLOCATABLE,TARGET,DIMENSION(:)::t_previous_orientations, t_current_orientations, t_next_orientations
  REAL*4,ALLOCATABLE,TARGET,DIMENSION(:)::t_previous_translational_velocities, t_current_translational_velocities
  REAL*4,ALLOCATABLE,TARGET,DIMENSION(:)::t_previous_rotational_velocities, t_current_rotational_velocities

  REAL*4,POINTER::tmp_pointer(:)
  REAL*4,POINTER::previous_positions(:),current_positions(:),next_positions(:)
  REAL*4,POINTER::previous_orientations(:),current_orientations(:),next_orientations(:)
  REAL*4,POINTER::previous_translational_velocities(:),current_translational_velocities(:)
  REAL*4,POINTER::previous_rotational_velocities(:),current_rotational_velocities(:)

  REAL*4,ALLOCATABLE,DIMENSION(:,:)::a_matrix
  REAL*4,ALLOCATABLE,DIMENSION(:)::b_vector

  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_points
  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_weights
  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS,NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::legendre_polynomials

  REAL*4,SAVE,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::lambda
  REAL*4,SAVE,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::eigen

  REAL*4,SAVE,DIMENSION(DIMENSIONS)::external_force

  INTEGER::i,j,force_index,force_index_i, force_index_j,quadrature_index_i,quadrature_index_j
  REAL*4,DIMENSION(DIMENSIONS)::position_i, orientation_i
  REAL*4,DIMENSION(DIMENSIONS)::position_j, orientation_j
  REAL*4::quadrature_weight, legendre_polynomial

  REAL*4,DIMENSION(6)::T
  REAL*4,DIMENSION(3)::Q
  REAL*4,DIMENSION(3)::TF, TFA0, TFA1, TFA1_TMP
  REAL*4::QF

  REAL*4,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS, 6)::G
  REAL*4,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS, 3)::GF
  REAL*4,DIMENSION(DIMENSIONS)::position_on_fiber_i,position_on_fiber_j,difference,difference2, oriented_force
  REAL*4::invDistance,invDistance3,invDistance5
  REAL*4,DIMENSION(6)::K
  REAL*4,DIMENSION(DIMENSIONS)::force_on_fiber_j

#if defined(ANALYTICAL)
  INTEGER::f, n

  REAL*4,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3)::I1,I3,I5

  REAL*4,DIMENSION(DIMENSIONS)::R0

  REAL*4::b,l,invL,llimit,m,s_upper,s_lower,u_upper,u_lower,pow_s_upper,pow_s_lower
  REAL*4::i1n2,i1n1,i1n0,i3n2,i3n1,i3n0,i5n2,i5n1,i5n0
  REAL*4::L01,L03,L05,L13,L15,L23,L25
  REAL*4::G11,G22,G33,G12,G13,G23
#endif

#if defined(GMRES)
  INTEGER,DIMENSION(8)::icntl
  REAL*4,DIMENSION(5)::cntl

  REAL*4,DIMENSION(GMRES_LWORK)::work

  INTEGER,DIMENSION(5)::irc

  INTEGER::gmres_done
  INTEGER,DIMENSION(3)::gmres_info
  REAL*4,DIMENSION(2)::gmres_rinfo
#endif

#if defined(CONDITION_NUMBER)
  REAL*4,DIMENSION(4*TOTAL_NUMBER_OF_ROWS)::cond_work
  INTEGER,DIMENSION(4*TOTAL_NUMBER_OF_ROWS)::cond_iwork
  REAL*4,EXTERNAL::slange
  REAL*4::a_norm, rcond
  INTEGER::cond_info
  REAL*4,DIMENSION(TOTAL_NUMBER_OF_ROWS)::cond_ipvt
#endif

  INTEGER::x_row_index,y_row_index,z_row_index
  INTEGER::x_col_index,y_col_index,z_col_index
  REAL*4::c,d,e,cc,D1,gamma

  INTEGER,SAVE,DIMENSION(TOTAL_NUMBER_OF_ROWS)::IPIV
  INTEGER::INFO

  INTEGER::IDUMMY,ind

  INTEGER::current_timestep
  CHARACTER(len=8)::str_timestep

  ALLOCATE( t_previous_positions(DIMENSIONS*NUMBER_OF_FIBERS))
  ALLOCATE( t_current_positions(DIMENSIONS*NUMBER_OF_FIBERS))
  ALLOCATE( t_next_positions(DIMENSIONS*NUMBER_OF_FIBERS))

  ALLOCATE( t_previous_orientations(DIMENSIONS*NUMBER_OF_FIBERS))
  ALLOCATE( t_current_orientations(DIMENSIONS*NUMBER_OF_FIBERS))
  ALLOCATE( t_next_orientations(DIMENSIONS*NUMBER_OF_FIBERS))

  ALLOCATE( t_previous_translational_velocities(DIMENSIONS*NUMBER_OF_FIBERS))
  ALLOCATE( t_current_translational_velocities(DIMENSIONS*NUMBER_OF_FIBERS))

  ALLOCATE( t_previous_rotational_velocities(DIMENSIONS*NUMBER_OF_FIBERS))
  ALLOCATE( t_current_rotational_velocities(DIMENSIONS*NUMBER_OF_FIBERS))

  ALLOCATE( a_matrix(TOTAL_NUMBER_OF_ROWS,TOTAL_NUMBER_OF_ROWS))
  ALLOCATE( b_vector(TOTAL_NUMBER_OF_ROWS))

  previous_positions => t_previous_positions
  current_positions => t_current_positions
  next_positions => t_next_positions
  previous_orientations => t_previous_orientations
  current_orientations => t_current_orientations
  next_orientations => t_next_orientations
  previous_translational_velocities => t_previous_translational_velocities
  current_translational_velocities => t_current_translational_velocities
  previous_rotational_velocities => t_previous_rotational_velocities
  current_rotational_velocities => t_current_rotational_velocities

  c = LOG(SLENDERNESS * SLENDERNESS * EXP(1.0))
  d = -c
  e = 2.0
  cc = 1.0
  D1 = 0.75 / (d - 2.0 * cc)

  a_matrix = 0.0
  b_vector = 0.0
  external_force = (/ 0.0, 0.0, -0.5 /)

  !--------------------------------------------------
  ! Load positions and orientations
  !--------------------------------------------------
  !me:  Intializes the positons and @todo tVecs with data from the specified
  !     input file.
  CALL GETARG(0, program_name)
	CALL GETARG(1, data_name)

  PRINT *, TRIM(data_name)
  PRINT *, NUMBER_OF_FIBERS

  OPEN(10,file=TRIM(data_name));
  READ(10,*) IDUMMY
  DO i=1,NUMBER_OF_FIBERS
     ind=(i-1)*3
     READ(10,*) current_positions(ind+1),current_positions(ind+2),current_positions(ind+3),current_orientations(ind+1),current_orientations(ind+2),current_orientations(ind+3)
  END DO
  CLOSE(10)
  PRINT *,"Read initial data from file. "
  ! PRINT '(*(F16.8))', current_positions

  !--------------------------------------------------
  ! Precompute constants
  !--------------------------------------------------
  CALL precomputeLegendrePolynomials(quadrature_points, quadrature_weights, legendre_polynomials)
  CALL precomputeLambda(lambda, eigen)

  !--------------------------------------------------
  ! Simulation Step
  !--------------------------------------------------

  DO current_timestep = 0, NUMBER_OF_TIMESTEPS-1
#if defined(VALIDATE)
    WRITE(str_timestep, '(I1.1)') current_timestep
#endif
#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(total_count1, total_count_rate, total_count_max)
#endif

    !--------------------------------------------------
    ! 1. Assemble System
    !--------------------------------------------------
#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)
#endif

#include "assemble_matrix.incl"

#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:assemble_system:", CPU_p
#endif

#if defined(VALIDATE)
    OPEN(10,file=""//TRIM(str_timestep)//"_AMat.out");
    DO i=1,TOTAL_NUMBER_OF_ROWS
      WRITE(10,'(*(F16.8))') (a_matrix(i,j),j=1,TOTAL_NUMBER_OF_ROWS)
    END DO
    CLOSE(10)

    OPEN(10,file=""//TRIM(str_timestep)//"_BVec.out");
    DO i=1,TOTAL_NUMBER_OF_ROWS
      WRITE(10,'(*(F16.8))') (b_vector(i))
    END DO
    CLOSE(10)
#endif

    !--------------------------------------------------
    ! 2. Solve System
    !--------------------------------------------------
#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)
#endif

#include "solve_system.incl"

#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:solve_system:", CPU_p
#endif

#if defined(VALIDATE)
    OPEN(10,file=""//TRIM(str_timestep)//"_XVec.out");
    DO i=1,TOTAL_NUMBER_OF_ROWS
      WRITE(10,'(*(F16.8))') (b_vector(i))
    END DO
    CLOSE(10)
#endif

    !--------------------------------------------------
    ! 3. Update System
    !--------------------------------------------------
    !--------------------------------------------------
    ! 3.1. Update Velocity
    !--------------------------------------------------
#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)
#endif

#include "update_velocities.incl"

#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:update_velocities:", CPU_p
#endif

#if defined(VALIDATE)
    OPEN(10,file=""//TRIM(str_timestep)//"_TRANSVel.out");
    DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
     WRITE(10,'(*(F16.8))') (current_translational_velocities(i))
    END DO
    CLOSE(10)
    OPEN(10,file=""//TRIM(str_timestep)//"_ROTVel.out");
    DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
     WRITE(10,'(*(F16.8))') (current_rotational_velocities(i))
    END DO
    CLOSE(10)
#endif

    !--------------------------------------------------
    ! 3.2. Update Fibers
    !--------------------------------------------------
#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)
#endif

    IF (current_timestep == 0) THEN
#include "update_fibers_firststep.incl"
    ELSE
#include "update_fibers.incl"
    END IF

#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:update_fibers:", CPU_p
#endif

#if defined(VALIDATE)
    OPEN(10,file=""//TRIM(str_timestep)//"_POS.out");
    DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
      WRITE(10,'(*(F16.8))') (next_positions(i))
    END DO
    CLOSE(10)
    OPEN(10,file=""//TRIM(str_timestep)//"_ORIENT.out");
    DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
      WRITE(10,'(*(F16.8))') (next_orientations(i))
    END DO
    CLOSE(10)
#endif

#if !defined(BENCHMARK) && !defined(VALIDATE)
    IF (STATE_SAVE_INTERVAL > 0 .AND. mod(current_timestep,STATE_SAVE_INTERVAL) == 0) THEN
      WRITE(str_timestep, '(I0.5)') current_timestep
      OPEN(10,file=""//TRIM(str_timestep)//".state");
      WRITE(10,*) NUMBER_OF_FIBERS
      DO i=1,NUMBER_OF_FIBERS
       WRITE(10,'(*(F16.8))') current_positions((i-1)*DIMENSIONS+1),current_positions((i-1)*DIMENSIONS+2),current_positions((i-1)*DIMENSIONS+3)
       WRITE(10,'(*(F16.8))') current_orientations((i-1)*DIMENSIONS+1),current_orientations((i-1)*DIMENSIONS+2),current_orientations((i-1)*DIMENSIONS+3)
      END DO
      CLOSE(10)
    END IF
    IF (VELOCITY_SAVE_INTERVAL > 0 .AND. mod(current_timestep,VELOCITY_SAVE_INTERVAL) == 0) THEN
      OPEN(20,file=""//TRIM(str_timestep)//".velocity")
      DO i=1,NUMBER_OF_FIBERS
       WRITE(20,'(*(F16.8))') current_translational_velocities((i-1)*DIMENSIONS+1),current_translational_velocities((i-1)*DIMENSIONS+2),current_translational_velocities((i-1)*DIMENSIONS+3)
      END DO
      CLOSE(20)
    END IF
#endif

    tmp_pointer => previous_translational_velocities
    previous_translational_velocities => current_translational_velocities
    current_translational_velocities => tmp_pointer

    tmp_pointer => previous_rotational_velocities
    previous_rotational_velocities => current_rotational_velocities
    current_rotational_velocities => tmp_pointer

    tmp_pointer => previous_positions
    previous_positions => current_positions
    current_positions => next_positions
    next_positions => tmp_pointer

    tmp_pointer => previous_orientations
    previous_orientations => current_orientations
    current_orientations => next_orientations
    next_orientations => tmp_pointer

#if defined(BENCHMARK)
    CALL SYSTEM_CLOCK(total_count2, total_count_rate, total_count_max)
    total_CPU_p = real(total_count2-total_count1)/total_count_rate
    PRINT *,"BENCHMARK:total:", total_CPU_p
#endif

  END DO

END PROGRAM fibers

SUBROUTINE precomputeLegendrePolynomials(quadrature_points, quadrature_weights, legendre_polynomials)

  REAL*8::p0,p1,p2
  REAL*8::w0,w1,w2

  REAL*8::lower_bound
  REAL*8::interval_size

  REAL*4,INTENT(OUT),DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_points
  REAL*4,INTENT(OUT),DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_weights
  REAL*4,INTENT(OUT),DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS,NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::legendre_polynomials
  REAL*8,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::internal_quadrature_points
  REAL*8,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::internal_quadrature_weights
  REAL*8,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS,NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::internal_legendre_polynomials

  INTEGER::interval_index, force_index, point_index

  p0 = -SQRT(15.0d0) / 5.0d0
  p1 = 0.0d0
  p2 = SQRT(15.0d0) / 5.0d0

  w0 = 5.0d0 / 9.0d0
  w1 = 8.0d0 / 9.0d0
  w2 = 5.0d0 / 9.0d0

  lower_bound = -1.0d0

  interval_size = 2.0d0 / NUMBER_OF_QUADRATURE_INTERVALS

  DO interval_index = 0, NUMBER_OF_QUADRATURE_INTERVALS-1
    internal_quadrature_points(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 1) = &
      (2.0d0 * lower_bound + interval_size + p0 * interval_size) / 2.0d0
    internal_quadrature_points(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 2) = &
      (2.0d0 * lower_bound + interval_size + p1 * interval_size) / 2.0d0
    internal_quadrature_points(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 3) = &
      (2.0d0 * lower_bound + interval_size + p2 * interval_size) / 2.0d0

    internal_quadrature_weights(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 1) = &
      w0 / NUMBER_OF_QUADRATURE_INTERVALS
    internal_quadrature_weights(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 2) = &
      w1 / NUMBER_OF_QUADRATURE_INTERVALS
    internal_quadrature_weights(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 3) = &
      w2 / NUMBER_OF_QUADRATURE_INTERVALS

    lower_bound = lower_bound + interval_size
  END DO

  quadrature_points = internal_quadrature_points
  quadrature_weights = internal_quadrature_weights

  DO force_index = 0, NUMBER_OF_TERMS_IN_FORCE_EXPANSION-1
    DO point_index = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
      internal_legendre_polynomials(point_index + 1,force_index + 1) = &
        calculateLegendrePolynomial(internal_quadrature_points(point_index+1), force_index+1)
    END DO
  END DO

  legendre_polynomials = internal_legendre_polynomials

END SUBROUTINE precomputeLegendrePolynomials

SUBROUTINE precomputeLambda(lambda, eigen)

  REAL*4,INTENT(OUT),DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::lambda
  REAL*4,INTENT(OUT),DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::eigen
  REAL*8,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::internal_lambda
  REAL*8,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::internal_eigen

  INTEGER::force_index

  REAL*8::c, d, e, cc
  c = LOG(SLENDERNESS * SLENDERNESS * EXP(1.0d0))
  d = -c
  e = 2.0d0
  cc = 1.0d0

  internal_lambda(1) = 2.0d0
  internal_eigen(1) = ((d - e - cc * internal_lambda(1)) / 2.0d0) / (d - cc * internal_lambda(1))

  DO force_index = 2, NUMBER_OF_TERMS_IN_FORCE_EXPANSION
    internal_lambda(force_index) = internal_lambda(force_index - 1) + 2.0d0 / force_index;
    internal_eigen(force_index) = ((d - e - cc * internal_lambda(force_index)) / 2.0d0) / (d - cc * internal_lambda(force_index))
  END DO

  lambda = internal_lambda
  eigen = internal_eigen

END SUBROUTINE precomputeLambda

FUNCTION calculateLegendrePolynomial(x, n)

  REAL*8,INTENT(IN)::x
  INTEGER,INTENT(IN)::n

  IF (n == 0) THEN
    calculateLegendrePolynomial = 1.0d0
  ELSEIF (n == 1) THEN
    calculateLegendrePolynomial = x
  ELSEIF (n == 2) THEN
    calculateLegendrePolynomial = (1.0d0 / 2.0d0) * &
      (3.0d0 * x**2 - 1.0d0)
  ELSEIF (n == 3) THEN
    calculateLegendrePolynomial = (1.0d0 / 2.0d0) * &
      (5.0d0 * x**3 - 3.0d0 * x)
  ELSEIF (n == 4) THEN
    calculateLegendrePolynomial = (1.0d0 / 8.0d0) * &
      (35.0d0 * x**4 - 30.0d0 * x**2 + 3.0d0)
  ELSEIF (n == 5) THEN
    calculateLegendrePolynomial = (1.0d0 / 8.0d0) * &
      (63.0d0 * x**5 - 70.0d0 * x**3 + 15.0d0 * x)
  ELSEIF (n == 6) THEN
    calculateLegendrePolynomial = (1.0d0 / 16.0d0) * &
      (231.0d0 * x**6 - 315.0d0 * x**4 + 105.0d0 * x**2 - 5.0d0)
  ELSEIF (n == 7) THEN
    calculateLegendrePolynomial = (1.0d0 / 16.0d0) * &
      (429.0d0 * x**7 - 693.0d0 * x**5 + 315.0d0 * x**3 - 35.0d0 * x)
  ELSEIF (n == 8) THEN
    calculateLegendrePolynomial = (1.0d0 / 128.0d0) * &
      (6435.0d0 * x**8 - 12012.0d0 * x**6 + 6930.0d0 * x**4 - 1260.0d0 * x**2 + 35.0d0)
  ELSE
    PRINT *,"Could not precompute legendre polynomials - n not in range [1..8]: "
  ENDIF

END FUNCTION calculateLegendrePolynomial
