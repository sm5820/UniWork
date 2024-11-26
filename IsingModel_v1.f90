! Name: Swathi Mahashetti
! CID: u5665541
! Ising Model: Clamped Boundaries, no extensions 

MODULE IsingModel

    USE ISO_FORTRAN_ENV
    ! USE ascii_display
    USE write_netcdf

    IMPLICIT NONE 

    CONTAINS

    ! Subroutine returns 2D Ising lattice with all spin up (+1)
    ! Call this subroutine for initialising ferromagnet ('F')
    SUBROUTINE initialise_plus1(n_rows, n_cols, isinglattice)
        
        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        INTEGER :: i, j

            DO i = 1, n_rows
                DO j = 1, n_cols
                    isinglattice(i,j) = 1
                END DO
            END DO

        ! Model tested by display_3val_array subroutine of ascii_display module 
        ! CALL display_3val_array(isinglattice, .TRUE., 'Ferromagnetic Ising model: ')

        END SUBROUTINE initialise_plus1

    
    ! Subroutine  returns 2D Ising lattice with randomly assigned spin up (+1) and down (-1) states
    ! Call this subroutine for initialising random state ('R')
    SUBROUTINE initialise_random(n_rows, n_cols, isinglattice)
        
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER :: i, j
        REAL :: rand_val

        DO i = 1, n_rows
            DO j = 1, n_cols

                ! Call FORTRAN's intrinsic subroutine to generate a random number 
                CALL RANDOM_NUMBER(rand_val) 
                ! rand_val is a random number in the range [0, 1) creating a bias
                ! This is adjusted using FLOOR (gives the greatest integer less than or equal to arguments)
                ! 0 <= rand_val < 0.5: -1 states are initialised 
                ! 0.5 <= rand_val < 1: +1 states are initialised 
                isinglattice(i,j) = (FLOOR(2 * rand_val)*2) - 1

                ! Print to make sure random {-1,1} values are being assigned for all combinations of i, j
                ! print*, 'isinglattice(', i, ',', j, ') =', isinglattice(i, j)

            END DO
        END DO

        ! Model tested by display_3val_array subroutine of ascii_display module 
        ! CALL display_3val_array(isinglattice, .TRUE., 'Random Ising model: ')
    
    END SUBROUTINE initialise_random

    ! Subroutine  returns 2D Ising lattice with alternating spin up (+1) and down (-1) states
    ! Call this subroutine for initialising antiferromagnet state ('A')
    SUBROUTINE initialise_chessboard(n_rows, n_cols, isinglattice)
        
        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        INTEGER :: i, j

        DO i = 1, n_rows
            DO j = 1, n_cols

                ! If the remainder of i + j by 2 is not 0 (i.e., i+j is odd) the domain is spin down (-1)
                ! Else spin up (+1)
                IF (MODULO((i+j), 2) /= 0) THEN
                    isinglattice(i,j) = -1
                ELSE
                    isinglattice(i,j) = 1

                END IF
            END DO
        END DO

        ! Model tested by display_3val_array subroutine of ascii_display module 
        ! CALL display_3val_array(isinglattice, .TRUE., 'Antiferromagnet Ising model: ')

    END SUBROUTINE initialise_chessboard

    ! Subroutine clamps edges to spin up (+1) states 
    SUBROUTINE edges(n_rows, n_cols, isinglattice)

        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        INTEGER :: i, j

        DO i = 1, n_rows
            DO j = 1, n_cols

                ! Checks if state is in 1st or last row; 1st or last column 
                ! If any of the conditions are true spin up state is fixed 
                IF (i == 1 .OR. i == n_rows .OR. j == 1 .OR. j == n_cols) THEN 
                    isinglattice(i,j) = 1
                END IF
            END DO
        END DO

    END SUBROUTINE edges 

    ! Function - returns total magentisation of the lattice 
    FUNCTION magnetisation(n_rows, n_cols, isinglattice)  

        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        REAL :: magnetisation

        magnetisation = SUM(REAL(isinglattice))

    END FUNCTION magnetisation

    ! Function - returns magnetisation per state of the lattice
    FUNCTION avg_magnetisation(n_rows, n_cols, isinglattice)  

        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        REAL :: avg_magnetisation

        avg_magnetisation = SUM(REAL(isinglattice))/(n_cols*n_rows)

    END FUNCTION avg_magnetisation

    ! Function returns the energy of the lattice
    FUNCTION energy(n_rows, n_cols, isinglattice, J_value) RESULT(energy_val)

        INTEGER, INTENT(IN) :: n_rows, n_cols
        REAL, INTENT(IN) :: J_value
        INTEGER :: i, j
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        REAL :: energy_val, i_j, j_i, spin

        ! Initialise energy value 
        energy_val = 0.0

        ! Since all spins are clamped at +1 in the boundaries, only the energy of the internal spins are calculated here
        ! This function is used for further steps where spin flipping may/may not occur
        ! Spin flipping can only happen to the internal lattice 
        DO i = 2, n_rows-1
            DO j = 2, n_cols-1

                ! Energy of state with respect to its 4 nearest neighbours are calculated 
                ! Here indices are defined keeping the boundary condition and on the assumption that the 
                ! internal spins interact with the boundary. 
                ! For states i-1 and j-1 for i = 2 or j = 2 (first indices of internal spin lattice): up and left are +1
                ! For states i+1 and j+1 for i = n_rows-1 or j = n_cols-1 (last indices): down and left states are +1
                
                spin = REAL(isinglattice(i,j), KIND = REAL32)
                i_j  = (spin * REAL(isinglattice(i-1,j), KIND = REAL32)) + (spin * REAL(isinglattice(i+1,j), KIND = REAL32))
                j_i = (spin * REAL(isinglattice(i, j-1), KIND = REAL32)) + (spin * REAL(isinglattice(i, j+1), KIND = REAL32))
                energy_val = energy_val + (i_j + j_i)
                ! print*, energy_val

            END DO
        END DO
        
        ! Energy is multiplied by the coupling constant and divided by 2.0 to ensure each interaction is only counted once  
        energy_val =  -(J_value * energy_val) / 2.0
        
    END FUNCTION energy 

    SUBROUTINE montecarlostep(n_rows, n_cols, isinglattice, J_value, Beta)

        INTEGER, INTENT(IN) :: n_rows, n_cols
        REAL, INTENT(IN) :: J_value, Beta
        INTEGER :: random_i, random_j
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        REAL :: val_i, val_j, random_prob, energy_0, energy_1, delta_e, exp_val 

        ! A random number in the range [2, n_rows-1] is called for the random row index, boundaries are kept constant 
        CALL RANDOM_NUMBER(val_i)
        random_i = 2 + FLOOR(val_i * (n_rows-2))

        ! A random number in the range [2, n_cols-1] is called for the random column index, boundaries are kept constant 
        CALL RANDOM_NUMBER(val_j)
        random_j = 2 + FLOOR(val_j * (n_cols-2))

        ! A random probability of for the Monte-Carlo step in the range [1e-10, 1 - 1e-10] is called
        ! Since random_number calls numbers between [0, 1), and we need the probability to be in the range (0, 1] 
        ! the number is adjusted 
        CALL RANDOM_NUMBER(random_prob)
        random_prob = (random_prob * (1.0 - 1e-10)) + 1e-10

        ! PRINT*, random_prob, random_i, random_j

        ! Initial energy before any spin flip is calculated 
        energy_0 = energy(n_rows, n_cols, isinglattice, J_value)

        ! Random state is flipped 
        isinglattice(random_i, random_j) = -isinglattice(random_i, random_j)

        ! Energy after flip of random state 
        energy_1 = energy(n_rows, n_cols, isinglattice, J_value)

        ! Difference between senergies of the states 
        delta_e = energy_1 - energy_0

        ! print '(F8.3)', delta_e
        
        IF (delta_e < 0) THEN 
            RETURN ! Early exit - spin flip accepted because it is energetically favourable.  

        ELSE IF (delta_e > 0) THEN

            Exp_val = EXP(- delta_e * Beta) ! Probability of spin flip calculated to compare to a random probability 
            ! Print*, exp_val

            IF (random_prob <= Exp_val) THEN
                RETURN ! Early exit spin flip accepted, probability of the flip is high 

            ELSE IF (random_prob > Exp_val) THEN
                isinglattice(random_i, random_j) = -isinglattice(random_i, random_j) ! Spin flip rejected 
                ! Probability of the spin flip is low, initial state restored. 

            END IF
        END IF


    END SUBROUTINE montecarlostep

    SUBROUTINE montecarlo_iter(T, n_rows, n_cols, isinglattice, J_value, Beta, filename, run)

        INTEGER, INTENT(IN) :: n_rows, n_cols, T
        REAL, INTENT(IN) :: J_value, Beta
        INTEGER, DIMENSION(n_rows,n_cols), INTENT(INOUT) :: isinglattice
        CHARACTER(LEN=500), INTENT(IN) :: filename 

        TYPE(run_data), INTENT(IN) :: run


        INTEGER :: t_step
        REAL, DIMENSION(T) :: magnetisation_history
        REAL :: final_magnetisation, initial_magnetisation, magnetisation_chang
        REAL :: average_initial, average_final, change_average
        
        ! Edges clamped for all cases 
        CALL edges(n_rows, n_cols, isinglattice)


        initial_magnetisation = magnetisation(n_rows, n_cols, isinglattice)
        average_initial = avg_magnetisation(n_rows, n_cols, isinglattice)

        ! Monte-Carlo step iterated for T iterations 
        DO t_step = 1, T
            CALL montecarlostep(n_rows, n_cols, isinglattice, J_value, Beta)
            magnetisation_history(t_step) = avg_magnetisation(n_rows, n_cols, isinglattice)
        END DO 
    
        final_magnetisation = magnetisation(n_rows, n_cols, isinglattice)
        average_final = avg_magnetisation(n_rows, n_cols, isinglattice)

        CALL write_display(isinglattice, magnetisation_history, filename, run)
        ! CALL display_3val_array(isinglattice, .TRUE., 'After T trials: ')
        

        magnetisation_chang = ABS(final_magnetisation - initial_magnetisation)*100/initial_magnetisation
        change_average = ABS(average_final - average_initial)*100/average_initial


    END SUBROUTINE montecarlo_iter
    

END MODULE 

PROGRAM test 

    USE IsingModel
    USE command_line

    IMPLICIT NONE

    INTEGER:: n_rows, n_cols, N, T, status
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: ising 
    CHARACTER(LEN=10) :: init
    REAL :: beta, J_value
    LOGICAL :: success
    CHARACTER(LEN=500) :: filename 
    TYPE(run_data) :: run

    CALL parse_args ! Subroutine from module command_line called to turn command line prompts to inputs

    success = get_arg('N', N)
    success = get_arg('T', T)
    success = get_arg('init', init)
    success = get_arg('Beta', beta)
    success = get_arg('J', J_value)

    ! print*, N, T, init, Beta, J_value

    ! Checks if N>0
    IF (N .le. 0) THEN 
        status = 1
        CALL EXIT(status)
        PRINT*, 'Choose N > 0'
    END IF

    n_rows = N
    n_cols = N
    ALLOCATE(ising(n_rows,n_cols))  ! NxN lattice allocated  

    ! Init can take arguments 'A', 'F' or 'R', which initialise different models
    SELECT CASE (init)
    CASE ('A')
        CALL initialise_chessboard(n_rows, n_cols, ising)
    CASE ('R')
        CALL initialise_random(n_rows, n_cols, ising)
    CASE ('F')
        CALL initialise_plus1(n_rows, n_cols, ising)
    CASE DEFAULT
        STOP
    END SELECT

    filename = 'ising_fx.nc'
    ! Initialises run_data types 
    CALL initialise(run, N, T, init, beta, J_value) 
    ! Simulation 
    CALL montecarlo_iter(T, n_rows, n_cols, ising, J_value, Beta, filename, run)

END PROGRAM 
    