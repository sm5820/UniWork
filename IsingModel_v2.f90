! Name: Swathi Mahashetti
! CID: u5665541
! Ising Model: Periodic Boundaries, CONTAINS EXTENSION RELATED CODE: writes each lattice at every timestep into a folder ising_op.
!              Recommended not to use T > 1000 and recommended small lattice sizes

MODULE IsingModel

    USE ISO_FORTRAN_ENV
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
        
        INTEGER, INTENT(IN) :: n_rows, n_cols
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
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
        INTEGER :: i, j, up, down, left, right
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        REAL :: energy_val, i_j, j_i, spin

        ! Initialise energy value 
        energy_val = 0.0

        ! Since all spins are clamped at +1 in the boundaries, only the energy of the internal spins are calculated here
        ! This function is used for further steps where spin flipping may/may not occur
        ! Spin flipping can only happen to the internal lattice 
        DO i = 1, n_rows
            DO j = 1, n_cols

                up = MOD(i-2 + n_rows, n_rows) + 1     ! Up: If i = 1, it wraps to n_rows
                down = MOD(i, n_rows) + 1              ! Down: If i = n_rows, it wraps to 1
                left = MOD(j-2 + n_cols, n_cols) + 1  ! Left: If j = 1, it wraps to n_cols
                right = MOD(j, n_cols) + 1          ! Right: If j = n_cols, it wraps to 1

                spin = REAL(isinglattice(i,j))
                i_j  = (spin * REAL(isinglattice(up,j))) + (spin * REAL(isinglattice(down,j)))
                j_i = (spin * REAL(isinglattice(i, left))) + (spin * REAL(isinglattice(i, right)))
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

        ! Unlike the fixed boundary model, any state in the lattice can change so random_i/j are defined to return a random 
        ! number anywhere between 1 and n_rows/n_cols. 
        CALL RANDOM_NUMBER(val_i)
        random_i = 1 + FLOOR(val_i * n_rows)

        CALL RANDOM_NUMBER(val_j)
        random_j = 1 + FLOOR(val_j * n_cols)

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

        ! Difference between energies of the states 
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
        INTEGER :: t_step
        INTEGER, DIMENSION(n_rows,n_cols) :: isinglattice 
        REAL :: final_magnetisation, initial_magnetisation, magnetisation_chang
        REAL :: average_initial, average_final, change_average
        REAL, DIMENSION(T+1) :: magnetisation_history
        CHARACTER(LEN=500), INTENT(INOUT) :: filename
        TYPE(run_data) :: run

        ! Creates new directory to store all .nc files
        CALL SYSTEM('mkdir -p ising_op')  
        

        initial_magnetisation = magnetisation(n_rows, n_cols, isinglattice)
        average_initial = avg_magnetisation(n_rows, n_cols, isinglattice)

        ! Changes directory to newly created directory and Monte-Carlo step is iterated for T iterations
        ! while write_display writes lattice information about each step into the files in the folder
        CALL SYSTEM('cd ising_op') 

        magnetisation_history(1) = average_initial

        WRITE(filename, '(A,I0, A)') 'ising_op/File1.nc'
        ! File 1 has the initialised system

        CALL write_display(isinglattice, magnetisation_history,filename, run)

        DO t_step = 2, T+1
            CALL montecarlostep(n_rows, n_cols, isinglattice, J_value, Beta)
            magnetisation_history(t_step) = avg_magnetisation(n_rows, n_cols, isinglattice)
            WRITE(filename, '(A,I0, A)') "ising_op/File", t_step,".nc"
            CALL write_display(isinglattice, magnetisation_history, filename , run)
        END DO 
    
        final_magnetisation = magnetisation(n_rows, n_cols, isinglattice)
        average_final = avg_magnetisation(n_rows, n_cols, isinglattice)

        
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
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: isinglattice 
    CHARACTER(LEN=10) :: init
    CHARACTER(LEN=500) :: filename
    REAL :: beta, J_value
    LOGICAL :: success
    TYPE(run_data) :: run

    CALL parse_args ! Subroutine from module command_line called to turn command line prompts to inputs

    success = get_arg('N', N)
    success = get_arg('T', T)
    success = get_arg('init', init)
    success = get_arg('Beta', beta)
    success = get_arg('J', J_value)

    ! print*, N, T, init, Beta, J_value

    IF (N .le. 0) THEN 
        status = 1
        CALL EXIT(status)
        PRINT*, 'Choose N > 0'
    END IF

    n_rows = N
    n_cols = N
    ALLOCATE(isinglattice(n_rows,n_cols))  ! NxN lattice allocated  

    ! Init can take arguments 'A', 'F' or 'R', which initialise different models
    SELECT CASE (init)
    CASE ('A')
        CALL initialise_chessboard(n_rows, n_cols, isinglattice)
    CASE ('R')
        CALL initialise_random(n_rows, n_cols, isinglattice)
    CASE ('F')
        CALL initialise_plus1(n_rows, n_cols, isinglattice)
    CASE DEFAULT
        STOP
    END SELECT
    
    CALL initialise(run, N, T, init, beta, J_value)
    CALL montecarlo_iter(T, n_rows, n_cols, isinglattice, J_value, Beta, filename, run)

END PROGRAM 
    