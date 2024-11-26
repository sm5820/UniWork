MODULE write_netcdf

    USE ISO_FORTRAN_ENV
    USE netcdf

    IMPLICIT NONE 

    ! run_data is a type for storing the input arguments which are 
    TYPE :: run_data 
      CHARACTER(LEN=10) :: init_arg
      INTEGER(INT64) :: N_size, T_arg
      REAL(REAL64) :: J_arg, Beta_arg
    END TYPE run_data 

    CONTAINS
    
    SUBROUTINE initialise(run, size, iter, initial_lat, beta_val, j_val)
      
      TYPE(run_data), INTENT(OUT) :: run
      CHARACTER(LEN=10), INTENT(IN) :: initial_lat
      INTEGER, INTENT(IN) :: size, iter
      REAL, INTENT(IN) :: beta_val, j_val

      run%N_size = size
      run%T_arg = iter
      run%init_arg = initial_lat
      run%Beta_arg = beta_val
      run%J_arg = j_val

    END SUBROUTINE initialise
    
    SUBROUTINE write_display(isinglattice, magnetisation_history, filename, run)
      
      ! Define intent(in) variables
      INTEGER, DIMENSION(:,:), INTENT(IN) :: isinglattice
      REAL, DIMENSION(:), INTENT(IN) :: magnetisation_history
      CHARACTER(LEN=500), INTENT(IN) :: filename
      ! SUBROUTINE initialise needs to be called to define run before this subroutine is called 
      TYPE(run_data), INTENT(IN) :: run

      ! Define subroutine variables
      INTEGER, PARAMETER :: ndims = 3 ! For 'x', 'y', and 't'
      INTEGER, DIMENSION(ndims) :: sizes, dim_ids, dim_size
      CHARACTER(LEN=10), DIMENSION(ndims) :: dims=(/"x", "y", "t" /) 
      INTEGER :: file_id, var_id_ising, var_id_mag, i
      
      ! Define intent(out) variable
      INTEGER :: ierr

      ! EXTENSION: Logicals defined if dimensions have been correctly defined, if data is written into the file and if file exists 
      LOGICAL :: file_created, write_global_att, write_data, def_dims

      ! Set logicals that check writing as true 
      write_global_att = .TRUE.
      write_data = .TRUE.
      def_dims = .TRUE. 
          
      ! Collecting data for defining the dimensions of the NetCDF file 
      ! sizes has dimension = 3, first two of which are occupied by the lattice, given by SHAPE() which returns X and Y dimensions. 
      ! 3: is given by a 1D array magnetisation_history which is the magnetisation at each timestep. 
      sizes(1:2) = SHAPE(isinglattice) 
      sizes(3:) = SHAPE(magnetisation_history) 
      

      ! Create the netCDF file
      ! NF90_CLOBBER overwrites any existing dataset with the same file name, and buffer and cache accesses for efficiency.
      ierr = nf90_create(filename, NF90_CLOBBER, file_id)

      ! Check if there is an error in defining the file 
      IF (ierr /= nf90_noerr) THEN
        ! If ierr = nf90_noerr (which is 0 then - successful), else an error string is returned
        ! NF90_STRERROR returns a static reference to an error message string corresponding to an integer netCDF error status or to a system error number.
        ! TRIM(): Fortran initrinsic function that removes trailing spaces from character string
        PRINT*, 'Error creating file: ', TRIM(nf90_strerror(ierr))
        RETURN

      END IF
      
      ! Adds a new dimension to an open netCDF dataset in define mode. 
      ! It returns a dimension ID, given the netCDF ID, the dimension name, and the dimension length. 
      DO i = 1, ndims
        ierr = nf90_def_dim(file_id, dims(i), sizes(i), dim_ids(i))
        IF (ierr /= nf90_noerr) THEN
          PRINT*, 'Error defining dimensions: ', TRIM(nf90_strerror(ierr))
          RETURN
        END IF
      END DO 
    

      ! Adds a new variable to an open netCDF dataset in define mode. 
      ! It returns a variable ID, the variable name, the variable type, the number of dimensions, and a list of the dimension IDs 
      ierr = nf90_def_var(file_id, 'ISING LATTICE', NF90_INT, dim_ids(1:2), var_id_ising)
      ierr = nf90_def_var(file_id, 'MAGNETISATION', NF90_REAL, dim_ids(3), var_id_mag)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, 'Error defining variables: ', TRIM(nf90_strerror(ierr))
        RETURN 
      END IF

      ! Information about length of dimension is written, checks if dimension size for ids 1,2 = N and id 3 = T
      ierr = nf90_inquire_dimension(file_id, dim_ids(1),  LEN = dim_size(1))
      IF (ierr /= nf90_noerr.OR. dim_size(1) .NE. run%N_size) THEN
        def_dims = .FALSE.
      END IF
    
      ierr = nf90_inquire_dimension(file_id, dim_ids(2),  LEN= dim_size(2))
      IF (ierr /= nf90_noerr.OR. dim_size(2) .NE. run%N_size) THEN
        def_dims = .FALSE.
      END IF

      ierr = nf90_inquire_dimension(file_id, dim_ids(3),  LEN = dim_size(3))
      IF (ierr /= nf90_noerr.OR. dim_size(3) .NE. run%T_arg) THEN
        def_dims = .FALSE.
      END IF
    
      ! Global attributes are not associated with any variables and are identified by using NF90_GLOBAL 
      ! They are associated with the data as a whole, so command line arguments are global attributes here
      ! They are written in the define mode and not the data mode  
      ierr = nf90_put_att(file_id, NF90_GLOBAL, 'N', run%N_size)
      ierr = nf90_put_att(file_id, NF90_GLOBAL, 'T', run%T_arg)
      ierr = nf90_put_att(file_id, NF90_GLOBAL, 'init', run%init_arg)
      ierr = nf90_put_att(file_id, NF90_GLOBAL, 'Beta', run%Beta_arg)
      ierr = nf90_put_att(file_id, NF90_GLOBAL, 'J', run%J_arg)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, 'Error writing global attributes: ', TRIM(nf90_strerror(ierr))
        ! If global attribute hasn't been written, returns logical as false (ext) 
        write_global_att = .FALSE. 
        RETURN 
      END IF

      ierr = nf90_enddef(file_id)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF 
      ! Finish defining metadata, end define mode.   
      
      ! Adds or changes a variable attribute or global attribute of an open netCDF dataset. 
      ! If this attribute is new, or if the space required to store the attribute is greater than before, the netCDF dataset must be in define mode.
      ierr = nf90_put_var(file_id, var_id_ising, isinglattice)
      ierr = nf90_put_var(file_id, var_id_mag, magnetisation_history)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        ! Similarly logical becomes false here if variables haven't been written 
        write_data = .TRUE.
        RETURN
      END IF

      ierr = nf90_close(file_id)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    
    ! Only if the required data has been written in the file, code enquires if the file has been created. 
    IF (write_data .AND. write_global_att .AND. def_dims) THEN 
      INQUIRE(FILE=filename, EXIST=file_created)
    ELSE
      RETURN
    END IF

    ! Print statement to confirm creation of file 
    PRINT*, 'File successfully created: ', file_created

    END SUBROUTINE write_display 
    
END MODULE write_netcdf