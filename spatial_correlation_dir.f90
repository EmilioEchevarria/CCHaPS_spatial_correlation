PROGRAM spatial_corr

    USE netcdf
    USE mpi
    USE iso_fortran_env, only: real64, int64
    IMPLICIT NONE

    INTEGER                        :: the_node, the_node_local
    CHARACTER(LEN=200)             :: node_str, output_var, logfile, line, param_name, &
                                      param_value
    CHARACTER(LEN=250),ALLOCATABLE :: files_in(:)
    REAL                           :: the_lon, the_lat, max_dist, sum_x, sum_x2,       &
                                      sum_x_g, sum_x2_g
    INTEGER                        :: start_year, end_year, n_years, C, sum_N, sum_N_g
    INTEGER                        :: fid, status, dimID_time, dimID_node, n_nodes,    &
                                      n_times, time_ID, nodex_ID, nodey_ID, idx,       &
                                      count_close_nodes, idx_node0, idx_node1, dir_ID,  &
                                      n_close_nodes, I, J, N, F, Y, M, n_times_total,  &
                                      num_files, num_files_per_node, extra_files,      &
                                      err=0, eq_pos, ex_pos, EXTCDE
    INTEGER, DIMENSION(2)          :: start_nodes, count_nodes
    REAL,ALLOCATABLE,DIMENSION(:)  :: schism_lons, schism_lats, distances, correlations, &
                                      close_node_indices, nearest_schism_lons,           &
                                      nearest_schism_lats, sum_y, sum_y2, sum_xy,        &
                                      sum_y_g, sum_y2_g, sum_xy_g

    REAL(real64)                   :: sum_sin_x, sum_cos_x, sum_sin_x_g, sum_cos_x_g
    REAL(real64),ALLOCATABLE,DIMENSION(:)  :: sum_sin_y, sum_cos_y, sum_sin_y_g, sum_cos_y_g
    REAL(real64)                   :: den_x, den_x_g
    REAL(real64),ALLOCATABLE       :: num(:), den_y(:), num_g(:), den_y_g(:), mu_y(:)
    REAL(real64)                   :: sx, sy, mu_x
    INTEGER :: T
    real(real64), parameter :: pi = acos(-1.0_real64)
    real(real64), parameter :: two_pi = 2.0_real64*pi
    real(real64) :: scale, offset,  fillv, missv
    logical :: has_scale, has_offset, has_fill, has_miss

    INTEGER,ALLOCATABLE,DIMENSION(:) :: nearest_schism_idxs
    REAL(real64),ALLOCATABLE       :: dir_around_node(:,:)
    CHARACTER(LEN=100)             :: file_out
        
    INTEGER                        :: mpi_err, mpi_rank, mpi_size, start_idx_mpi, &
                                      end_idx_mpi
    DOUBLE PRECISION               :: mpi_start_time, mpi_end_time, mpi_elapsed_time, & 
                                      mpi_max_time

    CHARACTER(LEN=250)             :: base_dir
    CHARACTER(LEN=4)               :: y4
    
    ! Defaults
    scale = 1.0_real64 ; offset = 0.0_real64
    has_scale = .false. ; has_offset = .false.
    has_fill  = .false. ; has_miss  = .false.

    CALL MPI_INIT(mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_err)

    OPEN( UNIT=20, FILE='spatial_corr.inp', STATUS="old", FORM="formatted", iostat=err )
    IF( err .NE. 0 ) THEN
        PRINT*,"[ERROR] Failed to open file spatial_corr.inp"
        PRINT*,"[ERROR] IOERROR: ", err
        EXTCDE = 2
        CALL EXIT(EXTCDE)
    ENDIF

    DO
        READ(20, '(A)', iostat=err) line
        IF (err .NE. 0) EXIT

        eq_pos = INDEX(line, '=')
        param_name = TRIM(line(:eq_pos-1))
        param_value = TRIM(line(eq_pos+1:))

        ! Removing any trailing comments
        ex_pos = INDEX(param_value, '!')
        IF (ex_pos .NE. 0) THEN
            param_value = param_value(:ex_pos-1)
        END IF

        param_value = TRIM(param_value)

        SELECT CASE(param_name)
        CASE ('output_var')
            READ(param_value, '(A)') output_var
        CASE ('start_year')
            READ(param_value, *) start_year
        CASE ('end_year')
            READ(param_value, *) end_year
        CASE ('the_lon')
            READ(param_value, *) the_lon
        CASE ('the_lat')
            READ(param_value, *) the_lat
        CASE ('max_dist')
            READ(param_value, *) max_dist
        END SELECT

    END DO

    CLOSE( 20 )
    
    base_dir = '/datasets/work/ev-acs-wp3-cchaps/reference/release/WP3/CSIRO/hindcast/ERA5/BARRA-R2_ORAS5_TPXO_WHACS/BARRA-R2_WHACS/SCHISM-WWMIII-v5.9/national-mesh/1hr/'//TRIM(output_var)

    IF ( mpi_rank .EQ. 0 ) THEN
        WRITE(*,*) 'base_dir = ', base_dir
        WRITE(*,*) 'output_var = ', output_var
        WRITE(*,*) 'start_year = ', start_year
        WRITE(*,*) 'end_year = ', end_year
        WRITE(*,*) 'the_lon = ', the_lon
        WRITE(*,*) 'the_lat = ', the_lat
        WRITE(*,*) 'max_dist = ', max_dist
    END IF

    n_years = end_year-start_year+1
    ALLOCATE( files_in(n_years) )
    C = 0
    DO Y = start_year, end_year
        C = C + 1
        WRITE(y4,'(I4.4)') Y
        files_in(C) = TRIM(base_dir)//'/'// &
                      TRIM(output_var)//'_cchaps_hindcast_BARRA-R2_WHACS_ERA5_1hr_'// &
                      y4//'01010000-'//y4//'12312300.nc'
    END DO

    WRITE(*,*) TRIM(files_in(1))
    
    status = NF90_OPEN(TRIM(files_in(1)), 0, fid)
    CALL nc_error(status, .TRUE., "Error reading netcdf file")

    ! Read dimensions IDs
    status = NF90_INQ_DIMID(fid, "time", dimID_time)
    CALL nc_error(status, .TRUE., "Error inquiring time ID")
    status = NF90_INQ_DIMID(fid, "nSCHISM_hgrid_node", dimID_node)
    CALL nc_error(status, .TRUE., "Error inquiring nSCHISM_hgrid_node ID")
    ! Get length of times and nSCHISM_hgrid_node
    status = NF90_INQUIRE_DIMENSION(fid, dimID_time, len=n_times)
    CALL nc_error(status, .TRUE., "Error inquiring time dimensions")
    status = NF90_INQUIRE_DIMENSION(fid, dimID_node, len=n_nodes)
    CALL nc_error(status, .TRUE., "Error inquiring nSCHISM_hgrid_node dimensions")

    ! Read variable IDs
    status = NF90_INQ_VARID(fid, "time", time_ID)
    CALL nc_error(status, .TRUE., "Error inquiring time var ID")
    !
    status = NF90_INQ_VARID(fid, "SCHISM_hgrid_node_x", nodex_ID)
    CALL nc_error(status, .TRUE., "Error inquiring SCHISM_hgrid_node_x var ID")
    status = NF90_INQ_VARID(fid, "SCHISM_hgrid_node_y", nodey_ID)
    CALL nc_error(status, .TRUE., "Error inquiring SCHISM_hgrid_node_y var ID")

    ! Read variable values
    ALLOCATE( schism_lons(n_nodes) )
    ALLOCATE( schism_lats(n_nodes) )
    status = NF90_GET_VAR(fid, nodex_ID, schism_lons)
    CALL nc_error(status, .TRUE., "Error getting schism_lons values")
    status = NF90_GET_VAR(fid, nodey_ID, schism_lats)
    CALL nc_error(status, .TRUE., "Error getting schism_lats values")

    ! Calculate distances to the_node
    ALLOCATE( distances(n_nodes) )
    !ALLOCATE( close_node_indices(n_nodes) )

    count_close_nodes = 0 
    DO I = 1, n_nodes
        CALL haversine(the_lat, the_lon, &
                    schism_lats(I), schism_lons(I), distances(I))
        IF (distances(I) < max_dist) THEN
            count_close_nodes = count_close_nodes + 1
            !close_node_indices(count_close_nodes) = I
        END IF
    END DO

    WRITE(*,*) 'count_close_nodes = ', count_close_nodes

    ALLOCATE( close_node_indices(count_close_nodes) )
    C = 0
    DO I = 1, n_nodes
        IF (distances(I) < max_dist) THEN
            C = C + 1
            close_node_indices(C) = I
        END IF
    END DO

    the_node = MINLOC(distances, DIM=1)

    WRITE(node_str, '(I0)') the_node
    file_out = 'schism_corrs_node_' // TRIM(node_str) // '.txt'

    DEALLOCATE( distances )

    WRITE(*,*) 'SIZE(close_node_indices) = ', SIZE(close_node_indices)

    idx_node0 = MINVAL(close_node_indices(1:count_close_nodes)) 
    idx_node1 = MAXVAL(close_node_indices(1:count_close_nodes))
    n_close_nodes = idx_node1 - idx_node0 + 1
    the_node_local = the_node - idx_node0 + 1
    WRITE(*,*) the_node, idx_node0, the_node_local

    WRITE(*,*) 'idx_node0, idx_node1, n_close_nodes = ', idx_node0, idx_node1, n_close_nodes

    ALLOCATE( nearest_schism_lons(n_close_nodes) )
    ALLOCATE( nearest_schism_lats(n_close_nodes) )
    ALLOCATE( nearest_schism_idxs(n_close_nodes) )
    
    nearest_schism_lons = schism_lons(idx_node0:idx_node1)
    nearest_schism_lats = schism_lats(idx_node0:idx_node1)
    nearest_schism_idxs = [(i, i = idx_node0, idx_node1+1)]
    DEALLOCATE( schism_lons )
    DEALLOCATE( schism_lats )

    n_times_total = 0
    DO F=1,SIZE(files_in)
        WRITE(*,*) TRIM(files_in(F))
        status = NF90_OPEN(TRIM(files_in(F)), 0, fid)
        CALL nc_error(status, .TRUE., "Error reading netcdf file")
        ! Read dimensions IDs
        status = NF90_INQ_DIMID(fid, "time", dimID_time)
        CALL nc_error(status, .TRUE., "Error inquiring time ID")
        ! Get length of times and nSCHISM_hgrid_node
        status = NF90_INQUIRE_DIMENSION(fid, dimID_time, len=n_times)
        CALL nc_error(status, .TRUE., "Error inquiring time dimensions")
        WRITE(*,*) 'ntimes = ', n_times
        n_times_total = n_times_total + n_times
        
    END DO

    ALLOCATE( sum_sin_y(n_close_nodes) )
    ALLOCATE( sum_cos_y(n_close_nodes) )
    sum_sin_y = 0.0
    sum_cos_y = 0.0

    IF (mpi_rank .EQ. 0) THEN
        ALLOCATE( sum_sin_y_g(n_close_nodes) )
        ALLOCATE( sum_cos_y_g(n_close_nodes) )
    END IF

    sum_sin_x = 0.0
    sum_cos_x = 0.0
    sum_N     = 0

    mpi_start_time = MPI_WTIME()

    WRITE(logfile, "('output_',I0,'.log')") mpi_rank
    OPEN(UNIT=25, FILE=logfile, STATUS='replace', ACTION='write')

    IF ( mpi_rank .EQ. 0 ) WRITE(25,*) 'the_node = ', the_node
    IF ( mpi_rank .EQ. 0 ) WRITE(25,*) 'Loop through files, extract dir at nearest nodes:'

    num_files = SIZE(files_in)
    num_files_per_node = num_files / mpi_size
    extra_files = MOD(num_files, mpi_size)

    IF ( mpi_rank < extra_files ) THEN
        start_idx_mpi = mpi_rank * (num_files_per_node + 1) + 1
        end_idx_mpi = start_idx_mpi + num_files_per_node
    ELSE
        start_idx_mpi = mpi_rank * num_files_per_node + extra_files + 1
        end_idx_mpi = start_idx_mpi + num_files_per_node - 1
    END IF
    !start_idx_mpi = mpi_rank * (num_files / mpi_size) + 1
    !end_idx_mpi = (mpi_rank + 1) * (num_files / mpi_size)

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    WRITE(25,*) 'RANK = ', mpi_rank
    WRITE(25,*) 'end_idx = ', end_idx_mpi
    WRITE(*,*) ' len(files_in) = ', num_files
    DO F = start_idx_mpi, end_idx_mpi
        WRITE(25,*) 'Reading ', TRIM(files_in(F))
        status = NF90_OPEN(TRIM(files_in(F)), 0, fid)
        CALL nc_error(status, .TRUE., "Error reading netcdf file")
        ! Read dimensions IDs
        status = NF90_INQ_DIMID(fid, "time", dimID_time)
        CALL nc_error(status, .TRUE., "Error inquiring time ID")
        ! Get length of times and nSCHISM_hgrid_node
        status = NF90_INQUIRE_DIMENSION(fid, dimID_time, len=n_times)
        CALL nc_error(status, .TRUE., "Error inquiring time dimensions")

        start_nodes = (/idx_node0, 1/) ! First node idx, irst time step idx (1)
        count_nodes = (/n_close_nodes, n_times/) ! n_nodes_around_THE_node, n_times 

        ! Read dir around closest nodes:
        ALLOCATE( dir_around_node(n_close_nodes, n_times) )

        status = NF90_INQ_VARID(fid, TRIM(output_var), dir_ID)
        CALL nc_error(status, .TRUE., "Error inquiring output_var ID")
        status = NF90_GET_ATT(fid, dir_ID, "scale_factor", scale)
        IF (status == NF90_NOERR) has_scale = .TRUE.
        status = NF90_GET_ATT(fid, dir_ID, "add_offset", offset)
        IF (status == NF90_NOERR) has_offset = .TRUE.
        status = NF90_GET_ATT(fid, dir_ID, "_FillValue", fillv)
        IF (status == NF90_NOERR) has_fill = .TRUE.
        status = NF90_GET_ATT(fid, dir_ID, "missing_value", missv)
        IF (status == NF90_NOERR) has_miss = .TRUE.

        status = NF90_GET_VAR(fid, dir_ID, dir_around_node, start=start_nodes, count=count_nodes)
        CALL nc_error(status, .TRUE., "Error getting dir values")

        IF (has_scale .or. has_offset) THEN
            IF (has_fill .and. has_miss) THEN
                where (dir_around_node /= fillv .and. dir_around_node /= missv)
                    dir_around_node = dir_around_node * scale + offset
                end where
            ELSE IF (has_fill) THEN
                where (dir_around_node /= fillv)
                    dir_around_node = dir_around_node * scale + offset
                end where
            ELSE IF (has_miss) THEN
                where (dir_around_node /= missv)
                    dir_around_node = dir_around_node * scale + offset
                end where
            ELSE
                dir_around_node = dir_around_node * scale + offset
            END IF
        END IF

        ! Convert angles to radians:
        dir_around_node = dir_around_node * (3.14159265358979323846 / 180.0)

        sum_N = sum_N + n_times

        ! Reference node (buoy) circular sums
        sum_sin_x = sum_sin_x + SUM( SIN(dir_around_node(the_node_local,:)) )
        sum_cos_x = sum_cos_x + SUM( COS(dir_around_node(the_node_local,:)) )

        DO N = 1, n_close_nodes
            sum_sin_y(N) = sum_sin_y(N) + SUM( SIN(dir_around_node(N,:)) )
            sum_cos_y(N) = sum_cos_y(N) + SUM( COS(dir_around_node(N,:)) )

            IF (MOD(N,10000) == 0 .OR. N == n_close_nodes) THEN
                WRITE(25,'("Progress: ",I0,"/",I0," (",F5.1,"%)")') &
                    N, n_close_nodes, 100.0*N/n_close_nodes
                FLUSH(25)
            END IF

        END DO 

        DEALLOCATE( dir_around_node )
        status = NF90_CLOSE(fid)
        CALL nc_error(status, .TRUE., "Error closing netcdf file")

    END DO  

    mpi_end_time = MPI_WTIME()
    mpi_elapsed_time = mpi_end_time - mpi_start_time
    CALL MPI_REDUCE(mpi_elapsed_time, mpi_max_time, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)

    CALL MPI_REDUCE(sum_N,     sum_N_g,     1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_sin_x, sum_sin_x_g, 1, MPI_DOUBLE_PRECISION,    MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_cos_x, sum_cos_x_g, 1, MPI_DOUBLE_PRECISION,    MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_sin_y, sum_sin_y_g, n_close_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_cos_y, sum_cos_y_g, n_close_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    IF (mpi_rank .EQ. 0) THEN
        mu_x = ATAN2(sum_sin_x_g, sum_cos_x_g)

        WRITE(25,*) 'sum_sin_x_g, sum_cos_x_g, mu_x = ', sum_sin_x_g, sum_cos_x_g, mu_x*180/3.14
        ALLOCATE(mu_y(n_close_nodes))
        DO N = 1, n_close_nodes
            mu_y(N) = ATAN2(sum_sin_y_g(N), sum_cos_y_g(N))
        END DO
    END IF

    ! Make sure mu_y exists on all ranks
    IF (mpi_rank .NE. 0) THEN
        ALLOCATE(mu_y(n_close_nodes))
    END IF

    CALL MPI_BCAST(mu_x, 1,             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_BCAST(mu_y, n_close_nodes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    ALLOCATE( num(n_close_nodes) )
    ALLOCATE( den_y(n_close_nodes) )
    ALLOCATE( num_g(n_close_nodes) )
    ALLOCATE( den_y_g(n_close_nodes) )
    num   = 0.0
    den_y = 0.0
    num_g   = 0.0
    den_y_g = 0.0
    den_x = 0.0

    DO F = start_idx_mpi, end_idx_mpi
        WRITE(25,*) 'Reading ', TRIM(files_in(F))
        status = NF90_OPEN(TRIM(files_in(F)), 0, fid)
        CALL nc_error(status, .TRUE., "Error reading netcdf file")
        ! Read dimensions IDs
        status = NF90_INQ_DIMID(fid, "time", dimID_time)
        CALL nc_error(status, .TRUE., "Error inquiring time ID")
        ! Get length of times and nSCHISM_hgrid_node
        status = NF90_INQUIRE_DIMENSION(fid, dimID_time, len=n_times)
        CALL nc_error(status, .TRUE., "Error inquiring time dimensions")

        start_nodes = (/idx_node0, 1/) ! First node idx, irst time step idx (1)
        count_nodes = (/n_close_nodes, n_times/) ! n_nodes_around_THE_node, n_times 

        ! Read dir around closest nodes:
        ALLOCATE( dir_around_node(n_close_nodes, n_times) )

        status = NF90_INQ_VARID(fid, TRIM(output_var), dir_ID)
        CALL nc_error(status, .TRUE., "Error inquiring output_var ID")
        status = NF90_GET_ATT(fid, dir_ID, "scale_factor", scale)
        IF (status == NF90_NOERR) has_scale = .TRUE.
        status = NF90_GET_ATT(fid, dir_ID, "add_offset", offset)
        IF (status == NF90_NOERR) has_offset = .TRUE.
        status = NF90_GET_ATT(fid, dir_ID, "_FillValue", fillv)
        IF (status == NF90_NOERR) has_fill = .TRUE.
        status = NF90_GET_ATT(fid, dir_ID, "missing_value", missv)
        IF (status == NF90_NOERR) has_miss = .TRUE.

        status = NF90_GET_VAR(fid, dir_ID, dir_around_node, start=start_nodes, count=count_nodes)
        CALL nc_error(status, .TRUE., "Error getting hs values")

        IF (has_scale .or. has_offset) THEN
            IF (has_fill .and. has_miss) THEN
                where (dir_around_node /= fillv .and. dir_around_node /= missv)
                    dir_around_node = dir_around_node * scale + offset
                end where
            ELSE IF (has_fill) THEN
                where (dir_around_node /= fillv)
                    dir_around_node = dir_around_node * scale + offset
                end where
            ELSE IF (has_miss) THEN
                where (dir_around_node /= missv)
                    dir_around_node = dir_around_node * scale + offset
                end where
            ELSE
                dir_around_node = dir_around_node * scale + offset
            END IF
        END IF

        ! Convert angles to radians:
        dir_around_node = dir_around_node * (3.14159265358979323846 / 180.0)

        DO T = 1, n_times
            ! reference node
            sx = SIN(dir_around_node(the_node_local,T) - mu_x)
            den_x = den_x + sx*sx

            DO N = 1, n_close_nodes
                sy = SIN(dir_around_node(N,T) - mu_y(N))
                num(N)   = num(N)   + sx*sy
                den_y(N) = den_y(N) + sy*sy

                IF (MOD(N,10000) == 0 .OR. N == n_close_nodes) THEN
                    WRITE(25,'("Progress: ",I0,"/",I0," (",F5.1,"%)")') &
                        N, n_close_nodes, 100.0*N/n_close_nodes
                    FLUSH(25)
                END IF

            END DO
        END DO

        status = NF90_CLOSE(fid)
        CALL nc_error(status, .TRUE., "Error closing netcdf file")

        DEALLOCATE( dir_around_node )

    END DO  

    CALL MPI_REDUCE(den_x, den_x_g, 1,             MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(num,   num_g,   n_close_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(den_y, den_y_g, n_close_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    CALL MPI_FINALIZE(mpi_err)


    IF ( mpi_rank .EQ. 0 ) THEN 
        WRITE(25,*) 'Calculating correlations...'
        ALLOCATE( correlations(n_close_nodes) )

        DO N = 1, n_close_nodes
            IF ((den_x_g <= 1.0e-6) .OR. (den_y_g(N) <= 1.0e-6)) THEN
                correlations(N) = 0.0
            ELSE
                correlations(N) = num_g(N) / SQRT(den_x_g * den_y_g(N))
            END IF
        END DO

        WRITE(25,*) "Saving lon/lat/correlation arrays to file..."
        WRITE(*,*) SIZE(nearest_schism_lons), SIZE(nearest_schism_idxs), SIZE(correlations)
        OPEN(UNIT=10, FILE=TRIM(file_out), STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED')
        DO J = 1, n_close_nodes
            WRITE(10,'(F12.6,1X,F12.6,1X,I10,1X,F12.6)') nearest_schism_lons(J), nearest_schism_lats(J), &
                                                         nearest_schism_idxs(J), correlations(J)
        END DO
        CLOSE(10)

        DEALLOCATE( correlations )
        DEALLOCATE( nearest_schism_lons )
        DEALLOCATE( nearest_schism_lats )
        DEALLOCATE( nearest_schism_idxs )
        DEALLOCATE( close_node_indices )
        DEALLOCATE( num )
        DEALLOCATE( den_y )
        DEALLOCATE( num_g )
        DEALLOCATE( den_y_g )

    END IF 
    
    IF ( mpi_rank .EQ. 0 ) WRITE(25,*) 'Done.'

    CLOSE(UNIT=25)
    
END PROGRAM spatial_corr


SUBROUTINE nc_error(iret, lstop, func)

  INTEGER, INTENT(in) :: iret
  LOGICAL, INTENT(in) :: lstop
  CHARACTER(LEN=*), INTENT(in) :: func
  !
  IF (iret .NE. 0) THEN
        WRITE(*,*) 'ROUTINE: ', TRIM(func)
        WRITE(*,*) 'ERROR: ', iret
        WRITE(*,*) 'Message: See NetCDF documentation for error code details.'
        IF (lstop) STOP
  ENDIF

END SUBROUTINE nc_error

SUBROUTINE deg2rad(degree, rad)
      REAL, INTENT(IN) :: degree
      REAL, PARAMETER  :: deg_to_rad = atan(1.0)/45
      REAL             :: rad

      rad = degree*deg_to_rad
END SUBROUTINE deg2rad

SUBROUTINE haversine(deglat1, deglon1, deglat2, deglon2, dist)
      ! Great circle distance calculator
      REAL, INTENT(IN) :: deglat1, deglon1, deglat2, deglon2
      REAL             :: h_a, h_c, h_dist, h_dlat, h_dlon, h_lat1, h_lat2
      REAL, PARAMETER  :: radius = 6372.8

      CALL deg2rad(deglat2-deglat1, h_dlat)
      CALL deg2rad(deglon2-deglon1, h_dlon)
      CALL deg2rad(deglat1, h_lat1)
      CALL deg2rad(deglat2, h_lat2)

      h_a = (sin(h_dlat/2))**2 + cos(h_lat1)*cos(h_lat2)*(sin(h_dlon/2))**2
      h_c = 2*asin(sqrt(h_a))
      dist = radius*h_c

END SUBROUTINE haversine
