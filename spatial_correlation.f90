PROGRAM spatial_corr

    USE netcdf
    USE mpi
    IMPLICIT NONE

    INTEGER                        :: the_node, the_node_local
    CHARACTER(LEN=50)              :: node_str, output_var, logfile, line, param_name, &
                                      param_value
    CHARACTER(LEN=250),ALLOCATABLE :: files_in(:)
    REAL                           :: the_lon, the_lat, max_dist, sum_x, sum_x2,       &
                                      sum_x_g, sum_x2_g
    INTEGER                        :: start_year, end_year, n_years, C, sum_N, sum_N_g
    INTEGER                        :: fid, status, dimID_time, dimID_node, n_nodes,    &
                                      n_times, time_ID, nodex_ID, nodey_ID, idx,       &
                                      count_close_nodes, idx_node0, idx_node1, hs_ID,  &
                                      n_close_nodes, I, J, N, F, Y, M, n_times_total,  &
                                      num_files, num_files_per_node, extra_files,      &
                                      err=0, eq_pos, ex_pos, EXTCDE
    INTEGER, DIMENSION(2)          :: start_nodes, count_nodes
    REAL,ALLOCATABLE,DIMENSION(:)  :: schism_lons, schism_lats, distances, correlations, &
                                      close_node_indices, nearest_schism_lons,           &
                                      nearest_schism_lats, sum_y, sum_y2, sum_xy,        &
                                      sum_y_g, sum_y2_g, sum_xy_g
                                      
    INTEGER,ALLOCATABLE,DIMENSION(:) :: nearest_schism_idxs
    REAL,ALLOCATABLE               :: hs_around_node(:,:)
    CHARACTER(LEN=100)             :: file_out
        
    INTEGER                        :: mpi_err, mpi_rank, mpi_size, start_idx_mpi, &
                                      end_idx_mpi
    DOUBLE PRECISION               :: mpi_start_time, mpi_end_time, mpi_elapsed_time, & 
                                      mpi_max_time

    CHARACTER(LEN=250)             :: base_dir
    CHARACTER(LEN=4)               :: y4
    
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
    WRITE(*,*) 'nearest_schism_idxs(1:100) = ', nearest_schism_idxs(1:100)
    DEALLOCATE( schism_lons )
    DEALLOCATE( schism_lats )

    WRITE(*,*) 'SIZE(files_in) =', SIZE(files_in)

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

    ALLOCATE( sum_y(n_close_nodes) )
    ALLOCATE( sum_y2(n_close_nodes) )
    ALLOCATE( sum_xy(n_close_nodes) )
    sum_y = 0.0 ; sum_y2 = 0.0 ; sum_xy = 0.0

    IF ( mpi_rank .EQ. 0 ) THEN
        ALLOCATE( sum_y_g(n_close_nodes) )
        ALLOCATE( sum_y2_g(n_close_nodes) )
        ALLOCATE( sum_xy_g(n_close_nodes) )
    END IF
    
    sum_x = 0.0
    sum_x2 = 0.0
    sum_N = 0

    ! CALL MPI_INIT(mpi_err)
    ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)
    ! CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_err)

    mpi_start_time = MPI_WTIME()

    WRITE(logfile, "('output_',I0,'.log')") mpi_rank
    OPEN(UNIT=25, FILE=logfile, STATUS='replace', ACTION='write')

    IF ( mpi_rank .EQ. 0 ) WRITE(25,*) 'the_node = ', the_node
    IF ( mpi_rank .EQ. 0 ) WRITE(25,*) 'Loop through files, extract Hs at nearest nodes:'

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

        ! Read hs around closest nodes:
        ALLOCATE( hs_around_node(n_close_nodes, n_times) )

        status = NF90_INQ_VARID(fid, TRIM(output_var), hs_ID)
        CALL nc_error(status, .TRUE., "Error inquiring output_var ID")
        status = NF90_GET_VAR(fid, hs_ID, hs_around_node, start=start_nodes, count=count_nodes)
        CALL nc_error(status, .TRUE., "Error getting hs values")

        sum_N = sum_N + n_times
        sum_x = sum_x + SUM(hs_around_node(the_node_local,:))
        sum_x2 = sum_x2 + SUM(hs_around_node(the_node_local,:)**2)
        DO N = 1, n_close_nodes
            sum_y(N) = sum_y(N) + SUM(hs_around_node(N,:))
            sum_y2(N) = sum_y2(N) + SUM(hs_around_node(N,:)**2)
            sum_xy(N) = sum_xy(N) + SUM(hs_around_node(N,:)*hs_around_node(the_node_local,:))

            if (mod(N,1000) == 0 .or. N == n_close_nodes) then
                write(25,'("Progress: ",I0,"/",I0," (",F5.1,"%)")') &
                    N, n_close_nodes, 100.0*N/n_close_nodes
                flush(6)
            end if
        END DO 

        DEALLOCATE( hs_around_node )

    END DO  
    
    DEALLOCATE( files_in )

    mpi_end_time = MPI_WTIME()
    mpi_elapsed_time = mpi_end_time - mpi_start_time
    CALL MPI_REDUCE(mpi_elapsed_time, mpi_max_time, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)
    
    CALL MPI_REDUCE(sum_x,  sum_x_g,  1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_x2, sum_x2_g, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_N,  sum_N_g,  1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_y,  sum_y_g,  n_close_nodes, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_y2, sum_y2_g, n_close_nodes, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    CALL MPI_REDUCE(sum_xy, sum_xy_g, n_close_nodes, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

    IF (mpi_rank .EQ. 0) WRITE(25,*) '------------'
    IF (mpi_rank .EQ. 0) WRITE(25,*) 'Elapsed time = ', mpi_elapsed_time
    IF (mpi_rank .EQ. 0) WRITE(25,*) '------------'

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    CALL MPI_FINALIZE(mpi_err)

    IF ( mpi_rank .EQ. 0 ) THEN 
        WRITE(25,*) 'Calculating correlations...'
        ALLOCATE( correlations(n_close_nodes) )
        DO N = 1, n_close_nodes
            IF (sum_y_g(N) .LT. 0.1) THEN
                correlations(N) = 0
            ELSE
                correlations(N) = (sum_N_g * sum_xy_g(N) - sum_x_g * sum_y_g(N)) / &
                                SQRT((sum_N_g * sum_x2_g - sum_x_g**2) * (sum_N_g * sum_y2_g(N) - sum_y_g(N)**2))
            END IF
        END DO

        WRITE(25,*) "Saving lon/lat/correlation arrays to file..."
        WRITE(*,*) SIZE(nearest_schism_lons), SIZE(nearest_schism_idxs), SIZE(correlations)
        OPEN(UNIT=10, FILE=TRIM(file_out), STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED')
        DO J = 1, n_close_nodes
            !WRITE(10, '(*(F10.5))') nearest_schism_lons(J), nearest_schism_lats(J), close_node_indices(J), correlations(J)
            WRITE(10,'(F12.6,1X,F12.6,1X,I10,1X,F12.6)') nearest_schism_lons(J), nearest_schism_lats(J), &
                                                         nearest_schism_idxs(J), correlations(J)
        END DO
        CLOSE(10)

        DEALLOCATE( sum_y )
        DEALLOCATE( sum_y2 )
        DEALLOCATE( sum_xy )
        DEALLOCATE( correlations )
        DEALLOCATE( nearest_schism_lons )
        DEALLOCATE( nearest_schism_lats )
        DEALLOCATE( nearest_schism_idxs )
        DEALLOCATE( close_node_indices )

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
