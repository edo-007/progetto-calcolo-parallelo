PROGRAM CartesianGrid2D
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    ! -------------------------------------------------- !
    ! SETTING OF THE COMPUTATIONAL GRID
    INTEGER, PARAMETER :: NDIM=2            ! number of space dimensions
    INTEGER, PARAMETER :: IMAX=200          ! number of cells in x-dir
    INTEGER, PARAMETER :: JMAX=200          ! number of cells in y-dir
    REAL   , PARAMETER :: x0=0.0            ! min x-coord
    REAL   , PARAMETER :: x1=1.0            ! max x-coord
    REAL   , PARAMETER :: y0=0.0            ! min y-coord
    REAL   , PARAMETER :: y1=1.0            ! max y-coord
    ! -------------------------------------------------- !
    ! Local variable declaration
    TYPE tMPI
      INTEGER  :: myrank
      INTEGER  :: nCPU
      INTEGER  :: status(MPI_STATUS_SIZE)
      INTEGER  :: iStart, iEnd              ! idx of starting and ending cell in x-dir
      INTEGER  :: jStart, jEnd              ! idx of starting and ending cell in y-dir
      INTEGER  :: imax, jmax                ! number of cells within each rank
      INTEGER, ALLOCATABLE :: mycoords(:)   ! point coords of the subgrid
      INTEGER  :: iErr                      ! flag for errors in MPI functions    
      INTEGER  :: x_thread                  ! number of CPUs in x-dir
      INTEGER  :: y_thread                  ! number of CPUs in y-dir
    END TYPE tMPI
    TYPE(tMPI) :: MPI
    !
    LOGICAL, ALLOCATABLE :: periods(:)
    INTEGER, ALLOCATABLE :: dims(:)
    INTEGER              :: TCPU, BCPU, RCPU, LCPU  ! neighbor ranks of myrank
    INTEGER              :: i, j, idx, jdx, source
    INTEGER              :: COMM_CART               ! Cartesian MPI communicator
    REAL                 :: dx, dy
    REAL  , ALLOCATABLE  :: x(:), y(:)              ! grid coordinates
    ! -------------------------------------------------- !
    
    ! 1) MPI initialization
    
    CALL MPI_INIT(MPI%iErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPI%myrank,MPI%iErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,MPI%nCPU,MPI%iErr)
    ! check the number of CPUs
    IF(MOD(MPI%nCPU,2).NE.0) THEN
      PRINT *, ' ERROR. Number of CPU must be even!'
      CALL MPI_FINALIZE(MPI%iErr)
      STOP
    ENDIF    
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,MPI%iErr)
    
    ! 2) Create a Cartesian topology
    
    ! 2.1) check the number of cells
    IF(MOD(IMAX,2).NE.0) THEN
      PRINT *, ' ERROR. Number of x-cells must be even!'
      CALL MPI_FINALIZE(MPI%iErr)
      STOP
    ENDIF
    IF(MOD(JMAX,2).NE.0) THEN
      PRINT *, ' ERROR. Number of y-cells must be even!'
      CALL MPI_FINALIZE(MPI%iErr)
      STOP
    ENDIF
    
    ! 2.2) Domain decomposition
    MPI%x_thread = MPI%nCPU/2
    MPI%y_thread = MPI%nCPU - MPI%x_thread
    
    ALLOCATE( dims(NDIM), periods(NDIM), MPI%mycoords(NDIM) )
    
    dims    = (/ MPI%x_thread, MPI%y_thread /)
    periods = .FALSE.
    
    CALL MPI_CART_CREATE(MPI_COMM_WORLD,NDIM,dims,periods,.TRUE.,COMM_CART,MPI%iErr)
    CALL MPI_COMM_RANK(COMM_CART,MPI%myrank,MPI%iErr)
    
    ! 2.3) Find CPU neighbors
    CALL MPI_CART_SHIFT(COMM_CART,0, 1,source,RCPU,mpi%iErr)
    CALL MPI_CART_SHIFT(COMM_CART,0,-1,source,LCPU,mpi%iErr)
    CALL MPI_CART_SHIFT(COMM_CART,1, 1,source,TCPU,mpi%iErr)
    CALL MPI_CART_SHIFT(COMM_CART,1,-1,source,BCPU,mpi%iErr)
    
    ! "coordinates" (integer indexes!)
    CALL MPI_CART_COORDS(COMM_CART,MPI%myrank,NDIM,MPI%mycoords,MPI%iErr)
    MPI%imax = IMAX/MPI%x_thread
    MPI%jmax = JMAX/MPI%y_thread
    MPI%iStart = 1 + MPI%mycoords(1)*MPI%imax
    MPI%iEnd   = MPI%iStart + MPI%imax - 1 
    MPI%jStart = 1 + MPI%mycoords(2)*MPI%jmax
    MPI%jEnd   = MPI%jStart + MPI%jmax - 1
    
    ! 3) Compute the real mesh
    
    dx = (x1-x0)/REAL(IMAX-1)
    dy = (y1-y0)/REAL(JMAX-1)
    ALLOCATE( x(MPI%imax) )
    ALLOCATE( y(MPI%jmax) )
    
    idx = 0
    DO i = MPI%iStart, MPI%iEnd
        idx = idx + 1
        x(idx) = (x0-dx/2.) + (i-1)*dx
    ENDDO
    jdx = 0
    DO j = MPI%jStart, MPI%jEnd
        jdx = jdx + 1
        y(jdx) = (y0-dy/2.) + (j-1)*dy
    ENDDO
    
    ! 4) Plot the output and finalize the program
    
    CALL ASCII_Output(x,y,MPI%imax,MPI%jmax,MPI%myrank)
    CALL MPI_FINALIZE(MPI%iErr)
    
END PROGRAM CartesianGrid2D    