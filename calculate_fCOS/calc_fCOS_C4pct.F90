!short example of using the I/O API fortran library 
!
!performs these tasks:
!(1) open an existing I/O API file
!(2) read and display some brief information aboutthe variables it contains
!(3) create a new I/O API file with the same metadata as the original
!    except that the GPP variable is replaced by "dummy"
!
!Timothy W. Hilton, UC Merced, 15 January 2015

program calc_fCOS_C4pct

  USE ioapi_regrid_tools

  IMPLICIT NONE

  include 'PARMS3.EXT'	! i/o API
  include 'FDESC3.EXT'	! i/o API
  include 'IODECL3.EXT'	! i/o API
  include 'netcdf.inc'

  real, allocatable, dimension(:,:,:):: GPP, LRU, cos_co2_ratio
  character*80 vdesc_arg
  integer LOGDEV, status

  LOGDEV = INIT3()

  call open3_and_desc3('LRU_INPUT', FSREAD3, 'test_open')
  write(*,*) 'vars: ', VNAME3D(1), UNITS3D(1), VDESC3D(1)
  write(*,*) 'dimensions: ', NLAYS3D,NROWS3D,NCOLS3D

! ! test whether scripting with python, ecampbell300_data_paths will
! ! work, handle spaces in paths correctly -- YES!

  VNAME3D(1) = 'dummy'
  UNITS3D(1) = 'furlongs fortnight-1'
  VDESC3D(1) = 'this is a test :)'

  !call open3_and_desc3('OUTPUT', FSCREA3, 'python_test')
  call write_COS_CO2_ratio_ioapi(124, 124, 1.1)

  status = 0  !successful completion
  call M3EXIT('python_test', sdate3d, stime3d, 'program completed', status)

END program calc_fCOS_C4pct

!============================================================

subroutine calc_fcos(t_start, t_end, t_step, GPP_FILE, LRU_FILE, RATIO_FILE)

!t_start: starting date, YYYYDDD
!t_end: ending date, YYYYDDD
!t_step: time step, HHMMSS
!gpp_file: logical name of GPP I/O API file
!lru_file: logical name of LRU I/O API file
!ratio_file: logical name of COS/CO2 ratio I/O API file
!
! preconditions: GPP_FILE, LRU_FILE, and RATIO_FILE must exist and
! have identical horizontal dimensions.

  USE ioapi_regrid_tools

  IMPLICIT NONE

  include 'PARMS3.EXT'	! i/o API
  include 'FDESC3.EXT'	! i/o API
  include 'IODECL3.EXT'	! i/o API

  integer, intent(in) :: t_start, t_end, t_step
  character(len=*), intent(in) :: GPP_FILE, LRU_FILE, RATIO_FILE

  integer :: this_t, nrows, ncols
  real, dimension(:,:), allocatable :: this_t_GPP, this_t_LRU, this_t_ratio
  real, dimension(:,:,:), allocatable :: fCOS

  call open3_and_desc3('GPP_FILE', FSREAD3, 'calc_fCOS')
  ncols = NCOLS3D
  nrows = NROWS3D
  call open3_and_desc3('LRU_FILE', FSREAD3, 'calc_fCOS')
  IF (.not.((ncols.eq.NCOLS3D) .and. (nrows.eq.NROWS3D))) THEN
     write(*,*) 'LRU, GPP, and COS/CO2 ratio files must have same horizontal dimensions'
     stop
  ENDIF
  call open3_and_desc3('RATIO_FILE', FSREAD3, 'calc_fCOS')
  IF (.not.((ncols.eq.NCOLS3D) .and. (nrows.eq.NROWS3D))) THEN
     write(*,*) 'LRU, GPP, and COS/CO2 ratio files must have same horizontal dimensions'
     stop
  ENDIF

  write(*,*) 'dimensions: ', nrows, ncols

END SUBROUTINE calc_fcos

!============================================================

subroutine write_COS_CO2_ratio_ioapi(nrows, ncols, const_val_arg)

! create an I/O API file containing a constant avlue for COS/CO2
! ratio.  This is a sort of placeholder to allow for easy substitution
! of a spatially and temporally varying ratio later. 
! nrows: number of rows in the array
! ncols: number of columns in the array
! const_val_arg: the constant value to be placed in the array

    USE ioapi_regrid_tools

    IMPLICIT NONE

    include 'PARMS3.EXT'	! i/o API
    include 'FDESC3.EXT'	! i/o API
    include 'IODECL3.EXT'	! i/o API
    include 'netcdf.inc'

    integer, intent(in) :: nrows, ncols
    real, intent(in) :: const_val_arg
    real, allocatable, dimension(:,:,:) :: ratio_vals
    integer :: LOGDEV, ierr, jdate, jtime

    LOGDEV = INIT3()

    ! populate dimensions, date & time, file & variable descriptions, etc.
    FDESC3D(1) = "contains constant COS/CO2 ratio 1.1"
    nvars3d=1                 ! emission species number
    ftype3d=GRDDED3           ! file is in grided, Global dobson file header
    gdtyp3d=LATGRD3           ! lat-lon
    xcent3d=0.
    ycent3d=0.
    xorig3d=-180.
    yorig3d=-90.
    xcell3d=real(1)
    ycell3d=real(1)
    ncols3d=ncols
    nrows3d=nrows
    nlays3d=1
    vgtyp3d=VGSGPN3           !  non-hydrostatic sigma-p vertical coordinate
    vgtop3d=1.                ! domain top in meter
    vglvs3d(1)=1.             ! levels in meter
    vglvs3d(2)=0.             ! levels in meter
    !     do L=1,nvars3d
    units3d(1)='fraction'
    VNAME3D(1) = 'COS_CO2_ratio'
    VDESC3D(1) = 'COS/CO2 ratio'
    vtype3d=m3real
    !The ratios are time-independent.  from I/O API documentation,
    !section "Dates and Time Conventions": "by convention, time
    !increment dT=0 means that the data in the file is
    !time-independent, and routines like READ3() and WRITE3() which
    !deal with time-independent files ignore the date-and-time
    !arguments."
    sdate3d = 0
    stime3d = 0
    tstep3d = 0

    !--------------------------------------------------
    ! write the ratio value to the I/O API file
    !--------------------------------------------------
    !date/time arguments are ignored by write3 for time-independent data
    !(as per I/O API documentation, section "Dates and Time
    !Conventions").  Thus jdate and jtime are dummy variables here,
    !included to satisfy the compiler.
    jdate = 0 
    jtime = 0
    allocate(ratio_vals(nrows, ncols, 1), STAT=ierr)
    ratio_vals = const_val_arg

    print*, 'writing COS/CO2 ratios'
    call open3_and_desc3('RATIO_FILE', FSCREA3, 'write_COS_CO2_ratio')
    if (.not.write3('RATIO_FILE',vname3d(1), &
         jdate, jtime, ratio_vals(:,:,:))) then 
       write(*,*) 'unable to write COS/CO2 ratio to file'
       stop
    ENDIF

    print *, '-----------------------------------------------'

    return
  end subroutine write_COS_CO2_ratio_ioapi
