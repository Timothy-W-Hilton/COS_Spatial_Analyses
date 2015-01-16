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

  integer :: JULIAN

  real, allocatable, dimension(:,:,:):: GPP, LRU, cos_co2_ratio
  character*80 vdesc_arg
  integer LOGDEV, status, t_start, t_end, t_step

  LOGDEV = INIT3()

  call write_COS_CO2_ratio_ioapi(124, 124, 1.1)
  call write_LRU_ioapi_from_C4_pct('C4pct_INPUT', 'LRU_FILE', status)
  ! t_start = (2008 * 1000) + JULIAN(2008, 7, 1)
  ! t_end = (2008 * 1000) + JULIAN(2008, 9, 1)
  ! write(*,*) 'times: ', t_start, t_end
  ! ! t_start = 2008180
  ! ! t_end = 2008242
  ! t_step = 030000  !3 hours expressed as HHMMSS
  ! call calc_fcos_3D(t_start, t_end, t_step, 'GPP_INPUT', 'LRU_FILE', 'RATIO_FILE')

  status = 0  !successful completion
  call M3EXIT('python_test', sdate3d, stime3d, 'program completed', status)

END program calc_fCOS_C4pct

!============================================================

subroutine calc_fcos_3D(t_start, t_end, t_step, GPP_FILE, LRU_FILE, RATIO_FILE)

!t_start: starting date, YYYYDDD
!t_end: ending date, YYYYDDD
!t_step: time step, HHMMSS
!gpp_file: logical name of GPP I/O API file
!lru_file: logical name of LRU I/O API file
!ratio_file: logical name of COS/CO2 ratio I/O API file
!
! preconditions: GPP_FILE, LRU_FILE, and RATIO_FILE must exist and
! have identical horizontal dimensions.
!
! REFERENCES
  ! Campbell, J. E., Carmichael, G. R., Chai, T., Mena-Carrasco, M.,
  ! Tang, Y., Blake, D. R., Blake, N. J., Vay, S. A., Collatz,
  ! G. J., Baker, I., Berry, J. A., Montzka, S. A., Sweeney, C.,
  ! Schnoor, J. L., and Stanier, C. O.: Photosynthetic Control of
  ! Atmospheric Carbonyl Sulfide During the Growing Season, Science,
  ! 322, 1085–1088, doi:10.1126/science.1164015, 2008.
  
  USE ioapi_regrid_tools

  IMPLICIT NONE

  include 'PARMS3.EXT'	! i/o API
  include 'FDESC3.EXT'	! i/o API
  include 'IODECL3.EXT'	! i/o API

  integer, intent(in) :: t_start, t_end, t_step
  character(len=*), intent(in) :: GPP_FILE, LRU_FILE, RATIO_FILE

  integer :: currstep, currec
  integer :: nrows, ncols, ntimes, ierr, midnight, cdate, ctime, t
  integer :: this_yyyyddd, this_hhmmss, surface_layer
  real, dimension(:,:), allocatable :: this_t_GPP, this_t_LRU, this_t_ratio
  real, dimension(:,:,:,:), allocatable :: fCOS
  character*100 :: dim_err_msg

  ! determine horizontal directions
  dim_err_msg = 'LRU, GPP, and COS/CO2 ratio files must &
       & have same horizontal dimensions'
  call open3_and_desc3(GPP_FILE, FSREAD3, 'calc_fCOS_3D')
  ncols = NCOLS3D
  nrows = NROWS3D
  call open3_and_desc3(LRU_FILE, FSREAD3, 'calc_fCOS_3d')
  IF (.not.((ncols.eq.NCOLS3D) .and. (nrows.eq.NROWS3D))) THEN
     write(*,*) dim_err_msg
     stop
  ENDIF
  call open3_and_desc3(RATIO_FILE, FSREAD3, 'calc_fCOS_3D')
  IF (.not.((ncols.eq.NCOLS3D) .and. (nrows.eq.NROWS3D))) THEN
     write(*,*) dim_err_msg
     stop
  ENDIF
  write(*,*) 'dimensions: ', nrows, ncols

  ! determine number of timesteps
  midnight = 0  !midnight expressed as hour of day
  ntimes = currec(t_end, midnight, t_start, midnight, t_step, cdate, ctime)
  write(*,*) 'number of timesteps: ', ntimes

  allocate(this_t_GPP(ncols, nrows), STAT=ierr)
  allocate(this_t_LRU(ncols, nrows), STAT=ierr)
  allocate(this_t_ratio(ncols, nrows), STAT=ierr)
  allocate(fCOS(ntimes, ncols, nrows, 1), STAT=ierr)
  
  this_yyyyddd = t_start
  this_hhmmss = midnight
  surface_layer = 1
  DO t=1, ntimes
     ierr = READ3('GPP_FILE', 'GPP', surface_layer, &
          this_yyyyddd, this_hhmmss, this_t_GPP)
     ierr = READ3('LRU_FILE', 'LRU', surface_layer, &
          this_yyyyddd, this_hhmmss, this_t_LRU)
     ierr = READ3('RATIO_FILE', 'COS_CO2_ratio', surface_layer, &
          this_yyyyddd, this_hhmmss, this_t_ratio)
     ! equation 1, Campbell et al (2008)
     fCOS(t, surface_layer, :, :) = this_t_GPP * this_t_LRU * this_t_ratio
     call nextime(this_yyyyddd, this_hhmmss, t_step)
  ENDDO
END SUBROUTINE calc_fCOS_3D

!============================================================

SUBROUTINE write_COS_CO2_ratio_ioapi(nrows, ncols, const_val_arg)

! create an I/O API file containing a constant avlue for COS/CO2
! ratio.  This is a sort of placeholder to allow for easy substitution
! of a spatially and temporally varying ratio later. 
! nrows: number of rows in the array
! ncols: number of columns in the array
! const_val_arg: the constant value to be placed in the array

    USE ioapi_regrid_tools

    IMPLICIT NONE

    INCLUDE 'PARMS3.EXT'	! i/o API
    INCLUDE 'FDESC3.EXT'	! i/o API
    INCLUDE 'IODECL3.EXT'	! i/o API
    INCLUDE 'netcdf.inc'

    INTEGER, INTENT(in) :: nrows, ncols
    REAL, INTENT(in) :: const_val_arg
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ratio_vals
    INTEGER :: LOGDEV, ierr, jdate, jtime

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
    xcell3d=REAL(1)
    ycell3d=REAL(1)
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
    ALLOCATE(ratio_vals(nrows, ncols, 1), STAT=ierr)
    ratio_vals = const_val_arg

    PRINT*, 'writing COS/CO2 ratios'
    CALL open3_and_desc3('RATIO_FILE', FSCREA3, 'write_COS_CO2_ratio')
    IF (.NOT.write3('RATIO_FILE',vname3d(1), &
         jdate, jtime, ratio_vals(:,:,:))) THEN 
       WRITE(*,*) 'unable to write COS/CO2 ratio to file'
       STOP
    ENDIF

    PRINT *, '-----------------------------------------------'

  END SUBROUTINE write_COS_CO2_ratio_ioapi

!============================================================

SUBROUTINE write_LRU_ioapi_from_C4_pct(C4PCT_FNAME, LRU_FNAME, stat)
  
  !----------------------------------------------------------------------
  ! DESC: write I/O API file containing a weighted-average LRU value
  ! based on the C4 vegetation percentages of Still et al (2008) and
  ! the C3 and C4 LRU values of Stimler et al (2011).
  !----------------------------------------------------------------------
  ! INPUT:
  !     C4PCT_FNAME: string; logical name of the I/O API file
  !        containing C4 percentage data (see I/O API documentation
  !        for logical name discussion.  Must be 16 characters or
  !        less.
  !     LRU_FNAME: string; logical name of the LRU I/O API file to
  !        be written
  !----------------------------------------------------------------------
  ! OUTPUT:
  !     stat: integer.  0 on success, 1 on failure
  !----------------------------------------------------------------------
  ! REFERENCES:
  !
  ! Still, C.J. and J.A. Berry and G.J. Collatz and R.S. DeFries:
  !   ISLSCP II C4 Vegetation Percentage, in: ISLSCP Initiative II
  !   Collection, edited by Hall, Forrest G. and G. Collatz and
  !   B. Meeson and S. Los and E. Brown de Colstoun and D. Landis, Oak
  !   Ridge National Laboratory, Oak Ridge, Tennessee, U.S.A.,
  !   doi:10.3334/ORNLDAAC/932, available on-line
  !   [http://daac.ornl.gov/] from Oak Ridge National Laboratory
  !   Distributed Active Archive Center, 2009.
  !
  ! Stimler, K., Berry, J. A., Montzka, S. A., and Yakir, D.:
  !   Association between Carbonyl Sulfide Uptake and 18delta during
  !   Gas Exchange in C3 and C4 Leaves, Plant Physiology, 157,
  !   509–517, doi:10.1104/pp.111.176578, 2011.
  !----------------------------------------------------------------------
  ! HISTORY:
  ! By: Timothy W. Hilton
  ! On: 15 Jan 2015
  !----------------------------------------------------------------------

  USE ioapi_regrid_tools

  IMPLICIT NONE

  include 'PARMS3.EXT'	! i/o API
  include 'FDESC3.EXT'	! i/o API
  include 'IODECL3.EXT'	! i/o API

  character(len=*), intent(in) :: C4PCT_FNAME, LRU_FNAME
  integer, intent(out) :: stat
  integer :: i, j, nrows, ncols, ierr
  integer :: this_yyyyddd, this_hhmmss, sfc_layer, ntimes
  real :: C4_LRU_val, C3_LRU_val
  real, DIMENSION(:,:,:,:), allocatable :: LRU, C4frac
  call open3_and_desc3(C4PCT_FNAME, FSREAD3, 'write_LRU_ioapi')
  nrows = NROWS3D
  ncols = NCOLS3D

  stat = 1

  ! dummy time arguments.  Per the I/O API documentation, to be ignored
  ! by READ3 for time-independent data.
  this_yyyyddd = 0
  this_hhmmss = 0
  sfc_layer = 1
  ntimes = 1
  allocate(LRU(ntimes, sfc_layer, NCOLS3D, NROWS3D), STAT=ierr)
  allocate(C4frac(ntimes, sfc_layer, NCOLS3D, NROWS3D), STAT=ierr)
  ierr = read3(C4PCT_FNAME, 'C4pct', & 
       ALLAYS3, this_yyyyddd, this_hhmmss, C4frac)

  C4_LRU_val = 1.16 ! Stimler et al (2011)
  C3_LRU_val = 1.82 ! Stimler et al (2011)
  C4frac = C4frac / 100.0  !convert percentage to fraction
  LRU = (C4frac * C4_LRU_val) + ((1.0 - C4frac) * C3_LRU_val)

  VNAME3D(1) = "LRU"
  VDESC3D(1) = "LRU; weighted average of C3, C4 pct using &
       &Stimler et al (2011) LRU values"
  UNITS3D(1) = "normalized LRU"
  FDESC3D(1) = 'LRU calculated as weighted average using the C4 &
       &vegetation percentages of Still et al (2008) and the C3 &
       &and C4 LRU values of Stimler et al (2011).'
  call open3_and_desc3(LRU_FNAME, FSCREA3, 'LRU_from_C4pct')

  ! write new LRU values
  IF (.NOT.write3(LRU_FNAME,'LRU', &
       this_yyyyddd, this_hhmmss, LRU(:,:,:,:))) THEN 
     WRITE(*,*) 'unable to write LRU ratio to file'
     STOP
  ENDIF

  stat = 0 !success

end SUBROUTINE write_LRU_ioapi_from_C4_pct



