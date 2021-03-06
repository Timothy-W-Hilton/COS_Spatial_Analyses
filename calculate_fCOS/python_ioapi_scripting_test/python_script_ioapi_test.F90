!short example of using the I/O API fortran library 
!
!performs these tasks:
!(1) open an existing I/O API file
!(2) read and display some brief information aboutthe variables it contains
!(3) create a new I/O API file with the same metadata as the original
!    except that the GPP variable is replaced by "dummy"
!
!Timothy W. Hilton, UC Merced, 15 January 2015

program python_script_ioapi_test

  USE ioapi_regrid_tools

  IMPLICIT NONE

  include 'PARMS3.EXT'	! i/o API
  include 'FDESC3.EXT'	! i/o API
  include 'IODECL3.EXT'	! i/o API
  include 'netcdf.inc'

  real, allocatable, dimension(:,:,:):: LRU, cos_co2_ratio
  character*80 vdesc_arg
  integer LOGDEV, status

  LOGDEV = INIT3()

  call open3_and_fdesc3('GPP_INPUT', FSREAD3, 'test_open')

! ! test whether scripting with python, ecampbell300_data_paths will
! ! work, handle spaces in paths correctly -- YES!

  VNAME3D(1) = 'dummy'
  UNITS3D(1) = 'furlongs fortnight-1'
  VDESC3D(1) = 'this is a test :)'

  call open3_and_fdesc3('OUTPUT', FSCREA3, 'python_test')

  status = 0  !successful completion
  call M3EXIT('python_test', sdate3d, stime3d, 'program completed', status)

END program python_script_ioapi_test

!============================================================

