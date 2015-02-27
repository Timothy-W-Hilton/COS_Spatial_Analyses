!----------------------------------------------------------------------
! meld_whelan_kettle_fsoils.F90
!----------------------------------------------------------------------
! DESC: calculates a "melded" COS soil flux (fsoil) by combining the
! equation of Whelan et al (2015) with the fluxes of Kettle et al
! (2002).  The fluxes are calculated as a weighted average using the
! Whelan et al (2015) model is used for croplands and the Kettle et
! al (2002) fluxes are used for non-croplands.  The weights are
! assigned from the cropland percentage data of Ramankutty et al
! (2008).  The Whelan et al (2015) model is driven using soil
! moisture (volumetric water content, VWC) and soil temperature
! (Tsoil) from WRF (D'Allura et al, 2011).
! ----------------------------------------------------------------------
! PRECONDITIONS: Models-3 I/O API files containing the input data
! mentioned above must exist in the following logical names:
!
! - WRF_sfc_met: must contain Tsoil and VWC in the variables TSOIL and
!     SMOIS, respectively
! - kettle_fsoil: Kettle et al (2002) COS soil flux, in the variable cos
! - crop_pct: Cropland percentages of Ramankutty et al (2008) in the
!     variable crop_pct
! - fsoil_out: the file in which to write the calculated Fsoil.  Will
!     be created; an existing file of that logical name will be deleted.
! ----------------------------------------------------------------------
! REFERENCES
!
! D’Allura, A., S. Kulkarni, G. R. Carmichael, S. Finardi,
! B. Adhikary, C. Wei, D.  Streets, Q. Zhang, R. B. Pierce,
! J. A. Al-Saadi, G. Diskin, and P. Wennberg (2011), Meteorological
! and air quality forecasting using the WRF–STEM model during the 2008
! ARCTAS field campaign, Atmospheric Environment, 45(38), 6901–6910,
! doi:10.1016/j.atmosenv.2011.02.073.
!
! Kettle, A. J., U. Kuhn, M. von Hobe, J. Kesselmeier, and
! M. O. Andreae (2002), Global budget of atmospheric carbonyl sulfide:
! Temporal and spatial variations of the dominant sources and sinks,
! Journal of Geophysical Research: Atmospheres, 107(D22), ACH
! 25–1–ACH25–16, doi:10.1029/2002JD002187.
!
! Ramankutty, N., A. T. Evan, C. Monfreda, and J. A. Foley (2008),
! Farming the planet: 1.  Geographic distribution of global
! agricultural lands in the year 2000, Global Biogeochemical Cycles,
! 22(1), n/a–n/a, doi:10.1029/2007GB002952.
!
! Whelan, M. et al (2015).  In preparation.   Title TBD.
! ----------------------------------------------------------------------
! HISTORY: By: Timothy W. Hilton On: 26 Feb 2015
! ----------------------------------------------------------------------

MODULE HELPER_ROUTINES

  USE M3UTILIO
  IMPLICIT NONE

  ! INCLUDE 'PARMS3.EXT'	! i/o API
  ! INCLUDE 'FDESC3.EXT'	! i/o API
  ! INCLUDE 'IODECL3.EXT'	! i/o API

CONTAINS

  SUBROUTINE LASTTIME( SDATE,STIME,TSTEP,NRECS, EDATE,ETIME )

    INTEGER, INTENT(IN   ) :: SDATE, STIME, TSTEP, NRECS
    INTEGER, INTENT(  OUT) :: EDATE, ETIME

    INTEGER, PARAMETER :: S365 = 365 * 24 * 60 * 60     !!  seconds
    INTEGER, PARAMETER :: T365 = 365 * 24 * 10000       !!  time as H*MMSS

    !!.......   LOCAL VARIABLES:  day and year components of date

    INTEGER     ISEC, STEP
    INTEGER*8   IREC, SECS

    !!  Normalized copies of the arguments:

    EDATE = SDATE
    ETIME = STIME
    CALL NEXTIME( EDATE, ETIME, 0 )

    IF ( TSTEP .EQ. 0 ) THEN   !  time-independent case:
       RETURN
    ELSE IF ( NRECS .LE. 1 ) THEN   !  at most 1 record
       RETURN
    END  IF

    IREC = NRECS - 1
    SECS = IREC * TIME2SEC( ABS( TSTEP ) )

    DO
       IF ( SECS .LT. S365 ) EXIT
       SECS = SECS - S365
       CALL NEXTIME( EDATE, ETIME, T365 )
    END DO

    ISEC = SECS
    STEP = SEC2TIME( ISEC )
    CALL NEXTIME( EDATE, ETIME, STEP )

    RETURN

  END SUBROUTINE LASTTIME


  SUBROUTINE select_tstamp(sdate1, stime1, sdate2, stime2, &
       & sdate_out, stime_out, switch)
    ! ----------------------------------------------------------------------
    ! DESC: Compare two Models-3 I/O API timestamps, returning either
    ! the earlier or later as requested by the caller.
    ! ----------------------------------------------------------------------
    ! INPUTS
    !    sdate1, stime1, sdate2, stime2: integers; two Models-3 I/O
    !       API timestamps
    !    switch: string; either "earlier" or "later"
    ! OUTPUTS
    !    sdate_out, stime_out: integers; the earlier or later of the
    !      two inputs, as etermined by switch, in Models-3 I/O API
    !      timestamp format
    ! ----------------------------------------------------------------------
    ! author: Timothy W. Hilton, thilton@ucmerced.edu, Feb 2015
    ! ----------------------------------------------------------------------
    USE M3UTILIO
    IMPLICIT NONE

    INTEGER, INTENT(in) :: sdate1, stime1, sdate2, stime2
    CHARACTER(*), INTENT(in) :: switch
    INTEGER, INTENT(out) :: sdate_out, stime_out
    LOGICAL :: t1_earlier



    IF (.NOT.((switch .EQ. 'earlier') .OR. (switch .EQ. 'later'))) THEN
       STOP
    END IF

    sdate_out = sdate1
    stime_out = stime1

    t1_earlier = SECSDIFF(sdate1, stime1, sdate2, stime2) .GT. 0
    IF ((t1_earlier) .AND. (switch .EQ. 'later')) THEN
       sdate_out = sdate2
       stime_out = stime2
    ELSEIF (.NOT.(t1_earlier) .AND. (switch .EQ. 'earlier')) THEN
       sdate_out = sdate2
       stime_out = stime2
    END IF
  END SUBROUTINE select_tstamp

  SUBROUTINE get_melded_time_info(sdate, stime, tstep, nsteps)
    ! ----------------------------------------------------------------------
    ! DESC: The "melded" file should run from the period its two
    ! inputs overlap and use the shorter of the two timesteps.  Among
    ! WRF_sfc_met and kettle_fsoil this routine returns the later
    ! start timestamp, the earlier end timestamp, and the shorter
    ! timestep.
    ! ----------------------------------------------------------------------
    ! INPUTS
    !    none
    ! OUTPUTS
    !    sdate, stime, tstep, nsteps: integers; the start date, start
    !    time, time step, and number of time steps that the melded
    !    file will contain.
    ! ----------------------------------------------------------------------
    ! author: Timothy W. Hilton, thilton@ucmerced.edu, Feb 2015
    ! ----------------------------------------------------------------------
    USE M3UTILIO
    USE ioapi_regrid_tools
    IMPLICIT NONE

    INTEGER, INTENT(out) :: sdate, stime, tstep, nsteps
    INTEGER :: stime_w, sdate_w, tstep_w, nsteps_w, &
         & stime_k, sdate_k, tstep_k, nsteps_k, &
         & last_date_w, last_time_w, &
         & last_date_k, last_time_k, &
         & last_date, last_time, &
         & dummy1, dummy2

    CALL open3_and_desc3('WRF_sfc_met', FSREAD3, 'meld_fsoil')
    stime_w = STIME3D
    sdate_w = SDATE3D
    tstep_w = TSTEP3D
    nsteps_w = MXREC3D
    CALL open3_and_desc3('kettle_fsoil', FSREAD3, 'meld_fsoil')
    stime_k = STIME3D
    sdate_k = SDATE3D
    tstep_k = TSTEP3D
    nsteps_k = MXREC3D
    ! determine which *first* time stamp is later
    CALL select_tstamp(sdate_w, stime_w, sdate_k, &
         & stime_k, sdate, stime, "later")
    ! determine which *last* time stamp is earlier
    CALL LASTTIME(sdate_w, stime_w, tstep_w, nsteps_w, &
         & last_date_w, last_time_w)
    CALL LASTTIME(sdate_k, stime_k, tstep_k, nsteps_k, &
         & last_date_k, last_time_k)
    CALL select_tstamp(last_date_w, last_time_w, last_date_k, last_time_k, &
         & last_date, last_time, "earlier")

    ! determine which time step is smaller
    tstep = MIN(tstep_w, tstep_k)
    nsteps = currec(last_date, last_time, sdate, stime, tstep, &
         & dummy1, dummy2)

  END SUBROUTINE get_melded_time_info

  ! ----------------------------------------------------------------------
  SUBROUTINE calc_whelan_fsoil(sdate_req, stime_req, fsoil)
    ! ----------------------------------------------------------------------
    ! DESC: calculate the Whelan soil flux for a caller-specified
    ! timestamp from WRF soil temperature and soil moisture using the
    ! equation in Mary's email of 5 Feb 2015.
    ! ----------------------------------------------------------------------
    ! INPUTS
    !    sdate_req, stime_req; integers: the requested timestamp in
    !       models-3 I/O API timestamp format.
    ! OUTPUTS
    !    fsoil; real(:, :, :, :): the 2-D soil flux for the requested timestep.
    ! ----------------------------------------------------------------------
    ! author: Timothy W. Hilton, thilton@ucmerced.edu, Feb 2015
    ! ----------------------------------------------------------------------
    USE ioapi_regrid_tools
    USE M3UTILIO
    IMPLICIT NONE

    INTEGER, INTENT(in) :: sdate_req, stime_req
    REAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(out) :: fsoil
    INTEGER :: sdate, stime, layer, ierr
    REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: Tsoil, VWC

    layer = 1  ! these are 2-D variables, so only need layer 1


    CALL read_STEM_var(sdate_req, stime_req, 'WRF_sfc_met', 'TSOIL', Tsoil)
    CALL read_STEM_var(sdate_req, stime_req, 'WRF_sfc_met', 'SMOIS', VWC)

    ALLOCATE(fsoil(1, 1, SIZE(Tsoil, 3), SIZE(Tsoil, 4)), STAT=ierr)

    ! Mary's linear Fsoil model from her email of 5 Feb 2015
    fsoil(1, 1, :, :) = (VWC(1, 1, :, :) * -28.77873448) + &
         & (Tsoil(1, 1, :, :) * 0.88867741) - 252.76497309

    DEALLOCATE(Tsoil)
    DEALLOCATE(VWC)

  END SUBROUTINE calc_whelan_fsoil

  SUBROUTINE read_STEM_var(sdate_req, stime_req, fname, vname, arr)
    ! ----------------------------------------------------------------------
    ! DESC: read the Kettle soil flux for a caller-specified
    ! timestamp
    ! ----------------------------------------------------------------------
    ! INPUTS
    !    sdate_req, stime_req; integers: the requested timestamp in
    !       models-3 I/O API timestamp format.
    !    fname; character: logical name of the I/O API file
    !    vname; character: name of the variable to be read
    ! OUTPUTS
    !    arr; real(:, :, :, :): the values for the requested timestep.
    ! ----------------------------------------------------------------------
    ! author: Timothy W. Hilton, thilton@ucmerced.edu, Feb 2015
    ! ----------------------------------------------------------------------
    USE ioapi_regrid_tools
    USE M3UTILIO
    IMPLICIT NONE

    INTEGER, INTENT(in) :: sdate_req, stime_req
    CHARACTER(*), INTENT(in) :: fname, vname
    REAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(out) :: arr
    INTEGER :: sdate, stime, layer, ierr

    ALLOCATE(arr(1, 1, NROWS3D, NCOLS3D), STAT=ierr)

    CALL open3_and_desc3(fname, FSREAD3, vname)
    ierr = CURRSTEP(sdate_req, stime_req, SDATE3D, STIME3D, TSTEP3D, &
         & sdate, stime)
    layer = 1
    ierr = READ3(fname, vname, layer, sdate, stime, arr)

  END SUBROUTINE read_STEM_var

  ! ----------------------------------------------------------------------

  SUBROUTINE meld_fsoils(sdate, stime, tstep, nsteps, fsoil)
    ! ----------------------------------------------------------------------
    ! DESC: calculate a COS soil flux from a weighted average of
    ! Mary's soil flux (cropland areas) and Kettle's soil flux
    ! (non-cropland areas).  Weights are calculated from Ramankutty's
    ! cropland percentage dataset.
    ! ----------------------------------------------------------------------
    ! INPUTS
    ! sdate_req, stime_req; integers: the requested timestamp
    !   in models-3 I/O API timestamp format.
    ! tstep; integer: time step in models-3 I/O API format (HHMMSS)
    ! nsteps; integer: number of timesteps to calculate
    ! OUTPUTS
    ! fsoil; real(:, :, :, :): allocatable array into which hte
    !   calculated soil fluxes will be placed
    ! ----------------------------------------------------------------------
    ! author: Timothy W. Hilton, thilton@ucmerced.edu, Feb 2015
    ! ----------------------------------------------------------------------
    USE M3UTILIO
    IMPLICIT NONE

    INTEGER, INTENT(in) :: sdate, stime, tstep, nsteps
    REAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(out) :: fsoil
    INTEGER :: this_date, this_time, nrows, ncols, ierr, i
    REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: fsoil_w, fsoil_k, pct

    CALL read_STEM_var(sdate, stime, 'crop_pct', 'crop_pct', pct)
    nrows = SIZE(pct, 3)
    ncols = SIZE(pct, 4)
    PRINT *, 'allocating fsoil', nsteps, 1, nrows, ncols
    ALLOCATE(fsoil(nsteps, 1, nrows, ncols), STAT=ierr)

    this_date = sdate
    this_time = stime

    DO i = 1, nsteps
       CALL calc_whelan_fsoil(this_date, this_time, fsoil_w)
       CALL read_STEM_var(this_date, this_time, 'kettle_fsoil', 'cos', fsoil_k)

       fsoil(i, 1, :, :) = (pct(1, 1, :, :) * fsoil_w(1, 1, :, :)) + &
            & (1-pct(1, 1, :, :) * fsoil_k(1, 1, :, :))

       CALL NEXTIME( this_date, this_time, tstep)
    ENDDO

    DEALLOCATE(pct)

  END SUBROUTINE meld_fsoils

! ----------------------------------------------------------------------
  SUBROUTINE write_fsoil_to_ioapi(fsoil, sdate_arg, stime_arg, tstep, nsteps, desc)
    ! ----------------------------------------------------------------------
    ! DESC: write soil flux to a Models-3 I/O API file with a
    ! user-specified description string
    ! ----------------------------------------------------------------------
    ! INPUTS
    ! sdate_arg, stime_arg; integers: the starting timestamp in
    !   models-3 I/O API timestamp format.
    ! tstep; integer: time step in models-3 I/O API format (HHMMSS)
    ! nsteps; integer: number of timesteps in fsoil
    ! fsoil; real(nstep, 1, nrows, ncols): array containing the soil
    !   fluxes to be written
    ! desc; character string: description of the file contents to be
    !   placed in the Models-3 I/O API FDESC3D field.  80-character
    !   maximum
    ! ----------------------------------------------------------------------
    ! author: Timothy W. Hilton, thilton@ucmerced.edu, Feb 2015
    ! ----------------------------------------------------------------------
    USE M3UTILIO
    IMPLICIT NONE

    REAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(in) :: fsoil
    INTEGER, INTENT(in) :: sdate_arg, stime_arg, tstep, nsteps
    CHARACTER(*), INTENT(in) :: desc
    INTEGER :: ierr, i, this_date, this_time

    ierr = DSCGRID( 'ARCNAGRID', GDNAM3D, &
         GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, &
         XORIG3D, YORIG3D, XCELL3D, YCELL3D, NCOLS3D, NROWS3D, NTHIK3D )
    FDESC3D(1) = desc
    nvars3d=1                 ! emission species number
    ftype3d=GRDDED3           ! file is in gridded
    nlays3d=1
    vgtyp3d=VGSGPN3           ! non-hydrostatic sigma-p vertical coordinate
    vgtop3d=1.                ! domain top in meter
    vglvs3d(1)=1.             ! levels in meter
    vglvs3d(2)=0.             ! levels in meter
    units3d(1)='pmol COS m-2 s-1'
    vname3d(1) = 'fsoil'
    VDESC3D(1) = 'hybrid Whelan-Kettle COS Fsoil'
    vtype3d(1)=m3real
    SDATE3D = sdate_arg
    STIME3D = stime_arg
    TSTEP3D = tstep

    ierr = OPEN3('fsoil_out', FSNEW3, 'meld_fsoils')

    this_date = sdate_arg
    this_time = stime_arg
    DO i=1, nsteps
       ierr = write3('fsoil_out', 'fsoil', this_date, this_time, &
            & fsoil(i, 1, :, :))
       CALL NEXTIME(this_date, this_time, tstep)
    ENDDO
  END SUBROUTINE write_fsoil_to_ioapi


END MODULE HELPER_ROUTINES
! ----------------------------------------------------------------------

PROGRAM meld_whelan_kettle_soils
  ! ----------------------------------------------------------------------
  ! main program.  See description at the top of this file.
  ! ----------------------------------------------------------------------
  USE M3UTILIO
  USE HELPER_ROUTINES

  IMPLICIT NONE

  INTEGER :: i, j, i_stat, sdate_out, stime_out, tstep, nsteps, &
       & date_end, time_end
  REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: fsoil
  CHARACTER(len=200) :: desc

  CALL get_melded_time_info(sdate_out, stime_out, tstep, nsteps)
  CALL meld_fsoils(sdate_out, stime_out, tstep, nsteps, fsoil)
  desc = "hybrid soil COS flux from applying Whelan model (of 5 Feb 2015 email) to croplands, Kettle soil COS flux elsewhere"
  CALL write_fsoil_to_ioapi(fsoil, sdate_out, stime_out, tstep, nsteps, desc)
  DEALLOCATE(fsoil)

  i_stat = 0  ! success
  CALL LASTTIME(sdate_out, stime_out, tstep, nsteps, date_end, time_end)
  CALL M3EXIT('meld_fsoils', date_end, time_end, 'program completed', i_stat)
END PROGRAM meld_whelan_kettle_soils
