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
!   SMOIS, respectively
! - kettle_fsoil: Kettle et al (2002) COS soil flux, in the variable cos
! - crop_pct: Cropland percentages of Ramankutty et al (2008) in the
!   variable crop_pct
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

END MODULE HELPER_ROUTINES
! ----------------------------------------------------------------------

PROGRAM meld_whelan_kettle_soils

  USE M3UTILIO
  USE HELPER_ROUTINES

  IMPLICIT NONE

  INTEGER :: i, j, i_stat, sdate, stime, tstep, nsteps

  CALL get_melded_time_info(sdate, stime, tstep, nsteps)
  PRINT *, 'sdate, stime, tstep, nsteps', sdate, stime, tstep, nsteps

END PROGRAM meld_whelan_kettle_soils
