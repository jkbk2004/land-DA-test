module rrtmgp_lw_aerosol_optics
  use machine,                   only: kind_phys
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_optical_props,          only: ty_optical_props_1scl
  use radiation_tools,           only: check_error_msg
  use rrtmgp_sw_gas_optics,      only: sw_gas_props
  use rrtmgp_lw_gas_optics,      only: lw_gas_props
  use module_radiation_aerosols, only: &
       NF_AESW,                  & ! Number of optical-fields in SW output (3=tau+g+omega)
       NF_AELW,                  & ! Number of optical-fields in LW output (3=tau+g+omega)
       setaer,                   & ! Routine to compute aerosol radiative properties (tau,g,omega)
       NSPC1                       ! Number of species for vertically integrated aerosol optical-depth
  use netcdf

  implicit none

  public rrtmgp_lw_aerosol_optics_init, rrtmgp_lw_aerosol_optics_run, rrtmgp_lw_aerosol_optics_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_aerosol_optics_init()
  ! #########################################################################################
  subroutine rrtmgp_lw_aerosol_optics_init()
  end subroutine rrtmgp_lw_aerosol_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_aerosol_optics_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_aerosol_optics_run
!! \htmlinclude rrtmgp_lw_aerosol_optics.html
!!
  subroutine rrtmgp_lw_aerosol_optics_run(doLWrad, nCol, nLev, nspc, nTracer, nTracerAer,   &
       p_lev, p_lay, p_lk, tv_lay, relhum, lsmask, tracer, aerfld, lon, lat,                &
       lw_optical_props_aerosol, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad                  ! Logical flag for longwave radiation call
    integer, intent(in) :: &
         nCol,                  & ! Number of horizontal grid points
         nLev,                  & ! Number of vertical layers
         nspc,                  & ! Number of aerosol optical-depths
         nTracer,               & ! Number of tracers
         nTracerAer               ! Number of aerosol tracers
    real(kind_phys), dimension(:), intent(in) :: &
         lon,                   & ! Longitude
         lat,                   & ! Latitude
         lsmask                   ! Land/sea/sea-ice mask
    real(kind_phys), dimension(:,:),intent(in) :: &
         p_lay,                 & ! Pressure @ layer-centers (Pa)
         tv_lay,                & ! Virtual-temperature @ layer-centers (K)
         relhum,                & ! Relative-humidity @ layer-centers
         p_lk                     ! Exner function @ layer-centers (1)
    real(kind_phys), dimension(:,:, :),intent(in) :: &
         tracer                   ! trace gas concentrations
    real(kind_phys), dimension(:,:, :),intent(in) :: &
         aerfld                   ! aerosol input concentrations
    real(kind_phys), dimension(:,:),intent(in) :: &
         p_lev                    ! Pressure @ layer-interfaces (Pa)

    ! Outputs
    type(ty_optical_props_1scl),intent(inout) :: &
         lw_optical_props_aerosol ! RRTMGP DDT: Longwave aerosol optical properties (tau)
    integer, intent(out) :: &
         errflg                   ! CCPP error flag
    character(len=*), intent(out) :: &
         errmsg                   ! CCPP error message

    ! Local variables
    real(kind_phys), dimension(nCol, nLev, lw_gas_props%get_nband(), NF_AELW) :: &
         aerosolslw            !
    real(kind_phys), dimension(nCol, nLev, sw_gas_props%get_nband(), NF_AESW) :: &
         aerosolssw
    real(kind_phys), dimension(nCol,nspc) :: aerodp
    integer :: iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile
    call setaer(p_lev/100., p_lay/100., p_lk, tv_lay, relhum, lsmask, tracer, aerfld, lon, lat, ncol, nLev, &
         nLev+1, .true., .true., aerosolssw, aerosolslw, aerodp)

    ! Copy aerosol optical information to RRTMGP DDT
    lw_optical_props_aerosol%tau = aerosolslw(:,:,:,1) * (1. - aerosolslw(:,:,:,2))

    lw_optical_props_aerosol%band_lims_wvn = lw_gas_props%get_band_lims_wavenumber()
    do iBand=1,lw_gas_props%get_nband()
       lw_optical_props_aerosol%band2gpt(1:2,iBand) = iBand
       lw_optical_props_aerosol%gpt2band(iBand)     = iBand
    end do

  end subroutine rrtmgp_lw_aerosol_optics_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_aerosol_optics_finalize()
  ! #########################################################################################
  subroutine rrtmgp_lw_aerosol_optics_finalize()
  end subroutine rrtmgp_lw_aerosol_optics_finalize
end module rrtmgp_lw_aerosol_optics
