module rrtmgp_sw_cloud_sampling
  use machine,                  only: kind_phys
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_optical_props,         only: ty_optical_props_2str
  use rrtmgp_sampling,          only: sampled_mask, draw_samples
  use mersenne_twister,         only: random_setseed, random_number, random_stat
  use radiation_tools,          only: check_error_msg
  use rrtmgp_sw_gas_optics,     only: sw_gas_props
  use netcdf

  implicit none

contains

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_sampling_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_sampling_run
!! \htmlinclude rrtmgp_sw_cloud_sampling.html
!!
  subroutine rrtmgp_sw_cloud_sampling_run(doSWrad, nCol, nDay, nLev, idxday, iovr,          &
       iovr_max, iovr_maxrand, iovr_rand, iovr_dcorr, iovr_exp, iovr_exprand, isubc_sw,     &
       icseed_sw, cld_frac, precip_frac, cloud_overlap_param, precip_overlap_param,         &
       sw_optical_props_cloudsByBand, sw_optical_props_precipByBand,                        &
       sw_optical_props_clouds, sw_optical_props_precip, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doSWrad                            ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                            & ! Number of horizontal gridpoints
         nDay,                            & ! Number of daylit points.
         nLev,                            & ! Number of vertical layers
         iovr,                            & ! Choice of cloud-overlap method
         iovr_max,                        & ! Flag for maximum cloud overlap method
         iovr_maxrand,                    & ! Flag for maximum-random cloud overlap method
         iovr_rand,                       & ! Flag for random cloud overlap method
         iovr_dcorr,                      & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,                        & ! Flag for exponential cloud overlap method
         iovr_exprand,                    & ! Flag for exponential-random cloud overlap method
         isubc_sw
    integer,intent(in),dimension(ncol) :: &
         idxday                             ! Indices for daylit points.
    integer,intent(in),dimension(ncol) :: &
         icseed_sw                          ! auxiliary special cloud related array when module 
                                            ! variable isubc_sw=2, it provides permutation seed 
                                            ! for each column profile that are used for generating 
                                            ! random numbers. when isubc_sw /=2, it will not be used.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac,                        & ! Total cloud fraction by layer
         precip_frac                        ! Precipitation fraction by layer
    real(kind_phys), dimension(ncol,nLev), intent(in)  :: &
         cloud_overlap_param,             & ! Cloud overlap parameter
         precip_overlap_param               ! Precipitation overlap parameter
    type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_cloudsByBand,   & ! RRTMGP DDT: Shortwave optical properties in each band (clouds)
         sw_optical_props_precipByBand      ! RRTMGP DDT: Shortwave optical properties in each band (precipitation)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                             ! Error message
    integer,          intent(out) :: &
         errflg                             ! Error flag
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_clouds,         & ! RRTMGP DDT: Shortwave optical properties at each spectral point (clouds) 
         sw_optical_props_precip            ! RRTMGP DDT: Shortwave optical properties at each spectral point (precipitation) 

    ! Local variables
    integer :: iday,iLay,iGpt
    integer,dimension(nday) :: ipseed_sw
    type(random_stat) :: rng_stat
    real(kind_phys) :: tauloc,asyloc,ssaloc
    real(kind_phys), dimension(sw_gas_props%get_ngpt(),nLev,nday) :: rng3D,rng3D2
    real(kind_phys), dimension(sw_gas_props%get_ngpt()*nLev) :: rng2D
    real(kind_phys), dimension(sw_gas_props%get_ngpt()) :: rng1D
    logical, dimension(nday,nLev,sw_gas_props%get_ngpt()) :: cldfracMCICA,precipfracSAMP

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. doSWrad) return
    if (nDay .gt. 0) then
       ! #################################################################################
       ! First sample the clouds...
       ! #################################################################################

       ! Allocate space RRTMGP DDTs [nday,nLev,nGpt]
       call check_error_msg('rrtmgp_sw_cloud_sampling_run', & 
            sw_optical_props_clouds%alloc_2str(nday, nLev, sw_gas_props))
 
       ! Change random number seed value for each radiation invocation (isubc_sw =1 or 2).
       if(isubc_sw == 1) then      ! advance prescribed permutation seed
          do iday = 1, nday
             ipseed_sw(iday) = sw_gas_props%get_ngpt() + iday
          enddo
       elseif (isubc_sw == 2) then ! use input array of permutaion seeds
          do iday = 1, nday
             ipseed_sw(iday) = icseed_sw(idxday(iday))
          enddo
       endif

       ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
       ! and layers. ([nGpts,nLev,nDayumn]-> [nGpts*nLev]*nDayumn)
       do iday=1,nday
          call random_setseed(ipseed_sw(iday),rng_stat)
          ! Use same rng for each layer
          if (iovr == iovr_max) then
             call random_number(rng1D,rng_stat)
             do iLay=1,nLev
                rng3D(:,iLay,iday) = rng1D
             enddo
          else
             do iLay=1,nLev
                call random_number(rng1D,rng_stat)
                rng3D(:,iLay,iday) = rng1D
             enddo
          endif
       enddo

       ! Cloud overlap.
       ! Maximum-random, random, or maximum cloud overlap
       if (iovr == iovr_maxrand .or. iovr == iovr_max .or. iovr == iovr_rand) then
          call sampled_mask(rng3D, cld_frac(idxday(1:nDay),:), cldfracMCICA)  
       endif
       ! Decorrelation-length overlap
       if (iovr == iovr_dcorr) then
          do iday=1,nday
             call random_setseed(ipseed_sw(iday),rng_stat)
             call random_number(rng2D,rng_stat)
             rng3D2(:,:,iday) = reshape(source = rng2D,shape=[sw_gas_props%get_ngpt(),nLev])
          enddo
          call sampled_mask(rng3D, cld_frac(idxday(1:nDay),:), cldfracMCICA,             &
	                        overlap_param = cloud_overlap_param(idxday(1:nDay),1:nLev-1),&
	                        randoms2      = rng3D2)
       endif 
       ! Exponential or exponential-random cloud overlap
       if (iovr == iovr_exp .or. iovr == iovr_exprand) then
          call sampled_mask(rng3D, cld_frac(idxday(1:nDay),:), cldfracMCICA, &
                            overlap_param = cloud_overlap_param(idxday(1:nDay),1:nLev-1))
       endif

       !
       ! Sampling. Map band optical depth to each g-point using McICA
       !
       call check_error_msg('rrtmgp_sw_cloud_sampling_run_draw_samples', & 
            draw_samples(cldfracMCICA, .true.,              &
                         sw_optical_props_cloudsByBand,     &
                         sw_optical_props_clouds))
         
       ! #################################################################################       
       ! Next sample precipitation (same as clouds for now)
       ! #################################################################################

       ! Allocate space RRTMGP DDTs [nday,nLev,nGpt]
       call check_error_msg('rrtmgp_sw_cloud_sampling_run', &
           sw_optical_props_precip%alloc_2str( nday, nLev, sw_gas_props))

       ! Precipitation overlap
       ! Maximum-random, random or maximum precipitation overlap
       if (iovr == iovr_maxrand .or. iovr == iovr_max .or. iovr == iovr_rand) then
          call sampled_mask(rng3D, precip_frac(idxday(1:nDay),:), precipfracSAMP)       
       endif
       ! Exponential decorrelation length overlap
       if (iovr == iovr_dcorr) then
          call sampled_mask(rng3D, precip_frac(idxday(1:nDay),:), precipfracSAMP,         & 
                            overlap_param = precip_overlap_param(idxday(1:nDay),1:nLev-1),& 
                            randoms2 = rng3D2)
       endif
       if (iovr == iovr_exp .or. iovr == iovr_exprand) then
          call sampled_mask(rng3D, precip_frac(idxday(1:nDay),:),precipfracSAMP, &
                            overlap_param = precip_overlap_param(idxday(1:nDay),1:nLev-1))
       endif
 
       !
       ! Sampling. Map band optical depth to each g-point using McICA
       !
       call check_error_msg('rrtmgp_sw_precip_sampling_run_draw_samples', & 
            draw_samples(precipfracSAMP, .true.,            &
                         sw_optical_props_precipByBand,     &
                         sw_optical_props_precip))                  
         
       ! #################################################################################        
       ! Just add precipitation optics to cloud-optics
       ! #################################################################################        
       do iGpt=1,sw_gas_props%get_ngpt()
          do iday=1,nDay
             do iLay=1,nLev
                tauloc = sw_optical_props_clouds%tau(iday,iLay,iGpt) + &
                         sw_optical_props_precip%tau(iday,iLay,iGpt)
                if (sw_optical_props_precip%tau(iday,iLay,iGpt) > 0) then
                   ssaloc = (sw_optical_props_clouds%tau(iday,iLay,iGpt)  * &
                             sw_optical_props_clouds%ssa(iday,iLay,iGpt)  + &
                             sw_optical_props_precip%tau(iday,iLay,iGpt)  * &
                             sw_optical_props_precip%ssa(iday,iLay,iGpt)) / &
                             tauloc
                   if (ssaloc > 0) then
                      asyloc = (sw_optical_props_clouds%tau(iday,iLay,iGpt) * &
                                sw_optical_props_clouds%ssa(iday,iLay,iGpt) * &
                                sw_optical_props_clouds%g(iday,iLay,iGpt)   + &
                                sw_optical_props_precip%tau(iday,iLay,iGpt) * &
                                sw_optical_props_precip%ssa(iday,iLay,iGpt) * &
                                sw_optical_props_precip%g(iday,iLay,iGpt))  / &
                                (tauloc*ssaloc)
                   else
                      tauloc = sw_optical_props_clouds%tau(iday,iLay,iGpt) 
                      ssaloc = sw_optical_props_clouds%ssa(iday,iLay,iGpt)
                      asyloc = sw_optical_props_clouds%g(iday,iLay,iGpt)            
                   endif
                   sw_optical_props_clouds%tau(iday,iLay,iGpt) = tauloc	
                   sw_optical_props_clouds%ssa(iday,iLay,iGpt) = ssaloc   
                   sw_optical_props_clouds%g(iday,iLay,iGpt)   = asyloc
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine rrtmgp_sw_cloud_sampling_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_sampling_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_sw_cloud_sampling_finalize()
  end subroutine rrtmgp_sw_cloud_sampling_finalize 

end module rrtmgp_sw_cloud_sampling
