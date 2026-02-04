!> filepath: src/simulation/m_thermal_src.fpp
!! @file m_thermal_src.fpp
!! @brief Contains module m_thermal_src

#:include 'macros.fpp'

!> @brief The module contains subroutines for thermal source term calculations
module m_thermal_src

    use m_derived_types        !< Definitions of the derived types
    use m_global_parameters    !< Definitions of the global parameters
    use m_nvtx                 !< For profiling
    use m_mpi_proxy

    implicit none
    private; public :: s_initialize_thermal_src, &
                       s_compute_thermal_src

    ! Module-level variables

    integer, allocatable, dimension(:) :: frequency    !< Laser frequency
    $:GPU_DECLARE(create='[frequency]')

    real(wp), allocatable, dimension(:) :: E_pulse, radius, pulseDur, laserDur, TimeStart
    $:GPU_DECLARE(create='[E_pulse,radius,pulseDur,laserDur,TimeStart]')

    real(wp), allocatable, target, dimension(:,:) :: thermal_loc
    $:GPU_DECLARE(create='[thermal_loc]')

    ! Precompute Pulse Parameters
    real(wp), allocatable, dimension(:) :: pulse_period, sigma_t, spatial_int, temporal_int, A0
    $:GPU_DECLARE(create='[pulse_period,sigma_t,spatial_int,temporal_int,A0]')

    ! Time Keeping
    logical, allocatable, dimension(:) :: laserStart
    $:GPU_DECLARE(create='[laserStart]')

contains

    !> This subroutine initializes the thermal source module
    impure subroutine s_initialize_thermal_src

        integer :: i, j   !< Generic loop variables

        @:ALLOCATE(thermal_loc(1:num_source_th,1:3), E_pulse(1:num_source_th), radius(1:num_source_th), pulseDur(1:num_source_th), laserDur(1:num_source_th))
        @:ALLOCATE(pulse_period(1:num_source_th),sigma_t(1:num_source_th),spatial_int(1:num_source_th),temporal_int(1:num_source_th),A0(1:num_source_th),laserStart(1:num_source_th),TimeStart(1:num_source_th))
        @:ALLOCATE(frequency(1:num_source_th))

        do i = 1, num_source_th
            do j = 1, 3
                thermal_loc(i, j) = thermal_s(i)%loc(j)
            end do
            
            ! ------ INPUTS ------- !
            frequency(i)    = thermal_s(i)%frequency
            E_pulse(i)      = thermal_s(i)%E_pulse
            radius(i)       = thermal_s(i)%radius
            pulseDur(i)     = thermal_s(i)%pulse_duration
            laserDur(i)     = thermal_s(i)%laser_duration

            ! ------ PULSE PARAMS ------ !
            pulse_period(i) = 1 / frequency(i)            !< Time between each pulse
            sigma_t(i)      = pulseDur(i) / 2.355         !< FWHM definition of st. dev

            if (p > 0) then
                ! 3D: Volume integral of 3D Gaussian
                spatial_int(i) = (2.0_wp * pi)**(1.5_wp) * radius(i)**3
            else
                ! 2D: Area integral of 2D Gaussian
                spatial_int(i) = 2.0_wp * pi * radius(i)**2
            end if

            temporal_int(i) = sigma_t(i) * sqrt(2.0_wp * pi)  !< Temporal Integral for Amplitude
            A0(i)           = E_pulse(i) / (spatial_int(i) * temporal_int(i))  !< Amplitude

            ! ----- TIME KEEPING ----- !
            laserStart(i)   = .true.
            TimeStart(i)    = 0.0_wp
        end do

        $:GPU_UPDATE(device='[thermal_loc,frequency,E_pulse,radius,pulseDur,laserDur,pulse_period,sigma_t,spatial_int,temporal_int,A0,laserStart,TimeStart]')



    end subroutine s_initialize_thermal_src

    !> This subroutine computes thermal source contributions to RHS
    !! @param rhs_vf RHS variables to be updated
    !! @param t_step Current time step
    impure subroutine s_compute_thermal_src(rhs_vf, t_step)
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer, intent(in) :: t_step

        ! --- Thermal Source local variable definitions --- !
        integer :: i, j, k, ti       !< Loop variable for number of thermal sources
        real(wp) :: dist_sq, source_val, spatial, temporal, t_mod, t_center !< Calculation variables
        real(wp) :: loc_x, loc_y, loc_z, rad, amp  !< Local copies for GPU
        integer :: th_i0, th_i1, th_j0, th_j1, th_k0, th_k1    ! < Loop boundary definitions
  

        do ti = 1, num_source_th

          ! Save time laser was turned on 
          if (laserStart(ti) .eqv. .true.) then
              TimeStart(ti) = mytime
              laserStart(ti) = .false.
          end if

          if (mytime > TimeStart(ti) + laserDur(ti)) cycle

          ! -------- TEMPORAL SECTION -------- !

          t_mod = mod(mytime,pulse_period(ti))
          t_center = 0.5 * pulse_period(ti)

          ! Skip if temporally far from pulse (beyond 4 sigma, contribution < 0.03%)
          if (abs(t_mod - t_center) > 4.0_wp * sigma_t(ti)) cycle

          temporal = exp(-((t_mod - t_center)**2) / (2.0_wp * sigma_t(ti)**2))

          ! -------- TEMPORAL SECTION -------- !



          ! -------- SPATIAL SECTION -------- !

          ! Non-halo cells of the core
          th_i0 = idwint(1)%beg
          th_i1 = idwint(1)%end
          th_j0 = idwint(2)%beg
          th_j1 = idwint(2)%end
          th_k0 = idwint(3)%beg
          th_k1 = idwint(3)%end

          ! Copy to local scalars for GPU (avoids array indexing inside parallel region)
            loc_x = thermal_loc(ti, 1)
            loc_y = thermal_loc(ti, 2)
            loc_z = thermal_loc(ti, 3)
            rad   = radius(ti)
            amp   = A0(ti)

          $:GPU_PARALLEL_LOOP(collapse=3, private='[dist_sq,spatial,source_val]') 
          do k = th_k0, th_k1
              do j = th_j0, th_j1
                  do i = th_i0, th_i1
                      source_val = 0.0_wp   !< Initial set to 0 to avoid NaNs

                      if (p > 0) then
                            ! 3D case
                            dist_sq = (x_cc(i) - loc_x)**2 + &
                                      (y_cc(j) - loc_y)**2 + &
                                      (z_cc(k) - loc_z)**2
                        else
                            ! 2D case
                            dist_sq = (x_cc(i) - loc_x)**2 + &
                                      (y_cc(j) - loc_y)**2
                        end if
                              
                      spatial = exp(-dist_sq / (2.0_wp * rad**2) )

                      source_val = amp * spatial * temporal

                      ! After computing source_val:
                      if (.not. (source_val == source_val)) then  ! NaN check
                          print *, "NaN detected: ti=", ti, "i,j,k=", i,j,k, &
                                  "dist_sq=", dist_sq, "spatial=", spatial, "temporal=", temporal
                          call s_mpi_abort("NaN in thermal source")
                      end if

                      rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) + source_val
                  end do
              end do
          end do

          ! -------- SPATIAL SECTION -------- !

        end do

        
    end subroutine s_compute_thermal_src


end module m_thermal_src
