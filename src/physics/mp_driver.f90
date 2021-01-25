module module_mp_driver
    use domain_interface,   only: domain_t
    use module_mp_thompson, only: thompson_init, mp_gt_driver
    implicit none
    logical :: initialized = .false.
  contains

    subroutine mp_init(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        call thompson_init()

        initialized = .true.

    end subroutine mp_init

    subroutine process_subdomain(domain, dt, its,ite, jts,jte, kts,kte)
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt
        integer,        intent(in)    :: its,ite, jts,jte, kts,kte

        call mp_gt_driver(  qv=domain%water_vapor%local,            &
                            th=domain%potential_temperature%local,  &
                            qc=domain%cloud_water_mass%local,       &
                            qi=domain%cloud_ice_mass%local,         &
                            ni=domain%cloud_ice_number%local,       &
                            qr=domain%rain_mass%local,              &
                            nr=domain%rain_number%local,            &
                            qs=domain%snow_mass%local,              &
                            qg=domain%graupel_mass%local,           &
                            pii=domain%exner,                       &
                            p=domain%pressure,                      &
                            w=domain%w%local,                       &
                            dz=domain%dz_mass,                      &
                            dt_in=dt,                               &
                            RAINNC=domain%accumulated_precipitation,&
                            SNOWNC=domain%accumulated_snowfall,     &
                            has_reqc=0, has_reqi=0, has_reqs=0,     &
                            ids=domain%ids,ide=domain%ide,          & ! domain dims
                            jds=domain%jds,jde=domain%jde,          &
                            kds=domain%kds,kde=domain%kde,          &
                            ims=domain%ims,ime=domain%ime,          & ! memory dims
                            jms=domain%jms,jme=domain%jme,          &
                            kms=domain%kms,kme=domain%kme,          &
                            its=its,ite=ite,          & ! tile dims
                            jts=jts,jte=jte,          &
                            kts=kts,kte=kte)


    end subroutine process_subdomain

    subroutine process_halo(domain, dt, halo)
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt
        integer,        intent(in)    :: halo

        integer :: its,ite, jts,jte, kts,kte

        its = domain%its
        ite = domain%ite
        jts = domain%jts
        jte = domain%jte
        kts = domain%kts
        kte = domain%kte

        ! process the western halo
        ite = its+halo-1
        call process_subdomain(domain, dt, its,ite, jts,jte, kts,kte)
        ite = domain%ite

        ! process the eastern halo
        its = ite-halo+1
        call process_subdomain(domain, dt, its,ite, jts,jte, kts,kte)
        its = domain%its

        its = its+halo
        ite = ite-halo

        ! process the southern halo
        jte = jts+halo-1
        call process_subdomain(domain, dt, its,ite, jts,jte, kts,kte)
        jte = domain%jte

        ! process the northern halo
        jts = jte-halo+1
        call process_subdomain(domain, dt, its,ite, jts,jte, kts,kte)
        ! jts = domain%jts

    end subroutine process_halo

    subroutine microphysics(domain, dt, halo, subset, t, report_convection, &
         total_num_particles)
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt
        integer,        intent(in),   optional :: halo, subset
        integer,        intent(in),   optional :: t, total_num_particles
        logical,        intent(in),   optional :: report_convection
        integer,        intent(in),   optional :: total_num_particles

        integer                       :: dz_lb(3), i, n_particles

        if (.not. initialized) call mp_init(domain)

        if (present(total_num_particles)) then
           n_particles = total_num_particles
        else
           n_particles = 0
        end if


        if (present(subset)) then

          dz_lb = lbound(domain%dz_interface)

          if (n_particles .gt. 0) then
          !  if (domain%convection_obj%do_replacement() .eqv. .false.) then
          !     print *, "---Artless: unsure why this is needed now---"
          !     stop
          !     call domain%convection_obj%process( &
          !         domain%nx_global, domain%ny_global, &
          !         domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
          !         dt, domain%dz_interface(dz_lb(1),dz_lb(2),dz_lb(3)), &
          !         domain%temperature, domain%z_interface(:,1,:), &
          !         domain%its, domain%ite, domain%kts, domain%kte, domain%jts, &
          !         domain%jte)
             ! else

             if (1 .eq. 1 ) then  ! take every time step
             do i=1,int(dt)
                if ((i .gt. 1)) then
                   call domain%convection_obj%retrieve(no_sync=.false.)
                end if
                call domain%convection_obj%process( &
                  domain%nx_global, domain%ny_global, &
                  domain%ims, domain%ime, domain%kms, &
                  domain%kme, domain%jms, domain%jme, &
                  1.0, domain%dz_interface(dz_lb(1),dz_lb(2),dz_lb(3)), &
                  domain%temperature, domain%z_interface(:,1,:), &
                  domain%its, domain%ite, domain%kts, domain%kte, domain%jts, &
                  domain%jte, domain%z, domain%potential_temperature, &
                  domain%pressure, domain%u, domain%v, domain%w, &
                  int((t-1)*dt + i) )
                if (report_convection .eqv. .true.) &
                     call domain%report_convection(int((t-1)*dt + i))
             end do
             end if
             ! end if
             if (0 .eq. 1 ) then  ! no extra time steps taken
                call domain%convection_obj%process( &
                  domain%nx_global, domain%ny_global, &
                  domain%ims, domain%ime, domain%kms, &
                  domain%kme, domain%jms, domain%jme, &
                  dt, domain%dz_interface(dz_lb(1),dz_lb(2),dz_lb(3)), &
                  domain%temperature, domain%z_interface(:,1,:), &
                  domain%its, domain%ite, domain%kts, domain%kte, domain%jts, &
                  domain%jte, domain%z, domain%potential_temperature, &
                  domain%pressure, domain%u, domain%v, domain%w)
             end if

          end if

            call process_subdomain(domain, dt,                               &
                                   domain%its + subset, domain%ite - subset, &
                                   domain%jts + subset, domain%jte - subset, &
                                   domain%kts,          domain%kte)
        endif

        if (present(halo)) then
            call process_halo(domain, dt, halo)
        endif

        if ((.not.present(halo)).and.(.not.present(subset))) then

            call mp_gt_driver(  qv=domain%water_vapor%local,            &
                                th=domain%potential_temperature%local,  &
                                qc=domain%cloud_water_mass%local,       &
                                qi=domain%cloud_ice_mass%local,         &
                                ni=domain%cloud_ice_number%local,       &
                                qr=domain%rain_mass%local,              &
                                nr=domain%rain_number%local,            &
                                qs=domain%snow_mass%local,              &
                                qg=domain%graupel_mass%local,           &
                                pii=domain%exner,                       &
                                p=domain%pressure,                      &
                                w=domain%w%local,                       &
                                dz=domain%dz_mass,                      &
                                dt_in=dt,                               &
                                RAINNC=domain%accumulated_precipitation,&
                                SNOWNC=domain%accumulated_snowfall,     &
                                has_reqc=0, has_reqi=0, has_reqs=0,     &
                                ids=domain%ids,ide=domain%ide,          & ! domain dims
                                jds=domain%jds,jde=domain%jde,          &
                                kds=domain%kds,kde=domain%kde,          &
                                ims=domain%ims,ime=domain%ime,          & ! memory dims
                                jms=domain%jms,jme=domain%jme,          &
                                kms=domain%kms,kme=domain%kme,          &
                                its=domain%its,ite=domain%ite,          & ! tile dims
                                jts=domain%jts,jte=domain%jte,          &
                                kts=domain%kts,kte=domain%kte)

        endif

    end subroutine microphysics

end module module_mp_driver
