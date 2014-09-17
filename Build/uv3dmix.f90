      MODULE uv3dmix_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes horizontal viscosity of momentum.             !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC uv3dmix2
      CONTAINS
      SUBROUTINE uv3dmix2 (ng, tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes harmonic mixing of momentum, rotated along    !
!  geopotentials,  from the  horizontal  divergence  of the  stress    !
!  tensor.  A transverse  isotropy is assumed so the stress tensor     !
!  is split into vertical and horizontal subtensors.                   !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Wajsowicz, R.C, 1993: A consistent formulation of the           !
!         anisotropic stress tensor for use in models of the           !
!         large-scale ocean circulation, JCP, 105, 333-338.            !
!                                                                      !
!      Sadourny, R. and K. Maynard, 1997: Formulations of              !
!         lateral diffusion in geophysical fluid dynamics              !
!         models, In Numerical Methods of Atmospheric and              !
!         Oceanic Modelling. Lin, Laprise, and Ritchie,                !
!         Eds., NRC Research Press, 547-556.                           !
!                                                                      !
!      Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic              !
!         friction with a Smagorinsky-like viscosity for               !
!         use in large-scale eddy-permitting ocean models,             !
!         Monthly Weather Rev., 128, 8, 2935-2946.                     !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 31)
      CALL uv3dmix2_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nnew(ng),                           &
     &                    GRID(ng) % om_p,                              &
     &                    GRID(ng) % om_r,                              &
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_p,                              &
     &                    GRID(ng) % on_r,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    MIXING(ng) % visc2_p,                         &
     &                    MIXING(ng) % visc2_r,                         &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    COUPLING(ng) % rufrc,                         &
     &                    COUPLING(ng) % rvfrc)
      CALL wclock_off (ng, iNLM, 31)
      RETURN
      END SUBROUTINE uv3dmix2
!
!***********************************************************************
      SUBROUTINE uv3dmix2_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nnew,                             &
     &                          om_p, om_r, om_u, om_v,                 &
     &                          on_p, on_r, on_u, on_v,                 &
     &                          pm, pn,                                 &
     &                          Hz, z_r,                                &
     &                          visc2_p, visc2_r,                       &
     &                          u, v,                                   &
     &                          rufrc, rvfrc)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nnew
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)
      real(r8), intent(inout) :: rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k, k1, k2
      real(r8) :: cff, fac1, fac2, pm_p, pn_p
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: cff5, cff6, cff7, cff8
      real(r8) :: dmUdz, dnUdz, dmVdz, dnVdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: UFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: UFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: VFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: VFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dmUde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dmVde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dnUdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dnVdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dUdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dVdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_r
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic viscosity along geopotential surfaces.
!-----------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients.  Notice the recursive
!  blocking sequence.  The vertical placement of the gradients is:
!
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k1) k      rho-points
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k2) k+1    rho-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k1) k      psi-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k2) k+1    psi-points
!                UFse, UFsx, dUdz(:,:,k1) k-1/2  WU-points
!                UFse, UFsx, dUdz(:,:,k2) k+1/2  WU-points
!                VFse, VFsx, dVdz(:,:,k1) k-1/2  WV-points
!                VFse, VFsx, dVdz(:,:,k2) k+1/2  WV-points
!
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=0.5_r8*(pm(i-1,j)+pm(i,j))
              UFx(i,j)=cff*(z_r(i  ,j,k+1)-                             &
     &                      z_r(i-1,j,k+1))
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=0.5_r8*(pn(i,j-1)+pn(i,j))
              VFe(i,j)=cff*(z_r(i,j  ,k+1)-                             &
     &                      z_r(i,j-1,k+1))
            END DO
          END DO
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              dZdx_p(i,j,k2)=0.5_r8*(UFx(i,j-1)+                        &
     &                               UFx(i,j  ))
              dZde_p(i,j,k2)=0.5_r8*(VFe(i-1,j)+                        &
     &                               VFe(i  ,j))
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              dZdx_r(i,j,k2)=0.5_r8*(UFx(i  ,j)+                        &
     &                               UFx(i+1,j))
              dZde_r(i,j,k2)=0.5_r8*(VFe(i,j  )+                        &
     &                               VFe(i,j+1))
            END DO
          END DO
!
!  Compute momentum horizontal (1/m/s) and vertical (1/s) gradients.
!
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=0.5_r8*pm(i,j)
              dnUdx(i,j,k2)=cff*((pn(i  ,j)+pn(i+1,j))*                 &
     &                           u(i+1,j,k+1,nrhs)-                     &
     &                           (pn(i-1,j)+pn(i  ,j))*                 &
     &                           u(i  ,j,k+1,nrhs))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                      &
     &                      pn(i-1,j-1)+pn(i,j-1))
              dmUde(i,j,k2)=cff*((pm(i-1,j  )+pm(i,j  ))*               &
     &                           u(i,j  ,k+1,nrhs)-                     &
     &                           (pm(i-1,j-1)+pm(i,j-1))*               &
     &                           u(i,j-1,k+1,nrhs))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                      &
     &                      pm(i-1,j-1)+pm(i,j-1))
              dnVdx(i,j,k2)=cff*((pn(i  ,j-1)+pn(i  ,j))*               &
     &                           v(i  ,j,k+1,nrhs)-                     &
     &                           (pn(i-1,j-1)+pn(i-1,j))*               &
     &                           v(i-1,j,k+1,nrhs))
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=0.5_r8*pn(i,j)
              dmVde(i,j,k2)=cff*((pm(i,j  )+pm(i,j+1))*                 &
     &                           v(i,j+1,k+1,nrhs)-                     &
     &                           (pm(i,j-1)+pm(i,j  ))*                 &
     &                           v(i,j  ,k+1,nrhs))
            END DO
          END DO
        END IF
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              dUdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              dVdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              UFsx(i,j,k2)=0.0_r8
              UFse(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              VFsx(i,j,k2)=0.0_r8
              VFse(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-z_r(i-1,j,k)+          &
     &                            z_r(i  ,j,k+1)-z_r(i  ,j,k)))
              dUdz(i,j,k2)=cff*(u(i,j,k+1,nrhs)-                        &
     &                          u(i,j,k  ,nrhs))
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-z_r(i,j-1,k)+          &
     &                            z_r(i,j  ,k+1)-z_r(i,j  ,k)))
              dVdz(i,j,k2)=cff*(v(i,j,k+1,nrhs)-                        &
     &                          v(i,j,k  ,nrhs))
            END DO
          END DO
        END IF
!
!  Compute components of the rotated viscous flux (m5/s2) along
!  geopotential surfaces in the XI- and ETA-directions.
!
        IF (k.gt.0) THEN
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
              cff3=MIN(dZde_r(i,j,k1),0.0_r8)
              cff4=MAX(dZde_r(i,j,k1),0.0_r8)
              cff=Hz(i,j,k)*                                            &
     &            (on_r(i,j)*(dnUdx(i,j,k1)-                            &
     &                        0.5_r8*pn(i,j)*                           &
     &                        (cff1*(dUdz(i  ,j,k1)+                    &
     &                               dUdz(i+1,j,k2))+                   &
     &                         cff2*(dUdz(i  ,j,k2)+                    &
     &                               dUdz(i+1,j,k1))))-                 &
     &             om_r(i,j)*(dmVde(i,j,k1)-                            &
     &                        0.5_r8*pm(i,j)*                           &
     &                        (cff3*(dVdz(i,j  ,k1)+                    &
     &                               dVdz(i,j+1,k2))+                   &
     &                         cff4*(dVdz(i,j  ,k2)+                    &
     &                               dVdz(i,j+1,k1)))))
              UFx(i,j)=on_r(i,j)*on_r(i,j)*visc2_r(i,j)*cff
              VFe(i,j)=om_r(i,j)*om_r(i,j)*visc2_r(i,j)*cff
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              pm_p=0.25_r8*(pm(i-1,j-1)+pm(i-1,j)+                      &
     &                      pm(i  ,j-1)+pm(i  ,j))
              pn_p=0.25_r8*(pn(i-1,j-1)+pn(i-1,j)+                      &
     &                      pn(i  ,j-1)+pn(i  ,j))
              cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
              cff3=MIN(dZde_p(i,j,k1),0.0_r8)
              cff4=MAX(dZde_p(i,j,k1),0.0_r8)
              cff=0.25_r8*                                              &
     &            (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                           &
     &             Hz(i-1,j-1,k)+Hz(i,j-1,k))*                          &
     &            (on_p(i,j)*(dnVdx(i,j,k1)-                            &
     &                        0.5_r8*pn_p*                              &
     &                        (cff1*(dVdz(i-1,j,k1)+                    &
     &                               dVdz(i  ,j,k2))+                   &
     &                         cff2*(dVdz(i-1,j,k2)+                    &
     &                               dVdz(i  ,j,k1))))+                 &
     &             om_p(i,j)*(dmUde(i,j,k1)-                            &
     &                        0.5_r8*pm_p*                              &
     &                        (cff3*(dUdz(i,j-1,k1)+                    &
     &                               dUdz(i,j  ,k2))+                   &
     &                         cff4*(dUdz(i,j-1,k2)+                    &
     &                               dUdz(i,j  ,k1)))))
              UFe(i,j)=om_p(i,j)*om_p(i,j)*visc2_p(i,j)*cff
              VFx(i,j)=on_p(i,j)*on_p(i,j)*visc2_p(i,j)*cff
            END DO
          END DO
!
!  Compute vertical flux (m2/s2) due to sloping terrain-following
!  surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                cff=0.25_r8*(visc2_r(i-1,j)+visc2_r(i,j))
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                dnUdz=cff*dUdz(i,j,k2)
                dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                dmUdz=cff*dUdz(i,j,k2)
                dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                UFsx(i,j,k2)=fac1*                                      &
     &                       (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+        &
     &                        cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+        &
     &                        cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+        &
     &                        cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                UFse(i,j,k2)=fac2*                                      &
     &                       (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+        &
     &                        cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+        &
     &                        cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+        &
     &                        cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff5=MIN(dZdx_p(i,j  ,k1),0.0_r8)
                cff6=MIN(dZdx_p(i,j+1,k2),0.0_r8)
                cff7=MAX(dZdx_p(i,j  ,k2),0.0_r8)
                cff8=MAX(dZdx_p(i,j+1,k1),0.0_r8)
                UFsx(i,j,k2)=UFsx(i,j,k2)+                              &
     &                       fac1*                                      &
     &                       (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+        &
     &                        cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+        &
     &                        cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+        &
     &                        cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff5=MIN(dZde_r(i-1,j,k1),0.0_r8)
                cff6=MIN(dZde_r(i  ,j,k2),0.0_r8)
                cff7=MAX(dZde_r(i-1,j,k2),0.0_r8)
                cff8=MAX(dZde_r(i  ,j,k1),0.0_r8)
                UFse(i,j,k2)=UFse(i,j,k2)-                              &
     &                       fac2*                                      &
     &                       (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+        &
     &                        cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+        &
     &                        cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+        &
     &                        cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
              END DO
            END DO
!
            DO j=JstrV,Jend
              DO i=Istr,Iend
                cff=0.25_r8*(visc2_r(i,j-1)+visc2_r(i,j))
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                dnUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dnVdz=cff*dVdz(i,j,k2)
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                dmUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dmVdz=cff*dVdz(i,j,k2)
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                VFsx(i,j,k2)=fac1*                                      &
     &                       (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+        &
     &                        cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+        &
     &                        cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+        &
     &                        cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                VFse(i,j,k2)=fac2*                                      &
     &                       (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+        &
     &                        cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+        &
     &                        cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+        &
     &                        cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff5=MIN(dZdx_r(i,j-1,k1),0.0_r8)
                cff6=MIN(dZdx_r(i,j  ,k2),0.0_r8)
                cff7=MAX(dZdx_r(i,j-1,k2),0.0_r8)
                cff8=MAX(dZdx_r(i,j  ,k1),0.0_r8)
                VFsx(i,j,k2)=VFsx(i,j,k2)-                              &
     &                       fac1*                                      &
     &                       (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+        &
     &                        cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+        &
     &                        cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+        &
     &                        cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff5=MIN(dZde_p(i  ,j,k1),0.0_r8)
                cff6=MIN(dZde_p(i+1,j,k2),0.0_r8)
                cff7=MAX(dZde_p(i  ,j,k2),0.0_r8)
                cff8=MAX(dZde_p(i+1,j,k1),0.0_r8)
                VFse(i,j,k2)=VFse(i,j,k2)+                              &
     &                       fac2*                                      &
     &                       (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+        &
     &                        cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+        &
     &                        cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+        &
     &                        cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
              END DO
            END DO
          END IF
!
!  Time-step harmonic, geopotential viscosity term.  Notice that
!  momentum at this stage is HzU and HzV and has m2/s units.  Add
!  contribution for barotropic forcing terms.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=dt(ng)*0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
              cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))
              cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j))
              cff3=UFsx(i,j,k2)-UFsx(i,j,k1)
              cff4=UFse(i,j,k2)-UFse(i,j,k1)
              cff5=cff*(cff1+cff2)
              cff6=dt(ng)*(cff3+cff4)
              rufrc(i,j)=rufrc(i,j)+cff1+cff2+cff3+cff4
              u(i,j,k,nnew)=u(i,j,k,nnew)+cff5+cff6
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=dt(ng)*0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
              cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
              cff3=VFsx(i,j,k2)-VFsx(i,j,k1)
              cff4=VFse(i,j,k2)-VFse(i,j,k1)
              cff5=cff*(cff1-cff2)
              cff6=dt(ng)*(cff3+cff4)
              rvfrc(i,j)=rvfrc(i,j)+cff1-cff2+cff3+cff4
              v(i,j,k,nnew)=v(i,j,k,nnew)+cff5+cff6
            END DO
          END DO
        END IF
      END DO K_LOOP
      RETURN
      END SUBROUTINE uv3dmix2_tile
      END MODULE uv3dmix_mod