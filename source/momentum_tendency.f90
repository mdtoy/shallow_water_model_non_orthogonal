module momentum_tendency

!-----------------------------------------------------------------------
! PURPOSE: Calculates the tendencies of u and v.
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters


implicit none
save



contains



!======================================================================
! BEGINNING OF GET_U1COVF_U2COVF
!======================================================================

subroutine get_u1covf_u2covf ( u1_cov, u2_cov, u1_cont, u2_cont,     &
                 F_u1, F_u2, h, h_star, sqrt_G_h_star,               &
                 pv, zeta, ke_horiz, u1_cov_f, u2_cov_f )

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of the x and y components of velocity,
!   i.e. u and v respectively.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm), intent(in) ::            &
          u1_cov, u2_cov,  &    ! covariant velocity components (m/s)
          u1_cont, u2_cont      ! contravariant velocity components (m/s)

real (kind=dbl_kind), dimension(im,jm), intent(in) ::            &
          F_u1, F_u2,  & ! Mass fluxes in u1 and u2 directions
          h,           & ! height of free surface (m)
          h_star,      & ! thickness of fluid (m) (h-hs)
          sqrt_G_h_star  ! sqrt_G times h_star (m)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm), intent(out) ::           &
          pv,      &     ! vert. component of pot. vorticity (m-1 s-1)
          zeta,    &     ! relative vorticity (s^-1)
          ke_horiz       ! contribution to kinetic energy
                         ! from horizontal velocity (J/kg)

real (kind=dbl_kind), dimension(im,jm), intent(out) ::           &
          u1_cov_f, u2_cov_f      ! tendency of u1 and u2 (m/s^2)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j

real (kind=dbl_kind), parameter ::                                   &
          inv24 = c1/24.00000_dbl_kind

real (kind=dbl_kind), dimension(im,jm) ::                            &
          abs_vort,        &                  ! absolute vorticity (s^-1)
          h_q,             &                  ! SW height interpolated to q-pts (m)
          alfa, beta, gamm, delt, epsln, fi   ! linear combinations of pv





!---------------------------------------------------------------------------
! Calculate the (PV k) cross (m times velocity) term.       (term 1 of 5)
! Potential enstrophy and energy conserving scheme of Arakawa and 
! Lamb (1981) is used.
!---------------------------------------------------------------------------

! calculate potential vorticity
zeta(:,:) = inv_sqrt_G_q(:,:) *                                        &
                  ( invdx2 * ( u1_cov(:,jm1(:)) - u1_cov(:,:) )  +     &
                    invdx1 * ( u2_cov(:,:) - u2_cov(im1(:),:) ) )
abs_vort(:,:) = f_cor(:,:) + zeta(:,:)
h_q(:,:) = inv_sqrt_G_q(:,:) * p25 * ( sqrt_G_h_star(:,:) +            &
              sqrt_G_h_star(im1(:),:) + sqrt_G_h_star(im1(:),jm1(:)) + &
              sqrt_G_h_star(:,jm1(:)) )
pv(:,:) = abs_vort(:,:) / h_q(:,:)


! calculate linear combinations of pv ( eqn. (3.34) of AL (1981) )
alfa(:,:)  = inv24 * ( c2*pv(ip1(:),jp1(:)) + pv(:,jp1(:)) +  &
                       c2*pv(:,:) + pv(ip1(:),:) )
beta(:,:)  = inv24 * ( pv(:,jp1(:)) + c2*pv(im1(:),jp1(:)) +  &
                       pv(im1(:),:) + c2*pv(:,:) )
gamm(:,:)  = inv24 * ( c2*pv(:,jp1(:)) + pv(im1(:),jp1(:)) +  &
                       c2*pv(im1(:),:) + pv(:,:) )
delt(:,:)  = inv24 * ( pv(ip1(:),jp1(:)) + c2*pv(:,jp1(:)) +  &
                       pv(:,:) + c2*pv(ip1(:),:) )
epsln(:,:) = inv24 * ( pv(ip1(:),jp1(:)) + pv(:,jp1(:)) -     &
                       pv(:,:) - pv(ip1(:),:) )
fi(:,:)    = inv24 * (-pv(ip1(:),jp1(:)) + pv(:,jp1(:)) +     &
                       pv(:,:) - pv(ip1(:),:) ) 


!
! calculate u and v tendencies -- term 1 of 5 ( see eqns. (3.5) and (3.6)
!                                 of AL (1981) )

u1_cov_f(:,:) = alfa(:,:)*F_u2(:,jp1(:)) + beta(:,:)*F_u2(im1(:),jp1(:)) + &
           gamm(:,:)*F_u2(im1(:),:) + delt(:,:)*F_u2(:,:) -                &
           epsln(:,:)*F_u1(ip1(:),:) + epsln(im1(:),:)*F_u1(im1(:),:)
u2_cov_f(:,:) = -gamm(ip1(:),:)*F_u1(ip1(:),:) - delt(:,:)*F_u1(:,:) -     &
           alfa(:,jm1(:))*F_u1(:,jm1(:)) -                                 &
           beta(ip1(:),jm1(:))*F_u1(ip1(:),jm1(:)) -                       &
           fi(:,:)*F_u2(:,jp1(:)) + fi(:,jm1(:))*F_u2(:,jm1(:))



!---------------------------------------------------------------------------
! Add contribution of the horiz. gradient of kinetic energy.   (term 2 of 5)
! And also divergence damping!!!
!---------------------------------------------------------------------------


! Calculate contribution to kinetic energy
! from the horizontal velocity  ( see eqn (3.41) of AL (1981) )
! Note:  expect SICK to result from use of this K.E.
ke_horiz(:,:) = p25 * ( u1_cov(:,:)*u1_cont(:,:) +                   &
                        u1_cov(ip1(:),:)*u1_cont(ip1(:),:) ) +       &
                p25 * ( u2_cov(:,:)*u2_cont(:,:) +                   &
                        u2_cov(:,jp1(:))*u2_cont(:,jp1(:)) )

! Add contribution of horiz. gradient of K.E.
u1_cov_f(:,:) = u1_cov_f(:,:) - invdx1*(ke_horiz(:,:)-ke_horiz(im1(:),:))
u2_cov_f(:,:) = u2_cov_f(:,:) - invdx2*(ke_horiz(:,:)-ke_horiz(:,jm1(:)))



!---------------------------------------------------------------------------
! Add contributions due to vert. advection of horiz. momentum and
! subgrid-scale turbulent momentum flux.            (terms 3 & 4 of 5)
! XXXXXXXXXXXXXXX NOT APPLICABLE TO SHALLOW WATER EQUATIONS XXXXXXXXXXXX
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
! Add contribution of the horizontal pressure gradient force.  (term 5 of 5)
!---------------------------------------------------------------------------

u1_cov_f(:,:) = u1_cov_f(:,:) - invdx1 * grav * ( h(:,:) - h(im1(:),:) )
u2_cov_f(:,:) = u2_cov_f(:,:) - invdx2 * grav * ( h(:,:) - h(:,jm1(:)) )




end subroutine get_u1covf_u2covf

!======================================================================
! END OF GET_UF_VF
!======================================================================



end module momentum_tendency
