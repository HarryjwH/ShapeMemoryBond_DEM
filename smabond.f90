subroutine bonds_sma (ia_bnd)

use grcom
use flcom
use adh
use ncblk1
use ncblk
use ioproc
use fxconst
use frd
use smadata

implicit none
integer, intent (in) :: ia_bnd
integer ::  i, ibody_fiber
integer, dimension(2) :: ipath
real (kind = kind(0.d0) ), dimension(3) :: tem_ft, tem_mn,     &
  & tem_mt, tem_ft1, tem_mn1, tem_mt1, np1p2
real (kind = kind(0.d0) ), dimension(3) :: x_rel_global, vel1, vel2,   &
  &   wrot1, wrot2, r1_cen2con, r2_cen2con, vc1_global, vc2_global, vn_rel, vt_rel, &
    &   w_rel_global, wn_rel_local, wt_rel_local, v1_w, v2_w, torq1,  &
    &   torq2,  ftot, x_p1_global, x_p2_global, ftan, mom_nor, &
    &   mom_tan, vc_rel_global, w_rel_local, vc_rel_local, &
    &   vcn_rel_local, vct_rel_local, x_rel_local, xn_rel_local, &
    &   xt_rel_local, vcn_rel_global, vct_rel_global
real (kind = kind(0.d0) ) :: nom, vn12, wn12, ovlp, ovlp_n, ovlp_t, ovlp2
real (kind = kind(0.d0) ) :: len_ft, len_mn, len_mt, len_ft1,   &
  &   len_mn1, len_mt1, norm_ft, norm_mn, norm_mt, dis_p1p2,  &
&   delta_n
real :: stiff_norm,  stiff_tan, a_bnd, j_bnd, mstiff_norm, mstiff_tan
real :: mstar, moi_star, damp_n, damp_t, damp_tor, damp_bend, fbnd_norm
real :: stiff_com, sigma_bnd, tau_bnd, sigma_lmt, tau_lmt
real :: fn_bnd, ft_bnd, mn_bnd, mt_bnd, thres1, thres2, thres3
real, dimension (3) :: fdamp_n_global, fdamp_t_global, mdamp_n_local, mdamp_t_local
real :: dens_temp, moi_star_bend, moi_star_twist
real (kind = kind(0.d0) ) :: r11, r12, r13, r21, r22, r23, r31, r32, r33
real (kind = kind(0.d0) ), dimension(3) :: fnorm_global, ftan_global, mom_twist_local, &
   &     mom_bend_local, ftot_local,  mom_twist_tot_local, mom_bend_tot_local, &
&     ftot_global, mom_twist_tot_global, mom_bend_tot_global, moment_temp,  &
&     unitVectorAxis
real :: kelastic, kep, kcp, melast, mplast, theta_elastic, theta_plastic
real :: mom_magnitude, theta_tot, omega_bend, dtheta, moment_axis, vcn_global_scaler
real :: ubx_local, uby_local, ubz_local, wn_mag
real:: mms, mmf, dmt, mthe, theta_temp, smafrac, sma_crimoment, SigmaMs, SigmaMf, Mup, Mdown
real:: theta0, theta1, theta2, theta3, theta4, theta5, theta6

!!--check whether the bond is already broken
if (ia(ia_bnd+19) .eq. 0) return

iab1 = ia(ia_bnd+1)
iab2 = ia(ia_bnd+2)

ibody_fiber = ia(iab1+22)
lbond = clen (ibody_fiber,ia(ia_bnd+67))

thres = r(ibody_fiber)
thres3 = r(ibody_fiber)*0.0001

stiff_norm = a(ia_bnd+3)
stiff_tan  = a(ia_bnd+4)
mstiff_norm = a(ia_bnd+17)

mstiff_tan  = a(ia_bnd+18)

sigma_lmt = a(ia_bnd+20)
tau_lmt = a(ia_bnd+21)
!------unit vector of the bond direcion in the template
ubx_local = a(ia_bnd+64)
uby_local = a(ia_bnd+65)
ubz_local = a(ia_bnd+66)

!---read normal vectors of particles---------
do i = 1, 3
wrot1(i) = a(iab1+12+i)
wrot2(i) = a(iab2+12+i)
w_rel_global(i) = wrot2(i)-wrot1(i)
enddo

do i=1, 3
x_p1_global(i) =a(iab1+i-1)+a(iab1+i+2)
x_p2_global(i) =a(iab2+i-1)+a(iab2+i+2)
x_rel_global(i) = x_p2_global(i) - x_p1_global(i)
r1_cen2con (i) = 0.5*(x_p1_global(i)+x_p2_global(i))-x_p1_global(i)
r2_cen2con (i) = 0.5*(x_p1_global(i)+x_p2_global(i))-x_p2_global(i)
 enddo

dis_p1p2 = dsqrt (x_rel_global(1)**2.0+x_rel_global(2)**2.0+x_rel_global(3)**2.0)

do i=1, 3
 np1p2 (i) = x_rel_global(i)/dis_p1p2
enddo

! delta_n for normal bond force
delta_n = dis_p1p2 - lbond

 if (delta_n .gt. del(1)) then
  print*, 'The normal displacement in ', &
    &  ' a bonded particle is greater than the cell size!'
  stop
 endif
vcn_global_scaler =  vc_rel_global(1)*np1p2(1) + vc_rel_global(2)*np1p2(2) &
&                    + vc_rel_global(3)*np1p2(3)

 do i=1, 3
vcn_rel_global (i) = vcn_global_scaler*np1p2(i)
vct_rel_global (i) = vc_rel_global (i) -  vcn_rel_global (i)
 enddo

a_bnd=pi*r(ibody_fiber)**2.0;
j_bnd=0.5*pi*r(ibody_fiber)**4.0;

thres1=thres * (mstiff_tan*j_bnd/lbond)
thres2=thres * (mstiff_norm*(0.5*j_bnd)/lbond)

do i = 1, 3

a(ia_bnd+4+i) = (stiff_norm*a_bnd/lbond)*delta_n*np1p2 (i)
a(ia_bnd+7+i) = a(ia_bnd+7+i) + (stiff_tan*a_bnd/lbond)*vct_rel_global (i)*tdel
a(ia_bnd+21+i) = a(ia_bnd+21+i) + (mstiff_tan*j_bnd/lbond)*wn_rel_local(i)*tdel
if (abs(a(ia_bnd+21+i)) .gt. thres1) then
a(ia_bnd+10+i) = a(ia_bnd+10+i) + a(ia_bnd+21+i)
a(ia_bnd+21+i) = 0.0
endif

enddo

! Shape memory alloy related parameters changes with temperature. 
SigmaMs=Crisigma1(Temperature)
SigmaMf=Crisigma2(Temperature)

kelastic = (Ea+MartenFrac(Temperature)*(Em-Ea))*(0.5*j_bnd)/lbond
kep = (Ea+MartenFrac(Temperature)*(Em-Ea))*(SigmaMf-SigmaMs) &
& /( (Ea+MartenFrac(Temperature)*(Em-Ea)) * Hmax *MartenFrac(Temperature) + (SigmaMf-SigmaMs) ) &
& * (0.5*j_bnd)/lbond
kcp = 0.0

mms = SigmaMs*(0.5*j_bnd)/r(ibody_fiber)
mmf = SigmaMf*(0.5*j_bnd)/r(ibody_fiber)

theta_elastic = mms/kelastic
theta_plastic = theta_elastic + (mmf-mms)/kep

smafrac=0.005
sma_crimoment=(Ea+Em)/2.0*(0.5*j_bnd)/lbond*(Hmax*lbond/r(ibody_fiber))*smafrac

ipath (1) = ia(ia_bnd+38)
ipath (2) = ia(ia_bnd+39)


if( (Temperature .lt. As) .or. (Temperature .gt. Af) .or. (Temperature .ge. As .and. Temperature .le. Af .and. Rotc .eq. 0 ) )then
a(ia_bnd+70) = 0
a(ia_bnd+71) = 0
end if

! Internal bending moment is updated and the loading and unloading path is checked for switching
! when the temperature rises within the phase transition temperature range.
if( Temperature .ge. As .and. Temperature .lt. Af .and. (Rotc .gt. 0.0)   )then
! in the elastic range
if ((ipath (1).eq.1) .and. (ipath (2).eq.1)) then

if( abs( a(ia_bnd+48) ) .le.theta_elastic ) then

a(ia_bnd+70) = kelastic * a(ia_bnd+48)
dmt=a(ia_bnd+70)-a(ia_bnd+46)

else
a(ia_bnd+70) = kelastic*theta_elastic +kep*( a(ia_bnd+48)-theta_elastic )
dmt=a(ia_bnd+70)-a(ia_bnd+46)
ipath(1) = 2

end if

a(ia_bnd+46)=a(ia_bnd+46)+dmt
 a(ia_bnd+14) = a(ia_bnd+46)*a(ia_bnd+40)
 a(ia_bnd+15) = a(ia_bnd+46)*a(ia_bnd+41)
 a(ia_bnd+16) = a(ia_bnd+46)*a(ia_bnd+42)

else

!Beyond the elastic range, update the bending moment and path along the two main axes.
do  i =1, 2

if (ipath(i) .eq. 1) then
if( abs(a(ia_bnd+47+i)) .le. theta_elastic ) then
  mthe = kelastic*a(ia_bnd+47+i)
  dmt=mthe-a(ia_bnd+45+i)
else
    dmt=kelastic*theta_elastic+(a(ia_bnd+47+i)-theta_elastic)*kep-a(ia_bnd+45+i)
end if

elseif (ipath(i) .eq. 2) then

    mthe=kelastic*(a(ia_bnd+47+i)-theta_plastic)+mmf

elseif ( (ipath(i) .eq. 3) .or. (ipath(i).eq.5) ) then
a(ia_bnd+71+i)=Hmax*lbond/r(ibody_fiber)*MartenFrac(Temperature)

    mthe=kelastic*(a(ia_bnd+47+i)-theta_plastic)+mmf

elseif (ipath(i) .eq. 4) then

a(ia_bnd+77+i)= MartenFrac(Temperature)  * a(ia_bnd+71+i)
a(ia_bnd+73+i) = ( mmf-kep*theta_plastic +kelastic*a(ia_bnd+77+i) )/(kelastic-kep)
a(ia_bnd+75+i) = a(ia_bnd+73+i)-2*theta_elastic
Mup = kelastic*(a(ia_bnd+73+i)-a(ia_bnd+77+i))
Mdown = kelastic*(a(ia_bnd+75+i)-a(ia_bnd+77+i))

  mthe=kelastic * ( a(ia_bnd+47+i)-a(ia_bnd+77+i) )
  dmt=mthe-a(ia_bnd+45+i)

elseif (ipath(i) .eq. 6) then
dmt=0
else
print*, 'ipath = ', ipath(i), '   is out of range!'
stop
endif

enddo
endif

!Update bond moments and path types when temperature decreases

elseif (Temperature .ge. As .and. Temperature .lt. Af .and. (Rotc .lt. 0.0)   )then

if ((ipath (1).eq.1) .and. (ipath (2).eq.1)) then

if( abs( a(ia_bnd+48) ) .le.theta_elastic ) then

a(ia_bnd+70) = kelastic * a(ia_bnd+48)
dmt=a(ia_bnd+70)-a(ia_bnd+46)

else
a(ia_bnd+70) = kelastic*theta_elastic +kep*( a(ia_bnd+48)-theta_elastic )
dmt=a(ia_bnd+70)-a(ia_bnd+46)

end if

 a(ia_bnd+14) = a(ia_bnd+46)*a(ia_bnd+40)
 a(ia_bnd+15) = a(ia_bnd+46)*a(ia_bnd+41)
 a(ia_bnd+16) = a(ia_bnd+46)*a(ia_bnd+42)

else
do i=1,2
if (ipath(i) .eq. 1) then
mthe = kelastic*a(ia_bnd+47+i)
dmt=mthe-a(ia_bnd+45+i)
elseif(ipath(i) .eq.2)then
if(a(ia_bnd+45+i).gt.0)then
    mthe = kelastic* (a(ia_bnd+47+i)+theta_plastic) - mmf
else
    mthe = kelastic* (a(ia_bnd+47+i)-theta_plastic) + mmf
endif

dmt = mthe-a(ia_bnd+45+i)

elseif((ipath(i) .eq. 3) .or. (ipath(i).eq.5))then


theta_temp= a(ia_bnd+47+i) + dtheta

theta0 = a(ia_bnd+47+i)-a(ia_bnd+45+i)/kelastic
theta1 = (kelastic*theta0 + kep*theta_plastic - mmf)/(kelastic - kep)
theta2 = (kelastic*theta0 - kep*theta_plastic + mmf)/(kelastic - kep)
theta3 = theta_plastic
theta4 = theta2
theta5 = theta1
theta6 = -theta_plastic

mthe=kep*(a(ia_bnd+47+i)+theta_plastic)-mmf
end if
dmt=mthe-a(ia_bnd+45+i)
elseif(ipath(i).eq.4)then
theta_temp= a(ia_bnd+47+i) + dtheta


theta0 = a(ia_bnd+47+i)-a(ia_bnd+45+i)/kelastic
theta1 = (kelastic*theta0 + kep*theta_plastic - mmf)/(kelastic - kep)
theta2 = (kelastic*theta0 - kep*theta_plastic + mmf)/(kelastic - kep)

dmt=mthe-a(ia_bnd+45+i)

elseif (ipath(i) .eq. 6)then

if(a(ia_bnd+45+i) .gt. 0)then
    mthe=kelastic*a(ia_bnd+47+i)

dmt=mthe-a(ia_bnd+45+i)

else

print*, 'ipath = ', ipath(i), '   is out of range!'
stop

endif

a(ia_bnd+69+i) = dmt
a(ia_bnd+45+i)=a(ia_bnd+45+i)+dmt

end do
endif

end if

end subroutine bonds_sma
