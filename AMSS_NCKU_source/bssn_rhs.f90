

#include "macrodef.fh"

  function compute_rhs_bssn(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0               
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Aupxx,Aupxy,Aupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Aupyy,Aupyz,Aupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: d2bx_xx, d2bx_xy, d2bx_xz, d2bx_yy, d2bx_yz, d2bx_zz
  real*8, dimension(ex(1),ex(2),ex(3)) :: d2by_xx, d2by_xy, d2by_xz, d2by_yy, d2by_yz, d2by_zz
  real*8, dimension(ex(1),ex(2),ex(3)) :: d2bz_xx, d2bz_xy, d2bz_xz, d2bz_yy, d2bz_yz, d2bz_zz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chixx, chixy, chixz, chiyy, chiyz, chizz
  real*8, dimension(ex(1),ex(2),ex(3)) :: lapxx, lapxy, lapxz, lapyy, lapyz, lapzz

  integer :: i, j, k
  real*8 :: L_gxx, L_gxy, L_gxz, L_gyy, L_gyz, L_gzz
  real*8 :: L_Axx, L_Axy, L_Axz, L_Ayy, L_Ayz, L_Azz
  real*8 :: L_gupxx, L_gupxy, L_gupxz, L_gupyy, L_gupyz, L_gupzz
  real*8 :: L_Gamxa, L_Gamya, L_Gamza
  real*8 :: L_gxxx, L_gxyx, L_gxzx, L_gyyx, L_gyzx, L_gzzx
  real*8 :: L_gxxy, L_gxyy, L_gxzy, L_gyyy, L_gyzy, L_gzzy
  real*8 :: L_gxxz, L_gxyz, L_gxzz, L_gyyz, L_gyzz, L_gzzz
  real*8 :: L_Rxxt, L_Rxyt, L_Rxzt, L_Ryyt, L_Ryzt, L_Rzzt
  real*8 :: L_Rxx, L_Rxy, L_Rxz, L_Ryy, L_Ryz, L_Rzz
  real*8 :: L_chix, L_chiy, L_chiz
  real*8 :: L_Dxx_chi, L_Dxy_chi, L_Dxz_chi, L_Dyy_chi, L_Dyz_chi, L_Dzz_chi
  real*8 :: L_gxxx_p, L_gxxy_p, L_gxxz_p
  real*8 :: L_Gamxxx_p, L_Gamyxx_p, L_Gamzxx_p, L_Gamxyy_p, L_Gamyyy_p, L_Gamzyy_p
  real*8 :: L_Gamxzz_p, L_Gamyzz_p, L_Gamzzz_p, L_Gamxxy_p, L_Gamyxy_p, L_Gamzxy_p
  real*8 :: L_Gamxxz_p, L_Gamyxz_p, L_Gamzxz_p, L_Gamxyz_p, L_Gamyyz_p, L_Gamzyz_p
  real*8 :: L_lapx, L_lapy, L_lapz
  real*8 :: L_lapxx, L_lapxy, L_lapxz, L_lapyy, L_lapyz, L_lapzz
  real*8 :: L_lap_lap, L_S, L_f_trace, L_tf
  real*8 :: L_tmp_xx, L_tmp_xy, L_tmp_xz, L_tmp_yy, L_tmp_yz, L_tmp_zz
  real*8 :: L_div_beta, L_alpn1, L_chin1, L_f

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta
#endif

#if (GAUGE == 6 || GAUGE == 7)
  integer :: BHN,i,j,k
  real*8, dimension(9) :: Porg
  real*8, dimension(3) :: Mass
  real*8 :: r1,r2,M,A,w1,w2,C1,C2
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta

  call getpbh(BHN,Porg,Mass)
#endif

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

  call fdderivs(ex, betax, d2bx_xx, d2bx_xy, d2bx_xz, d2bx_yy, d2bx_yz, d2bx_zz, &
              X, Y, Z, ANTI, SYM, SYM, Symmetry, Lev)
  call fdderivs(ex, betay, d2by_xx, d2by_xy, d2by_xz, d2by_yy, d2by_yz, d2by_zz, &
              X, Y, Z, SYM, ANTI, SYM, Symmetry, Lev)
  call fdderivs(ex, betaz, d2bz_xx, d2bz_xy, d2bz_xz, d2bz_yy, d2bz_yz, d2bz_zz, &
              X, Y, Z, SYM, SYM, ANTI, Symmetry, Lev)
  call fdderivs(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fdderivs(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fdderivs(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fdderivs(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,X,Y,Z,ANTI,ANTI,SYM,Symmetry,Lev)
  call fdderivs(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,X,Y,Z,ANTI,SYM,ANTI,Symmetry,Lev)
  call fdderivs(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,X,Y,Z,SYM,ANTI,ANTI,Symmetry,Lev)
  call fdderivs(ex,chi,chixx,chixy,chixz,chiyy,chiyz,chizz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fdderivs(ex,Lap,lapxx,lapxy,lapxz,lapyy,lapyz,lapzz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  alpn1 = Lap + ONE
  chin1 = chi + ONE

  div_beta = betaxx + betayy + betazz

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta )

  gxx_rhs = - TWO * alpn1 * Axx    -  F2o3 * gxx * div_beta          + &
              TWO *(  gxx * betaxx +   gxy * betayx +   gxz * betazx)
  gyy_rhs = - TWO * alpn1 * Ayy    -  F2o3 * gyy * div_beta          + &
              TWO *(  gxy * betaxy +   gyy * betayy +   gyz * betazy)
  gzz_rhs = - TWO * alpn1 * Azz    -  F2o3 * gzz * div_beta          + &
              TWO *(  gxz * betaxz +   gyz * betayz +   gzz * betazz)
  gxy_rhs = - TWO * alpn1 * Axy    +  F1o3 * gxy    * div_beta       + &
                      gxx * betaxy                  +   gxz * betazy + &
                                       gyy * betayx +   gyz * betazx   &
                                                    -   gxy * betazz
  gyz_rhs = - TWO * alpn1 * Ayz    +  F1o3 * gyz    * div_beta       + &
                      gxy * betaxz +   gyy * betayz                  + &
                      gxz * betaxy                  +   gzz * betazy   &
                                                    -   gyz * betaxx
  gxz_rhs = - TWO * alpn1 * Axz    +  F1o3 * gxz    * div_beta       + &
                      gxx * betaxz +   gxy * betayz                  + &
                                       gyz * betayx +   gzz * betazx   &
                                                    -   gxz * betayy

  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gamx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                   +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                   +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                   +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gamy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                   +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                   +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                   +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gamz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                   +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                   +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                   +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                   +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                   +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                   +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  endif

  Gamxxx =HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx*(TWO*gxzz - gzzx ) + gupxy*(TWO*gyzz - gzzy ) + gupxz*gzzz)
  Gamyzz =HALF*( gupxy*(TWO*gxzz - gzzx ) + gupyy*(TWO*gyzz - gzzy ) + gupyz*gzzz)
  Gamzzz =HALF*( gupxz*(TWO*gxzz - gzzx ) + gupyz*(TWO*gyzz - gzzy ) + gupzz*gzzz)

  Gamxxy =HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx*gxxz + gupxy*( gxyz + gyzx - gxzy ) + gupxz*gzzx )
  Gamyxz =HALF*( gupxy*gxxz + gupyy*( gxyz + gyzx - gxzy ) + gupyz*gzzx )
  Gamzxz =HALF*( gupxz*gxxz + gupyz*( gxyz + gyzx - gxzy ) + gupzz*gzzx )

  Gamxyz =HALF*( gupxx*( gxyz + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy )
  Gamyyz =HALF*( gupxy*( gxyz + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy )
  Gamzyz =HALF*( gupxz*( gxyz + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy )

  Aupxx =    gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

  Aupyy =    gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

  Aupzz =    gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

  Aupxy =    gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
          (gupxx * gupyy       + gupxy * gupxy)* Axy                       + &
          (gupxx * gupyz       + gupxz * gupxy)* Axz                       + &
          (gupxy * gupyz       + gupxz * gupyy)* Ayz

  Aupxz =    gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
          (gupxx * gupyz       + gupxy * gupxz)* Axy                       + &
          (gupxx * gupzz       + gupxz * gupxz)* Axz                       + &
          (gupxy * gupzz       + gupxz * gupyz)* Ayz

  Aupyz =    gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
          (gupxy * gupyz       + gupyy * gupxz)* Axy                       + &
          (gupxy * gupzz       + gupyz * gupxz)* Axz                       + &
          (gupyy * gupzz       + gupyz * gupyz)* Ayz
        
  Gamx_rhs = - TWO * (   Lapx * Aupxx +   Lapy * Aupxy +   Lapz * Aupxz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Aupxx +   chiy * Aupxy +   chiz * Aupxz ) - &
              gupxx * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupxy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupxz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamxxx * Aupxx + Gamxyy * Aupyy + Gamxzz * Aupzz   + &
                TWO * ( Gamxxy * Aupxy + Gamxxz * Aupxz + Gamxyz * Aupyz ) )

   Gamy_rhs = - TWO * (   Lapx * Aupxy +   Lapy * Aupyy +   Lapz * Aupyz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Aupxy +  chiy * Aupyy +    chiz * Aupyz ) - &
              gupxy * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupyz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamyxx * Aupxx + Gamyyy * Aupyy + Gamyzz * Aupzz   + &
                TWO * ( Gamyxy * Aupxy + Gamyxz * Aupxz + Gamyyz * Aupyz ) )

   Gamz_rhs = - TWO * (   Lapx * Aupxz +   Lapy * Aupyz +   Lapz * Aupzz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Aupxz +  chiy * Aupyz +    chiz * Aupzz ) - &
              gupxz * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyz * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupzz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamzxx * Aupxx + Gamzyy * Aupyy + Gamzzz * Aupzz   + &
                TWO * ( Gamzxy * Aupxy + Gamzxz * Aupxz + Gamzyz * Aupyz ) )

  fxx = d2bx_xx + d2by_xy + d2bz_xz
  fxy = d2bx_xy + d2by_yy + d2bz_yz
  fxz = d2bx_xz + d2by_yz + d2bz_zz

  Gamxa =       gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz + &
          TWO*( gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz )
  Gamya =       gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz + &
          TWO*( gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz )
  Gamza =       gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz + &
          TWO*( gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz )

  Gamx_rhs =               Gamx_rhs +  F2o3 *  Gamxa * div_beta        - &
                     Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz  + &
             F1o3 * (gupxx * fxx    + gupxy * fxy    + gupxz * fxz    ) + &
                     gupxx * d2bx_xx   + gupyy * d2bx_yy   + gupzz * d2bx_zz    + &
              TWO * (gupxy * d2bx_xy   + gupxz * d2bx_xz   + gupyz * d2bx_yz  )

  Gamy_rhs =               Gamy_rhs +  F2o3 *  Gamya * div_beta        - &
                     Gamxa * betayx - Gamya * betayy - Gamza * betayz  + &
             F1o3 * (gupxy * fxx    + gupyy * fxy    + gupyz * fxz    ) + &
                     gupxx * d2by_xx   + gupyy * d2by_yy   + gupzz * d2by_zz    + &
              TWO * (gupxy * d2by_xy   + gupxz * d2by_xz   + gupyz * d2by_yz  )

  Gamz_rhs =               Gamz_rhs +  F2o3 *  Gamza * div_beta        - &
                     Gamxa * betazx - Gamya * betazy - Gamza * betazz  + &
             F1o3 * (gupxz * fxx    + gupyz * fxy    + gupzz * fxz    ) + &
                     gupxx * d2bz_xx   + gupyy * d2bz_yy   + gupzz * d2bz_zz    + &
              TWO * (gupxy * d2bz_xy   + gupxz * d2bz_xz   + gupyz * d2bz_yz  )

  !$OMP PARALLEL DO PRIVATE(i, j, k, &
  !$OMP& L_gxx, L_gxy, L_gxz, L_gyy, L_gyz, L_gzz, &
  !$OMP& L_Axx, L_Axy, L_Axz, L_Ayy, L_Ayz, L_Azz, &
  !$OMP& L_gupxx, L_gupxy, L_gupxz, L_gupyy, L_gupyz, L_gupzz, &
  !$OMP& L_Gamxa, L_Gamya, L_Gamza, &
  !$OMP& L_gxxx, L_gxyx, L_gxzx, L_gyyx, L_gyzx, L_gzzx, &
  !$OMP& L_gxxy, L_gxyy, L_gxzy, L_gyyy, L_gyzy, L_gzzy, &
  !$OMP& L_gxxz, L_gxyz, L_gxzz, L_gyyz, L_gyzz, L_gzzz, &
  !$OMP& L_Rxxt, L_Rxyt, L_Rxzt, L_Ryyt, L_Ryzt, L_Rzzt, &
  !$OMP& L_Rxx, L_Rxy, L_Rxz, L_Ryy, L_Ryz, L_Rzz, &
  !$OMP& L_chix, L_chiy, L_chiz, &
  !$OMP& L_Dxx_chi, L_Dxy_chi, L_Dxz_chi, L_Dyy_chi, L_Dyz_chi, L_Dzz_chi, &
  !$OMP& L_gxxx_p, L_gxxy_p, L_gxxz_p, &
  !$OMP& L_Gamxxx_p, L_Gamyxx_p, L_Gamzxx_p, L_Gamxyy_p, L_Gamyyy_p, L_Gamzyy_p, &
  !$OMP& L_Gamxzz_p, L_Gamyzz_p, L_Gamzzz_p, L_Gamxxy_p, L_Gamyxy_p, L_Gamzxy_p, &
  !$OMP& L_Gamxxz_p, L_Gamyxz_p, L_Gamzxz_p, L_Gamxyz_p, L_Gamyyz_p, L_Gamzyz_p, &
  !$OMP& L_lapx, L_lapy, L_lapz, &
  !$OMP& L_lapxx, L_lapxy, L_lapxz, L_lapyy, L_lapyz, L_lapzz, &
  !$OMP& L_lap_lap, L_S, L_f_trace, L_tf, &
  !$OMP& L_tmp_xx, L_tmp_xy, L_tmp_xz, L_tmp_yy, L_tmp_yz, L_tmp_zz, &
  !$OMP& L_div_beta, L_alpn1, L_chin1, L_f) COLLAPSE(2)
  do k = 1, ex(3)
    do j = 1, ex(2)
      !$OMP SIMD
      do i = 1, ex(1)
        L_chin1 = chin1(i,j,k)
        L_alpn1 = alpn1(i,j,k)
        L_div_beta = div_beta(i,j,k)

        L_gxx = gxx(i,j,k); L_gxy = gxy(i,j,k); L_gxz = gxz(i,j,k)
        L_gyy = gyy(i,j,k); L_gyz = gyz(i,j,k); L_gzz = gzz(i,j,k)

        L_Axx = Axx(i,j,k); L_Axy = Axy(i,j,k); L_Axz = Axz(i,j,k)
        L_Ayy = Ayy(i,j,k); L_Ayz = Ayz(i,j,k); L_Azz = Azz(i,j,k)

        L_gupxx = gupxx(i,j,k)
        L_gupxy = gupxy(i,j,k)
        L_gupxz = gupxz(i,j,k)
        L_gupyy = gupyy(i,j,k)
        L_gupyz = gupyz(i,j,k)
        L_gupzz = gupzz(i,j,k)

        L_Gamxa = Gamxa(i,j,k)
        L_Gamya = Gamya(i,j,k)
        L_Gamza = Gamza(i,j,k)

        L_gxxx = L_gxx * Gamxxx(i,j,k) + L_gxy * Gamyxx(i,j,k) + L_gxz * Gamzxx(i,j,k)
        L_gxyx = L_gxx * Gamxxy(i,j,k) + L_gxy * Gamyxy(i,j,k) + L_gxz * Gamzxy(i,j,k)
        L_gxzx = L_gxx * Gamxxz(i,j,k) + L_gxy * Gamyxz(i,j,k) + L_gxz * Gamzxz(i,j,k)
        L_gyyx = L_gxx * Gamxyy(i,j,k) + L_gxy * Gamyyy(i,j,k) + L_gxz * Gamzyy(i,j,k)
        L_gyzx = L_gxx * Gamxyz(i,j,k) + L_gxy * Gamyyz(i,j,k) + L_gxz * Gamzyz(i,j,k)
        L_gzzx = L_gxx * Gamxzz(i,j,k) + L_gxy * Gamyzz(i,j,k) + L_gxz * Gamzzz(i,j,k)

        L_gxxy = L_gxy * Gamxxx(i,j,k) + L_gyy * Gamyxx(i,j,k) + L_gyz * Gamzxx(i,j,k)
        L_gxyy = L_gxy * Gamxxy(i,j,k) + L_gyy * Gamyxy(i,j,k) + L_gyz * Gamzxy(i,j,k)
        L_gxzy = L_gxy * Gamxxz(i,j,k) + L_gyy * Gamyxz(i,j,k) + L_gyz * Gamzxz(i,j,k)
        L_gyyy = L_gxy * Gamxyy(i,j,k) + L_gyy * Gamyyy(i,j,k) + L_gyz * Gamzyy(i,j,k)
        L_gyzy = L_gxy * Gamxyz(i,j,k) + L_gyy * Gamyyz(i,j,k) + L_gyz * Gamzyz(i,j,k)
        L_gzzy = L_gxy * Gamxzz(i,j,k) + L_gyy * Gamyzz(i,j,k) + L_gyz * Gamzzz(i,j,k)

        L_gxxz = L_gxz * Gamxxx(i,j,k) + L_gyz * Gamyxx(i,j,k) + L_gzz * Gamzxx(i,j,k)
        L_gxyz = L_gxz * Gamxxy(i,j,k) + L_gyz * Gamyxy(i,j,k) + L_gzz * Gamzxy(i,j,k)
        L_gxzz = L_gxz * Gamxxz(i,j,k) + L_gyz * Gamyxz(i,j,k) + L_gzz * Gamzxz(i,j,k)
        L_gyyz = L_gxz * Gamxyy(i,j,k) + L_gyz * Gamyyy(i,j,k) + L_gzz * Gamzyy(i,j,k)
        L_gyzz = L_gxz * Gamxyz(i,j,k) + L_gyz * Gamyyz(i,j,k) + L_gzz * Gamzyz(i,j,k)
        L_gzzz = L_gxz * Gamxzz(i,j,k) + L_gyz * Gamyzz(i,j,k) + L_gzz * Gamzzz(i,j,k)

        L_Rxxt = L_gupxx*gxxxx(i,j,k) + L_gupyy*gxxyy(i,j,k) + L_gupzz*gxxzz(i,j,k) + &
                TWO*(L_gupxy*gxxxy(i,j,k) + L_gupxz*gxxxz(i,j,k) + L_gupyz*gxxyz(i,j,k))
        L_Ryyt = L_gupxx*gyyxx(i,j,k) + L_gupyy*gyyyy(i,j,k) + L_gupzz*gyyzz(i,j,k) + &
                TWO*(L_gupxy*gyyxy(i,j,k) + L_gupxz*gyyxz(i,j,k) + L_gupyz*gyyyz(i,j,k))
        L_Rzzt = L_gupxx*gzzxx(i,j,k) + L_gupyy*gzzyy(i,j,k) + L_gupzz*gzzzz(i,j,k) + &
                TWO*(L_gupxy*gzzxy(i,j,k) + L_gupxz*gzzxz(i,j,k) + L_gupyz*gzzyz(i,j,k))
        L_Rxyt = L_gupxx*gxyxx(i,j,k) + L_gupyy*gxyyy(i,j,k) + L_gupzz*gxyzz(i,j,k) + &
                TWO*(L_gupxy*gxyxy(i,j,k) + L_gupxz*gxyxz(i,j,k) + L_gupyz*gxyyz(i,j,k))
        L_Rxzt = L_gupxx*gxzxx(i,j,k) + L_gupyy*gxzyy(i,j,k) + L_gupzz*gxzzz(i,j,k) + &
                TWO*(L_gupxy*gxzxy(i,j,k) + L_gupxz*gxzxz(i,j,k) + L_gupyz*gxzyz(i,j,k))
        L_Ryzt = L_gupxx*gyzxx(i,j,k) + L_gupyy*gyzyy(i,j,k) + L_gupzz*gyzzz(i,j,k) + &
                TWO*(L_gupxy*gyzxy(i,j,k) + L_gupxz*gyzxz(i,j,k) + L_gupyz*gyzyz(i,j,k))

        L_Rxx = - HALF * L_Rxxt + &
                L_gxx * Gamxx(i,j,k) + L_gxy * Gamyx(i,j,k) + L_gxz * Gamzx(i,j,k) + &
                L_Gamxa * L_gxxx + L_Gamya * L_gxyx + L_Gamza * L_gxzx + &
                L_gupxx * ( &
                  TWO * (Gamxxx(i,j,k) * L_gxxx + Gamyxx(i,j,k) * L_gxyx + Gamzxx(i,j,k) * L_gxzx) + &
                  Gamxxx(i,j,k) * L_gxxx + Gamyxx(i,j,k) * L_gxxy + Gamzxx(i,j,k) * L_gxxz) + &
                L_gupyy * ( &
                  TWO * (Gamxxy(i,j,k) * L_gxyx + Gamyxy(i,j,k) * L_gyyx + Gamzxy(i,j,k) * L_gyzx) + &
                  Gamxxy(i,j,k) * L_gxyx + Gamyxy(i,j,k) * L_gxyy + Gamzxy(i,j,k) * L_gxyz) + &
                L_gupzz * ( &
                  TWO * (Gamxxz(i,j,k) * L_gxzx + Gamyxz(i,j,k) * L_gyzx + Gamzxz(i,j,k) * L_gzzx) + &
                  Gamxxz(i,j,k) * L_gxzx + Gamyxz(i,j,k) * L_gxzy + Gamzxz(i,j,k) * L_gxzz) + &
                L_gupxy * ( &
                  TWO * (Gamxxx(i,j,k) * L_gxyx + Gamyxx(i,j,k) * L_gyyx + Gamzxx(i,j,k) * L_gyzx + &
                  Gamxxy(i,j,k) * L_gxxx + Gamyxy(i,j,k) * L_gxyx + Gamzxy(i,j,k) * L_gxzx) + &
                  Gamxxy(i,j,k) * L_gxxx + Gamyxy(i,j,k) * L_gxxy + Gamzxy(i,j,k) * L_gxxz + &
                  Gamxxx(i,j,k) * L_gxyx + Gamyxx(i,j,k) * L_gxyy + Gamzxx(i,j,k) * L_gxyz) + &
                L_gupxz * ( &
                  TWO * (Gamxxx(i,j,k) * L_gxzx + Gamyxx(i,j,k) * L_gyzx + Gamzxx(i,j,k) * L_gzzx + &
                  Gamxxz(i,j,k) * L_gxxx + Gamyxz(i,j,k) * L_gxyx + Gamzxz(i,j,k) * L_gxzx) + &
                  Gamxxz(i,j,k) * L_gxxx + Gamyxz(i,j,k) * L_gxxy + Gamzxz(i,j,k) * L_gxxz + &
                  Gamxxx(i,j,k) * L_gxzx + Gamyxx(i,j,k) * L_gxzy + Gamzxx(i,j,k) * L_gxzz) + &
                L_gupyz * ( &
                  TWO * (Gamxxy(i,j,k) * L_gxzx + Gamyxy(i,j,k) * L_gyzx + Gamzxy(i,j,k) * L_gzzx + &
                  Gamxxz(i,j,k) * L_gxyx + Gamyxz(i,j,k) * L_gyyx + Gamzxz(i,j,k) * L_gyzx) + &
                  Gamxxz(i,j,k) * L_gxyx + Gamyxz(i,j,k) * L_gxyy + Gamzxz(i,j,k) * L_gxyz + &
                  Gamxxy(i,j,k) * L_gxzx + Gamyxy(i,j,k) * L_gxzy + Gamzxy(i,j,k) * L_gxzz)

        L_Ryy = - HALF * L_Ryyt + &
                L_gxy * Gamxy(i,j,k) + L_gyy * Gamyy(i,j,k) + L_gyz * Gamzy(i,j,k) + &
                L_Gamxa * L_gxyy + L_Gamya * L_gyyy + L_Gamza * L_gyzy + &
                L_gupxx * ( &
                  TWO * (Gamxxy(i,j,k) * L_gxxy + Gamyxy(i,j,k) * L_gxyy + Gamzxy(i,j,k) * L_gxzy) + &
                  Gamxxy(i,j,k) * L_gxyx + Gamyxy(i,j,k) * L_gxyy + Gamzxy(i,j,k) * L_gxyz) + &
                L_gupyy * ( &
                  TWO * (Gamxyy(i,j,k) * L_gxyy + Gamyyy(i,j,k) * L_gyyy + Gamzyy(i,j,k) * L_gyzy) + &
                  Gamxyy(i,j,k) * L_gyyx + Gamyyy(i,j,k) * L_gyyy + Gamzyy(i,j,k) * L_gyyz) + &
                L_gupzz * ( &
                  TWO * (Gamxyz(i,j,k) * L_gxzy + Gamyyz(i,j,k) * L_gyzy + Gamzyz(i,j,k) * L_gzzy) + &
                  Gamxyz(i,j,k) * L_gyzx + Gamyyz(i,j,k) * L_gyzy + Gamzyz(i,j,k) * L_gyzz) + &
                L_gupxy * ( &
                  TWO * (Gamxxy(i,j,k) * L_gxyy + Gamyxy(i,j,k) * L_gyyy + Gamzxy(i,j,k) * L_gyzy + &
                        Gamxyy(i,j,k) * L_gxxy + Gamyyy(i,j,k) * L_gxyy + Gamzyy(i,j,k) * L_gxzy) + &
                  Gamxyy(i,j,k) * L_gxyx + Gamyyy(i,j,k) * L_gxyy + Gamzyy(i,j,k) * L_gxyz + &
                  Gamxxy(i,j,k) * L_gyyx + Gamyxy(i,j,k) * L_gyyy + Gamzxy(i,j,k) * L_gyyz) + &
                L_gupxz * ( &
                  TWO * (Gamxxy(i,j,k) * L_gxzy + Gamyxy(i,j,k) * L_gyzy + Gamzxy(i,j,k) * L_gzzy + &
                        Gamxyz(i,j,k) * L_gxxy + Gamyyz(i,j,k) * L_gxyy + Gamzyz(i,j,k) * L_gxzy) + &
                  Gamxyz(i,j,k) * L_gxyx + Gamyyz(i,j,k) * L_gxyy + Gamzyz(i,j,k) * L_gxyz + &
                  Gamxxy(i,j,k) * L_gyzx + Gamyxy(i,j,k) * L_gyzy + Gamzxy(i,j,k) * L_gyzz) + &
                L_gupyz * ( &
                  TWO * (Gamxyy(i,j,k) * L_gxzy + Gamyyy(i,j,k) * L_gyzy + Gamzyy(i,j,k) * L_gzzy + &
                        Gamxyz(i,j,k) * L_gxyy + Gamyyz(i,j,k) * L_gyyy + Gamzyz(i,j,k) * L_gyzy) + &
                  Gamxyz(i,j,k) * L_gyyx + Gamyyz(i,j,k) * L_gyyy + Gamzyz(i,j,k) * L_gyyz + &
                  Gamxyy(i,j,k) * L_gyzx + Gamyyy(i,j,k) * L_gyzy + Gamzyy(i,j,k) * L_gyzz)

        L_Rzz = - HALF * L_Rzzt + &
                L_gxz * Gamxz(i,j,k) + L_gyz * Gamyz(i,j,k) + L_gzz * Gamzz(i,j,k) + &
                L_Gamxa * L_gxzz + L_Gamya * L_gyzz + L_Gamza * L_gzzz + &
                L_gupxx * ( &
                  TWO * (Gamxxz(i,j,k) * L_gxxz + Gamyxz(i,j,k) * L_gxyz + Gamzxz(i,j,k) * L_gxzz) + &
                  Gamxxz(i,j,k) * L_gxzx + Gamyxz(i,j,k) * L_gxzy + Gamzxz(i,j,k) * L_gxzz) + &
                L_gupyy * ( &
                  TWO * (Gamxyz(i,j,k) * L_gxyz + Gamyyz(i,j,k) * L_gyyz + Gamzyz(i,j,k) * L_gyzz) + &
                  Gamxyz(i,j,k) * L_gyzx + Gamyyz(i,j,k) * L_gyzy + Gamzyz(i,j,k) * L_gyzz) + &
                L_gupzz * ( &
                  TWO * (Gamxzz(i,j,k) * L_gxzz + Gamyzz(i,j,k) * L_gyzz + Gamzzz(i,j,k) * L_gzzz) + &
                  Gamxzz(i,j,k) * L_gzzx + Gamyzz(i,j,k) * L_gzzy + Gamzzz(i,j,k) * L_gzzz) + &
                L_gupxy * ( &
                  TWO * (Gamxxz(i,j,k) * L_gxyz + Gamyxz(i,j,k) * L_gyyz + Gamzxz(i,j,k) * L_gyzz + &
                        Gamxyz(i,j,k) * L_gxxz + Gamyyz(i,j,k) * L_gxyz + Gamzyz(i,j,k) * L_gxzz) + &
                  Gamxyz(i,j,k) * L_gxzx + Gamyyz(i,j,k) * L_gxzy + Gamzyz(i,j,k) * L_gxzz + &
                  Gamxxz(i,j,k) * L_gyzx + Gamyxz(i,j,k) * L_gyzy + Gamzxz(i,j,k) * L_gyzz) + &
                L_gupxz * ( &
                  TWO * (Gamxxz(i,j,k) * L_gxzz + Gamyxz(i,j,k) * L_gyzz + Gamzxz(i,j,k) * L_gzzz + &
                        Gamxzz(i,j,k) * L_gxxz + Gamyzz(i,j,k) * L_gxyz + Gamzzz(i,j,k) * L_gxzz) + &
                  Gamxzz(i,j,k) * L_gxzx + Gamyzz(i,j,k) * L_gxzy + Gamzzz(i,j,k) * L_gxzz + &
                  Gamxxz(i,j,k) * L_gzzx + Gamyxz(i,j,k) * L_gzzy + Gamzxz(i,j,k) * L_gzzz) + &
                L_gupyz * ( &
                  TWO * (Gamxyz(i,j,k) * L_gxzz + Gamyyz(i,j,k) * L_gyzz + Gamzyz(i,j,k) * L_gzzz + &
                        Gamxzz(i,j,k) * L_gxyz + Gamyzz(i,j,k) * L_gyyz + Gamzzz(i,j,k) * L_gyzz) + &
                  Gamxzz(i,j,k) * L_gyzx + Gamyzz(i,j,k) * L_gyzy + Gamzzz(i,j,k) * L_gyzz + &
                  Gamxyz(i,j,k) * L_gzzx + Gamyyz(i,j,k) * L_gzzy + Gamzyz(i,j,k) * L_gzzz)
        
        L_Rxy = HALF * ( - L_Rxyt + &
                L_gxx * Gamxy(i,j,k) + L_gxy * Gamyy(i,j,k) + L_gxz * Gamzy(i,j,k) + &
                L_gxy * Gamxx(i,j,k) + L_gyy * Gamyx(i,j,k) + L_gyz * Gamzx(i,j,k) + &
                L_Gamxa * L_gxyx + L_Gamya * L_gyyx + L_Gamza * L_gyzx + &
                L_Gamxa * L_gxxy + L_Gamya * L_gxyy + L_Gamza * L_gxzy ) + &
           L_gupxx * ( &
                Gamxxx(i,j,k) * L_gxxy + Gamyxx(i,j,k) * L_gxyy + Gamzxx(i,j,k) * L_gxzy + &
                Gamxxy(i,j,k) * L_gxxx + Gamyxy(i,j,k) * L_gxyx + Gamzxy(i,j,k) * L_gxzx + &
                Gamxxx(i,j,k) * L_gxyx + Gamyxx(i,j,k) * L_gxyy + Gamzxx(i,j,k) * L_gxyz ) + &
           L_gupxy * ( &
                Gamxxx(i,j,k) * L_gxyy + Gamyxx(i,j,k) * L_gyyy + Gamzxx(i,j,k) * L_gyzy + &
                Gamxxy(i,j,k) * L_gxyx + Gamyxy(i,j,k) * L_gyyx + Gamzxy(i,j,k) * L_gyzx + &
                Gamxxy(i,j,k) * L_gxyx + Gamyxy(i,j,k) * L_gxyy + Gamzxy(i,j,k) * L_gxyz + &
                Gamxxy(i,j,k) * L_gxxy + Gamyxy(i,j,k) * L_gxyy + Gamzxy(i,j,k) * L_gxzy + &
                Gamxyy(i,j,k) * L_gxxx + Gamyyy(i,j,k) * L_gxyx + Gamzyy(i,j,k) * L_gxzx + &
                Gamxxx(i,j,k) * L_gyyx + Gamyxx(i,j,k) * L_gyyy + Gamzxx(i,j,k) * L_gyyz ) + &
           L_gupxz * ( &
                Gamxxx(i,j,k) * L_gxzy + Gamyxx(i,j,k) * L_gyzy + Gamzxx(i,j,k) * L_gzzy + &
                Gamxxy(i,j,k) * L_gxzx + Gamyxy(i,j,k) * L_gyzx + Gamzxy(i,j,k) * L_gzzx + &
                Gamxxz(i,j,k) * L_gxyx + Gamyxz(i,j,k) * L_gxyy + Gamzxz(i,j,k) * L_gxyz + &
                Gamxxz(i,j,k) * L_gxxy + Gamyxz(i,j,k) * L_gxyy + Gamzxz(i,j,k) * L_gxzy + &
                Gamxyz(i,j,k) * L_gxxx + Gamyyz(i,j,k) * L_gxyx + Gamzyz(i,j,k) * L_gxzx + &
                Gamxxx(i,j,k) * L_gyzx + Gamyxx(i,j,k) * L_gyzy + Gamzxx(i,j,k) * L_gyzz ) + &
           L_gupyy * ( &
                Gamxxy(i,j,k) * L_gxyy + Gamyxy(i,j,k) * L_gyyy + Gamzxy(i,j,k) * L_gyzy + &
                Gamxyy(i,j,k) * L_gxyx + Gamyyy(i,j,k) * L_gyyx + Gamzyy(i,j,k) * L_gyzx + &
                Gamxxy(i,j,k) * L_gyyx + Gamyxy(i,j,k) * L_gyyy + Gamzxy(i,j,k) * L_gyyz ) + &
           L_gupyz * ( &
                Gamxxy(i,j,k) * L_gxzy + Gamyxy(i,j,k) * L_gyzy + Gamzxy(i,j,k) * L_gzzy + &
                Gamxyy(i,j,k) * L_gxzx + Gamyyy(i,j,k) * L_gyzx + Gamzyy(i,j,k) * L_gzzx + &
                Gamxxz(i,j,k) * L_gyyx + Gamyxz(i,j,k) * L_gyyy + Gamzxz(i,j,k) * L_gyyz + &
                Gamxxz(i,j,k) * L_gxyy + Gamyxz(i,j,k) * L_gyyy + Gamzxz(i,j,k) * L_gyzy + &
                Gamxyz(i,j,k) * L_gxyx + Gamyyz(i,j,k) * L_gyyx + Gamzyz(i,j,k) * L_gyzx + &
                Gamxxy(i,j,k) * L_gyzx + Gamyxy(i,j,k) * L_gyzy + Gamzxy(i,j,k) * L_gyzz ) + &
           L_gupzz * ( &
                Gamxxz(i,j,k) * L_gxzy + Gamyxz(i,j,k) * L_gyzy + Gamzxz(i,j,k) * L_gzzy + &
                Gamxyz(i,j,k) * L_gxzx + Gamyyz(i,j,k) * L_gyzx + Gamzyz(i,j,k) * L_gzzx + &
                Gamxxz(i,j,k) * L_gyzx + Gamyxz(i,j,k) * L_gyzy + Gamzxz(i,j,k) * L_gyzz )
          
        L_Rxz = HALF * ( - L_Rxzt + &
                L_gxx * Gamxz(i,j,k) + L_gxy * Gamyz(i,j,k) + L_gxz * Gamzz(i,j,k) + &
                L_gxz * Gamxx(i,j,k) + L_gyz * Gamyx(i,j,k) + L_gzz * Gamzx(i,j,k) + &
                L_Gamxa * L_gxzx + L_Gamya * L_gyzx + L_Gamza * L_gzzx + &
                L_Gamxa * L_gxxz + L_Gamya * L_gxyz + L_Gamza * L_gxzz ) + &
           L_gupxx * ( &
                Gamxxx(i,j,k) * L_gxxz + Gamyxx(i,j,k) * L_gxyz + Gamzxx(i,j,k) * L_gxzz + &
                Gamxxz(i,j,k) * L_gxxx + Gamyxz(i,j,k) * L_gxyx + Gamzxz(i,j,k) * L_gxzx + &
                Gamxxx(i,j,k) * L_gxzx + Gamyxx(i,j,k) * L_gxzy + Gamzxx(i,j,k) * L_gxzz ) + &
           L_gupxy * ( &
                Gamxxx(i,j,k) * L_gxyz + Gamyxx(i,j,k) * L_gyyz + Gamzxx(i,j,k) * L_gyzz + &
                Gamxxz(i,j,k) * L_gxyx + Gamyxz(i,j,k) * L_gyyx + Gamzxz(i,j,k) * L_gyzx + &
                Gamxxy(i,j,k) * L_gxzx + Gamyxy(i,j,k) * L_gxzy + Gamzxy(i,j,k) * L_gxzz + &
                Gamxxy(i,j,k) * L_gxxz + Gamyxy(i,j,k) * L_gxyz + Gamzxy(i,j,k) * L_gxzz + &
                Gamxyz(i,j,k) * L_gxxx + Gamyyz(i,j,k) * L_gxyx + Gamzyz(i,j,k) * L_gxzx + &
                Gamxxx(i,j,k) * L_gyzx + Gamyxx(i,j,k) * L_gyzy + Gamzxx(i,j,k) * L_gyzz ) + &
           L_gupxz * ( &
                Gamxxx(i,j,k) * L_gxzz + Gamyxx(i,j,k) * L_gyzz + Gamzxx(i,j,k) * L_gzzz + &
                Gamxxz(i,j,k) * L_gxzx + Gamyxz(i,j,k) * L_gyzx + Gamzxz(i,j,k) * L_gzzx + &
                Gamxxz(i,j,k) * L_gxzx + Gamyxz(i,j,k) * L_gxzy + Gamzxz(i,j,k) * L_gxzz + &
                Gamxxz(i,j,k) * L_gxxz + Gamyxz(i,j,k) * L_gxyz + Gamzxz(i,j,k) * L_gxzz + &
                Gamxzz(i,j,k) * L_gxxx + Gamyzz(i,j,k) * L_gxyx + Gamzzz(i,j,k) * L_gxzx + &
                Gamxxx(i,j,k) * L_gzzx + Gamyxx(i,j,k) * L_gzzy + Gamzxx(i,j,k) * L_gzzz ) + &
           L_gupyy * ( &
                Gamxxy(i,j,k) * L_gxyz + Gamyxy(i,j,k) * L_gyyz + Gamzxy(i,j,k) * L_gyzz + &
                Gamxyz(i,j,k) * L_gxyx + Gamyyz(i,j,k) * L_gyyx + Gamzyz(i,j,k) * L_gyzx + &
                Gamxxy(i,j,k) * L_gyzx + Gamyxy(i,j,k) * L_gyzy + Gamzxy(i,j,k) * L_gyzz ) + &
           L_gupyz * ( &
                Gamxxy(i,j,k) * L_gxzz + Gamyxy(i,j,k) * L_gyzz + Gamzxy(i,j,k) * L_gzzz + &
                Gamxyz(i,j,k) * L_gxzx + Gamyyz(i,j,k) * L_gyzx + Gamzyz(i,j,k) * L_gzzx + &
                Gamxxz(i,j,k) * L_gyzx + Gamyxz(i,j,k) * L_gyzy + Gamzxz(i,j,k) * L_gyzz + &
                Gamxxz(i,j,k) * L_gxyz + Gamyxz(i,j,k) * L_gyyz + Gamzxz(i,j,k) * L_gyzz + &
                Gamxzz(i,j,k) * L_gxyx + Gamyzz(i,j,k) * L_gyyx + Gamzzz(i,j,k) * L_gyzx + &
                Gamxxy(i,j,k) * L_gzzx + Gamyxy(i,j,k) * L_gzzy + Gamzxy(i,j,k) * L_gzzz ) + &
           L_gupzz * ( &
                Gamxxz(i,j,k) * L_gxzz + Gamyxz(i,j,k) * L_gyzz + Gamzxz(i,j,k) * L_gzzz + &
                Gamxzz(i,j,k) * L_gxzx + Gamyzz(i,j,k) * L_gyzx + Gamzzz(i,j,k) * L_gzzx + &
                Gamxxz(i,j,k) * L_gzzx + Gamyxz(i,j,k) * L_gzzy + Gamzxz(i,j,k) * L_gzzz )
        
        L_Ryz = HALF * ( - L_Ryzt + &
                L_gxy * Gamxz(i,j,k) + L_gyy * Gamyz(i,j,k) + L_gyz * Gamzz(i,j,k) + &
                L_gxz * Gamxy(i,j,k) + L_gyz * Gamyy(i,j,k) + L_gzz * Gamzy(i,j,k) + &
                L_Gamxa * L_gxzy + L_Gamya * L_gyzy + L_Gamza * L_gzzy + &
                L_Gamxa * L_gxyz + L_Gamya * L_gyyz + L_Gamza * L_gyzz ) + &
           L_gupxx * ( &
                Gamxxy(i,j,k) * L_gxxz + Gamyxy(i,j,k) * L_gxyz + Gamzxy(i,j,k) * L_gxzz + &
                Gamxxz(i,j,k) * L_gxxy + Gamyxz(i,j,k) * L_gxyy + Gamzxz(i,j,k) * L_gxzy + &
                Gamxxy(i,j,k) * L_gxzx + Gamyxy(i,j,k) * L_gxzy + Gamzxy(i,j,k) * L_gxzz ) + &
           L_gupxy * ( &
                Gamxxy(i,j,k) * L_gxyz + Gamyxy(i,j,k) * L_gyyz + Gamzxy(i,j,k) * L_gyzz + &
                Gamxxz(i,j,k) * L_gxyy + Gamyxz(i,j,k) * L_gyyy + Gamzxz(i,j,k) * L_gyzy + &
                Gamxyy(i,j,k) * L_gxzx + Gamyyy(i,j,k) * L_gxzy + Gamzyy(i,j,k) * L_gxzz + &
                Gamxyy(i,j,k) * L_gxxz + Gamyyy(i,j,k) * L_gxyz + Gamzyy(i,j,k) * L_gxzz + &
                Gamxyz(i,j,k) * L_gxxy + Gamyyz(i,j,k) * L_gxyy + Gamzyz(i,j,k) * L_gxzy + &
                Gamxxy(i,j,k) * L_gyzx + Gamyxy(i,j,k) * L_gyzy + Gamzxy(i,j,k) * L_gyzz ) + &
           L_gupxz * ( &
                Gamxxy(i,j,k) * L_gxzz + Gamyxy(i,j,k) * L_gyzz + Gamzxy(i,j,k) * L_gzzz + &
                Gamxxz(i,j,k) * L_gxzy + Gamyxz(i,j,k) * L_gyzy + Gamzxz(i,j,k) * L_gzzy + &
                Gamxyz(i,j,k) * L_gxzx + Gamyyz(i,j,k) * L_gxzy + Gamzyz(i,j,k) * L_gxzz + &
                Gamxyz(i,j,k) * L_gxxz + Gamyyz(i,j,k) * L_gxyz + Gamzyz(i,j,k) * L_gxzz + &
                Gamxzz(i,j,k) * L_gxxy + Gamyzz(i,j,k) * L_gxyy + Gamzzz(i,j,k) * L_gxzy + &
                Gamxxy(i,j,k) * L_gzzx + Gamyxy(i,j,k) * L_gzzy + Gamzxy(i,j,k) * L_gzzz ) + &
           L_gupyy * ( &
                Gamxyy(i,j,k) * L_gxyz + Gamyyy(i,j,k) * L_gyyz + Gamzyy(i,j,k) * L_gyzz + &
                Gamxyz(i,j,k) * L_gxyy + Gamyyz(i,j,k) * L_gyyy + Gamzyz(i,j,k) * L_gyzy + &
                Gamxyy(i,j,k) * L_gyzx + Gamyyy(i,j,k) * L_gyzy + Gamzyy(i,j,k) * L_gyzz ) + &
           L_gupyz * ( &
                Gamxyy(i,j,k) * L_gxzz + Gamyyy(i,j,k) * L_gyzz + Gamzyy(i,j,k) * L_gzzz + &
                Gamxyz(i,j,k) * L_gxzy + Gamyyz(i,j,k) * L_gyzy + Gamzyz(i,j,k) * L_gzzy + &
                Gamxyz(i,j,k) * L_gyzx + Gamyyz(i,j,k) * L_gyzy + Gamzyz(i,j,k) * L_gyzz + &
                Gamxyz(i,j,k) * L_gxyz + Gamyyz(i,j,k) * L_gyyz + Gamzyz(i,j,k) * L_gyzz + &
                Gamxzz(i,j,k) * L_gxyy + Gamyzz(i,j,k) * L_gyyy + Gamzzz(i,j,k) * L_gyzy + &
                Gamxyy(i,j,k) * L_gzzx + Gamyyy(i,j,k) * L_gzzy + Gamzyy(i,j,k) * L_gzzz ) + &
           L_gupzz * ( &
                Gamxyz(i,j,k) * L_gxzz + Gamyyz(i,j,k) * L_gyzz + Gamzyz(i,j,k) * L_gzzz + &
                Gamxzz(i,j,k) * L_gxzy + Gamyzz(i,j,k) * L_gyzy + Gamzzz(i,j,k) * L_gzzy + &
                Gamxyz(i,j,k) * L_gzzx + Gamyyz(i,j,k) * L_gzzy + Gamzyz(i,j,k) * L_gzzz )
            
        L_chix = chix(i,j,k); L_chiy = chiy(i,j,k); L_chiz = chiz(i,j,k)
        L_Dxx_chi = chixx(i,j,k) - Gamxxx(i,j,k) * L_chix - Gamyxx(i,j,k) * L_chiy - Gamzxx(i,j,k) * L_chiz
        L_Dxy_chi = chixy(i,j,k) - Gamxxy(i,j,k) * L_chix - Gamyxy(i,j,k) * L_chiy - Gamzxy(i,j,k) * L_chiz
        L_Dxz_chi = chixz(i,j,k) - Gamxxz(i,j,k) * L_chix - Gamyxz(i,j,k) * L_chiy - Gamzxz(i,j,k) * L_chiz
        L_Dyy_chi = chiyy(i,j,k) - Gamxyy(i,j,k) * L_chix - Gamyyy(i,j,k) * L_chiy - Gamzyy(i,j,k) * L_chiz
        L_Dyz_chi = chiyz(i,j,k) - Gamxyz(i,j,k) * L_chix - Gamyyz(i,j,k) * L_chiy - Gamzyz(i,j,k) * L_chiz
        L_Dzz_chi = chizz(i,j,k) - Gamxzz(i,j,k) * L_chix - Gamyzz(i,j,k) * L_chiy - Gamzzz(i,j,k) * L_chiz

        L_f = L_gupxx * (L_Dxx_chi - F3o2/L_chin1 * L_chix * L_chix) + &
              L_gupyy * (L_Dyy_chi - F3o2/L_chin1 * L_chiy * L_chiy) + &
              L_gupzz * (L_Dzz_chi - F3o2/L_chin1 * L_chiz * L_chiz) + &
              TWO * L_gupxy * (L_Dxy_chi - F3o2/L_chin1 * L_chix * L_chiy) + &
              TWO * L_gupxz * (L_Dxz_chi - F3o2/L_chin1 * L_chix * L_chiz) + &
              TWO * L_gupyz * (L_Dyz_chi - F3o2/L_chin1 * L_chiy * L_chiz)

        Rxx(i,j,k) = L_Rxx + (L_Dxx_chi - L_chix*L_chix/L_chin1/TWO + L_gxx * L_f)/L_chin1/TWO
        Ryy(i,j,k) = L_Ryy + (L_Dyy_chi - L_chiy*L_chiy/L_chin1/TWO + L_gyy * L_f)/L_chin1/TWO
        Rzz(i,j,k) = L_Rzz + (L_Dzz_chi - L_chiz*L_chiz/L_chin1/TWO + L_gzz * L_f)/L_chin1/TWO
        Rxy(i,j,k) = L_Rxy + (L_Dxy_chi - L_chix*L_chiy/L_chin1/TWO + L_gxy * L_f)/L_chin1/TWO
        Rxz(i,j,k) = L_Rxz + (L_Dxz_chi - L_chix*L_chiz/L_chin1/TWO + L_gxz * L_f)/L_chin1/TWO
        Ryz(i,j,k) = L_Ryz + (L_Dyz_chi - L_chiy*L_chiz/L_chin1/TWO + L_gyz * L_f)/L_chin1/TWO

        L_gxxx_p = (L_gupxx * L_chix + L_gupxy * L_chiy + L_gupxz * L_chiz) / L_chin1
        L_gxxy_p = (L_gupxy * L_chix + L_gupyy * L_chiy + L_gupyz * L_chiz) / L_chin1
        L_gxxz_p = (L_gupxz * L_chix + L_gupyz * L_chiy + L_gupzz * L_chiz) / L_chin1

        L_Gamxxx_p = Gamxxx(i,j,k) - ( (L_chix + L_chix)/L_chin1 - L_gxx * L_gxxx_p ) * HALF
        L_Gamyxx_p = Gamyxx(i,j,k) - (                         - L_gxx * L_gxxy_p ) * HALF
        L_Gamzxx_p = Gamzxx(i,j,k) - (                         - L_gxx * L_gxxz_p ) * HALF
        L_Gamxyy_p = Gamxyy(i,j,k) - (                         - L_gyy * L_gxxx_p ) * HALF
        L_Gamyyy_p = Gamyyy(i,j,k) - ( (L_chiy + L_chiy)/L_chin1 - L_gyy * L_gxxy_p ) * HALF
        L_Gamzyy_p = Gamzyy(i,j,k) - (                         - L_gyy * L_gxxz_p ) * HALF
        L_Gamxzz_p = Gamxzz(i,j,k) - (                         - L_gzz * L_gxxx_p ) * HALF
        L_Gamyzz_p = Gamyzz(i,j,k) - (                         - L_gzz * L_gxxy_p ) * HALF
        L_Gamzzz_p = Gamzzz(i,j,k) - ( (L_chiz + L_chiz)/L_chin1 - L_gzz * L_gxxz_p ) * HALF
        L_Gamxxy_p = Gamxxy(i,j,k) - (  L_chiy         /L_chin1 - L_gxy * L_gxxx_p ) * HALF
        L_Gamyxy_p = Gamyxy(i,j,k) - (  L_chix         /L_chin1 - L_gxy * L_gxxy_p ) * HALF
        L_Gamzxy_p = Gamzxy(i,j,k) - (                         - L_gxy * L_gxxz_p ) * HALF
        L_Gamxxz_p = Gamxxz(i,j,k) - (  L_chiz         /L_chin1 - L_gxz * L_gxxx_p ) * HALF
        L_Gamyxz_p = Gamyxz(i,j,k) - (                         - L_gxz * L_gxxy_p ) * HALF
        L_Gamzxz_p = Gamzxz(i,j,k) - (  L_chix         /L_chin1 - L_gxz * L_gxxz_p ) * HALF
        L_Gamxyz_p = Gamxyz(i,j,k) - (                         - L_gyz * L_gxxx_p ) * HALF
        L_Gamyyz_p = Gamyyz(i,j,k) - (  L_chiz         /L_chin1 - L_gyz * L_gxxy_p ) * HALF
        L_Gamzyz_p = Gamzyz(i,j,k) - (  L_chiy         /L_chin1 - L_gyz * L_gxxz_p ) * HALF

        Gamxxx(i,j,k) = L_Gamxxx_p
        Gamyxx(i,j,k) = L_Gamyxx_p
        Gamzxx(i,j,k) = L_Gamzxx_p
        Gamxyy(i,j,k) = L_Gamxyy_p
        Gamyyy(i,j,k) = L_Gamyyy_p
        Gamzyy(i,j,k) = L_Gamzyy_p
        Gamxzz(i,j,k) = L_Gamxzz_p
        Gamyzz(i,j,k) = L_Gamyzz_p
        Gamzzz(i,j,k) = L_Gamzzz_p
        Gamxxy(i,j,k) = L_Gamxxy_p
        Gamyxy(i,j,k) = L_Gamyxy_p
        Gamzxy(i,j,k) = L_Gamzxy_p
        Gamxxz(i,j,k) = L_Gamxxz_p
        Gamyxz(i,j,k) = L_Gamyxz_p
        Gamzxz(i,j,k) = L_Gamzxz_p
        Gamxyz(i,j,k) = L_Gamxyz_p
        Gamyyz(i,j,k) = L_Gamyyz_p
        Gamzyz(i,j,k) = L_Gamzyz_p

        L_lapx = Lapx(i,j,k); L_lapy = Lapy(i,j,k); L_lapz = Lapz(i,j,k)

        L_lapxx = lapxx(i,j,k) - L_Gamxxx_p * L_lapx - L_Gamyxx_p * L_lapy - L_Gamzxx_p * L_lapz
        L_lapyy = lapyy(i,j,k) - L_Gamxyy_p * L_lapx - L_Gamyyy_p * L_lapy - L_Gamzyy_p * L_lapz
        L_lapzz = lapzz(i,j,k) - L_Gamxzz_p * L_lapx - L_Gamyzz_p * L_lapy - L_Gamzzz_p * L_lapz
        L_lapxy = lapxy(i,j,k) - L_Gamxxy_p * L_lapx - L_Gamyxy_p * L_lapy - L_Gamzxy_p * L_lapz
        L_lapxz = lapxz(i,j,k) - L_Gamxxz_p * L_lapx - L_Gamyxz_p * L_lapy - L_Gamzxz_p * L_lapz
        L_lapyz = lapyz(i,j,k) - L_Gamxyz_p * L_lapx - L_Gamyyz_p * L_lapy - L_Gamzyz_p * L_lapz

        L_lap_lap = L_gupxx * L_lapxx + L_gupyy * L_lapyy + L_gupzz * L_lapzz + &
                      TWO* ( L_gupxy * L_lapxy + L_gupxz * L_lapxz + L_gupyz * L_lapyz )
        
        L_S = L_chin1 * ( L_gupxx * Sxx(i,j,k) + L_gupyy * Syy(i,j,k) + L_gupzz * Szz(i,j,k) + &
              TWO * ( L_gupxy * Sxy(i,j,k) + L_gupxz * Sxz(i,j,k) + L_gupyz * Syz(i,j,k) ) )

        L_f_trace = F2o3 * trK(i,j,k) * trK(i,j,k) - ( &
                  L_gupxx * ( &
                  L_gupxx * L_Axx * L_Axx + L_gupyy * L_Axy * L_Axy + L_gupzz * L_Axz * L_Axz + &
                  TWO * (L_gupxy * L_Axx * L_Axy + L_gupxz * L_Axx * L_Axz + L_gupyz * L_Axy * L_Axz) ) + &
                  L_gupyy * ( &
                  L_gupxx * L_Axy * L_Axy + L_gupyy * L_Ayy * L_Ayy + L_gupzz * L_Ayz * L_Ayz + &
                  TWO * (L_gupxy * L_Axy * L_Ayy + L_gupxz * L_Axy * L_Ayz + L_gupyz * L_Ayy * L_Ayz) ) + &
                  L_gupzz * ( &
                  L_gupxx * L_Axz * L_Axz + L_gupyy * L_Ayz * L_Ayz + L_gupzz * L_Azz * L_Azz + &
                  TWO * (L_gupxy * L_Axz * L_Ayz + L_gupxz * L_Axz * L_Azz + L_gupyz * L_Ayz * L_Azz) ) + &
                  TWO * ( &
                  L_gupxy * ( &
                  L_gupxx * L_Axx * L_Axy + L_gupyy * L_Axy * L_Ayy + L_gupzz * L_Axz * L_Ayz + &
                  L_gupxy * (L_Axx * L_Ayy + L_Axy * L_Axy) + &
                  L_gupxz * (L_Axx * L_Ayz + L_Axz * L_Axy) + &
                  L_gupyz * (L_Axy * L_Ayz + L_Axz * L_Ayy) ) + &
                  L_gupxz * ( &
                  L_gupxx * L_Axx * L_Axz + L_gupyy * L_Axy * L_Ayz + L_gupzz * L_Axz * L_Azz + &
                  L_gupxy * (L_Axx * L_Ayz + L_Axy * L_Axz) + &
                  L_gupxz * (L_Axx * L_Azz + L_Axz * L_Axz) + &
                  L_gupyz * (L_Axy * L_Azz + L_Axz * L_Ayz) ) + &
                  L_gupyz * ( &
                  L_gupxx * L_Axy * L_Axz + L_gupyy * L_Ayy * L_Ayz + L_gupzz * L_Ayz * L_Azz + &
                  L_gupxy * (L_Axy * L_Ayz + L_Ayy * L_Axz) + &
                  L_gupxz * (L_Axy * L_Azz + L_Ayz * L_Axz) + &
                  L_gupyz * (L_Ayy * L_Azz + L_Ayz * L_Ayz) ) )) - 1.6d1 * PI * rho(i,j,k) + EIGHT * PI * L_S
        
        L_tf = -F1o3 * ( L_gupxx * L_lapxx + L_gupyy * L_lapyy + L_gupzz * L_lapzz + &
             TWO*(L_gupxy * L_lapxy + L_gupxz * L_lapxz + L_gupyz * L_lapyz) + L_alpn1 / L_chin1 * L_f_trace )

        L_lapxx = L_alpn1 * (Rxx(i,j,k) - EIGHT * PI * Sxx(i,j,k)) - L_lapxx
        L_lapxy = L_alpn1 * (Rxy(i,j,k) - EIGHT * PI * Sxy(i,j,k)) - L_lapxy
        L_lapxz = L_alpn1 * (Rxz(i,j,k) - EIGHT * PI * Sxz(i,j,k)) - L_lapxz
        L_lapyy = L_alpn1 * (Ryy(i,j,k) - EIGHT * PI * Syy(i,j,k)) - L_lapyy
        L_lapyz = L_alpn1 * (Ryz(i,j,k) - EIGHT * PI * Syz(i,j,k)) - L_lapyz
        L_lapzz = L_alpn1 * (Rzz(i,j,k) - EIGHT * PI * Szz(i,j,k)) - L_lapzz

        L_tmp_xx = L_lapxx - L_gxx * L_tf
        L_tmp_yy = L_lapyy - L_gyy * L_tf
        L_tmp_zz = L_lapzz - L_gzz * L_tf
        L_tmp_xy = L_lapxy - L_gxy * L_tf
        L_tmp_xz = L_lapxz - L_gxz * L_tf
        L_tmp_yz = L_lapyz - L_gyz * L_tf

        L_lapxx = L_gupxx * L_Axx * L_Axx + L_gupyy * L_Axy * L_Axy + L_gupzz * L_Axz * L_Axz + &
                TWO * (L_gupxy * L_Axx * L_Axy + L_gupxz * L_Axx * L_Axz + L_gupyz * L_Axy * L_Axz)
        L_lapyy = L_gupxx * L_Axy * L_Axy + L_gupyy * L_Ayy * L_Ayy + L_gupzz * L_Ayz * L_Ayz + &
                TWO * (L_gupxy * L_Axy * L_Ayy + L_gupxz * L_Axy * L_Ayz + L_gupyz * L_Ayy * L_Ayz)
        L_lapzz = L_gupxx * L_Axz * L_Axz + L_gupyy * L_Ayz * L_Ayz + L_gupzz * L_Azz * L_Azz + &
                TWO * (L_gupxy * L_Axz * L_Ayz + L_gupxz * L_Axz * L_Azz + L_gupyz * L_Ayz * L_Azz)
        L_lapxy = L_gupxx * L_Axx * L_Axy + L_gupyy * L_Axy * L_Ayy + L_gupzz * L_Axz * L_Ayz + &
                  L_gupxy *(L_Axx * L_Ayy + L_Axy * L_Axy) + &
                  L_gupxz *(L_Axx * L_Ayz + L_Axz * L_Axy) + &
                  L_gupyz *(L_Axy * L_Ayz + L_Axz * L_Ayy)
        L_lapxz = L_gupxx * L_Axx * L_Axz + L_gupyy * L_Axy * L_Ayz + L_gupzz * L_Axz * L_Azz + &
                  L_gupxy *(L_Axx * L_Ayz + L_Axy * L_Axz) + &
                  L_gupxz *(L_Axx * L_Azz + L_Axz * L_Axz) + &
                  L_gupyz *(L_Axy * L_Azz + L_Axz * L_Ayz)
        L_lapyz = L_gupxx * L_Axy * L_Axz + L_gupyy * L_Ayy * L_Ayz + L_gupzz * L_Ayz * L_Azz + &
                  L_gupxy *(L_Axy * L_Ayz + L_Ayy * L_Axz) + &
                  L_gupxz *(L_Axy * L_Azz + L_Ayz * L_Axz) + &
                  L_gupyz *(L_Ayy * L_Azz + L_Ayz * L_Ayz)

        Axx_rhs(i,j,k) = L_chin1 * L_tmp_xx+ L_alpn1 * (trK(i,j,k) * L_Axx - TWO * L_lapxx)  + &
                         TWO * (L_Axx * betaxx(i,j,k) + L_Axy * betayx(i,j,k) + L_Axz * betazx(i,j,k)) - &
                         F2o3 * L_Axx * L_div_beta
        Ayy_rhs(i,j,k) = L_chin1 * L_tmp_yy + L_alpn1 * (trK(i,j,k) * L_Ayy - TWO * L_lapyy) + &
                         TWO * (L_Axy * betaxy(i,j,k) + L_Ayy * betayy(i,j,k) + L_Ayz * betazy(i,j,k)) - &
                         F2o3 * L_Ayy * L_div_beta
        Azz_rhs(i,j,k) = L_chin1 * L_tmp_zz + L_alpn1 * (trK(i,j,k) * L_Azz - TWO * L_lapzz) + &
                         TWO * (L_Axz * betaxz(i,j,k) + L_Ayz * betayz(i,j,k) + L_Azz * betazz(i,j,k)) - &
                         F2o3 * L_Azz * L_div_beta
        Axy_rhs(i,j,k) = L_chin1 * L_tmp_xy + L_alpn1 * (trK(i,j,k) * L_Axy - TWO * L_lapxy) + &
                         L_Axx * betaxy(i,j,k) + L_Axz * betazy(i,j,k) + &
                         L_Ayy * betayx(i,j,k) + L_Ayz * betazx(i,j,k) + &
                         F1o3 * L_Axy * L_div_beta - L_Axy * betazz(i,j,k)

        Ayz_rhs(i,j,k) = L_chin1 * L_tmp_yz + L_alpn1 * (trK(i,j,k) * L_Ayz - TWO * L_lapyz) + &
                         L_Axy * betaxz(i,j,k) + L_Ayy * betayz(i,j,k) + &
                         L_Axz * betaxy(i,j,k) + L_Azz * betazy(i,j,k) + &
                         F1o3 * L_Ayz * L_div_beta - L_Ayz * betaxx(i,j,k)

        Axz_rhs(i,j,k) = L_chin1 * L_tmp_xz + L_alpn1 * (trK(i,j,k) * L_Axz - TWO * L_lapxz) + &
                         L_Axx * betaxz(i,j,k) + L_Axy * betayz(i,j,k) + &
                         L_Ayz * betayx(i,j,k) + L_Azz * betazx(i,j,k) + &
                         F1o3 * L_Axz * L_div_beta - L_Axz * betayy(i,j,k)

        L_S = L_chin1 * (L_gupxx * Sxx(i,j,k) + L_gupyy * Syy(i,j,k)+ L_gupzz * Szz(i,j,k) + &
              TWO * (L_gupxy * Sxy(i,j,k) + L_gupxz * Sxz(i,j,k) + L_gupyz * Syz(i,j,k)))

        trK_rhs(i,j,k) = - L_chin1 * L_lap_lap + L_alpn1 *(F1o3 * trK(i,j,k) * trK(i,j,k) + &
                        L_gupxx * L_lapxx + L_gupyy * L_lapyy + L_gupzz * L_lapzz + &
                        TWO * (L_gupxy * L_lapxy + L_gupxz * L_lapxz + L_gupyz * L_lapyz) + &
                        FOUR * PI * (rho(i,j,k) + L_S))
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  
!!!! gauge variable part

  Lap_rhs = -TWO*alpn1*trK
#if (GAUGE == 0)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz
#elif (GAUGE == 1)
  betax_rhs = Gamx - eta*betax
  betay_rhs = Gamy - eta*betay
  betaz_rhs = Gamz - eta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 2)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 3)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 4)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 5)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 6)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1/(ONE+w1*r1) + C2/(ONE+w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 1",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 7)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1*dexp(-w1*r1) + C2*dexp(-w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 2",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#endif  

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

!!!!!!!!!advection term part

  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "before",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "after",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1 
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
#endif

  endif

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

movx_Res = movx_Res - F2o3*Kx - F8*PI*sx
movy_Res = movy_Res - F2o3*Ky - F8*PI*sy
movz_Res = movz_Res - F2o3*Kz - F8*PI*sz
  endif

#if (ABV == 1)
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res,movx_Res,movy_Res,movz_Res,Gmx_Res,Gmy_Res,Gmz_Res, &
               Symmetry)
#endif 
#if 0
#define i 2
if(Lev == 1)then
write(*,*) X(i),Y(i),Z(i)
write(*,*) Axx(i,i,i),Axy(i,i,i),Axz(i,i,i),Ayy(i,i,i),Ayz(i,i,i),Azz(i,i,i)
write(*,*) 1+Lap(i,i,i),dtSfx(i,i,i),dtSfy(i,i,i),dtSfz(i,i,i)
write(*,*) betax(i,i,i),betay(i,i,i),betaz(i,i,i)
write(*,*) 1+chi(i,i,i),Gamx(i,i,i),Gamy(i,i,i),Gamz(i,i,i)
write(*,*) gxx(i,i,i),gxy(i,i,i),gxz(i,i,i),gyy(i,i,i),gyz(i,i,i),gzz(i,i,i)
write(*,*) trK(i,i,i)
write(*,*) "====="
write(*,*) Axx_rhs(i,i,i),Axy_rhs(i,i,i),Axz_rhs(i,i,i),Ayy_rhs(i,i,i),Ayz_rhs(i,i,i),Azz_rhs(i,i,i)
write(*,*) Lap_rhs(i,i,i),dtSfx_rhs(i,i,i),dtSfy_rhs(i,i,i),dtSfz_rhs(i,i,i)
write(*,*) betax_rhs(i,i,i),betay_rhs(i,i,i),betaz_rhs(i,i,i)
write(*,*) chi_rhs(i,i,i),Gamx_rhs(i,i,i),Gamy_rhs(i,i,i),Gamz_rhs(i,i,i)
write(*,*) gxx_rhs(i,i,i),gxy_rhs(i,i,i),gxz_rhs(i,i,i),gyy_rhs(i,i,i),gyz_rhs(i,i,i),gzz_rhs(i,i,i)
write(*,*) trK_rhs(i,i,i)
endif
#undef i
!!stop
#endif

  gont = 0

  return

  end function compute_rhs_bssn
