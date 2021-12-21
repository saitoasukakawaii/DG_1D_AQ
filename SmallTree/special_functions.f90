subroutine cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1 )

!*****************************************************************************80
!
!! CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, 
!    CDY1, the values of J0(z), J0'(z), J1(z), J1'(z), Y0(z), Y0'(z), 
!    Y1(z), Y1'(z).
!
  implicit none

  real ( kind = 8 ), save, dimension ( 12 ) :: a = (/ &
    -0.703125D-01,0.112152099609375D+00, &
    -0.5725014209747314D+00,0.6074042001273483D+01, &
    -0.1100171402692467D+03,0.3038090510922384D+04, &
    -0.1188384262567832D+06,0.6252951493434797D+07, &
    -0.4259392165047669D+09,0.3646840080706556D+11, &
    -0.3833534661393944D+13,0.4854014686852901D+15 /)
  real ( kind = 8 ) a0
  real ( kind = 8 ), save, dimension ( 12 ) :: a1 = (/ &
    0.1171875D+00,-0.144195556640625D+00, &
    0.6765925884246826D+00,-0.6883914268109947D+01, &
    0.1215978918765359D+03,-0.3302272294480852D+04, &
    0.1276412726461746D+06,-0.6656367718817688D+07, &
    0.4502786003050393D+09,-0.3833857520742790D+11, &
    0.4011838599133198D+13,-0.5060568503314727D+15 /)
  real ( kind = 8 ), save, dimension ( 12 ) :: b = (/ &
    0.732421875D-01,-0.2271080017089844D+00, &
    0.1727727502584457D+01,-0.2438052969955606D+02, &
    0.5513358961220206D+03,-0.1825775547429318D+05, &
    0.8328593040162893D+06,-0.5006958953198893D+08, &
    0.3836255180230433D+10,-0.3649010818849833D+12, &
    0.4218971570284096D+14,-0.5827244631566907D+16 /)
  real ( kind = 8 ), save, dimension ( 12 ) :: b1 = (/ &
    -0.1025390625D+00,0.2775764465332031D+00, &
    -0.1993531733751297D+01,0.2724882731126854D+02, &
    -0.6038440767050702D+03,0.1971837591223663D+05, &
    -0.8902978767070678D+06,0.5310411010968522D+08, &
    -0.4043620325107754D+10,0.3827011346598605D+12, &
    -0.4406481417852278D+14,0.6065091351222699D+16 /)
  complex ( kind = 8 ) cbj0
  complex ( kind = 8 ) cbj1
  complex ( kind = 8 ) cby0
  complex ( kind = 8 ) cby1
  complex ( kind = 8 ) cdj0
  complex ( kind = 8 ) cdj1
  complex ( kind = 8 ) cdy0
  complex ( kind = 8 ) cdy1
  complex ( kind = 8 ) ci
  complex ( kind = 8 ) cp
  complex ( kind = 8 ) cp0
  complex ( kind = 8 ) cp1
  complex ( kind = 8 ) cq0
  complex ( kind = 8 ) cq1
  complex ( kind = 8 ) cr
  complex ( kind = 8 ) cs
  complex ( kind = 8 ) ct1
  complex ( kind = 8 ) ct2
  complex ( kind = 8 ) cu
  real ( kind = 8 ) el
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  real ( kind = 8 ) pi
  real ( kind = 8 ) rp2
  real ( kind = 8 ) w0
  real ( kind = 8 ) w1
  complex ( kind = 8 ) z
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  rp2 = 2.0D+00 / pi
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0D+00 ) then
    cbj0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cbj1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdj0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdj1 = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
    cby0 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cby1 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdy0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdy1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    return
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    z1 = -z
  end if

  if ( a0 <= 12.0D+00 ) then

    cbj0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 40
      cr = -0.25D+00 * cr * z2 / ( k * k )
      cbj0 = cbj0 + cr
      if ( abs ( cr ) < abs ( cbj0 ) * 1.0D-15 ) then
        exit
      end if
    end do

    cbj1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 40
      cr = -0.25D+00 * cr * z2 / ( k * ( k + 1.0D+00 ) )
      cbj1 = cbj1 + cr
      if ( abs ( cr ) < abs ( cbj1 ) * 1.0D-15 ) then
        exit
      end if
    end do

    cbj1 = 0.5D+00 * z1 * cbj1
    w0 = 0.0D+00
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 40
      w0 = w0 + 1.0D+00 / k
      cr = -0.25D+00 * cr / ( k * k ) * z2
      cp = cr * w0
      cs = cs + cp
      if ( abs ( cp ) < abs ( cs ) * 1.0D-15 ) then
        exit
      end if
    end do

    cby0 = rp2 * ( log ( z1 / 2.0D+00 ) + el ) * cbj0 - rp2 * cs
    w1 = 0.0D+00
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 40
      w1 = w1 + 1.0D+00 / k
      cr = -0.25D+00 * cr / ( k * ( k + 1 ) ) * z2
      cp = cr * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
      cs = cs + cp
      if ( abs ( cp ) < abs ( cs ) * 1.0D-15 ) then
        exit
      end if
    end do

    cby1 = rp2 * ( ( log ( z1 / 2.0D+00 ) + el ) * cbj1 &
      - 1.0D+00 / z1 - 0.25D+00 * z1 * cs )

  else

    if ( a0 < 35.0D+00 ) then
      k0 = 12
    else if ( a0 < 50.0D+00 ) then
      k0 = 10
    else
      k0 = 8
    end if

    ct1 = z1 - 0.25D+00 * pi

    cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, k0
      cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
    end do

    cq0 = -0.125D+00 / z1
    do k = 1, k0
      cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
    end do

    cu = sqrt ( rp2 / z1 )
    cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
    cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
    ct2 = z1 - 0.75D+00 * pi

    cp1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, k0
      cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
    end do

    cq1 = 0.375D+00 / z1
    do k = 1, k0
      cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
    end do

    cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
    cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    if ( imag ( z ) < 0.0D+00 ) then
      cby0 =     cby0 - 2.0D+00 * ci * cbj0
      cby1 = - ( cby1 - 2.0D+00 * ci * cbj1 )
    else
      cby0 =     cby0 + 2.0D+00 * ci * cbj0
      cby1 = - ( cby1 + 2.0D+00 * ci * cbj1 )
    end if
    cbj1 = -cbj1
  end if

  cdj0 = -cbj1
  cdj1 = cbj0 - 1.0D+00 / z * cbj1
  cdy0 = -cby1
  cdy1 = cby0 - 1.0D+00 / z * cby1

  return
end
