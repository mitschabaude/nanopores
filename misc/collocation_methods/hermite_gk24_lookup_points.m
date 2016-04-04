function x = hermite_gk24_lookup_points ( n )

%*****************************************************************************80
%
%% HERMITE_GK24_LOOKUP_POINTS: abscissas of a Genz-Keister 24 Hermite rule.
%
%  Discussion:
%
%    The integral:
%
%      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
%
%    The quadrature rule:
%
%      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
%
%    A nested family of rules for the Hermite integration problem
%    was produced by Genz and Keister.  The structure of the nested
%    family was denoted by 1+2+6+10+24, that is, it comprised rules
%    of successive orders O = 1, 3, 9, 19, and 43.
%
%    The precisions of these rules are P = 1, 5, 15, 29, and 67.
%
%    Some of the data in this function was kindly supplied directly by
%    Alan Genz on 24 April 2011.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 April 2011
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Alan Genz, Bradley Keister,
%    Fully symmetric interpolatory rules for multiple integrals
%    over infinite regions with Gaussian weight,
%    Journal of Computational and Applied Mathematics,
%    Volume 71, 1996, pages 299-309
%
%    Thomas Patterson,
%    The Optimal Addition of Points to Quadrature Formulae,
%    Mathematics of Computation,
%    Volume 22, Number 104, October 1968, pages 847-856.
%
%  Parameters:
%
%    Input, integer N, the order.
%    Legal values are 1, 3, 9, 19, and 43.
%
%    Output, real X(N,1), the abscissas.
%
  x = zeros ( n, 1 );

  if ( n == 1 )

    x( 1) =   0.0000000000000000;

  elseif ( n == 3 )

    x( 1) =  -1.2247448713915889;
    x( 2) =   0.0000000000000000;
    x( 3) =   1.2247448713915889;

  elseif ( n == 9 )

    x( 1) =  -2.9592107790638380;
    x( 2) =  -2.0232301911005157;
    x( 3) =  -1.2247448713915889;
    x( 4) =  -0.52403354748695763;
    x( 5) =   0.0000000000000000;
    x( 6) =   0.52403354748695763;
    x( 7) =   1.2247448713915889;
    x( 8) =   2.0232301911005157;
    x( 9) =   2.9592107790638380;

  elseif ( n == 19 )

    x( 1) =  -4.4995993983103881;
    x( 2) =  -3.6677742159463378;
    x( 3) =  -2.9592107790638380;
    x( 4) =  -2.2665132620567876;
    x( 5) =  -2.0232301911005157;
    x( 6) =  -1.8357079751751868;
    x( 7) =  -1.2247448713915889;
    x( 8) =  -0.87004089535290285;
    x( 9) =  -0.52403354748695763;
    x(10) =   0.0000000000000000;
    x(11) =   0.52403354748695763;
    x(12) =   0.87004089535290285;
    x(13) =   1.2247448713915889;
    x(14) =   1.8357079751751868;
    x(15) =   2.0232301911005157;
    x(16) =   2.2665132620567876;
    x(17) =   2.9592107790638380;
    x(18) =   3.6677742159463378;
    x(19) =   4.4995993983103881;

  elseif ( n == 43 )

    x( 1) = -10.167574994881873;
    x( 2) =  -7.231746029072501;
    x( 3) =  -6.535398426382995;
    x( 4) =  -5.954781975039809;
    x( 5) =  -5.434053000365068;
    x( 6) =  -4.952329763008589;
    x( 7) =  -4.4995993983103881;
    x( 8) =  -4.071335874253583;
    x( 9) =  -3.6677742159463378;
    x(10) =  -3.295265921534226;
    x(11) =  -2.9592107790638380;
    x(12) =  -2.633356763661946;
    x(13) =  -2.2665132620567876;
    x(14) =  -2.089340389294661;
    x(15) =  -2.0232301911005157;
    x(16) =  -1.8357079751751868;
    x(17) =  -1.583643465293944;
    x(18) =  -1.2247448713915889;
    x(19) =  -0.87004089535290285;
    x(20) =  -0.52403354748695763;
    x(21) =  -0.196029453662011;
    x(22) =   0.0000000000000000;
    x(23) =   0.196029453662011;
    x(24) =   0.52403354748695763;
    x(25) =   0.87004089535290285;
    x(26) =   1.2247448713915889;
    x(27) =   1.583643465293944;
    x(28) =   1.8357079751751868;
    x(29) =   2.0232301911005157;
    x(30) =   2.089340389294661;
    x(31) =   2.2665132620567876;
    x(32) =   2.633356763661946;
    x(33) =   2.9592107790638380;
    x(34) =   3.295265921534226;
    x(35) =   3.6677742159463378;
    x(36) =   4.071335874253583;
    x(37) =   4.4995993983103881;
    x(38) =   4.952329763008589;
    x(39) =   5.434053000365068;
    x(40) =   5.954781975039809;
    x(41) =   6.535398426382995;
    x(42) =   7.231746029072501;
    x(43) =  10.167574994881873;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'HERMITE_GK24_LOOKUP_POINTS - Fatal error!\n' );
    fprintf ( 1, '  Illegal input value of N.\n' );
    fprintf ( 1, '  N must be 1, 3, 9, 19, or 43.\n' );
    error ( 'HERMITE_GK24_LOOKUP_POINTS - Fatal error!' );

  end

  return
end
