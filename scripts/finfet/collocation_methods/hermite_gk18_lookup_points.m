function x = hermite_gk18_lookup_points ( n )

%*****************************************************************************80
%
%% HERMITE_GK18_LOOKUP_POINTS: abscissas of a Genz-Keister 18 Hermite rule.
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
%    family was denoted by 1+2+6+10+18, that is, it comprised rules 
%    of successive orders O = 1, 3, 9, 19, and 37.
%
%    The precisions of these rules are P = 1, 5, 15, 29, and 55.
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
%    30 April 2011
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
%    Florian Heiss, Viktor Winschel,
%    Likelihood approximation by numerical integration on sparse grids,
%    Journal of Econometrics,
%    Volume 144, 2008, pages 62-80.
%
%    Thomas Patterson,
%    The Optimal Addition of Points to Quadrature Formulae,
%    Mathematics of Computation,
%    Volume 22, Number 104, October 1968, pages 847-856.
%
%  Parameters:
%
%    Input, integer N, the order.
%    N must be 1, 3, 9, 19, or 37.
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

  elseif ( n == 35 )

    x( 1) =  -6.853200069757519;
    x( 2) =  -6.124527854622158;
    x( 3) =  -5.521865209868350;
    x( 4) =  -4.986551454150765;
    x( 5) =  -4.499599398310388;
    x( 6) =  -4.057956316089741;
    x( 7) =  -3.667774215946338;
    x( 8) =  -3.315584617593290;
    x( 9) =  -2.959210779063838;
    x(10) =  -2.597288631188366;
    x(11) =  -2.266513262056788;
    x(12) =  -2.023230191100516;
    x(13) =  -1.835707975175187;
    x(14) =  -1.561553427651873;
    x(15) =  -1.224744871391589;
    x(16) =  -0.870040895352903;
    x(17) =  -0.524033547486958;
    x(18) =  -0.214618180588171;
    x(19) =   0.000000000000000;
    x(20) =   0.214618180588171;
    x(21) =   0.524033547486958;
    x(22) =   0.870040895352903;
    x(23) =   1.224744871391589;
    x(24) =   1.561553427651873;
    x(25) =   1.835707975175187;
    x(26) =   2.023230191100516;
    x(27) =   2.266513262056788;
    x(28) =   2.597288631188366;
    x(29) =   2.959210779063838;
    x(30) =   3.315584617593290;
    x(31) =   3.667774215946338;
    x(32) =   4.057956316089741;
    x(33) =   4.499599398310388;
    x(34) =   4.986551454150765;
    x(35) =   5.521865209868350;
    x(36) =   6.124527854622158;
    x(37) =   6.853200069757519;

  else

    fprintf ( stderr, '\n' );
    fprintf ( stderr, 'HERMITE_GK18_LOOKUP_POINTS - Fatal error!\n' );
    fprintf ( stderr, '  Illegal input value of N.\n' );
    fprintf ( stderr, '  N must be 1, 3, 9, 19, or 37.\n' );
    error ( 'HERMITE_GK18_LOOKUP_POINTS - Fatal error!' );

  end

  return
end
