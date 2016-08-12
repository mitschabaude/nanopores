function w = hc_compute_weights_from_points ( m, x )

%*****************************************************************************80
%
%% HC_COMPUTE_WEIGHTS_FROM_POINTS: Hermite-Cubic weights, user-supplied points.
%
%  Discussion:
%
%    An interval [A,B] has been divided by NHALF points X; at each
%    point both function and derivative information is available.
%
%    The piecewise cubic Hermite interpolant is constructed for this data.
%
%    A quadrature rule is determined for the interpolant.
%
%    There will be N=2*M weights.  If the quadrature rule is to be written 
%    out, one would normally list each point twice, so that the number of points
%    and weights are equal.  The listing of the same point value twice is an
%    implicit indication that both function and derivative values should be
%    used.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    31 March 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the number of points, not counting 
%    repetitions.
%
%    Input, real X(M), the points, without repetition.
%
%    Output, real W(2,M), the weights.  The first two weights are 
%    associated with the first point, and so on.
%
  w = zeros ( 2, m );

  w(1,1) = 0.5 * ( x(2) - x(1) );
  w(2,1) = ( x(2) - x(1) )^2 / 12.0;

  for j = 2 : m - 1
    w(1,j) = 0.5 * ( x(j+1) - x(j-1) );
    w(2,j) = ( x(j+1) - x(j-1) ) * ( x(j+1) - 2.0 * x(j) + x(j-1) ) / 12.0;
  end

  w(1,m) = 0.5 * ( x(m) - x(m-1) );
  w(2,m) = - ( x(m-1) - x(m) )^2 / 12.0;

  return
end
