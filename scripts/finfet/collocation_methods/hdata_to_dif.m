function [ xd, yd ] = hdata_to_dif ( n, x, y, yp )

%*****************************************************************************80
%
%% HDATA_TO_DIF sets up a divided difference table from Hermite data.
%
%  Discussion:
%
%    The polynomial represented by the divided difference table can be
%    evaluated by calling DIF_VALS.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 May 2011
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Carl deBoor,
%    A Practical Guide to Splines,
%    Springer, 2001,
%    ISBN: 0387953663,
%    LC: QA1.A647.v27.
%
%  Parameters:
%
%    Input, integer N, the number of items of data
%    ( X(I), Y(I), YP(I) ).
%
%    Input, real X(N), the abscissas.
%    These values must be distinct.
%
%    Input, real Y(N), YP(N), the function and
%    derivative values.
%
%    Output, real XD(2*N), YD(2*N), the divided
%    difference coefficients.
%
  x = x ( : );
  y = y ( : );
  yp = yp ( : );
%
%  Copy the data:
%
  nd = 2 * n;
  xd = zeros ( nd, 1 );
  xd(1:2:nd-1) = x(1:n);
  xd(2:2:nd  ) = x(1:n);
%
%  Carry out the first step of differencing.
%
  yd = zeros ( nd, 1 );
  yd(1) = y(1);
  yd(3:2:nd-1) = ( y(2:n) - y(1:n-1) ) ./ ( x(2:n) - x(1:n-1) );
  yd(2:2:nd  ) = yp(1:n);
%
%  Carry out the remaining steps in the usual way.
%
  for i = 3 : nd
    for j = nd : -1 : i
      yd(j) = ( yd(j) - yd(j-1) ) / ( xd(j) - xd(j+1-i) );
    end
  end

  return
end
