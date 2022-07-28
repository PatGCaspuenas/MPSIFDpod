function sfun = CubicSpline(xi,yi)
%%                           CubicSpline.m
%--------------------------------------------------------------------------
%
% Creates a third-order polynomial interpolant with 1st-order continuity at
% the nodes
%
% INPUTS
%
%   xi   : nodes
%   yi   : points at which the function will be passing through ( f(xi) )
%
% OUTPUT
%
%   sfun : 3rd-order polynomial as a function of x
%
%--------------------------------------------------------------------------
%

    x1 = xi(1); x2 = xi(2); x3 = xi(3);
    B1 = [yi(1);yi(2);yi(2);yi(3);0;0;0;0];
    %
    %    a1    b1       c1       d1       a2     b2     c2       d2
    % -------------------------------------------------------------------
    A1 = [1     0        0        0        0     0       0        0;
          1  (x2-x1)  (x2-x1)^2 (x2-x1)^3  0     0       0        0;
          0     0        0        0        1     0       0        0;
          0     0        0        0        1  (x2-x1) (x2-x1)^2 (x2-x1)^3;
          0     1        0        0        0     0       0        0;
          0     0        0        0        0     1   2*(x3-x2) 3*(x3-x2)^2;
          0     1   2*(x2-x1) 3*(x2-x1)^2  0    -1       0        0;
          0     0        2    6*(x2-x1)    0    0       -2        0];

%

    coeff1 = A1\B1;

%
% Coefficients
%

    a1 = coeff1(1); b1 = coeff1(2); c1 = coeff1(3); d1 = coeff1(4);
    a2 = coeff1(5); b2 = coeff1(6); c2 = coeff1(7); d2 = coeff1(8); 

%
% Piece-wise interpolators
%

    sfun = @(x) (a1 + b1*(x-x1) + c1*(x-x1).^2 + d1*(x-x1).^3).*((x >= x1) & (x < x2)) + (a2 + b2*(x-x2) + c2*(x-x2).^2 + d2*(x-x2).^3).*((x >= x2) & (x <= x3));

%
end