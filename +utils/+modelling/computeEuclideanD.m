function d = computeEuclideanD(m,n,x0,y0)
%%                       computeEuclideanD.m
%--------------------------------------------------------------------------
%
% Calculates euclidean distance from point (x0,y0) to straight line defined
% by slope m and y-interceipt n. We will define our straight line as
%
%                            y = m*x + n
%         ax + by + c = 0        ->         by = -ax -c
%
%
% INPUTS
%
%   m       : line slope
%   n       : y-interceipt
%   (x0,y0) : external point to line
%
% OUTPUT
%
%   d       : euclidean distance
%
%--------------------------------------------------------------------------
%
% Curve coefficients a,b,c
%

    a = -m; 
    c = -n; 
    b = 1;

%
% Projected point (x0,y0) on line
%

    x1 = (b*(b*x0-a*y0)-a*c)/(a^2+b^2);
    y1 = (a*(-b*x0+a*y0)-b*c)/(a^2+b^2);

%
% Distance
%

    d = abs(a*x0+b*y0+c)/sqrt(a^2+b^2);
    d = d - 2*d.*(x1<x0);

%