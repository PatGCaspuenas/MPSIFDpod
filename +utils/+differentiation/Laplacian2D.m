function [D,fx,fy] = Laplacian2D(f, dx, dy)
% 2D Laplacian operator on f
% calculates the values on the edges by linearly extrapolating the second
% differences from the interior
%
% Woii junwei.chen@alumnos.uc3m.es 20210922 v1.0
% Experimental Aerodynamics and Propulsion Lab
% of Univeridad Carlos III de Madrid
%
fx = 0.*f;
fx(:,2:end-1) = (f(:,1:end-2)+f(:,3:end)-2*f(:,2:end-1))./dx^2;
fx(:,1) = (2*f(:,1) - 5*f(:,2) + 4*f(:,3) - f(:,4))./dx^2;
fx(:,end) = (2*f(:,end) - 5*f(:,end-1) + 4*f(:,end-2) - f(:,end-3))./dx^2;
%
% Another option of doing 2nd order derivatives at extreme points
%
% fx(:,1) = 2*fx(:,2) - fx(:,3);
% fx(:,end) = 2*fx(:,end-1) - fx(:,end-2);
%
fy = 0.*f;
fy(2:end-1,:) = (f(1:end-2,:)+f(3:end,:)-2*f(2:end-1,:))./dy^2;
%
fy(1,:) = (2*f(1,:) - 5*f(2,:) + 4*f(3,:) - f(4,:))./dy^2;
fy(end,:) = (2*f(end,:) - 5*f(end-1,:) + 4*f(end-2,:) - f(end-3,:))./dy^2;
%
% Another option of doing 2nd order derivatives at extreme points
%
% fy(1,:) = 2*fy(2,:) - fy(3,:);
% fy(end,:) = 2*fy(end-1,:) - fy(end-2,:);
%
D = fx + fy;
end