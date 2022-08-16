function dY = odeSINDyGlobal(t,Y,S)
%%                         odeSINDyGlobal.m
%--------------------------------------------------------------------------
%
% Function handle for integration purposes (odefun of ode45, e.g.)
%
% INPUTS
%
%   t        : time instant of integration
%   Y        : states of the system at time instant t. Dimensions of 1 x Nr
%   S        : structure containing matrix of coefficients Chi (Nf x Nr)
%              and polynomial order (vector from 0 up to chosen order)
%
% OUTPUT
%
%   dY       : derivative of states of the system, of dimension Nr x 1
%
% UTILS
%
%   utils.SINDy.poolpolyData
%
%-------------------------------------------------------------------------- 
%
Chi = S.Chi;                     % Coefficient matrix
PO = S.PO;                       % Polynomial order
%
Y = Y';                          % States of the system
Nx = size(Chi,2);                % Number of sttes
%
Theta = utils.SINDy.poolpolyData(Y,Nx,PO);   % Library of functions in terms of states at time instant t
%
dY = (Theta*Chi)';               % Derivatives of the states of the system
%
end