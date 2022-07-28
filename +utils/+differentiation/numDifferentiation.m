function du = numDifferentiation(u,t)
%%                        numDifferentiation.m
%--------------------------------------------------------------------------
%
% Retrieves temporal derivative of a given matrix, whose columns advance
% irregularly in time. Performs a central/forward/backward finite
% difference method with second order errors O(Delta t ^2) 
%
% INPUTS
%
%   u  : matrix with dimensions Np x Nt
%   t  : time vector of irregular or regular spacing
%
% OUTPUT
%
%   du : temporal derivative of dimensions Np x Nt
%
%-------------------------------------------------------------------------- 
%

for i = 1:size(u,2)
    %
    if i == 1                                             % at start point
        %
        Dti = t(i+1) - t(i);                              % time spacing between instant i+1 and i
        Dtii = t(i+2) - t(i+1);                           % time spacing between instant i+2 and i+1
        %
        b = 1/Dti + 1/Dtii;                               % differentiation coefficient scheme 
        c = -Dti/Dtii/(Dti + Dtii);
        a = -b - c;
        %
        du(:,i) = c*u(:,i+2) + b*u(:,i+1) + a*u(:,i); 
        %
    elseif i == size(u,2)                                 % at end point
        %
        Dti = t(i-1) - t(i-2); 
        Dtii = t(i) - t(i-1);
        %
        b = -1/Dti - 1/Dtii;
        c = Dtii/Dti/(Dti + Dtii);
        a = -b - c;
        %
        du(:,i) = c*u(:,i-2) + b*u(:,i-1) + a*u(:,i);
        %        
    else                                                  % at middle points
        %
        Dt = t(i) - t(i-1); 
        Dti = t(i+1) - t(i);
        %
        du(:,i) = 1/(Dt + Dti)*(...
            (u(:,i) - u(:,i-1))*Dti/Dt + ...
            (u(:,i+1) - u(:,i))*Dt/Dti);
        %
    end
    %
end