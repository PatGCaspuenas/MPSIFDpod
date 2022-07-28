function E = getEnergyFluctuations(a)
%%                       getEnergyFluctuations.m
%--------------------------------------------------------------------------
%
% Computes energy fluctuations temporal history from temporal mode set.
% Definition retrieved from Rubini et al (2022) in "A priori
% sparsification of Galerkin models"
%
% INPUT
%
%   a : temporal mode set of dimensions Nt x Nr
%
% OUTPUT
%
%   E : energy fluctuations of dimensions Nt x 1
%
%--------------------------------------------------------------------------
%
% Parameters
%

    Nr = size(a,2);                    % number of modes

%
% Kinetic energy history of each mode
%

    for i = 1:Nr
        %
        Ei(:,i) = 1/2*a(:,i).^2;
        %
    end

%
% Total energy fluctuations in time
%

    E = sum(Ei,2);

%
end