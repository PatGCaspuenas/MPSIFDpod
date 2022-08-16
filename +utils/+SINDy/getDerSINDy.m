function dX = getDerSINDy(SINDytr,X)
%%                          getDerSINDy.m
%--------------------------------------------------------------------------
%
% Retrieves derivatives of states of the system when known the dynamical
% system and their coefficients
%
% INPUTS
%
%   SINDytr  : structure containing matrix of ceofficients of the dynamical
%              system, as well as polynomial order of the library
%   X        : states of the system, with dimensions Nt x Nr
%
% OUTPUT
%
%   dX       : derivatives of the states of the system, with dimensions Nt x Nr   
%
% UTILS
%
%   utils.SINDy.poolpolyData
%
%-------------------------------------------------------------------------- 
%

    dX = X.*0;
    %
    for i = 1:SINDytr.Nx
        %
        Di = strcat('D',num2str(i));
        %
        Xi = X(:,1:SINDytr.Nx);
        Theta = utils.SINDy.poolpolyData(Xi,SINDytr.Nx,SINDytr.(Di).PO);
        %
        dX(:,i) = Theta*SINDytr.(Di).Chi;
        %
    end

%
end