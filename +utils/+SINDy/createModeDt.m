function SINDy = createModeDt(SINDy,SVD)
%%                          createModeDt.m
%--------------------------------------------------------------------------
%
% Arranges structure for SINDy optimization process of each system state
%
% INPUTS
%
%   SINDy     : SINDy structure. No relevant information is contained, it
%               is introduced not to loose the existing variables
%   SVD       : structure containing spatial and temporal modes, with their
%               singular values
%
% OUTPUT
%
%   SINDy     : update structure with derivative and values of each state
%               of the system 
%
%--------------------------------------------------------------------------
%

Nx = SVD.r;
SINDy.Nx = Nx;                     % Number of states of the system
%

for i = 1:SINDy.Nx
    %
    Di = strcat('D',num2str(i));
    SINDy.(Di).dX = SVD.dar(:,i);  % Derivative of the state
    SINDy.(Di).X = SVD.ar(:,i);    % State values
    %  
end