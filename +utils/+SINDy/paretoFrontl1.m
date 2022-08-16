function [Alpha,Delta,R2,Nact] = paretoFrontl1(Theta,dX,alpha,delta,ftol,tol)
%%                          paretoFrontl1.m
%--------------------------------------------------------------------------
%
% For a given state of the system, retrieves accuracy and sparsity of the
% fit for each possible combination of the ALASSO tunning parameters.
%
% INPUTS
%
%   Theta  : library of functions expressed in terms of the states of the system
%   dX     : derivative of one state of the system 
%   alpha  : l1-norm constraint parameter
%   delta  : Adaptive LASSO updating parameter
%   ftol   : inner-loop tolerance
%   tol    : sparsity tolerance
%
% OUTPUTS
%
%   Alpha  : grid of all tried alphas
%   Delta  : grid of all tried deltas
%   R2     : corresponding R-square determination coefficient for each
%            (alpha,delta) combination
%   Nact   : corresponding number of active terms for each
%            (alpha,delta) combination
%
% UTILS
%
%   utils.error_eval.getR2factor
%   utils.SINDy.ALASSO
%--------------------------------------------------------------------------
%
% Parameters
%

    Nx = length(alpha);                            % alpha vector length
    Ny = length(delta);                            % delta vector length
    Nt = size(Theta,1);                            % time samples   

%
% Allocate parameters and results grid
%

    [Alpha,Delta] = meshgrid(alpha,delta);
    R2(:,:) = zeros(Ny,Nx);
    Nact(:,:) = zeros(Ny,Nx);

%
% Separate 80% training 20% testing (to avoid overfitting)
%

    Thetaa = Theta(1:floor(Nt*0.8),:);                  % Training
    dXa = dX(1:floor(Nt*0.8),:);
    Thetab = Theta(floor(Nt*0.8)+1:end,:);              % Testing
    dXb = dX(floor(Nt*0.8)+1:end,:);

%
% Retrieve results for each (alpha,delta) combination
%

 for i = 1:Ny
    for j = 1:Nx
        %
        Chi = zeros(size(Theta,2),1);                   % allocate coefficient vector
        Chi = utils.SINDy.ALASSO(Thetaa, dXa,...
        Alpha(i,j), Delta(i,j),ftol,tol);
        %
        R2(i,j) = utils.error_eval.getR2factor('c',dXb,Thetab*Chi);
        Nact(i,j) = nnz(Chi);                           % number of nonzero terms in Chi
        %
    end
 end

%   
