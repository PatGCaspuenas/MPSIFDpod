function [dudx,dudy] = getGradient(u,SNPM)
%%                        getGradient.m
%--------------------------------------------------------------------------
%
% Retrieves gradient of snapshot matrix with 2nd order errors
%
% INPUTS
%
%   u        : snapshot matrix with dimensions Np x Nt
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%
% OUTPUT
%
%   dudx     : gradient in x of u, with dimensions Np x Nt
%   dudy     : gradient in y of u, with dimensions Np x Nt
%
% NOTE - Computational time comparison
%
% Given a snapshot matrix of 101 x 201 x 3001, gradient alone takes 0.6s,
% getGradient 0.73s (for all Nt). Another implementation was tried in
% getGradient, such that the whole Nt snapshots were introduced in the
% function and differentiated with the same scheme + a manual central diff
% scheme. This took 1.3s, thus it was discarded as an option
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    m = size(SNPM.Y,1);            % columns of snapshot (length in x)
    n = size(SNPM.X,2);            % rows of snapshot (length in y)
    Dx = SNPM.Dx;                  % regular x spacing
    Dy = SNPM.Dy;                  % regular y spacing
    Nt = size(u,2);                % time samples

%
% Reshape snapshot matrix and compute gradient at middle points
%

    u = reshape(u,[m n Nt]);
    [dudx,dudy] = gradient(u);     % function of 2nd order errors in the middle
    dudx = dudx/Dx;
    dudy = dudy/Dy;

% 
% Compute gradient in x/y at start and end columns/rows
%

    dudx(:,1,:) = (-3*u(:,1,:) + 4*u(:,2,:) - u(:,3,:))/(2*Dx); % start region
    dudy(1,:,:) = (-3*u(1,:,:) + 4*u(2,:,:) - u(3,:,:))/(2*Dy);
    %
    dudx(:,end,:) = (3*u(:,end,:) - 4*u(:,end-1,:) + u(:,end-2,:))/(2*Dx); % end region
    dudy(end,:,:) = (3*u(end,:,:) - 4*u(end-1,:,:) + u(end-2,:,:))/(2*Dy);

%
% Reshape gradient into original dimension of snapshot matrix
%

    dudx = reshape(dudx,[m*n Nt]);
    dudy = reshape(dudy,[m*n Nt]);