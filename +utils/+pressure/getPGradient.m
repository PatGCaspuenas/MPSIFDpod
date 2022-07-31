function M = getPGradient(M,SNPM,param)
%%                        getPGradient.m
%--------------------------------------------------------------------------
%
% Retrieves gradient of snapshot matrix with 2nd order errors
%
% INPUTS
%
%   M        : structure containing, among others, velocity and
%              acceleration fields
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   M        : structure adding pressure gradient fields
%
% UTILS
%
%   utils.differentiation.Laplacian2D
%   utils.differentiation.getGradient
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    if isfield(param,'SaveLaplacianFields')                 % chooses to save or not laplacian
        SaveLaplacian = param.SaveLaplacianFields;
    else
        SaveLaplacian = 0;
    end
    %
    if isfield(param,'SaveGradientVelFields')               % chooses to save or not velocity gradients
        SaveGradientVel = param.SaveGradientVelFields;
    else
        SaveGradientVel = 0;
    end
    %
    m = SNPM.m;                                             % columns of snapshot (length in x)
    n = SNPM.n;                                             % rows of snapshot (length in y)
    Dx = SNPM.Dx;                                           % regular x spacing
    Dy = SNPM.Dy;                                           % regular y spacing
    Nt = size(M.u,2);                                       % time samples
    nu = SNPM.nu;                                           % non-dimensional kinematic viscosity, equivalent to 1/Re

%
% Laplacian fields
%

    D2u = []; D2v = [];
    %
    for it = 1:Nt
        %
        [ d2u,~,~ ] = utils.differentiation.Laplacian2D(reshape(M.u(:,it), [ m n ]), Dx, Dy);
        [ d2v,~,~ ] = utils.differentiation.Laplacian2D(reshape(M.v(:,it), [ m n ]), Dx, Dy);
        D2u(:,it) = reshape(d2u, [ m*n 1 ]);
        D2v(:,it) = reshape(d2v, [ m*n 1 ]);
        %
    end
    %
    if SaveLaplacian
        M.D2u = D2u;
        M.D2v = D2v;
    end

%
% Spatial derivatives
%

    var_str = {'dudx','dudy','dvdx','dvdy'};                % check if gradients are already available
    %
    for i = 1:length(var_str)
        %
        if ~isfield(M,var_str{i})
            %
            dudx = []; dudy = []; dvdx = []; dvdy = []; 
            %
            for it = 1:size(M.u,2)
                %
                [ Dudx,Dudy ] = utils.differentiation.getGradient(M.u(:,it),SNPM);
                %
                dudx(:,it) = reshape(Dudx, [ m*n 1 ]);
                dudy(:,it) = reshape(Dudy, [ m*n 1 ]);
                %
                [ Dvdx, Dvdy] = utils.differentiation.getGradient(M.v(:,it),SNPM);
                %
                dvdx(:,it) = reshape(Dvdx, [ m*n 1 ]);
                dvdy(:,it) = reshape(Dvdy, [ m*n 1 ]);
                %
            end
            %
            if SaveGradientVel
                %
                M.dudx = dudx;
                M.dudy = dudy;
                M.dvdx = dvdx;
                M.dvdy = dvdy;
                %
            end
            %
        else
            %
            dudx = M.dudx;
            dudy = M.dudy;
            dvdx = M.dvdx;
            dvdy = M.dvdy;
            %
        end
        %
    end

%
% Pressure gradient
%

    M.dpdx = nu*D2u - M.du ...
                        - M.u.*dudx - M.v.*dudy;
    M.dpdy = nu*D2v - M.dv ...
                        - M.u.*dvdx - M.v.*dvdy;

%
end