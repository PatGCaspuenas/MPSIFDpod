function p = pIntegratorSNPM(SNP,SNPM,param)
%%                        pIntegratorSNPM.m
%--------------------------------------------------------------------------
%
% Integrates pressure gradient snapshot matrix and imposes a single point
% or region Boundary Condition
%
% INPUTS
%
%   SNP      : structure containing, among others, pressure gradient fields
%              with or without temporal resolution
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   p        : integrated pressure fields from gradient
%
% UTILS
%
%   utils.pressure.pIntegrator
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    ImposeBC = param.PressureBoundaryCondition{1};        % Impose Boundary Condition (1) or not (0)
    p0 = param.PressureBoundaryCondition{2};              % BC history of selected region of dimensions Nx x Nt
    ix = param.PressureBoundaryCondition{3};              % x region of dimensions 1 x Nx (indexes of X grid)
    iy = param.PressureBoundaryCondition{4};              % y region of dimensions Nx x 1 (indexes of Y grid)
    Mask = SNPM.BD.M;                                     % body mask
    Dx = SNPM.Dx;                                         % regular x spacing
    Dy = SNPM.Dy;                                         % regular y spacing

%
% Boundary Condition region
%

    for i = 1:length(ix)
        %
        ireg(i) = sub2ind(size(SNPM.X),iy(i),ix(i));
        %
    end

%
% Integration of each time instant
%

    p = SNP.dpdx.*0;
    %
    for it = 1:size(SNP.dpdx,2)
        %
        px = reshape(SNP.dpdx(:,it),size(Mask));
        py = reshape(SNP.dpdy(:,it),size(Mask));
        %
        if it == 1                                       % initial guess is a null matrix
            %
            p(:,it) = utils.pressure.pIntegrator(px,py,0*px,Mask,Dx,Dy);
            %
        else                                             % initial guess is the previous solution
            %
            p(:,it) = utils.pressure.pIntegrator(px,py,reshape(p(:,it-1),size(Mask)),Mask,Dx,Dy);
            %
        end
        %
        % Correction of BC
        %
    
            if ImposeBC
                %
                p(:,it) = p(:,it) + ones(length(ireg),1)\( -p(ireg,it) + p0(:,it) );
                %
            end
        
        %
    end

%     
end