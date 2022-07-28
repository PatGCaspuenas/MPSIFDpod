function SI = reconstructSI(SNP,SI,SVDtr,SNPM,param)
%%                         reconstructSI.m
%--------------------------------------------------------------------------
%
% From SI states with temporal resolution, retrieves the original basis
% back, that is, the velocity fields. Integration of pressure is given as
% an option
%
% INPUTS
%
%   SNP      : structure containing, among others, mean velocity fields
%   SI       : structure containing states of the system with temporal
%              resolution from SINDy-BFI process
%   SVDtr    : structure containing spatial truncated modes
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   SI       : structure adding reconstructed velocity and pressure fields,
%              if selected
%
% UTILS
%
%   utils.pressure.pIntegratorSNPM
%   utils.pressure.getPGradient
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    if isfield(param,'IntegratePressure')
        %
        IntPressure = param.IntegratePressure;                   % choose to retrieve pressure (1) or not (0)
        %
    else
        %
        IntPressure = 0;
        %
    end
    % param.PressureboundaryCondition from utils.pressure.pIntegratorSNPM
    % param.SaveLaplacianFields from utils.pressure.getPGradient
    % param.SaveGradientVelFields from utils.pressure.getPGradient
    m = SNPM.m;                                                   % number of columns in snapshot
    n = SNPM.n;                                                   % number of rows in snapshot

%
% Rerconstructed velocity fields
%

    if isfield(SNP,'vdt') && isfield(SNP,'udt')
        %
        uvdt = SVDtr.Ur*SI.X';
        %
        SI.udt = uvdt(1:m*n,:);
        SI.vdt = uvdt((1+m*n):end,:);
        %
        SI.u = SI.udt + SNP.u_m;
        SI.v = SI.vdt + SNP.v_m;
        %
        duv = SVDtr.Ur*SI.dX';
        SI.du = duv(1:m*n,:);
        SI.dv = duv((1+m*n):end,:);
        %
    elseif isfield(SNP,'vdt')
        %
        SI.vdt = SVDtr.Ur*SI.X';
        SI.v = SI.vdt + SNP.v_m;
        SI.dv = SVDtr.Ur*SI.dX';
        %
    else
        %
        SI.udt = SVDtr.Ur*SI.X';
        SI.u = SI.udt + SNP.u_m;
        SI.du = SVDtr.Ur*SI.dX';
        %
    end

%
% Reconstructed pressure field (if selected)
%

    if IntPressure
        %
        SI = utils.pressure.getPGradient(SI,SNPM,param);
        SI.p = utils.pressure.pIntegratorSNPM(SI,SNPM,param);
        %
    end

%
end