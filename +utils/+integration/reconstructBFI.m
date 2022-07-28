function BFI = reconstructBFI(SNP,BFI,SVDtr,SINDytr,SNPM,param)
%%                         reconstructBFI.m
%--------------------------------------------------------------------------
%
% From BFI states with temporal resolution, retrieves the original basis
% back, that is, the velocity fields. Integration of pressure is given as
% an option
%
% INPUTS
%
%   SNP      : structure containing, among others, mean velocity fields
%   BFI      : structure containing states of the system with temporal
%              resolution from SINDy-BFI process
%   SVDtr    : structure containing spatial truncated modes
%   SINDytr  : structure containing matrix of ceofficients of the dynamical
%              system
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   BFI      : structure adding reconstructed velocity and pressure fields,
%              if selected
%
% UTILS
%
%   utils.SINDy.getDerSINDy
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
% Derivatives of the reconstructed states
%

    BFI.dX = getDerSINDy(SINDytr,BFI.X);

%
% Rerconstructed velocity fields
%

    if isfield(SNP,'vdt') && isfield(SNP,'udt')
        %
        duv = SVDtr.Ur*BFI.dX';
        uvdt = SVDtr.Ur*BFI.X';
        %
        BFI.udt = uvdt(1:m*n,:);
        BFI.vdt = uvdt((1+m*n):end,:);
        %
        BFI.u = BFI.udt + SNP.u_m;
        BFI.v = BFI.vdt + SNP.v_m;
        %
        BFI.du = duv(1:m*n,:);
        BFI.dv = duv((1+m*n):end,:);
        %
    elseif isfield(SNP,'vdt')
        %
        BFI.vdt = SVDtr.Ur*BFI.X';
        BFI.v = BFI.vdt + SNP.v_m;
        BFI.dv = SVDtr.Ur*BFI.dX';
        %
    else
        %
        BFI.udt = SVDtr.Ur*BFI.X';
        BFI.u = BFI.udt + SNP.u_m;
        BFI.du = SVDtr.Ur*BFI.dX';
        %
    end

%
% Reconstructed pressure field (if selected)
%

    if IntPressure
        %
        BFI = utils.pressure.getPGradient(BFI,SNPM,param);
        BFI.p = utils.pressure.pIntegratorSNPM(BFI,SNPM,param);
        %
    end

%
end