function TH = reconstructTH(SNP,TH,SVDtr,SNPM,param)
%%                         reconstructTH.m
%--------------------------------------------------------------------------
%
% From TH propagated horizontal and vertical velocity fields with temporal
% resolution, extends reconstruction to velocity fluctuations, states of
% the system, acceleration fields and pressure
%
% INPUTS
%
%   SNP      : structure containing, among others, mean velocity fields
%   TH       : structure containing horizontal and vertical velocity fields
%              with temporal resolution, of dimensions Np x Nt
%   SVDtr    : structure containing spatial truncated modes
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   TH       : structure adding velocity fluctuations, states of system,
%              acceleration and pressure fields
%
% UTILS
%
%   utils.pressure.pIntegratorSNPM
%   utils.pressure.getPGradient
%   utils.differentiation.numDifferentiation
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

%
% Reconstructed velocity fields and truncated temporal mode set
%

    TH.udt = TH.u - SNP.u_m;     
    TH.vdt = TH.v - SNP.v_m;
    TH.ar = (SVDtr.Ur'*[ TH.udt; TH.vdt ])';

%
% Reconstructed pressure field (if selected)
%

    if IntPressure
        %
        TH.du = utils.differentiation.numDifferentiation(TH.u,TH.t);
        TH.dv = utils.differentiation.numDifferentiation(TH.v,TH.t);
        %
        TH = utils.pressure.getPGradient(TH,SNPM,param);
        TH.p = utils.pressure.pIntegratorSNPM(TH,SNPM,param);
        %
    end

%     
end