function  SI = getSI(SNP,SVDtr,SNPM)
%%                              getSI.m
%--------------------------------------------------------------------------
%
% Interpolates states of the system in between available ICs using a cubic
% spline function
%
% INPUTS
%
%   SNP      : structure containing, among others, velocity and
%              temporal fields with no temporal resolution
%   SVDtr    : structure containing spatial truncated modes
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%
% OUTPUT
%
%   SI       : structure containing states of the system with temporal
%              resolution
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    ts = SNP.Dt;                                               % Time separation of NTR set 
    SI.Dt = SNPM.Dt;                                           % Time separation objective: TR

%
% Initial conditions for that time separation ts
%

    if isfield(SNP,'vdt') && isfield(SNP,'udt')
        %
        FLOW = [SNP.udt; SNP.vdt];
        dFLOW = [SNP.du; SNP.dv];
        %
    elseif isfield(SNP,'vdt')
        %
        FLOW = SNP.vdt;
        dFLOW = SNP.dv;
        %
    else
        %
        FLOW = SNP.udt;
        dFLOW = SNP.du;
        %
    end
    %
    SI.Xi = (SVDtr.Ur'*FLOW)';                                 % States used as initial conditions
    SI.dXi = (SVDtr.Ur'*dFLOW)';                               % Derivative of the states used as initial conditions
    SI.ti = SNP.t;                                             % Time instants at ICs
    SI.t = SI.ti(1):SI.Dt:SI.ti(end);                          % Time vector with temporal resolution

%
% Interpolation process
%

    for i = 1:size(SI.Xi,2)
        %
        SI.X(:,i) = spline(SI.ti, SI.Xi(:,i),SI.t);
        SI.dX(:,i) = spline(SI.ti, SI.dXi(:,i),SI.t);
        %
    end

%
end


