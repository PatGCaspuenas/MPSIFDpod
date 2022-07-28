function BFI = getBFI(SNP,SINDytr,SVDtr,SNPM,param)
%%                              getBFI.m
%--------------------------------------------------------------------------
%
% Performs a physically-informed integration, that is, a weighted
% backward-forward integration process in between available snapshots.
% Leverages the retrieved dynamical system converged from the SINDy
% training process
%
% INPUTS
%
%   SNP      : structure containing, among others, velocity and
%              temporal fields with no temporal resolution
%   SINDytr  : structure containing matrix of ceofficients of the dynamical
%              system
%   SVDtr    : structure containing spatial truncated modes
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   BFI      : structure containing states of the system with temporal
%              resolution
%
% UTILS
%
%   utils.integration.(param.odeIntegratorOption)
%   utils.fits.CubicSpline
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    odefun = eval(strcat('@',param.odeIntegratorOption));      % Integrator function handle. Should  
                                                               % always have as inputs (t,y,s)
                                                               % with s being a structure
    ts = SNP.Dt;                                               % Time separation of NTR set 
    BFI.Dt = SNPM.Dt;                                          % Time separation objective: TR

%
% Initial conditions (snapshots and time instants)
%
    
    if isfield(SNP,'vdt') && isfield(SNP,'udt')
        %
        FLOW = [ SNP.udt; SNP.vdt ];
        %
    elseif isfield(SNP,'vdt')
        %
        FLOW = SNP.vdt;
        %
    else
        %
        FLOW = SNP.udt;
        %
    end
    %
    BFI.Xi = (SVDtr.Ur'*FLOW)';                                % States of system used as ICs
    BFI.ti = SNP.t;                                            % Time instants of ICs

%
% Weighting function
%
    xi = [0; ts/2; ts];
    yi = [0; 0.5; 1];
    CS = utils.fits.CubicSpline(xi,yi);                        % Spline interpolation for weights in BFI
    w = CS(0:BFI.Dt:ts);                                       % Vector of weights such that  
                                                               % w = 0 at ti
                                                               % w = 1 at ti+1 
                                                               % w = 0.5 at 1/2*(ti+t+1)

%
% Start piece-wise integration for FI, BI and BFI
%

    options = odeset('RelTol',1e-9,'AbsTol',1e-9*ones(1,SINDytr.Nx));
    %
    BFI.f.t = [];
    BFI.b.t = [BFI.ti(1)];
    BFI.bf.t = [];
    BFI.f.X = [];
    BFI.b.X = [BFI.Xi(1,:)];
    BFI.bf.X = [];
    %
    for i = 1:(size(BFI.Xi,1) - 1)
        %
        % Forward integration
        %
    
            [ tf, xf ]= ode45(@(t,Y) odefun(t,Y,SINDytr),...
                BFI.ti(i):BFI.Dt:BFI.ti(i+1),BFI.Xi(i,:),options);
    
        %
        % Backward integration
        %
    
            [ tb, xb ]= ode45(@(t,Y) odefun(t,Y,SINDytr),...
                BFI.ti(i+1):-BFI.Dt:BFI.ti(i),BFI.Xi(i+1,:),options);
    
        %
        % If integration diverges and stops before, set states as null values
        %
    
            if length(tf) < length(w)
                %
                xf((length(tf)+1):length(w),:) = 0;
                tf = (BFI.ti(i):BFI.Dt:BFI.ti(i+1))';
                %
            end
            if length(tb) < length(w)
                %
                xb((length(tb)+1):length(w),:) = 0;
                tb = (BFI.ti(i+1):-BFI.Dt:BFI.ti(i))';
                %
            end
    
        %
        % Weight both contributions : Gaussian distribution
        %
    
            tbf = tf;
            xbf = flip(xb).*w' + (1 - w').*xf;
    
        %
        % Save results advanced in time
        %
    
            BFI.f.t  = [BFI.f.t ; tf(1:end-1,:)];
            BFI.b.t  = [BFI.b.t ; flip(tb(1:end-1))];
            BFI.bf.t = [BFI.bf.t; tbf(1:end-1,:)];
            %
            BFI.f.X  = [BFI.f.X ; xf(1:end-1,:)];
            BFI.b.X  = [BFI.b.X ; flip(xb(1:end-1,:))];
            BFI.bf.X = [BFI.bf.X; xbf(1:end-1,:)];
    
        %
    end

%
% Get final values
%

    BFI.f.X  = [BFI.f.X; BFI.Xi(end,:)];
    BFI.f.t  = [BFI.f.t; BFI.ti(end)];
    BFI.bf.X = [BFI.bf.X; BFI.Xi(end,:)];
    BFI.bf.t = [BFI.bf.t; BFI.ti(end)];

%
end