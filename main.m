%% Multi-Pulse System Identification of Fluid Dynamics with POD basis
%
% Implementation of methodology proposed in the Bachelor Thesis "Non-linear
% system identification for aerodynamic flows" by Patricia
% García-Caspueñas, supervised by Stefano Discetti.
%
% This code serves as an example of the data preparation, training and
% testing processes of the Fluidic Pinball configuration. Can be easily
% generalized for the PIV0018 dataset.
%
% This code requires in input a structure of parameters for some given
% functions. These parameters are explained inside each function.
% Nevertheless, for a better understanding there exists an interactive
% table containing all the required parameters along the process, found in:
%
% https://perfect-vibraphone-20c.notion.site/Code-8998ca572e30440f851aedbde6312d18
%
% Enjoy! =)
%
clc; clear all; close all;
warning('off','all');
%
%% I. DATA PREPARATION
%



%
% General settings
%

    utils.gen_settings.plot_config;                               % default plot settings

%
% Data loading
%

    TrS.TR = utils.data_settings.createTRSet('.\+data\FP_14k_24k.mat');    % create time-resolved training set
    TS.TR  = utils.data_settings.createTRSet('.\+data\FP_10k_13k.mat');    % create time-resolved testing set
    SNPM   = utils.data_settings.prepareSNPM('FP');                % general grid settings

%
% Create NTR data for training
%

    param.ts = 1;                                                 % average time separation between snapshots
    param.FlagNTRTimeSpacing = 'irr';                             % flag to determine regular/irregular spacing
    %
    TrS.NTR = utils.data_settings.createNTRSet(TrS.TR,SNPM,param);% non-time resolved training set
    TrS.TR = [];                                                  % free unnecessary space from time-resolved training

%
% Create TR testing set and compute velocity fluctuations wrt training mean
%

    param.u_m = TrS.NTR.SNP.u_m;                                       % mean horizontal velocity field from training
    param.v_m = TrS.NTR.SNP.v_m;                                       % mean vertical velocity field from training
    param.FlagNTRTimeSpacing = 'r';                                        % flag to determine regular/irregular spacing
    param.ts = SNPM.Dt;                                                % average time separation between snapshots (TR separation)
    %
    TS.TR =utils.data_settings.createNTRSet(TS.TR,SNPM,param);         % time resolved testing set
    
%



%
%% II. TRANSFORMATION INTO REDUCED-ORDER MODEL
%



%
% Perform POD on NTR training set and truncate it
%

    TrS.NTR.SVD = utils.POD.performPOD(TrS.NTR.SNP);                   % POD
    %
    param.TruncationMethodPOD = {'elbow',''};                          % elbow of the cumulative energy curve to choose Nr
    TrS.NTR.SVD = utils.POD.truncatePOD(TrS.NTR.SNP,TrS.NTR.SVD,param); % truncated POD

%



%
%% III. SYSTEM IDENTIFICATION OF POD TEMPORAL MODES DYNAMICS
%



%

param.OptimizationMethodSINDy = {'ALASSO','FP'};
TrS.NTR.SINDy = utils.SINDy.SINDy_config(TrS.NTR.SVD,param);           % SINDy settings

%
% Training to retrieve each temporal mode dynamics
%

    poly_order = 0:2;                                                  % polynomial order of library matrix
    TrS.NTR.SINDy.PO = poly_order;
    %
    TrS.NTR.SINDy = utils.SINDy.createModeDt(TrS.NTR.SINDy,TrS.NTR.SVD); % arrange structure to include each mode setting
    %
    % Look for optimal coefficients of each mode derivative function
    %

        for j = 1:TrS.NTR.SVD.r
            %
            j
            Dj = strcat('D',num2str(j));
            %
            TrS.NTR.SINDy.(Dj).PO = poly_order;
            %
            % Create library
            %
    
                TrS.NTR.SINDy.(Dj).Theta = ...
                    utils.SINDy.poolpolyData(TrS.NTR.SINDy.X,...
                    TrS.NTR.SINDy.Nx,...
                    TrS.NTR.SINDy.(Dj).PO);
    
            %
            % Find parsimonious model 
            %
    
                tic
                [TrS.NTR.SINDy.(Dj).Chi,...
                    TrS.NTR.SINDy.(Dj).alpha,...
                    TrS.NTR.SINDy.(Dj).delta,...
                    TrS.NTR.SINDy.(Dj).R2,...
                    TrS.NTR.SINDy.(Dj).Nact] = ...
                    utils.SINDy.paretoElbowl1(TrS.NTR.SINDy.(Dj).Theta,...
                    TrS.NTR.SINDy.(Dj).dX,...
                    TrS.NTR.SINDy.alpha,...
                    TrS.NTR.SINDy.delta, ...
                    TrS.NTR.SINDy.ftol,...
                    TrS.NTR.SINDy.tol);
                toc
            %
            % Create global matrix of coefficients
            %

                TrS.NTR.SINDy.Chi(:,j) = TrS.NTR.SINDy.(Dj).Chi;

        end

    %

%



%
%% IV. TEMPORAL MODE SET INTEGRATION AND VELOCITY AND PRESSURE RECONSTRUCTION
%



%
% Create ICs snapshots
%

    % param.u_m & param.v_m from before
    param.ts = 5.3;                                                     % Time separation between available snapshots
    param.FlagNTRTimeSpacing = 'r';                                   % Regular time separation
    %
    TS.NTR = utils.data_settings.createNTRSet(TS.TR.SNP,SNPM,param);  % non-time-resolved testing set

%
% Integrate ICs and reconstruct
%

    param.odeIntegratorOption = 'utils.integration.odeSINDyGlobal';                     % choose the integrator
    param.IntegratePressure = 1;                                      % integrates pressure gradient in space or not
    param.PressureBoundaryCondition = {1,TS.TR.SNP.p(1,:),1,1};       % boundary conditions options
    %
    TS.NTR.BFI = utils.integration.getBFI(...
        TS.NTR.SNP,TrS.NTR.SINDy,TrS.NTR.SVD,SNPM,param);
    TS.NTR.BFI.bf = utils.integration.reconstructBFI(...
        TS.TR.SNP,TS.NTR.BFI.bf,TrS.NTR.SVD,TrS.NTR.SINDy,SNPM,param);

%



%
%% V. COMPARISON MODELS
%



%
% Time-resolved truncates testing
%

    % param.IntegratePressure & param.PressureBoundaryCondition from before
    param.TruncateFlowAcceleration = 1;                               % obtain truncated acceleration
    param.TruncateFlowAccelerationNumDiff = 0;                        % obtain truncated acceleration from numerical differentiation
    param.TruncateFlowAccelerationSINDy = 0;                          % obtain truncated acceleration from SINDy's dynamical system
    %
    [TS.TR.SNPr,~] = utils.POD.truncateTRorNTR(TS.TR.SNP,TrS.NTR.SVD,SNPM,0,param);
    TS.TR.SVD.ar = (TrS.NTR.SVD.Ur'*[TS.TR.SNP.udt;TS.TR.SNP.vdt])';  % truncated time-resolved temporal mode set (ground truth)

%
% Spline Interpolated model
%

    % param.IntegratePressure & param.PressureBoundaryCondition from before
    %
    TS.NTR.SI = utils.comparison.getSI(TS.NTR.SNP,TrS.NTR.SVD,SNPM);  % cubic spline interpolated temporal modes
    TS.NTR.SI = utils.comparison.reconstructSI(TS.TR.SNP,TS.NTR.SI,TrS.NTR.SVD,SNPM,param); % cubic spline interpolated velocity and pressure fields

%
% Taylor Hypothesis model
%

    % param.IntegratePressure & param.PressureBoundaryCondition from before  
    param.TaylorHypothesisMask = ones(size(SNPM.X));                   % mask to avoid advection of body BCs
    param.TaylorHypothesisMask( SNPM.X > 0.5 ) = 0;
    %
    TS.NTR.TH = utils.comparison.getTH(TS.NTR.SNP,TrS.NTR.SVD,SNPM,param);
    TS.NTR.TH = utils.comparison.reconstructTH(TS.TR.SNP,TS.NTR.TH,TrS.NTR.SVD,SNPM,param);

%
% Make reference DNS and retrieved models of equal length
%
    
    endt = find( TS.TR.SNP.t == TS.NTR.TH.t(end) ,1,'last');
    TS.TR.SNPr = utils.data_settings.removeFinalInstants(TS.TR.SNPr,TS.NTR.TH);
    TS.TR.SNP = utils.data_settings.removeFinalInstants(TS.TR.SNP,TS.NTR.TH);
    TS.TR.SVD.ar = TS.TR.SVD.ar(1:endt,:);
%



%
%% VI. ERROR DETERMINATION (section included to show error functions performances)
%



%

    % param.ts from before
    param.ErrorTypeRMSE = 'w';                                         % single-value error
    param.ComputeRMSEatMiddle = 1;                                     % compute error at middle point between available ICs snapshots

%
% Root Mean Square Error between truncated DNS reference and BFI/SI/TH
%

    err_ptrBFI = utils.error_eval.computeRMSE(...
                 TS.TR.SNPr.p,TS.NTR.BFI.bf.p,SNPM,param); 
    err_ptrSI = utils.error_eval.computeRMSE(...
                 TS.TR.SNPr.p,TS.NTR.SI.p,SNPM,param);
    err_ptrTH = utils.error_eval.computeRMSE(...
                 TS.TR.SNPr.p,TS.NTR.TH.p,SNPM,param);

%
% R-squared coefficients between truncated DNS reference temporal mode set
% and integrated temporal modes with BFI/SI/TH
%

    R2_BFI  = utils.error_eval.getR2factor(...
              'c',TS.TR.SVD.ar,TS.NTR.BFI.bf.X);
    R2_SI  = utils.error_eval.getR2factor(...
              'c',TS.TR.SVD.ar,TS.NTR.SI.X);
    R2_TH  = utils.error_eval.getR2factor(...
              'c',TS.TR.SVD.ar,TS.NTR.TH.ar);

%
