%% INSERT TITLE, DESCRIPTION OF THE METHOD, COPYRIGHT ETC ETC


%



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

    TrS.TR = utils.data_setings.createTRSet('FP_14k_24k.mat');    % create time-resolved training set
    TS.TR  = utils.data_setings.createTRSet('FP_10k_13k.mat');    % create time-resolved testing set
    SNPM   = utils.data_settings.prepareSNPM('FP');               % general grid settings

%
% Create NTR data for training
%

    param.ts = 1;                                                 % average time separation between snapshots
    param.FlagNTRTimeSpacing = 'irr';                             % flag to determine regular/irregular spacing
    %
    TrS.NTR = utils.data_setings.createNTRSet(TrS.TR,SNPM,param); % non-time resolved training set
    TrS.TR = [];                                                  % free unnecessary space from time-resolved training

%
% Create TR testing set and compute velocity fluctuations wrt training mean
%

    param.u_m = TrS.NTR.SNP.u_m;                                       % mean horizontal velocity field from training
    param.v_m = TrS.NTR.SNP.v_m;                                       % mean vertical velocity field from training
    param.NTRTimeSpacing = 'r';                                        % flag to determine regular/irregular spacing
    param.ts = SNPM.Dt;                                                % average time separation between snapshots (TR separation)
    %
    TS.TR = createNTRSet(TS.TR,SNPM,param);                            % time resolved testing set
    
%



%
%% II. TRANSFORMATION INTO REDUCED-ORDER MODEL
%



%
% Perform POD on NTR training set and truncate it
%

    TrS.NTR.SVD = utils.POD.performPOD(TrS.NTR.SNP);                   % POD
    %
    param.TruncationMethodPOD = {'elbow',''}; 
    TrS.NTR.SVD = utils.POD.truncatePOD(TrS.NTR.SNP,TrS.NTR.SVD,SNPM,param); % truncated POD

%



%
%% III. SYSTEM IDENTIFICATION OF POD TEMPORAL MODES DYNAMICS
%



%

param.methodSINDy = 'ALASSO';
param.IntegratePressure = 1; param.TruncateFlowAccelerationNumDiff = 0; param.TruncateFlowAcceleration = 1;
param.ErrorTypeRMSE = 'w'; param.mask = SNPM.BD.M(:);
param.ComputeRMSEatMiddle = 1;
param.mean = 0;
param.TruncateFlowAccelerationSINDy = 1;
%


%

TrS.NTR.SINDy = SINDyconfig(TrS.NTR.SVD,param); % SINDy settings
%
%% Training
%
poly_n = 0:2; % polynomial order of library matrix
%

TrS.NTR.SINDy = createModeDt(TrS.NTR.SINDy,TrS.NTR.SVD,param);
%
% Look for optimal coefficients of each mode derivative function
%
for j = 1:TrS.NTR.SVD.r
    %
    j
    Dj = strcat('D',num2str(j));
    %
    TrS.NTR.SINDy.(Dj).PO = poly_n;
    %
    % Create library
    %
    TrS.NTR.SINDy.(Dj).Theta = ...
        poolpolyData(TrS.NTR.SINDy.(Dj).X,...
        TrS.NTR.SINDy.(Dj).Nx,...
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
        paretoElbowl1(TrS.NTR.SINDy.(Dj).Theta,...
        TrS.NTR.SINDy.(Dj).dX,...
        TrS.NTR.SINDy.alpha,...
        TrS.NTR.SINDy.delta,...
        TrS.NTR.SINDy.tol);
    toc
end
filename = strcat('MPIV_FP_t1_O2_I.mat');
% save(fullfile(fname, filename),'TrS','SNPM','param','-v7.3');
%
TrS.TR = [];
%
%% Testing
%
ts_lin = 50; %[ 1 2 4 5 6 8 10 ];
param.IntegratePressure = 1;
param.TruncateFlowAccelerationSINDy = 0;
param.PressureBC = TS.TR.SNP.p(1,:);
param.NeglectIntermediateResults = 1;
param.TaylorHypothesisMask = ones(size(SNPM.X));
param.TaylorHypothesisMask( SNPM.X > 0.5 ) = 0;
%
for i = 1:param.Nr
    %
    Dj = strcat('D',num2str(i));
    TrS.NTR.SINDy.Chi(:,i) = TrS.NTR.SINDy.(Dj).Chi;
    %
end
% [TS.TR.SNPr,~] = truncateTRorNTR(TS.TR.SNP,TrS.NTR.SVD,SNPM,0,param);
%
for k = 1:length(ts_lin)
    %
      param.NTRTimeSpacing = 'r'; param.ts =ts_lin(k);
    %
    TS.NTR = createNTRSet(TS.TR.SNP,SNPM,param); % non-time resolved testing set 
    %
%     TS.NTR.SNPr.Dt = TS.NTR.SNP.Dt;
%     TS.NTR.SNPr.t = TS.NTR.SNP.t;
%     uvr = TrS.NTR.SVD.Ur*(TrS.NTR.SVD.Ur'*TS.NTR.SNP.uvdt) + [TS.NTR.SNP.u_m;TS.NTR.SNP.v_m];
%     TS.NTR.SNPr.u = uvr(1:SNPM.m*SNPM.n,:);
%     TS.NTR.SNPr.v = uvr((1+SNPM.m*SNPM.n):end,:);
%     TS.NTR.SNPr.uvdt = TrS.NTR.SVD.Ur*(TrS.NTR.SVD.Ur'*TS.NTR.SNP.uvdt);
    %
    param.odeOption = 'odeSINDyGlobal';
    tic
    TS.NTR.BFI = getBFI(TS.NTR.SNP,TrS.NTR.SINDy,TrS.NTR.SVD,param);
    toc
%     TS.NTR.BFI.bf = reconstructBFI(TS.TR.SNP,TS.NTR.BFI.bf,TrS.NTR.SVD,TrS.NTR.SINDy,SNPM,param);
%     %
%     TS.NTR.SI = getSI(TS.NTR.SNP,TrS.NTR.SVD,param);
%     TS.NTR.SI = reconstructSI(TS.TR.SNP,TS.NTR.SI,TrS.NTR.SVD,SNPM,param);
%     %
%     TS.NTR.TH = getTHBFI3(TS.NTR.SNPr,SNPM,param);
%     TS.NTR.TH.u = TS.NTR.TH.U(1:SNPM.m*SNPM.n,:);   TS.NTR.TH.v = TS.NTR.TH.U((1+SNPM.m*SNPM.n):end,:);
%     TS.NTR.TH.udt = TS.NTR.TH.u - param.u_m;     TS.NTR.TH.vdt = TS.NTR.TH.v - param.v_m;
%     TS.NTR.TH.uvdt = [TS.NTR.TH.udt; TS.NTR.TH.vdt];
%     %
%     TS.NTR.TH.du = numDifferentiation(TS.NTR.TH.u,TS.NTR.TH.t);
%     TS.NTR.TH.dv = numDifferentiation(TS.NTR.TH.v,TS.NTR.TH.t);
%     %
%     TS.NTR.TH = getPGradient(TS.NTR.TH,SNPM);
%     TS.NTR.TH.p = pIntegratorSNPM(TS.NTR.TH,SNPM,param);
    %
    TS.TR.SVD.ar = (TrS.NTR.SVD.Ur'*TS.TR.SNP.uvdt)';
%     TS.NTR.TH.ar = (TrS.NTR.SVD.Ur'*[ TS.NTR.TH.udt; TS.NTR.TH.vdt ])';
    %
    % RMSE pressure fields
    %
    Nt = size(TS.NTR.BFI.bf.p,2);
    %
    RMSE.ptr(k) = computeRMSE(TS.TR.SNP.p,TS.TR.SNPr.p,SNPM,param);
    RMSE.pBFI(k) = computeRMSE(TS.TR.SNP.p(:,1:Nt),TS.NTR.BFI.bf.p,SNPM,param);
    RMSE.ptrBFI(k) = computeRMSE(TS.TR.SNPr.p(:,1:Nt),TS.NTR.BFI.bf.p,SNPM,param);
    RMSE.pSI(k) = computeRMSE(TS.TR.SNP.p(:,1:Nt),TS.NTR.SI.p,SNPM,param);
    RMSE.ptrSI(k) = computeRMSE(TS.TR.SNPr.p(:,1:Nt),TS.NTR.SI.p,SNPM,param);
    RMSE.ptrTH(k) = computeRMSE(TS.TR.SNPr.p(:,1:Nt),TS.NTR.TH.p,SNPM,param);
    RMSE.pTH(k) = computeRMSE(TS.TR.SNP.p(:,1:Nt),TS.NTR.TH.p,SNPM,param);
    %
    % RMSE velocity fields
    %
    RMSE.uvtr(k) = computeRMSE(TS.TR.SNP.uvdt,TS.TR.SNPr.uvdt,SNPM,param);
    RMSE.uvBFI(k) = computeRMSE(TS.TR.SNP.uvdt(:,1:Nt),TS.NTR.BFI.bf.uvdt,SNPM,param); 
    RMSE.uvtrBFI(k) = computeRMSE(TS.TR.SNPr.uvdt(:,1:Nt),TS.NTR.BFI.bf.uvdt,SNPM,param);
    RMSE.uvSI(k) = computeRMSE(TS.TR.SNP.uvdt(:,1:Nt),TS.NTR.SI.uvdt,SNPM,param);
    RMSE.uvtrSI(k) = computeRMSE(TS.TR.SNPr.uvdt(:,1:Nt),TS.NTR.SI.uvdt,SNPM,param);
    RMSE.uvTH(k) = computeRMSE(TS.TR.SNP.uvdt(:,1:Nt),TS.NTR.TH.uvdt,SNPM,param);
    RMSE.uvtrTH(k) = computeRMSE(TS.TR.SNPr.uvdt(:,1:Nt),TS.NTR.TH.uvdt,SNPM,param);
    %
    % Correlation factor R2
    %
    RMSE.R2.BFI_X(k,:) = getR2factor('c',TS.TR.SVD.ar,TS.NTR.BFI.bf.X,param);
    RMSE.R2.SI_X(k,:) = getR2factor('c',TS.TR.SVD.ar,TS.NTR.SI.X,param);
    RMSE.R2.TH_X(k,:) = getR2factor('c',TS.TR.SVD.ar,TS.NTR.TH.ar,param);
    %
    filename = strcat('MPIV_FP_TS_ts',num2str(ts_lin(k)),'.mat');
    %
    TS.NTR.BFI.bf.dpdx = []; TS.NTR.BFI.bf.dpdy = [];
    TS.NTR.BFI.bf.duv = []; TS.NTR.BFI.bf.uv = [];    TS.NTR.BFI.bf.uvdt = [];
    %
    TS.NTR.SI.dpdx = []; TS.NTR.SI.dpdy = [];
    TS.NTR.SI.duv = []; TS.NTR.SI.uv = [];    TS.NTR.SI.uvdt = [];
    %
    TS_NTR_BFI = TS.NTR.BFI;
    TS_NTR_SI = TS.NTR.SI;
    TS_NTR_TH = TS.NTR.TH;
    save(fullfile(fname, filename),'TS_NTR_BFI','TS_NTR_SI','TS_NTR_TH','RMSE','-v7.3')
    %
    TS.NTR = []; clear TS_NTR_SI; clear TS_NTR_BFI; clear TS_NTR_TH;
end