function TH = getTH(SNP,SNPM,param)
%%                            getTH.m
%--------------------------------------------------------------------------
%
% Propagates velocity fields backward and forward in time
% (interpolated/integration process). Follows procedure of Scarano&Moore
% (2012) in "An advecytion-based model to increase the temporal resolution
% of PIV time series" with some modifications. Valid only for available
% horizontal and vertical velocity fields
%
% INPUTS
%
%   SNP      : structure containing, among others, velocity and
%              temporal fields with no temporal resolution
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   TH       : structure adding temporal resolyution to velocity fields
%              (horizontal + vertical)
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    if isfield(param,'TaylorHypothesisMask')                   % Mask to include in case
                                                               % body is inside the snapshot
                                                               % region (avoid advection of
                                                               % body null values)
        %
        Mask = param.TaylorHypothesisMask;
        %
        j = find( Mask(1,:) == 1 ,1,'last');                   % Smooth transition between 
                                                               % mask and no-mask regions
        Nc = 8;
        Mask(:,j-Nc:j) = (0.9:-0.1:0.1).*ones(m,Nc+1);
        %
    else
        %
        Mask = zeros(size(X));
        %
    end
    %
    ts = SNP.Dt;                                               % Time separation of NTR set 
    Dts = SNPM.Dt;                                             % Time separation of TR set 
    Nti = length(SNP.t);                                       % Time samples
    X = SNPM.X;                                                % X grid
    Y = SNPM.Y;                                                % Y grid
    xmin = min(X,[],"all");                                    % X/Y grid minimum and maximum values
    xmax = max(X,[],"all");
    ymin = min(Y,[],"all"); 
    ymax = max(Y,[],"all");

%
% Initial conditions for that time separation ts
%

    TH.Dt = Dts;                                               % Time separation objective: TR
    TH.ti = SNP.t;                                             % Time instants at ICs
    %
    TH.ui = [ SNP.u ];                                         % Snapshots used as initial conditions
    TH.vi = [ SNP.v ];  
    %
    Ui = reshape(TH.ui,[size(X) Nti]);
    Vi = reshape(TH.vi,[size(X) Nti]);

%
% Start propagation/integration process
%

    c = 0;                                                    % Counting variable
    %
    for i = 1:(size(TH.ui,2) - 1)
        %
        t = TH.ti(i):TH.Dt:TH.ti(i+1);                        % temporal vector with time 
                                                              % resolution in between ICs
        %
        Uc1 = imgaussfilt(Ui(:,:,i),7);                       % advection velocity for time instant i corresponding to horizontal velocity
                                                              % a local filter is applied each space in between ICs
                                                              % to get the largest fluctuations, responsible for advection
        Uc3 = imgaussfilt(Ui(:,:,i+1),7);                     % advection velocity for time instant i+1 corresponding to horizontal velocity
        %
        Vc1 = imgaussfilt(Ui(:,:,i),7);                       % advection velocity for time instant i corresponding to vertical velocity     
        Vc3 = imgaussfilt(Ui(:,:,i+1),7);                     % advection velocity for time instant i+1 corresponding to vertical velocity  
        %
        for j = 1:length(t)-1                                 % propagation process
            %
            c = c + 1;
    
    
            X1 = X + (t(1)-t(j))*Uc1;                         % points in the grid where the values at (X,Y) have propagated
            Y1 = Y + (t(1)-t(j))*Vc1;  
            %
            X3 = X + (t(end)-t(j))*Uc3;     
            Y3 = Y + (t(end)-t(j))*Vc3;
            %
            X1(X1 < xmin) = xmin;           X1(X1 > xmax) = xmax;  % if propagated points are out of limits, resterict them to the limits
            X3(X3 < xmin) = xmin;           X3(X3 > xmax) = xmax;
            Y1(Y1 < ymin) = ymin;           Y1(Y1 > ymax) = ymax;
            Y3(Y3 < ymin) = ymin;           Y3(Y3 > ymax) = ymax;
    
    
            
            U1 = interp2(X,Y,Ui(:,:,i),X1,Y,'spline',0);       % Interpolate velocity at propagated points at time instants i or i+1
            U1 = U1 + (Ui(:,:,i) - U1).*Mask;                  % Maintain region of mask as available time instants
            %
            U3 = interp2(X,Y,Ui(:,:,i+1),X3,Y,'spline',0);    
            U3 = U3 + (Ui(:,:,i+1) - U3).*Mask; 
            %
            V1 = interp2(X,Y,Vi(:,:,i),X1,Y,'spline',0);      
            V1 = V1 + (Vi(:,:,i) - V1).*Mask; 
            %
            V3 = interp2(X,Y,Vi(:,:,i+1),X3,Y,'spline',0);    
            V3 = V3 + (Vi(:,:,i+1) - V3).*Mask;
            
    
            
            U = (t(end) - t(j))/ts*U1 + (t(j) - t(1))/ts*U3;   % Weighted forward and backward predictions
            V = (t(end) - t(j))/ts*V1 + (t(j) - t(1))/ts*V3;
            
    
            
            TH.u(:,c) = U(:);
            TH.v(:,c) = V(:);
            TH.t(c) = t(j);
            
    

        end
    end

%
% Final values
%

    TH.u(:,c+1) = TH.ui(:,end); 
    TH.v(:,c+1) = TH.vi(:,end);
    TH.t(c+1) = TH.ti(end);

end