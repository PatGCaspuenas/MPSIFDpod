function  SNPM = prepareSNPM(flow)
%%                           prepareSNPM.m
%--------------------------------------------------------------------------
%
% Generates structure with geometric and general information about the flow
% configuration
%
% INPUTS
%
%   flow     : string that indicates the type of flow chosen. It can take
%              values
%              -> 'SC'        : Single Cylinder
%              -> 'FP'        : Fluidic Pinball
%              -> 'PIV0018'   : stalled 2D wing (Junwei Chen & Marco Raiola)
%              -> 'PIV0018_2' : stalled 2D wing (Stefano Discetti & Patricia GCaspueñas)
%
% OUTPUTS
%
%   SNPM     : structure containing geometrical and general info about the
%              flow configuration
%
%--------------------------------------------------------------------------
%
%% Main grid parameters
%

    if strcmp(flow,'SC') == 1                   
        %
        load('.\+data\SC_grid.mat');
        %
        SNPM.Dt = 1/30;                         % Non-dimensional time spacing
        SNPM.nu = 1/100;                        % Non-dimensional kinematic viscosity, equivalent to 1/Re
        %
    elseif strcmp(flow,'FP') == 1
        %
        load('.\+data\FP_grid.mat');
        %
        SNPM.Dt = 1/10;                         % Non-dimensional time spacing
        SNPM.nu = 1/130;                        % Non-dimensional kinematic viscosity, equivalent to 1/Re
        %
    elseif strcmp(flow,'PIV0018')
        %
        load('.\+data\PIV0018_grid.mat');
        %
        % Dimensional parameters
        %

        SNPM.D.Rs = 0.1197*1e-3;                % resolution, in m/px
        SNPM.D.Vs = 10;                         % vector spacing, in px
        SNPM.D.c = 0.08;                        % airfoil chord, in m
        SNPM.D.Uinf = 0.06;                     % free stream velocity, in m/s
        SNPM.D.nu = 1e-6;                       % kinematic viscosity of water, in m2/s
        SNPM.D.Dt = 1/30;                       % time separation between snapshots, in s

        %
        % Non-dimensional parameters
        %

        SNPM.nu = SNPM.D.nu/(SNPM.D.Uinf*SNPM.D.c); % Non-dimensional kinematic viscosity, equivalent to 1/Re
        SNPM.Dt = SNPM.D.Dt*SNPM.D.Uinf/SNPM.D.c;   % Non-dimensional time spacing

        %
        % Region of interest
        %

        m = size(X,1);             
        n = size(X,2);
        %
        ROI = [ 1 m 1 n ];
        % ROI = [ 24 124 41 152 ];
        %
        X = X(ROI(1):ROI(2),ROI(3):ROI(4));
        Y = Y(ROI(1):ROI(2),ROI(3):ROI(4));
        XD = XD(ROI(1):ROI(2),ROI(3):ROI(4));
        YD = YD(ROI(1):ROI(2),ROI(3):ROI(4));

        %
    elseif strcmp(flow,'PIV0018_2')
        %
        load('.\+data\PIV0018_grid2.mat');
        %
        % Dimensional parameters
        %

        SNPM.D.Rs = 1/6.9*1e-3;                 % resolution, in m/px
        SNPM.D.Vs = 10;                         % vector spacing, in px
        SNPM.D.c = 0.08;                        % airfoil chord, in m
        SNPM.D.Uinf = 11.3*SNPM.D.Rs/0.0357;    % free stream velocity, in m/s
        SNPM.D.nu = 1e-6;                       % kinematic viscosity of water, in m2/s
        SNPM.D.Dt = 0.0357;                     % time separation between snapshots, in s

        %
        % Non-dimensional parameters
        %

        SNPM.nu = SNPM.D.nu/(SNPM.D.Uinf*SNPM.D.c); % Non-dimensional kinematic viscosity, equivalent to 1/Re
        SNPM.Dt = SNPM.D.Dt*SNPM.D.Uinf/SNPM.D.c;   % Non-dimensional time spacing

        %
        % Region of interest
        %

        m = size(X,1);             
        n = size(X,2);
        %
        Nr = 3;                                 % Take first and last Nr rows and columns to avoid noisy data
        ROI = [ 1+Nr m-Nr 1+Nr n-Nr ];
        % ROI = [ 21 m 5 n ];
        %
        X = X(ROI(1):ROI(2),ROI(3):ROI(4));
        Y = Y(ROI(1):ROI(2),ROI(3):ROI(4));
        XD = XD(ROI(1):ROI(2),ROI(3):ROI(4));
        YD = YD(ROI(1):ROI(2),ROI(3):ROI(4));

    end
    %
    % Load grid parameters into SNPM
    %

        SNPM.X = X;                     SNPM.Y = Y;
        SNPM.m = size(X,1);             SNPM.n = size(X,2);
        SNPM.Dx = abs(X(1,1) - X(1,2)); SNPM.Dy = abs(Y(1,1) - Y(2,1));
        %
        SNPM.xmin = min(min(X));        SNPM.ymin = min(min(Y));
        SNPM.xmax = max(max(X));        SNPM.ymax = max(max(Y));

    % 

%
%% Body (Boundary Condition)
%

    if strcmp(flow,'SC') == 1
        %
        % Configuration of one cylinder
        %

        R = 0.5; xC = 0; yC = 0;

        %
        % Mask: 0 if body, 1 otherwise
        %

        SNPM.BD.M = ones(size(SNPM.X));
        SNPM.BD.M( ( SNPM.X - xC ).^2 + ( SNPM.Y - yC ).^2 <= R^2 ) = 0;

        % 
        % Mask function: 0 if body, NaN otherwise
        %

        SNPM.BD.Mf = @(X,Y) 0.*( (( X - xC ).^2 + ( Y - yC ).^2) <= R^2 ) + ...
                       NaN.*( (( X - xC ).^2 + ( Y - yC ).^2) > R^2 );

        %
        % Perimeter of body
        %

        SNPM.BD.x = [ linspace(-R+xC,R+xC,100) ]';
        SNPM.BD.y = [ yC + sqrt( R^2 - ( SNPM.BD.x - xC ).^2 ); yC - sqrt( R^2 - ( flip(SNPM.BD.x,1) - xC ).^2 ) ];
        SNPM.BD.x = [ SNPM.BD.x; flip(SNPM.BD.x,1) ];

        %
        % Additional left grid
        %

        SNPM.BD.X = [];
        SNPM.BD.Y = [];

        %
    elseif strcmp(flow,'FP') == 1
        %
        % Configuration of three cylinders (forward F, bottom B, top T) of D = 1
        %

        R = 0.5;
        %
        xF = -3/2*cos(30*pi/180);   yF = 0;
        xB = 0;              yB = -3/4;
        xT = 0;              yT =  3/4;

        %
        % Mask: 0 if body, 1 otherwise
        %

        SNPM.BD.M = ones(size(SNPM.X));
        SNPM.BD.M( ( SNPM.X - xF ).^2 + ( SNPM.Y - yF ).^2 <= R^2 ) = 0;
        SNPM.BD.M( ( SNPM.X - xB ).^2 + ( SNPM.Y - yB ).^2 <= R^2 ) = 0;
        SNPM.BD.M( ( SNPM.X - xT ).^2 + ( SNPM.Y - yT ).^2 <= R^2 ) = 0;

        % 
        % Mask function: 0 if body, NaN otherwise
        %

        SNPM.BD.Mf = @(X,Y) [ 0.*( (( X - xF ).^2 + ( Y - yF ).^2) <= R^2 ) + ...
                              NaN.*( (( X - xF ).^2 + ( Y - yF ).^2) > R^2 );
                              0.*( (( X - xB ).^2 + ( Y - yB ).^2) <= R^2 ) + ...
                              NaN.*( (( X - xB ).^2 + ( Y - yB ).^2) > R^2 );
                              0.*( (( X - xT ).^2 + ( Y - yT ).^2) <= R^2 ) + ...
                              NaN.*( (( X - xT ).^2 + ( Y - yT ).^2) > R^2 ); ];

        %
        % Perimeter of body
        %

        SNPM.BD.x = [ linspace(-R+xF,R+xF,100)', linspace(-R+xB,R+xB,100)', linspace(-R+xT,R+xT,100)' ];
        SNPM.BD.y = [ yF + sqrt( R^2 - ( SNPM.BD.x(:,1) - xF ).^2 ),         yB + sqrt( R^2 - ( SNPM.BD.x(:,2) - xB ).^2 ),         yT + sqrt( R^2 - ( SNPM.BD.x(:,3) - xT ).^2 );
                      yF - sqrt( R^2 - ( flip(SNPM.BD.x(:,1),1) - xF ).^2 ), yB - sqrt( R^2 - ( flip(SNPM.BD.x(:,2),1) - xB ).^2 ), yT - sqrt( R^2 - ( flip(SNPM.BD.x(:,3),1) - xT ).^2 ) ];
        SNPM.BD.x = [ SNPM.BD.x; flip(SNPM.BD.x,1) ];

        %
        % Additional left grid
        %

        SNPM.BD.X = [];
        SNPM.BD.Y = [];

        %
    elseif (strcmp(flow,'PIV0018') == 1) || (strcmp(flow,'PIV0018_2') == 1)
        %
        load('D:\MPSIFDpod\+data\NACA0018.dat');
        %
        % Map NACA airfoil
        %

        x = NACA0018(:,1)'-1; y = NACA0018(:,2)'; z = NACA0018(:,1)'.*0;
        N =   [cosd(10)   sind(10) 0;
               -sind(10)  cosd(10) 0;
                     0         0  1]*[x;y;z];
        x = N(1,:); xT = x(1,1:ceil(length(x)/2))'; xB = x(1,ceil(length(x)/2):end)';
        y = N(2,:); yT = y(1,1:ceil(length(x)/2))'; yB = y(1,ceil(length(x)/2):end)';
        %
        yTf = @(X) interp1(xT,yT,X);
        yBf = @(X) interp1(xB,yB,X);
        xL = min(x); xR = max(x);

        %
        % Perimeter of bodies
        %

        SNPM.BD.x = [ linspace(xL,xR,200)' ];
        SNPM.BD.y = [ yTf(SNPM.BD.x); yTf(SNPM.BD.x(end,1)); yBf(flip(SNPM.BD.x(2:end-1,1),1)); yTf(SNPM.BD.x(1,1)) ];
        SNPM.BD.x = [ SNPM.BD.x; flip(SNPM.BD.x,1) ];

        %
        % Mask of velocity fields
        %
        
        SNPM.BD.M = ones(size(SNPM.X));

        %
        % Empty left grid with NACA airfoil
        %

        xmin = -1;        xmax = SNPM.xmin;
        ymin = SNPM.ymin; ymax = SNPM.ymax;
        %
        x = linspace(xmin,xmax,SNPM.n);
        y = linspace(ymin,ymax,SNPM.m);
        %
        [SNPM.BD.X,SNPM.BD.Y] = meshgrid(x,y);
        SNPM.BD.Y = flip(SNPM.BD.Y,1);
        
        %
    end
%
end
