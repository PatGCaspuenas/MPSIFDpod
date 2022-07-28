function [Xi,Alpha,Delta,R2,Nact] = paretoElbowl1(Theta,dX,alpha,delta,ftol,tol)
%%                          paretoElbowl1.m
%--------------------------------------------------------------------------
%
% For a given state of the system, retrieves the elbow point of the
% curve(Nact,R2) and their corresponding tuning parameters
%
% INPUTS
%
%   Theta  : library of functions expressed in terms of the states of the system
%   dX     : derivative of one state of the system 
%   alpha  : l1-norm constraint parameter
%   delta  : Adaptive LASSO updating parameter
%   ftol   : inner-loop tolerance
%   tol    : sparsity tolerance
%
% OUTPUTS
%
%   Xi     : converged coefficient vector for given state
%   Alpha  : converged alpha value
%   Delta  : converged delta value
%   R2     : corresponding R-square determination coefficient 
%   Nact   : corresponding number of active terms 
%
% UTILS
%
%   utils.SINDy.paretoFrontl1
%   utils.fits.elbowFit
%   utils.SINDy.ALASSO
%
%--------------------------------------------------------------------------
%
% Grid parameters and results
%

    [ALPHA,DELTA,RSQUARED,NACT] = utils.SINDy.paretoFrontl1(...
                                  Theta,dX,alpha,delta,ftol,tol);

%
% Define straight line between minimum and maximum points of scatter curve
%

    [imax,jmax] = find(NACT == max(NACT,[],'all'),1,'first');      % maximum Nact point
    [imin,jmin] = find(NACT == min(NACT,[],'all'),1,'first');      % minimum Nact point
    %
    m = (RSQUARED(imax,jmax) - RSQUARED(imin,jmin))/...            % line slope
        (max(NACT,[],'all') - min(NACT,[],'all'));
    n = RSQUARED(imax,jmax) - m*max(NACT,[],'all');                % line y-interceipt

%
% Compute Euclidean distance grid
%

    for i = 1:size(Alpha,1)
        for j = 1:size(Delta,2)
            D(i,j) = computeEuclideanD(m,n,NACT(i,j),RSQUARED(i,j));
        end
    end

%    
% Reshape grid into vectors and sort them ordered in NACT
%

    NACT = reshape(NACT,[1 size(NACT,1)*size(NACT,2)]);
    RSQUARED = reshape(RSQUARED,[1 size(RSQUARED,1)*size(RSQUARED,2)]);
    ALPHA = reshape(ALPHA,[1 size(ALPHA,1)*size(ALPHA,2)]);
    DELTA = reshape(DELTA,[1 size(DELTA,1)*size(DELTA,2)]);
    %
    [NACT, iS] = sort(NACT);
    RSQUARED = RSQUARED(iS); 
    ALPHA = ALPHA(iS); 
    DELTA = DELTA(iS);

%
% Create a continuous curve of increasing R2
%

    N = NACT(1);                                                % start from first point in curve
    R = RSQUARED(1); 
    A = ALPHA(1); 
    d = DELTA(1);
    %
    c = 1;                                                      % counting variable
    %
    for j = 2:length(NACT)
        %
        if RSQUARED(j) > R(c)                                   % If R2 increases for an 
                                                                % increase of Nact, include 
                                                                % the point in the curve
            c = c + 1;
            N(c) = NACT(j); 
            R(c) = RSQUARED(j); 
            A(c) = ALPHA(j); 
            d(c) = DELTA(j);
            %
        end
    end

%
% Look for elbow in the curve
%

    ir = elbowFit(N,R);

%
% If there are points with lower Nact but still a high R2, make that our
% elbow
%

    stop = 0;
    %
    while (~stop) && (ir~=1)
        if R(ir-1) > 0.95                                       % selected R2 threshold for an accurate fit
            ir = ir-1;
        else
            stop = 1;
        end
    end
    %
    R2 = R(ir);
    Nact = N(ir);
    Alpha = A(ir);
    Delta = d(ir);

%
% With the converged (alpha,delta), recalculate the coefficient vector Xi
%

    Xi = zeros(size(Theta,2),1);
    Xi = ALASSO(Theta, dX,Alpha, Delta, ftol,tol);

%
end