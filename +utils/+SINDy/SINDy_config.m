function SINDy = SINDy_config(SVD,param)
%%                           SINDy_config.m
%--------------------------------------------------------------------------
%
% Loads main configuration parameters, such as the tolerances and tuning
% parameter grid
%
% INPUTS
%
%   SVD       : structure containing, among others, temporal vector of states
%   param     : structure of required parameters for the given function
%
% OUTPUT
%
%   SINDy     : structure with main configuration parameters
%
%--------------------------------------------------------------------------
%
% Parameters
%

    OptimizationMethod = param.OptimizationMethodSINDy{1}; % Optimization method: 'ALASSO'/'STR'
    FlowConfig = param.OptimizationMethodSINDy{2};         % Flow chosen: 'SC'/'FP'/'PIV0018'
    %
    if isfield(param,'alpha')                              % alpha tuning parameter
        alpha = param.alpha;
    else
        alpha = 0;
    end
    %
    if isfield(param,'delta')                              % delta tuning parameter
        delta = param.delta;
    else
        delta = 0;
    end
    %
    if isfield(param,'Ng')                                 % number of swept tuning parameter values
        Ng = param.Ng;
    else
        Ng = 0;
    end
    %
    if isfield(param,'lambda')                             % lambda tuning parameter
        lambda = param.lambda;
    else
        lambda = 0;
    end
    %
    SINDy.t = SVD.t;                                       % temporal vector
    SINDy.tol = 1e-3;                                      % sparsity tolerance
    SINDy.ftol = 1e-5;                                     % convergence tolerance

%
% Choose tuning parameters
%

    if strcmp(OptimizationMethod,'ALASSO')
        %
        if strcmp(FlowConfig,'PIV0018')
            %
            N = Ng.*(isfield(param,'Ng')) + 10.*(~isfield(param,'Ng'));
            SINDy.alpha = alpha.*(isfield(param,'alpha')) + logspace(-3.5,-1.3010,N).*(~isfield(param,'alpha')); % logspace(-3.5,-1,N);
            SINDy.delta = delta.*(isfield(param,'delta')) + logspace(-0.6990,-0.3979,N).*(~isfield(param,'delta'));  % logspace(-0.8239,-0.4559,N); 
            %
        elseif strcmp(FlowConfig,'FP')
            %
            N = Ng.*(isfield(param,'Ng')) + 15.*(~isfield(param,'Ng'));
            SINDy.alpha = alpha.*(isfield(param,'alpha')) + logspace(-3,-1.3,N).*(~isfield(param,'alpha'));
            SINDy.delta = delta.*(isfield(param,'delta')) + logspace(-1.3010,-0.6021,N).*(~isfield(param,'delta')); % logspace(-1.3010,-0.6021,N);
            %
        end
        %
        delta_diff = flip(cumsum(diff(SINDy.delta)));
        %
        SINDy.delta(1:N-1) = SINDy.delta(N) - delta_diff;
        %
    elseif strcmp(OptimizationMethod,'STR')
        %
        N = Ng.*(isfield(param,'Ng')) + 500.*(~isfield(param,'Ng'));
        SINDy.lambda = lambda.*(isfield(param,'lambda')) + logspace(-4,2,N).*(~isfield(param,'lambda')); 
        %
    end
    
%