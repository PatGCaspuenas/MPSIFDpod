function SVD = includeMeanSVD(SNP,SVD)
%%                         includeMeanSVD.m
%--------------------------------------------------------------------------
%
% Includes mean mode into POD basis
%
% INPUTS
%
%   SNP      : structure containing mean velocity fields
%   SVD      : structure containing spatial (U), temporal (V and a) modes,
%              along with their singular values (S) of the flow. Also
%              contains truncated basis and projected acceleration fields
%
% OUTPUT
%
%   SVD      : structure adding truncated sets with mean mode
%
%--------------------------------------------------------------------------
%
% Parameters
%

    Nt = size(SVD.V,1);                             % time samples

%  
% Mean spatial mode
%

    if isfield(SNP,'udt') && isfield(SNP,'vdt')
        %
        U0 = [ SNP.u_m; SNP.v_m];
        %
    elseif isfield(SNP,'vdt')
        %
        U0 = [ SNP.v_m ];
        %
    else
        %
        U0 = [ SNP.u_m ];
        %
    end

%
% Temporal mean mode
%

    a0 = ones(Nt,1); 
    sigma0 = norm(a0); 
    V0 = a0/norm(a0);

%
% Include mean mode into truncated spatial set
%

    SVD.Ur0 = [ U0 SVD.Ur ];

%
% Include mean mode into temporal sets
%

    SVD.Vr0 = [ V0 SVD.Vr ];
    SVD.ar0 = [ a0 SVD.ar ];
    SVD.dar0 = [ zeros(Nt,1) SVD.dar ];
    SVD.dVr0 = [ zeros(Nt,1) SVD.dVr ];

%
% Include mean mode into sigma matrix
%

    SVD.Sr0 = zeros(size(SVD.Sr,1)+1,size(SVD.Sr,2)+1);
    SVD.Sr0(1,1) = sigma0;
    SVD.Sr0(2:end,2:end) = SVD.Sr;
    
%
SVD.r0 = SVD.r + 1;
%