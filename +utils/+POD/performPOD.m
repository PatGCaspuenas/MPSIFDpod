function SVD = performPOD(SNP)
%%                            performPOD.m
%--------------------------------------------------------------------------
%
% Performs Singular Value Decomposition onto flow data
%
% INPUTS
%
%   SNP      : structure containing flow velocity fields and time vector
%
% OUTPUTS
%
%   SVD      : structure containing spatial (U), temporal (V and a) modes,
%              along with their singular values (S) of the flow
%
%--------------------------------------------------------------------------
%
if isfield(SNP,'udt') && isfield(SNP,'vdt')
    %
    [SVD.U,SVD.S,SVD.V] = svd([SNP.udt; SNP.vdt],'econ');
    %
elseif isfield(SNP,'vdt')
    %
    [SVD.U,SVD.S,SVD.V] = svd(SNP.vdt,'econ');
    %
else
    %
    [SVD.U,SVD.S,SVD.V] = svd(SNP.udt,'econ');
    %
end
%
SVD.t = SNP.t;
SVD.a = (SVD.S*SVD.V')';
%
