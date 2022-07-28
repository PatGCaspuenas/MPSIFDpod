function [M,n] = normalizeMatrixColumns(M)
%%                      normalizeMatrixColumn.m
%--------------------------------------------------------------------------
%
% Non-dimensionalizes the columns of a given matrix
%
% INPUTS
%
%   M : matrix
%
% OUTPUT
%
%   M : matrix with normalized columns
%   n : row vector with norm of each columns
%
%--------------------------------------------------------------------------
%

for j = 1:size(M,2)
    n(1,j) = norm(M(:,j),2);
    M(:,j) = M(:,j)/n(1,j);
end

%
end