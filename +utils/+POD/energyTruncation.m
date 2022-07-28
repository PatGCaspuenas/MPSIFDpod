function r = energyTruncation(S,E)
%%                        energyTruncation.m
%--------------------------------------------------------------------------
%
% Finds the threshold vector index for reaching a certain level of energy
%
% INPUTS
%
%   S : Singular value matrix
%   E : energy level
%
% OUTPUT
%
%   r : index truncation position 
%
%--------------------------------------------------------------------------
%

energy = cumsum(diag(S.^2))/sum(diag(S.^2)); % cumulative energy up to each mode
%
r = find(energy <= E,1,'last'); 

%
end