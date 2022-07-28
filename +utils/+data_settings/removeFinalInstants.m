function TR = removeFinalInstants(TR,NTR)
%%                       removeFinalInstants.m
%--------------------------------------------------------------------------
%
% If NTR set reaches time tj, remove time samples t > tj from TR set
%
% INPUT
%
%   TR   : structure containing time-resolved fields
%   NTR  : structure containing non-time-resolved fields (can also work for
%          time-resolved) with a lower t(end) than TR
%
% OUTPUT
%
%   TR   : structure with final time instant removed
%
%--------------------------------------------------------------------------
%
% Ending time instant
%

    endt = find( TR.t == NTR.t(end) ,1,'last');

%
% Possible data that can be contained in files
%

    var_str = {'u','v','du','dv','t','p','w','dudx','dudy','dvdx','dvdy','udt','vdt'}; 

%
for i = 1:length(var_str) % load existing data
    %
    if isfield(TR,var_str{i})
        %
        TR.(var_str{i}) = TR.(var_str{i})(:,1:endt);
        %
    end
    %
end
%
end