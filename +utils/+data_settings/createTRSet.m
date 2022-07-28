function TR = createTRSet(filename)
%%                            createTRset.m
%--------------------------------------------------------------------------
%
% Loads time-resolved data and arranges it into a struct
%
% INPUT
%
%   filename : string with name of file to be loaded
%
% OUTPUT
%
%   TR       : structure containing, if available, velocity, pressure,
%              vorticity, acceleration, time and gradient time-resolved fields
%
%--------------------------------------------------------------------------
%
load(filename);
%
var_str = {'u','v','du','dv','t','p','w','dudx','dudy','dvdx','dvdy'}; % possible data that can be contained in files
%
for i = 1:length(var_str) % load existing data
    %
    if exist(var_str{i},'var')
        %
        TR.(var_str{i}) = eval(var_str{i});
        %
    end
    %
end
%
end