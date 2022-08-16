function NTR = createNTRSet(TR,SNPM,param)
%%                            createNTRset.m
%--------------------------------------------------------------------------
%
% Generates non-time-resolved data from time-resolved dataset
%
% INPUTS
%
%   TR       : structure containing, if available, velocity, pressure,
%              vorticity, time and gradient time-resolved fields
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   param    : structure of flags and parameters used in the function
%
% OUTPUTS
%
%   NTR      : structure containing, if available, velocity, pressure,
%              vorticity, acceleration, time and gradient non-time-resolved fields
%
%--------------------------------------------------------------------------
%
% Parameters
%

    flagTimeSpacing = param.FlagNTRTimeSpacing; % If 'irr' : irregular time spacing
                                                %     else : regular time spacing
    ts = param.ts;                              % Average time separation between snapshots
    %
    if isfield(param,'u_m')                     % Mean horizontal velocity field
        u_m = param.u_m;
    else
        u_m = 0;
    end
    %
    if isfield(param,'v_m')                     % Mean vertical velocity field
        v_m = param.v_m;
    else
        v_m = 0;
    end
    %
    Dt = SNPM.Dt;                               % Time-resolved separation
    NtTR = length(TR.t);                        % Number of snapshots in TR dataset
    NtNTR = floor((NtTR-1)/floor(ts/Dt)) + 1;   % Number of snapshots in NTR dataset
    it = randperm(NtTR,NtNTR);
    it = sort(it);                              % Snapshot indexes for irregular time spacing (chosen randomnly)

%
% Create NTR Set
%

    var_str = {'u','v','du','dv','t','p','w','dudx','dudy','dvdx','dvdy'}; % possible data that can be contained in files
    %
    for i = 1:length(var_str)
        %
        if isfield(TR,var_str{i})
            %
            aux = TR.(var_str{i});
            %
            if strcmp(flagTimeSpacing,'irr')    % Irregular spacing
                %
                NTR.SNP.(var_str{i}) = aux(:,it);
                %
            else                                % Regular spacing
                % 
                NTR.SNP.(var_str{i}) = aux(:,1:(floor(ts/Dt)+1):end);
                %
            end
        end
        %
    end

%
% Include velocity fluctuations
%

    if isfield(TR,'u') && isfield(TR,'v')
        %
        NTR.SNP.u_m = u_m.*(isfield(param,'u_m')) + mean(NTR.SNP.u,2).*(~isfield(param,'u_m'));
        NTR.SNP.udt = NTR.SNP.u - NTR.SNP.u_m;
        %
        NTR.SNP.v_m = v_m.*(isfield(param,'v_m')) + mean(NTR.SNP.v,2).*(~isfield(param,'v_m'));
        NTR.SNP.vdt = NTR.SNP.v - NTR.SNP.v_m;
        %
    elseif isfield(TR,'v')
        %
        NTR.SNP.v_m = v_m.*(isfield(param,'v_m')) + mean(NTR.SNP.v,2).*(~isfield(param,'v_m'));
        NTR.SNP.vdt = NTR.SNP.v - NTR.SNP.v_m;
        %
    else
        %
        NTR.SNP.u_m = u_m.*(isfield(param,'u_m')) + mean(NTR.SNP.u,2).*(~isfield(param,'u_m'));
        NTR.SNP.udt = NTR.SNP.u - NTR.SNP.u_m;
        %
    end

%
% Undersample time array as the NTR Set
%

    NTR.SNP.Dt = ts;
    %
    if strcmp(flagTimeSpacing,'irr')
        NTR.SNP.t = TR.t(it);
    else
        NTR.SNP.t = TR.t(1:(floor(ts/Dt)+1):end);
    end