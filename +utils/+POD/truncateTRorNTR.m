function [SNPtr,SVDt ] = truncateTRorNTR(SNP,SVDtr,SNPM,SINDytr,param)
%%                         truncateTRorNTR.m
%--------------------------------------------------------------------------
%
% From BFI states with temporal resolution, retrieves the original basis
% back, that is, the velocity fields. Integration of pressure is given as
% an option
%
% INPUTS
%
%   SNP      : structure containing, among others, velocity and
%              acceleration fields
%   SVDtr    : structure containing spatial truncated modes
%   SNPM     : structure containing main relevant geometric parameters of
%              dataset
%   SINDytr  : structure containing matrix of ceofficients of the dynamical
%              system
%   param    : structure with parameters to be used for the given function
%
% OUTPUT
%
%   SNPtr    : structure containing truncated velocity and acceleration
%              fields
%   
%
% UTILS
%
%   utils.SINDy.getDerSINDy
%   utils.differentiation.numDifferentiation
%   utils.pressure.pIntegratorSNPM
%   utils.pressure.getPGradient
%
%-------------------------------------------------------------------------- 
%
% Parameters
%

    % param.PressureboundaryCondition from utils.pressure.pIntegratorSNPM
    % param.SaveLaplacianFields from utils.pressure.getPGradient
    % param.SaveGradientVelFields from utils.pressure.getPGradient
    var_flags = {'IntegratePressure','TruncateFlowAcceleration','TruncateFlowAccelerationNumDiff',...
        'TruncateFlowAccelerationSINDy'};
    %
    for k = 1:length(var_flags)
        %
        if isfield(param,var_flags{k})
           flag.(var_flags{k}) = param.(var_flags{k});
        else
           flag.(var_flags{k}) = 0;
        end
        %
    end
    %
%
% Truncated flow field
%

    SNPtr.t = SNP.t; 
    %
    if isfield(SNP,'vdt') && isfield(SNP,'udt')
        %
        m = SNPM.m; n = SNPM.n;
        SVDt.Vr = (inv(SVDtr.Sr)*SVDtr.Ur'*[SNP.udt; SNP.vdt])';
        SVDt.ar = (SVDtr.Ur'*[SNP.udt; SNP.vdt])';
        uvdt = SVDtr.Ur*SVDt.ar';
        %
        SNPtr.udt = uvdt(1:m*n,:);
        SNPtr.u = SNPtr.udt + SNP.u_m;
        SNPtr.vdt = uvdt((m*n+1):end,:);
        SNPtr.v = SNPtr.vdt + SNP.v_m;
        %
    elseif isfield(SNP,'vdt')
        %
        SVDt.Vr = (inv(SVDtr.Sr)*SVDtr.Ur'*SNP.vdt)';
        SVDt.ar = (SVDtr.Ur'*SNP.vdt)';
        SNPtr.vdt = SVDtr.Ur*SVDt.ar';
        SNPtr.v = SNPtr.vdt + SNP.v_m;
        %
    else
        %
        SVDt.Vr = (inv(SVDtr.Sr)*SVDtr.Ur'*SNP.udt)';
        SVDt.ar = (SVDtr.Ur'*SNP.udt)';
        SNPtr.udt = SVDtr.Ur*SVDt.ar';
        SNPtr.u = SNPtr.udt + SNP.u_m;
        %
    end

%
% Flow derivatives
%

    if flag.TruncateFlowAcceleration                       % choose to truncate flow acceleration (1) or not (1)
        %
        if flag.TruncateFlowAccelerationNumDiff            % truncated acceleration from differentiation of truncated velocity
            %
            if isfield(SNP,'udt') || isfield(SNP,'uvdt')
                %
                SNPtr.du = utils.differentiation.numDifferentiation(SNPtr.u,SNPtr.t);
                %
            end
            if isfield(SNP,'vdt') || isfield(SNP,'uvdt')
                %
                SNPtr.dv = utils.differentiation.numDifferentiation(SNPtr.v,SNPtr.t);
                %
            end
            %
        else
            %
            if isfield(SNP,'vdt') && isfield(SNP,'udt')
                %
                if flag.TruncateFlowAccelerationSINDy       % truncated acceleration from SINDy dynamical system
                    %
                    SVDt.dar = utils.SINDy.getDerSINDy(SINDytr,SVDt.ar);
                    %
                else                                        % truncated acceleration from projection of available non-truncated acceleration fields
                    %
                    SVDt.dVr = (inv(SVDtr.Sr)*SVDtr.Ur'*[SNP.du; SNP.dv])';
                    SVDt.dar = (SVDtr.Ur'*[SNP.du; SNP.dv])';
                    %
                end
                %
                duv = SVDtr.Ur*SVDt.dar';
                SNPtr.du = duv(1:m*n,:);
                SNPtr.dv = duv((m*n+1):end,:);
                %
            elseif isfield(SNP,'vdt')
                %
                if flag.TruncateFlowAccelerationSINDy
                    %
                    SVDt.dar = utils.SINDy.getDerSINDy(SINDytr,SVDt.ar);
                    %
                else
                    %
                    SVDt.dVr = (inv(SVDtr.Sr)*SVDtr.Ur'*SNP.dv)';
                    SVDt.dar = (SVDtr.Ur'*SNP.dv)';
                    %
                end
                %
                SNPtr.dv = SVDtr.Ur*SVDt.dar';
                %
            else
                %
                if flag.TruncateFlowAccelerationSINDy
                    %
                    SVDt.dar = utils.SINDy.getDerSINDy(SINDytr,SVDt.ar);
                    %
                else
                    %
                    SVDt.dVr = (inv(SVDtr.Sr)*SVDtr.Ur'*SNP.du)';
                    SVDt.dar = (SVDtr.Ur'*SNP.du)';
                end
                %
                SNPtr.du = SVDtr.Ur*SVDt.dar';
                %
            end
            %
        end
    end

%
% Reconstructed pressure field (if selected)
%

    if flag.IntegratePressure
        %
        % Truncate dudx,dudy,dvdx,dvdy directly, if possible
        %
            
            var_grad = {'dudx','dudy','dvdx','dvdy'};
            %
            if isfield(SNP,'udt') && isfield(SNP,'vdt')
                %
                if isfield(SNP,'dudx') && isfield(SNP,'dvdx')
                    %
                    duvdx = SVDtr.Ur*(SVDtr.Ur'*[ SNP.dudx; SNP.dvdx ]);
                    SNPtr.dudx = duvdx(1:m*n,:);
                    SNPtr.dvdx = duvdx((1+m*n):end,:);
                    %
                end
                %
                if isfield(SNP,'dudy') && isfield(SNP,'dvdy')
                    %
                    duvdy = SVDtr.Ur*(SVDtr.Ur'*[ SNP.dudy; SNP.dvdy ]);
                    SNPtr.dudy = duvdy(1:m*n,:);
                    SNPtr.dvdy = duvdy((1+m*n):end,:);
                    %
                end
                %
            else
                %
                for g = 1:length(var_grad)
                    %
                    aux = SVDtr.Ur*(SVDtr.Ur'*[ SNP.(var_grad{g}) ]);
                    SNPtr.(var_grad{g}) = aux;
                    %
                end
                %
            end

        %
        SNPtr = utils.pressure.getPGradient(SNPtr,SNPM,param);
        SNPtr.p = utils.pressure.pIntegratorSNPM(SNPtr,SNPM,param);
        %
    end


%
end
