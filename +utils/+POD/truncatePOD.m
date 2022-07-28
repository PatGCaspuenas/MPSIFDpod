function SVD = truncatePOD(SNP,SVD,param)
%%                           truncatePOD.m
%--------------------------------------------------------------------------
%
% Truncates POD basis depending on the selected criteria
%
% INPUTS
%
%   SNP      : structure containing flow velocity, acceleration fields and time vector
%   SVD      : structure containing spatial (U), temporal (V and a) modes,
%              along with their singular values (S) of the flowc
%   param    : structure of flags and parameters used in the function
%
% OUTPUT
%
%   SVD      : SVD structure adding truncated basis, projected acceleration
%              fields and (if needed) mean mode
%
% UTILS
%
%   utils.fits.elbowFit
%   utils.POD.energyTruncation
%   utils.POD.includeMeanSVD
%
%--------------------------------------------------------------------------
%
% Parameters
%

    methodPOD = param.TruncationMethodPOD{1};                    % Method of truncation chosen. Can be
                                                                 % 'manual'/'energy'/'elbow'
    inputPOD = param.TruncationMethodPOD{2};                     % Input depending on truncation method. Can be 
                                                                 % number of modes/empty/level of energy
    if isfield(param,'IncludeMeanPOD')                           % determines if mean is included in POD basis
        % 
        FlagMean = param.IncludeMeanPOD;                        
        %
    else
        %
        FlagMean = 0;
        %
    end                                                             
    Nt = size(SVD.V,1);                                          % time samples

%
% Use truncation method
%

    if strcmp(methodPOD,'elbow')                                 % chooses the elbow point of the cumulative energy curve
        %
        SVD.r = utils.fits.elbowFit(1:length(SVD.S),cumsum(diag(SVD.S.^2))/sum(diag(SVD.S.^2)));
        %
    elseif strcmp(methodPOD,'energy')                            % sets threshold at certain level of energy
        %
        SVD.r = utils.POD.energyTruncation(SVD.S,inputPOD);
        %
    elseif strcmp(methodPOD,'manual')                            % sets manual threshold
        %
        SVD.r = inputPOD;
        %
    end

%
% Truncate modes
%

    SVD.ar = SVD.a(:,1:SVD.r);
    SVD.Vr = SVD.V(:,1:SVD.r);
    SVD.Sr = SVD.S(1:SVD.r,1:SVD.r);
    SVD.Ur = SVD.U(:,1:SVD.r);

%
% Project acceleration fields onto SVD basis
%

    if  isfield(SNP,'udt') && isfield(SNP,'vdt')
        %
        aux = [SNP.du;SNP.dv];
        %
    elseif isfield(SNP,'vdt')
        %
        aux = SNP.dv;
        %
    else
        %
        aux = SNP.du;
        %
    end
    %
    SVD.dVr = []; SVD.dar = [];
    %
    for i = 1:Nt
        %
        SVD.dVr(i,:) = (pinv(SVD.Sr)*SVD.Ur'*aux(:,i))';
        SVD.dar(i,:) = (SVD.Ur'*aux(:,i))';
        %
    end

%
% Include mean in SVD
%

    if FlagMean
        SVD = utils.POD.includeMeanSVD(SNP,SVD);
    end
    
%