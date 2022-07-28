function RMSE = computeRMSE(Xtrue,X,SNPM,param)
%%                          computeRMSE.m
%--------------------------------------------------------------------------
%
% Computes the Root Mean Square Error of two given fields in space,time 
% or both, depending on the method chosen
%
% INPUTS
%
%   Xtrue   : true reference fields of dimension Np x Nt
%   X       : reconstructed fields of dimension Np x Nt
%   SNPM    : structure containing main relevant geometric parameters of
%              dataset
%   param   : structure of required parameters for the given function
%
% OUTPUT
%
%   RMSE    : error obtained as a punctual value, in time or in space
%
%--------------------------------------------------------------------------
%
% Parameters
%

    errorType = param.ErrorTypeRMSE;                     % Error type. Can either be
                                                         % 'w' : single RMSE value
                                                         % 't' : RMSE temporal history
                                                         % 's' : RMSE spatial distribution
                                                         % else: RMSE snapshot matrix                           
    %
    if isfield(param,'ComputeRMSEatMiddle')              % flag to evaluate error only in between available NTR snapshots
        if param.ComputeRMSEatMiddle
            Ts = param.ts;
            Xtrue = Xtrue(:,(1+Ts/2):Ts:end);
            X = X(:,(1+Ts/2):Ts:end);
        end
    end
    %
    mask = SNPM.BD.M(:);                                 % mask to avoid body BCs
    if size(Xtrue,1) > size(mask,1)
        %
        mask = [ mask; mask ];
        %
    end
    Xtrue = Xtrue.*mask;
    X = X.*mask;
    %
    Np = size(Xtrue,1);                                  % Grid points x data dimension
    Nt = size(Xtrue,2);                                  % time samples
    std = sqrt(  sum( ( Xtrue - ...                      % standard deviation of flow
          mean(Xtrue,2) ).^2 , 'all' )/...
          (Np*Nt)  );

%
% Compute error
%

    if errorType == 'w'
        %
        RMSE = sqrt( sum( ( Xtrue - X ).^2,'all' ) )/...
            sqrt(Np*Nt)/std;
        %
    elseif errorType == 't'
        %
        for i = 1:Nt
            %
            RMSE(i) = sqrt(sum( ( Xtrue(:,i) - X(:,i) ).^2,'all') )/...
             sqrt(Np)/std;
            %
        end
        %
    elseif errorType == 's'
        %
        for i = 1:Np
            %
            RMSE(i,1) = sqrt(sum( ( Xtrue(i,:) - X(i,:) ).^2,'all'))/...
                sqrt(Nt)/std;
            %
        end
        %
    else
        %   
        for i = 1:Np
            for j = 1:Nt
                %
                RMSE(i,j) = sqrt( ( Xtrue(i,j) - X(i,j) ).^2 )/...
                    std;
                %
            end
        end
        %
    end

%
end