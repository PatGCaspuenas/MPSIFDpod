function R2 = getR2factor(R2Method,Xtrue,X)
%%                          getR2factor.m
%--------------------------------------------------------------------------
%
% Gets the R square determination/correlation for each columns of two given
% matrices
%
% INPUTS
%
%   R2Method : char flag indicating the type of R2 coeff to compute   
%              If == 'd' -> coefficient of determination
%              If == 'c' -> correlation coefficient
%
%   Xtrue   : true reference matrix of dimension Nt x Nr
%   X       : reconstructed matrix of dimension Nt x Nr
%
% OUTPUT
%
%   R2       : R-squared coefficient as row vector
%
%--------------------------------------------------------------------------
%
% Parameters
%

    Nr = size(Xtrue,2);                                       % number of modes                       

%
% Get R2 coefficient for each mode
%

    for i = 1:Nr
        %
        if R2Method == 'd'                                    % determination coefficient 
            %
            R2(i) = 1 - sum( (Xtrue(:,i) - X(:,i)).^2 )/...
                sum( (Xtrue(:,i) - mean(Xtrue(:,i)) ).^2 );
            %
        elseif R2Method == 'c'                                % correlation coefficient
            %
            R2(i) = mean(Xtrue(:,i).*X(:,i))^2/...
                mean(Xtrue(:,i).^2)/mean(X(:,i).^2);
            %
            if isnan(R2(i))                                   % avoid singularities due to the denominator
                R2(i) = 1e-4;
            end
        end
        %
    end

%
end