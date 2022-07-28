function r = elbowFit(x,y)
%%                           elbowFit.m
%--------------------------------------------------------------------------
%
% Finds the elbow of the (x,y) curve in the x-position (the index value). 
% Follows procedure of Brindise&Vlachos (2017) in "Proper orthogonal 
% decomposition truncation method for data denoising and order reduction"
%
% INPUTS
%
%   x : array of x values
%   y : array of y values
%
% OUTPUT
%
%   r : elbow index position in x
%
%--------------------------------------------------------------------------
%
Ny = length(y);
%

    for i = 1:Ny
        %
        c1 = polyfit(x(1:i),y(1:i),1);
        c2 = polyfit(x(i:Ny),y(i:Ny),1);
        y1 = polyval(c1,x(1:i));                                     % Fit of first region of the curve                     
        y2 = polyval(c2,x(i+1:Ny));                                  % Fit of second region of the curve 
        yt = [ y1 y2 ];                                              % Fit of both regions of the curve
        %
        if size(yt,1) ~= size(y,1)
            yt = yt';
        end
        % 
        R2(i) = 1 - sum( ( y - yt ).^2)./sum( ( y - mean(y) ).^2);   % R-square determination coefficient of each fit
        err(i) = sqrt(sum( ( y - yt ).^2 )/Ny);                      % error of each fit
        %
    end
    % 
    [~,r] = max(R2./err);                                            % index position that minimizes the error among all fits
    %
    
end