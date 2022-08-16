function Xii = ALASSO(Theta, dX, alpha, delta, ftol,tol)
%%                              ALASSO.m
%--------------------------------------------------------------------------
%
% This function applies the Adaptive LASSO regression fit onto a set of data, 
% given a set of tuning parameters:
%
%       Theta*Chi = dy minimizing its l2-norm and alpha*sum(w*abs(Chi))
%
% Implemented method from Fukami et al (2020) : "Sparse identification of
% nonlinear dynamics with low-dimensionalized flow representations
%
% Only applies to a single state of the system
%
% INPUTS
%
%   Theta  : library of functions expressed in terms of the states of the system
%   dX     : derivative of one state of the system 
%   alpha  : l1-norm constraint parameter
%   delta  : Adaptive LASSO updating parameter
%   ftol   : inner-loop tolerance
%   tol    : sparsity tolerance
%
% OUTPUT
%
%   Xii    : sparse vector of coefficients. Solution of the regression
%            problem Theta*Chi = dy
%
% UTILS
%
%   utils.modelling.normalizeMatrixColumns
%--------------------------------------------------------------------------
%
% Normalization of columns of the library matrix and state derivative
% matrix. Necessary for using the same sparsification criteria in the
% optimization
%

    [Theta,normT] = utils.modelling.normalizeMatrixColumns(Theta); 
    [dX,normdX]   = utils.modelling.normalizeMatrixColumns(dX);

%
% Parameters
%

    W = ones(size(Theta,2),1);                  % weights of columns of library matrix
    err = 2*ftol;                               % error between i and i+1 Chi solutions
    c = 0;                                      % counting index
    Nmax = 100;                                 % maximum number of iterations for convergence

%
% Initial guess: non-regularized  regression to decrease computational time
% of cvx optimization
%

    Xi0 = Theta\dX;                             % solution from linear regression
    %
    zi = find(abs(Xi0(:,1)) <= max(abs(Xi0(:,1)))*tol); % indexes of negligible terms
    nzi = find(abs(Xi0(:,1)) > max(abs(Xi0(:,1)))*tol); % indexes of non-negligible terms
    %
    Thetai = Theta(:,nzi);                      % reduced library matrix containing only non-negligible functionals
    Xii = Xi0;                                  % allocation of Xii
    Nf = size(Thetai,2);                        % number of non-negligible functionals

%
% Optimization process
%

    while err > ftol
        %
        c = c + 1;
        %
        cvx_begin quiet                         % each iteration the convergence problem dimensions are reduced and gets faster
            variable Chi(Nf,1)
            LS_norm = sum_square( dX - Thetai*Chi );
            minimize( LS_norm +  alpha*sum(W(nzi)'*abs(Chi),1) )
        cvx_end
        %
        % Sparsify solution
        %
            
            Xii(nzi,1) = Chi./W(nzi);    % weighted solution
            %
            zi = find(abs(Xii(:,1)) <= max(abs(Xii(:,1)))*tol); 
            nzi = find(abs(Xii(:,1)) > max(abs(Xii(:,1)))*tol);
            %
            W = (abs(Xii)).^(-delta);           % update weighting vector
            %
            Xii(zi,1) = 0;                      % make zero negligible coefficients

        %
        err = (abs(norm(Xii) - norm(Xi0))/abs(norm(Xii))); % error between i and i+1 solution
        Xi0 = Xii;
        %
        % If max number of iterations is exceeded, relax the regularization
        % paramters and start over
        %

            if c >= Nmax
                %
                delta = delta - 0.05;
                W = ones(size(Theta,2),1);
                c = 0;
                %
                Xi0 = Theta\dX;                             % solution from linear regression
                %
                zi = find(abs(Xi0(:,1)) <= max(abs(Xi0(:,1)))*tol); % indexes of negligible terms
                nzi = find(abs(Xi0(:,1)) > max(abs(Xi0(:,1)))*tol); % indexes of non-negligible terms
                %
                Thetai = Theta(:,nzi);                      % reduced library matrix containing only non-negligible functionals
                Xii = Xi0;                                  % allocation of Xii
                Nf = size(Thetai,2);                        % number of non-negligible functionals
                %
            end

        %
        % Update library matrix and reduce its dimension
        %

            Thetai = Theta(:,nzi)./W(nzi)';
            Nf = size(Thetai,2);

        %
    end

%
% Dimensionalize back the coefficient matrix
%

    Xii = Xii./normT'*normdX;

%
end