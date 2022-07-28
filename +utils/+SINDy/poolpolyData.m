function Theta = poolpolyData(X,Nr,polyorder)
%%                          poolpolyData.m
%--------------------------------------------------------------------------
%
% This function obtains the library of functions Theta from the regression 
% fit dyin/dt = Theta(yin)*Xi, where yin denotes the states of the system. The
% library Theta will be composed of a set of polynomials of order polyorder
%
% INPUTS
%
%   X         : matrix containing states of the system for a certain
%               timespan of dimension Nt x Nr
%   Nr        : number of states of the system
%   polyorder : array containing polynomial orders span by the library
%
% OUTPUT
%
%   Theta     : library of functions of dimension Nt x Nf
%
%--------------------------------------------------------------------------
%
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

Nt = size(X,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);
ind = 1;
%
    if ~isempty(find(polyorder == 0))
        Theta(:,ind) = ones(Nt,1);
        ind = ind+1;
    end
%
    if ~isempty(find(polyorder == 1))
        % 
        for i=1:Nr
            Theta(:,ind) = X(:,i);
            ind = ind+1;
        end
    end
%
    if ~isempty(find(polyorder == 2))
        % poly order 2
        for i=1:Nr
            for j=i:Nr
                Theta(:,ind) = X(:,i).*X(:,j);
                ind = ind+1;
            end
        end
    end
%
    if ~isempty(find(polyorder == 3))
        % poly order 3
        for i=1:Nr
            for j=i:Nr
                for k=j:Nr
                    Theta(:,ind) = X(:,i).*X(:,j).*X(:,k);
                    ind = ind+1;
                end
            end
        end
    end
%
    if ~isempty(find(polyorder == 4))
        % poly order 4
        for i=1:Nr
            for j=i:Nr
                for k=j:Nr
                    for l=k:Nr
                        Theta(:,ind) = X(:,i).*X(:,j).*X(:,k).*X(:,l);
                        ind = ind+1;
                    end
                end
            end
        end
    end
%
    if ~isempty(find(polyorder == 5))
        % poly order 5
        for i=1:Nr
            for j=i:Nr
                for k=j:Nr
                    for l=k:Nr
                        for m=l:Nr
                            Theta(:,ind) = X(:,i).*X(:,j).*X(:,k).*X(:,l).*X(:,m);
                            ind = ind+1;
                        end
                    end
                end
            end
        end
    end
%
    if ~isempty(find(polyorder == 6))
        % poly order 5
        for i=1:Nr
            for j=i:Nr
                for k=j:Nr
                    for l=k:Nr
                        for m=l:Nr
                            for Nt = m:Nr
                                Theta(:,ind) = X(:,i).*X(:,j).*X(:,k).*X(:,l).*X(:,m).*X(:,Nt);
                                ind = ind+1;
                            end
                        end
                    end
                end
            end
        end
    end
%
    if ~isempty(find(polyorder == 7))
        % poly order 5
        for i=1:Nr
            for j=i:Nr
                for k=j:Nr
                    for l=k:Nr
                        for m=l:Nr
                            for Nt = m:Nr
                                for o = Nt:Nr
                                    Theta(:,ind) = X(:,i).*X(:,j).*X(:,k).*X(:,l).*X(:,m).*X(:,Nt).*X(:,o);
                                    ind = ind+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
%
end