%% precond_Ruiz - Ruiz equilibration algorithm
%
% This function implements the Ruiz equilibration algorithm.
% 
% Consider a matrix R = [P A'; A 0]. This method returns a diagonal
% matrix S = blkdiag(D, E) such that the condition number of S*R*S
% is small.
%
% INPUTS:
%   - P, A: Submatrices of matrix R
%   - tau: Exit tolerance of the Ruiz algorithm
% 
% OUTPUTS:
%   - D, E: Submatrices of matrix S
%   - P, A: Preconditioned matrices P and A
%   - R: Preconditioned matrix R
%   - delta: Final value of delta at the exit of the algorithm
%   - k: Number of iterations of the modified Ruiz algorithm
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [D, E, P, A, R, delta, k] = precond_Ruiz(P, A, tau)

    %% Initialize
    n = size(P, 1);
    m = size(A, 1);
    D = eye(n);
    E = eye(m);
    delta = ones(m+n, 1);
    R = [P, A'; A zeros(m)];
    k = 0;
    done = false;
    
    while ~done
        
        k = k+1;
        
        for i=1:(n+m)
            ri = max(abs(R(:, i)));
            if ri > tau
                delta(i) = 1/sqrt(ri);
            end
        end
        hatD = diag(delta(1:n));
        hatE = diag(delta(n+1:end));
        D = hatD*D;
        E = hatE*E;
        P = hatD*P*hatD;
        A = hatE*A*hatD;
        R = [P, A'; A zeros(m)];
        
        if max(abs(ones(m+n, 1) - delta)) <= tau
            done = true;
        end
        
    end
    
end
