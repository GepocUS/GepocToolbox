%% gen_rand_boxQP - generate random box-constrained QP problem
%
% INPTUS:
%   - n: Number of decision variables
%   - m: Number of equality constraints
%   - condH: Condition number of H
%   - q_mag: Magnitude of q
%   - seed: Seed of rng(), for reproducibility
%
% OUTPUTS:
%   - H, q: Cost function matrix and vector
%   - G, b: Equality constraints matrix and vector
%   - LB, UB: Lower and upper bound of the box constraints
%
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
%

function [H, q, G, b, LB, UB] = gen_rand_boxQP(n, m, condH, q_mag, seed)

    % Check arguments
    if nargin == 5
        rng(seed);
    elseif nargin == 3
        condH = 100;
    elseif nargin == 2
        condH = 100;
        q_mag = 1;
    end
    
    % Compute H and q
    S=diag(exp(-log(condH)/4:log(condH)/2/(n-1):log(condH)/4));
    [U,~] = qr((rand(n,n)-.5)*200);
    [V,~] = qr((rand(n,n)-.5)*200);
    
    H = U*S*V';
    H = H'*H;
    q = q_mag*randn(n,1);
    
    % Compute inequality constraints
    LB = -rand(n,1)-.1;
    UB = rand(n,1)+.1;
    
    % Compute equality constraints
    G = randn(m, n);
    b = rand(m, 1);

end

