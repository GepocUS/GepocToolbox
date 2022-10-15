%% gen_rand_Lasso - Generates random Lasso problem(s)
%
% This function returns the ingredients A \in R^{m \times n},
% b \in \R^m and w \in \R^n defining a random Lasso problem:
% 
%   min_z 1/(2n) || A*z - b ||^2_2 + || W*z ||_1,
%
% where W = diag(w) and n > m.
%
% INPUTS:
%   - n: Integer >0. Number of decision variables of the Lasso problem
%   - m: Integer >0. Number of rows of matrix A
% 
% OPTIONAL INPUTS (must be passed using 'name' 'value' pairs):
%   - num_prob: Number of Lasso problems generated
%   - seed: Sets Matlab's rng seed to this number
%   - alpha: Scalar >0. Upper bound for the elements of w. Each element
%            in w is selected using a uniform distribution in (0, alpha]
%
% OUTPUTS:
%   - A: Matrix of dimensions (m, n, num_prob)
%   - b: Matrix of dimensions (m, 1, num_prob)
%   - w: Matrix of dimensions (n, 1, num_prob)
%   Each (:, :, i) selection of A, b and w defines a Lasso problem
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [A, b, w] = gen_rand_Lasso(n, m, varargin)

    %% Default values
    def_num_prob = 1; % Number of Lasso problems generated
    def_seed = []; % Seed used for Matlab's rng (if seed=[] then rng is not set)
    def_alpha = 0.01; % Upper bound for the elements of w

    %% Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'gen_rand_Lasso';
    % Required
    addRequired(par, 'n' , @(x) isnumeric(x) && (x>0) && x==floor(x));
    addRequired(par, 'm' , @(x) isnumeric(x) && (x>0) && (x<n) && x==floor(x));
    % Name-value parameters
    addParameter(par, 'num_prob', def_num_prob, @(x) isnumeric(x) && (x>0) && x==floor(x));
    addParameter(par, 'seed', def_seed, @(x) isempty(x) || isnumeric(x) && (x>=0) && x==floor(x));
    addParameter(par, 'alpha', def_alpha, @(x) isnumeric(x) && (x>0));
    % Parse
    parse(par, n, m, varargin{:})
    % Rename
    n = par.Results.n;
    m = par.Results.m;
    num_prob = par.Results.num_prob;
    seed = par.Results.seed;
    alpha = par.Results.alpha;

    %% Set random number generator
    if ~isempty(seed)
        rng(seed);
    end

    %% Initialize
    A = zeros(m, n, num_prob);
    b = zeros(m, 1, num_prob);
    w = zeros(n, 1, num_prob);

    %% Compute matrices
    for i = 1:num_prob
        % Compute random matrices
        A_aux = randn(m, n);
        elim = rand(m, n) < 0.1;
        A_aux = A_aux.*elim;
        b_aux = randn(m,1);
        w_aux = alpha*rand(n,1);
        % Save matrices
        A(:,:,i) = A_aux;
        b(:,:,i) = b_aux;
        w(:,:,i) = w_aux;
    end

end
