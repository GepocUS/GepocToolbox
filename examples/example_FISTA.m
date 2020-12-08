%% This script is an example of use of the FISTA solver
% Example script for the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 
% We solve an unconstrained QP problem
%
clear; clc;

%% Define the optimization problem
% Our problem is to minimize
%
%   min 0.5*z'*H*z + q'*z,
%
% where H is a positive definite matrix.

n = 20; % Number of decision variables of the problem

% Construct a positive definite matrix
H = 10*rand(n);
H =  (H + H')/2 + n*eye(n);

% Construct vector q
q = 20*randn(n, 1);

%% Define function handlers for FISTA

% Fista requires a function handler that returns the gradient of the objetcive function. We define this bellow,
grad = @(z) H*z + q;

% Optionalyl, it can also be given a function handler that returns the value of the objective function.
func =  @(z) 0.5*z'*H*z + q'*z;

%% Define other ingredients of FISTA

% The only other required argument of FISTA is a metric R, which is a vector whose dimension must be the same as the number of decision variables
R = sum(abs(H), 2);

% There are other optional arguments. For example
k_max = 500;  % Maximum number of iterations
tol = 1e-4; % Exit tolerance of the algorithm
verbose = 3; % Determines how much information is displayed to the user (0, 1, 2 or 3)
gen_Hist = 2; % Determines what variables are stored during the algorithms execution
z_0 = zeros(n, 1); % Initial value of the decision variables

%% Solve the problem
% The second argument is the plus_op function handler. Please refer to the functions documentation for more information about it.
[z_opt, f_opt, e_flag, Hist] = FISTA(grad, R, [], z_0, 'k_max', k_max, 'tol', tol, 'verbose', verbose, 'genHist', gen_Hist, 'z_0', z_0, 'f_eval', func);

%% Compare the result to quadprog
z_qp = quadprog(H, q, [], []);
diff_with_quadprog = max(abs(z_opt - z_qp))

