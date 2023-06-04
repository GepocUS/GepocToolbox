%% This script is an example of use of the FISTA solver
% Example script for the GepocToolbox: https://github.com/GepocUS/GepocToolbox
%
% This example shows you how to use the RFISTA solver to solve a Lasso optimization problem
% 
clear; clc;
rng(1234); % Comment this line to get different random Lasso problems on each script execution

%% Define the optimization problem
n = 300; % Number of decision variables
m = 200; % Dimension of Lasso problem
alpha = 0.01; % Determines the weight of the one-norm term

[A, b, w] = gen_rand_Lasso(n, m, 'alpha', alpha); % Generate random Lasso problem

%% We now need to construct a few things required by RFISTA.
% In particular:
%   - Function handlers for evaluation of the objective function and gradient of the smooth part
%     of the objective function. We will use the Oracle class for this.
%   - The R metric used by the solver. You can use the Lipschitz constant L of the optimization
%     problem times the identity matrix, i.e., L*eye(n), but other metrics can be used.
%   - A function handler that returns the solution of the Composite Gradient Mapping operator
%     for the optimization problem. The one used by RFISTA by default is only valid for smooth
%     unconstrained optimization problems. We will need to use a different one for the Lasso
%     problem, since it is not smooth.

% Create an Oracle instance
paramLasso = struct('H', A'*A/n, 'q', -A'*b/n, 'offset', b'*b/n, 'A', A, 'b', b, 'w', w, 'n', n, 'm', m);

f_L = @(x, p) 0.5*(x'*p.H*x + p.offset) + p.q'*x + norm(p.w.*x, 1); % Function handler for objective function value
g_L = @(x, p) p.H*x + p.q; % Function handler for the gradient of the smooth part of the objective function

lasso = Oracle(paramLasso, n, 'f_eval', f_L, 'g_eval', g_L); % Oracle for the Lasso optimization problem

% Compute metric R
% In the case of a Lasso problem, a cheap way to get a metric R for RFISTA is to take
R = diag(sum(abs(paramLasso.H)));
% Another alternative would be to compute the Lipschitz constant L of the smooth part, which in this case is equal to 
% the maximum eigenvalue of paramLasso.H, i.e., to take
% R = max(eig(paramLasso.H))*eye(n);
% However, the point of using a restart scheme is to speed the FISTA algorithm in the presence of oscillations, which
% commonly appear if you overestimate the value of L, as we have done with the above computation of R.
% The selling point here is that for sufficiently large Lasso problems, computing the exact L and then solving FISTA is
% more time-consuming than computing our R and using RFISTA.
% In any case, in this example we use the above R to highlight the restart schemes.

% Function handler for the computation of the Composite Gradient Mapping operator associated with the Lasso problem
% The function CGM_Lasso() performs this operation
cgm_handle = @(x, g, r) gt_utils.CGM_Lasso(x, g, r, w);

%% Select the other options of the RFISTA solver
x0 = []; % Initial condition of the decision variables
k_max = 3000; % Maximum number of iterations
j_max = 20; % Maximum number of restarts
tol = 1e-9; % Exit tolerance of the solver
restart = 'l'; % This string determines the restart strategy used by the solver. See the RFISTA documentation for the available options
genHist = 2; % This option determines how much information of the iterates is computed and returned by the solver.
             % Small values make the solver go faster, but less information is returned.
verbose = 3; % Determines the verbose level
exit_at_tol = true; % If true, then the solver exits as soon as it finds a point satisfying the suboptimality condition
                    % Otherwise, it only checks it at restart points

%% Call the RFISTA solver
[x_opt, f_opt, e_flag, Hist] = RFISTA(lasso.get_f_hand(), lasso.get_g_hand(), R, cgm_handle, x0, 'k_max', k_max,...
    'j_max', j_max, 'tol', tol, 'restart', restart, 'genHist', genHist, 'verbose', verbose, 'exit', exit_at_tol);

% You should see the report of the solver on the Command window, indicating the progress of the solver during the restart points
% x_opt is the (sub)optimal solution of the Lasso problem
% f_opt is (sub)optimal value
% e_flag indicates the exit flag of the solver (see its documentation for possible values)
% Hist is a structure that contains information of the solver iterates. Its fields depend on the value of genHist

%% Lets compare against the non-restarted version
% To do so, we change
restart = 'No';
% which disables the restart (it simply runs FISTA).
[x_No, f_No, e_flag_No, Hist_No] = RFISTA(lasso.get_f_hand(), lasso.get_g_hand(), R, cgm_handle, x0, 'k_max', k_max,...
    'j_max', j_max, 'tol', tol, 'restart', restart, 'genHist', genHist, 'verbose', verbose, 'exit', exit_at_tol);

% Lets plot the distance the evolution of both calls to the solver in terms of the distance to the optimal solution, i.e.,
% we plot f(x_k) - f_opt, in logarithmic scale, using the restart scheme 'l' and no restart scheme ('restart' = 'No').
% The plot is made using the Fig class
f1 = Fig(1);
f1.y_scale("log");
f1.plot(0:Hist_No.k, Hist_No.fx(:) - f_opt); % No restart, plotted in blue
f1.plot(0:Hist.k, Hist.fx(:) - f_opt); % Selected restart scheme, plotted in red

% In general, the blue line in the plot should have a significant amount of oscillations, whereas the red line, which is the 
% one corresponding to the selected restart scheme, should have reached the optimal solution in less time (and have less oscillations).

