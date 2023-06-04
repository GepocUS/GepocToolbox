%% FISTA - Fast Iterative Shrinking Threshold Algorithmn
% [z_opt, f_opt, e_flag, Hist] = FISTA(g_eval, R, plus_op, z_0) computes the
% minimizer x_opt of the optimization problem,
% 
%   z_opt = arg min F(z)
%              s.t. z in Z
% 
% with F(z) = f(z) + d(z), f(z) a smooth convex function and d(z) a closed convex function.
% 
% OUTPUTS
%  z_opt: Solution of the optimization problem.
%  f_opt: Value of F at the optimum, i.e. f_opt = F(z_opt)
%  e_flag: Integer that indicates the exit condition of the algorithm.
%           1: Algorithm converged successfuly to a solution within the given tolerance.
%          -1: Algorithm did not converge before reaching the maximum allowed iterations (k_max).
%  Hist: Structure containing additional data. See documentation for further details.
% 
% INPUTS
%  grad: Handler to a function that returns the gradient of F at point z. g(z) = grad(z)
%  R: Metric used by the algorithm. Diagonal matrix. R > 0
%  plus_op: Handler to a function that returns the solution of the plus operator, which is the 
%           solution of the following optimization problem for a given value of y, a function 
%           handler to grad and the metric R ( y+(y) = plus_op(y, grad, R) ),
%               y_plus(y) = arg min_z d(x) + grad(y)'*(x - y) + 0.5*(x-y)'*R*(x-y)
%                               s.t. x in Z
%           If no value is provided, a default function plus_op is used, which assumes
%           that the is no non-smooth function d(z), i.e. F(z) = f(z), and that z is not constrained.
%  z_0: Initial value of z (at iteration 0).
%       If no value is provided it defaults to a vector of zeros.
% 
% OPTIONAL INPUTS (must be passed using 'name' 'value' pairs)
%  f_eval: Handler to a function that returns the value of F at a point z. f(z) = f_eval(z)
%  k_min: Initial value of the minimum number of iterations of FISTA algorithm
%         Integer >= 0. Defaults to 0.
%  k_max: Maximum number of iterations. Integer > 0. Defaults to 10000
%  tol: Exit tolerance of the algorithm. Must be > 0. Default to 1e-9
%  genHist: Determines how much information is generated for output Hist. Can be 0, 1 or 2.
%           Higher values increase the amount of information returned in Hist. Defaults to 1.
%           A value of 2 can increase the computation time of the function significantly.
%  verbose: Integer that can take value 0, 1, 2 or 3. It indicates how much information is
%           displayed to the user. No information is displayed if set to 0. Defaults to 1
%  ISTA: If true, then FISTA is turned into ISTA (acceleration is disabled)
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [z_opt, f_opt, e_flag, Hist] = FISTA(g_eval, R, varargin)
    timer_total = tic;

    %% Default values
    def_z0 = zeros(length(R), 1); % Default value of the initial condition
    def_plus_op = @default_plus_operator; % Default value of the plus operator
    def_f_eval = @(x) []; % Default value of f_eval
    def_k_min = 0; % Default value of k_min
    def_k_max = 10000; % Default value of k_max
    def_tol = 1e-5; % Default value of the exit tolerance
    def_genHist = 0; % Default amount of data generated for Hist: Only stores fx, gy and basic data
    def_verbose = 1; % Default amount of information displayed: Displayes end of algorithm info only
    def_t0 = 1; % Default value of t0
    def_initialize_Hist = true; % Default value of initialize_Hist
    def_print_counter_step = 10; % Number of iterations between verbose to user (only if verbose >= 3)
    def_ISTA = false; % If true, implement ISTA instead of FISTA, i.e., disable acceleration
    
    %% Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'FISTA';
    % Required
    addRequired(par, 'g_eval', @(x) isa(x, 'function_handle'));
    addRequired(par, 'R', @(x) isnumeric(x) && ((size(x,1)==size(x,2)) || (min(size(x))==1)));
    % Optional
    addOptional(par, 'plus_op', def_plus_op, @(x) isa(x, 'function_handle') || isempty(x));
    addOptional(par, 'z_0', def_z0, @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    % Name-value parameters
    addParameter(par, 'f_eval', def_f_eval, @(x) isa(x, 'function_handle'));
    addParameter(par, 'k_min', def_k_min, @(x) isnumeric(x) && (x>=0) && x==floor(x));
    addParameter(par, 'k_max', def_k_max, @(x) isnumeric(x) && (x>0) && x==floor(x));
    addParameter(par, 'tol', def_tol, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'genHist', def_genHist, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'verbose', def_verbose, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 't_0', def_t0, @(x) isnumeric(x));
    addParameter(par, 'initialize_Hist', def_initialize_Hist, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'print_counter_step', def_print_counter_step, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'ISTA', def_ISTA, @(x) islogical(x) || x==1 || x==0);
    % Parse
    parse(par, g_eval, R, varargin{:})
    % Rename
    g_eval = par.Results.g_eval;
    R = par.Results.R;
    plus_op = par.Results.plus_op;
    f_eval = par.Results.f_eval;
    z_0 = par.Results.z_0;
    k_min = par.Results.k_min;
    k_max = par.Results.k_max;
    tol = par.Results.tol;
    genHist = par.Results.genHist;
    verbose = par.Results.verbose;
    t0 = par.Results.t_0;
    initialize_Hist = par.Results.initialize_Hist;
    print_counter_step = par.Results.print_counter_step;
    % Check arguments
    if isempty(z_0); z_0 = def_z0; end % Default x_0 if an empty array is provided
    if size(z_0, 2)>1; z_0 = z_0'; end % Make x_0 a column vector
    if k_min >= k_max
        if verbose>=3
            warning('FISTA:InputWarning', 'Value of k_min >= k_max. Algorithm will exit at k_max iterations allways');
        end
    end
    if tol==0; warning('FISTA:InputWarning', 'Value of tol=0. Algorithm will never converge to this tolerance'); end
    if genHist>2
        genHist = 2;
    end
    if verbose>3
        verbose = 3;
    end
    if isempty(plus_op)
        plus_op = def_plus_op;
    end
    
    %% Initialization
    done = false; % Flag that indicates the end of the algorithm
    k = 0; % Counter for the total number of FISTA iterations performed
    yk1 = plus_op(z_0, g_eval, R); % Value of y at iteration k-1. This is the value of y_0
    xk = zeros(length(z_0), 1); % Value of x at iteration k
    xk1 = yk1; % Value of x at iteration k-1. This is the value of x_0
    tk1 = t0; % Value of t at iteration k-1
    print_counter = 0; % Counter used to determine when to show info to user
    print_title = 0; % Conter used to determine when to reprint titles line
    
    % Declare historics
    if initialize_Hist
        if genHist >= 1
            h_functional = zeros(1, k_max); % Historic for the functional
            h_dual_comp_grad_norm = zeros(1, k_max); % Historic for the norm of the dual of the composite gradient mapping of y_k (used as the exit condition)
        end
        if genHist >= 2
            x = zeros(length(R), k_max); % Historic of x_k
            y = zeros(length(R), k_max); % Historic of y_k
            t = zeros(1, k_max); % Historic of t_k
            h_gradient = zeros(length(R), k_max); % Historic for the gradient
        end
    end
    
    % Save historics
    if genHist >= 1
        h_functional(1) = f_eval(xk1);
    end
    if genHist >= 2
        x(:,1) = xk1;
        y(:,1) = yk1;
        t(1) = tk1;
    end
    
    % Initial verbose to user
    if verbose >= 2
        fprintf('\tStarting FISTA algorithm\n');
    end
    if verbose >= 3
        fprintf('\tParamenets:\tk_min = %d\tk_max = %d\ttol = %g\n', k_min, k_max, tol);
        fprintf('\t\tNumber of decision variables: %d\n', length(R));
        fprintf('\n');
        fprintf('\t   Iter ||g(y_k)||_*  f(x_k)');
        fprintf('\n');
    end
    
    %% Algorithhm
    while ~done
        % Increate k and k_total counters
        k = k + 1;
        
        % Main steps
        xk = plus_op(yk1, g_eval, R);
        tk = 0.5*(1 + sqrt(1+4*tk1^2));
        yk = xk + (tk1 - 1)*(xk - xk1)/tk;
        if par.Results.ISTA
            yk = xk;
        end
        
        % Compute variables
        gy = g_dual_norm(yk1, xk, R); % Compute || g(y_{k-1})||_*. It is computed here to take advantage of x_k = y_{k-1}+
        
        % Save historics
        if genHist >= 1
            h_functional(k+1) = f_eval(xk);
            h_dual_comp_grad_norm(k) = gy;
        end
        if genHist >= 2
            x(:,k+1) = xk;
            y(:,k+1) = yk;
            t(k+1) = tk;
        end
        
        % Exit condition
        if (gy <= tol) && (k >= k_min)
            done = true;
            e_flag = 1;
        elseif k >= k_max
            done = true;
            e_flag = -1;
        end
        
        % Update variables
        xk1 = xk;
        yk1 = yk;
        tk1 = tk;
        
        % Verbose info to user
        if verbose >= 3
            if print_title == 10
                print_title = 0;
                fprintf('\t   Iter ||g(y_k)||_*  f(x_k)');
                fprintf('\n');
            end
            if genHist == 0
                fx = f_eval(xk);
            else
                fx = h_functional(k+1);
            end
            print_counter = print_counter + 1;
            if print_counter == print_counter_step
                print_title = print_title + 1;
                print_counter = 0;
                fprintf('\t%6d   %4.5e  %4.5e ', k, gy, fx);
                fprintf('\n');
            end
        end
        
    end
    
    % Compute the last value of gy
    gy = g_dual_norm_handler(yk, plus_op, g_eval, R);
    
    % Create Historic and return variables
    Hist.k = k; % Number of iterations
    Hist.exit_tol = gy; % Exit tolerance achieved. Value of the dual norm of the composite gradient mapping of y_{k-1} at the exit of the algoritm
    if genHist >= 1
        Hist.functional = h_functional(1:k);
        h_dual_comp_grad_norm(k+1) = gy;
        Hist.dual_comp_grad_norm = h_dual_comp_grad_norm;
    end
    if genHist >= 2
        Hist.x = x;
        Hist.y = y;
        Hist.t = t;
    end
    z_opt = xk;
    f_opt = f_eval(xk);
    
    % Timers
    Hist.time.total = toc(timer_total);
    
    % Final verbose to user
    if verbose >= 3
        fprintf('\n');
    end
    if e_flag <= 0
        if verbose >= 1
            if e_flag == -1
                fprintf('\tWARNING: Algorithm did not converge within the allowed number of iterations %d\n', k_max);
            end
            fprintf('\tFISTA exit flag is %d', e_flag);
            fprintf('\tElapsed time between start and end of algorithm: %fs\n', Hist.time.total);
        end
    else
        if verbose >=2
            fprintf('\tAlgorithm exited successfully with exit flag %d and the following results\n', e_flag);
            fprintf('\t\tNumber of iterations (k): %d (of a maximum of %d)\n', k, k_max);
            fprintf('\t\tf(x_opt) = %g\t\t||g(x_opt)||_* = %g\n', f_opt, gy);
            fprintf('\t\tElapsed time between start and end of algorithm: %fs\n', Hist.time.total);
        end
    end
        
end

%% Auxiliary functions

% Computes ||g(y)||_* for the given value y, function handlers plus_op and g, and the metric R
% R must be a vector
function gy = g_dual_norm_handler(y, plus_op, g, R) % TODO: can't this be simplified
	y_plus = plus_op(y, g, R);
    gy = sqrt( (y - y_plus)'*R*(y - y_plus) );
end

% Computes ||g(y)||_* for the given y, y+ and metric R
% R must be a vector
function gy = g_dual_norm(y, y_plus, R)
    gy = sqrt( (y - y_plus)'*R*(y - y_plus) );
end

% Computes the composite gradient mapping g(y) for the given y, function handlers plus_op and g, and metric R
% R must be a vector
function gy = comp_g_handler(y, plus_op, g, R)
    gy = R*(y - plus_op(y, g, R));
end

% Computes the composite gradient maping g(y) for the given y, y+ and metrix R
function gy = comp_g(y, y_plus, R)
    gy = R*(y - y_plus);
end

% Computes the plus operator for phi=0 and no restrictions
% R must be a vector
% This function is used as the default function for the plus_op function handler if non is provided by the user
function y_plus = default_plus_operator(y, g, R)
    y_plus = y - R\g(y);
end

