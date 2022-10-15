%% RFISTA - Restart Fast Iterative Shrinking Threshold Algorithmn
% 
% [x_opt, f_opt, e_flag, Hist] = RFISTA(f_eval, grad, R, plus_op, x_0) computes the
% minimizer x_opt of the optimization problem,
% 
%   x_opt = arg min F(x)
%              s.t. x in X
% 
% with F(x) = f(x) + d(x), f(x) a smooth convex function and d(x) a closed convex function.
% 
% OUTPUTS:
%  - x_opt: Solution of the optimization problem.
%  - f_opt: Value of F at the optimum, i.e. f_opt = F(x_opt)
%  - e_flag: Integer that indicates the exit condition of the algorithm.
%         1: Algorithm converged successfully to a solution within the given tolerance.
%        -1: Algorithm did not converge before reaching the maximum allowed iterations (k_max).
%        -2: Algorithm did not converge before reaching the maximum allowed restarts (j_max).
%  - Hist: Structure containing additional data. See documentation for further details.
% 
% INPUTS:
%  - f_eval: Handler to a function that returns the value of F at a point x. f(x) = f_eval(x)
%  - grad: Handler to a function that returns the gradient of F at point x. g(x) = grad(x)
%  - R: Metric used by the algorithm. R > 0
%  - plus_op: Handler to a function that returns the solution of the plus operator, which is the 
%             solution of the following optimization problem for a given value of y, a function 
%             handler to grad and the metric R ( y+(y) = plus_op(y, grad, R) ),
%               y_plus(y) = arg min_x d(x) + grad(y)'*(x - y) + 0.5*(x-y)'*R*(x-y)
%                               s.t. x in X
%             If no value is provided, a default function plus_op is used, which assumes
%             that the is no non-smooth function d(x), i.e. F(x) = f(x), and that x is not constrained.
%  - x_0: Initial value of x (at iteration 0).
%         If no value is provided it defaults to a vector of zeros.
% 
% OPTIONAL INPUTS (must be passed using 'name' 'value' pairs):
%  - k_min: Initial value of the minimum number of iterations of FISTA algorithm
%         Integer >= 0. Defaults to 0.
%  - k_max: Maximum number of iterations. Integer > 0. Defaults to 10000
%  - j_max: Maximum number of restarts. Integer > 0. Defaults to Inf
%  - tol: Exit tolerance of the algorithm. Must be > 0. Default to 1e-9
%  - restart: Restart method used by the algorithm. Must be a number between 0 and 6 or
%             one of the following strings of characters. Defaults to 1 (or 'l').
%             restart = {'No', 'l', 'f', 'g', 'f*', '1/e', 'k+1'}
%             Numbers 0 to 6 correspond to one of the 7 previous strings.
%             Can also be a function handler to a restart condition function. See documentation
%  - genHist: Determines how much information is generated for output Hist. Can be 0, 1 or 2.
%             Higher values increase the amount of information returned in Hist. Defaults to 1.
%             A value of 2 can increase the computation time of the function significantly.
%  - verbose: Integer that can take value 0, 1, 2 or 3. It indicates how much information is
%             displayed to the user. No information is displayed if set to 0. Defaults to 1
%  - f_opt: Optimal value of the objective function. It is used in restart scheme number 4 ('f*').
%           This parameter is only used if restart 4 is selected. It defaults to 0.
%  - exit_at_tol: If 1, the algorithm exits as soon as a value of x_k is found that satisfies the
%                 exit condition. If 0, the exit condition is only checked at points r_j.
%                 It defaults to 1.
%  - disp_j_step: Determines how many j iterations information is displayed to the user. Only used
%                 if verbose = 3. It defaults to 1 (display information every iteration).
%                 This option has not been fully tested. Use at your own discretion.
%  - use_MFISTA: If 1, the algorithm used is the monotone version of FISTA (MFISTA), which
%                can be found in "Fast Gradient-Based Algorithms for Constrained Total Variation
%                Image Denoising and Deblurring Problems", by Beck A. and Teboulle M. in 
%                IEEE Transaction on Image Processing. If 0, the standard FISTA algorithm is used.
%                Defaults to 0.
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [x_opt, f_opt, e_flag, Hist] = RFISTA(f_eval, grad, R, varargin)

    %% Default values
    def_x0 = zeros(length(R), 1); % Default value of the initial condition
    def_plus_op = @default_plus_operator; % Default value of the plus operator
    def_k_min = 0; % Default value of k_min
    def_k_max = 10000; % Default value of k_max
    def_j_max = Inf; % Default alue of j_max
    def_tol = 1e-9; % Default value of the exit tolerance
    def_restart = 1; % Default exit condition: Restart with 1/sqrt(mu) convergence
    max_restart = 8; % (Not a default) Value of the maximum allowed value for restart
    def_genHist = 0; % Default amount of data generated for Hist: Only stores fx, gy and basic data
    def_verbose = 1; % Default amount of information displayed: Displays end of algorithm info only
    def_fopt_ext = 0; % Default value of the optimal value of the objective function (it is used in one of the restart schemes)
    def_t0 = 1; % Default value of t0
    def_exit_at_tol = 1; % Default value of variable exit_at_tol. Determines if the algorithm terminates as soon as g(y_k)<=tol
    def_disp_j_step = 1; % Default value of variable disp_j_step. Determines how many j iterations information is displayed to the user
    def_use_MFISTA = 0; % Default value of use_MFISTA. Determines is subFISTA or subMFISTA are used as the FISTA algorithm
    % Restart methods
    restart_chars = {'No', 'l', 'f', 'g', 'fopt', 'e', 'k+1', 'delay'};
    
    %% Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'RFISTA';
    % Required
    addRequired(par, 'f_eval' , @(x) isa(x, 'function_handle'));
    addRequired(par, 'grad'   , @(x) isa(x, 'function_handle'));
    addRequired(par, 'R', @(x) isnumeric(x) && ((size(x,1)==size(x,2)) || (min(size(x))==1)));
    % Optional
    addOptional(par, 'plus_op', def_plus_op, @(x) isa(x, 'function_handle') || isempty(x));
    addOptional(par, 'x_0', def_x0, @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    % Name-value parameters
    addParameter(par, 'k_min', def_k_min, @(x) isnumeric(x) && (x>=0) && x==floor(x));
    addParameter(par, 'k_max', def_k_max, @(x) isnumeric(x) && (x>0) && x==floor(x));
    addParameter(par, 'j_max', def_j_max, @(x) isnumeric(x) && (x>0) && x==floor(x));    
    addParameter(par, 'tol', def_tol, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'restart', def_restart,...
                 @(x) ischar(x) || (x==floor(x) && (x>=0) && (x<=max_restart)) || isa(x, 'function_handle'));
    addParameter(par, 'genHist', def_genHist, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'verbose', def_verbose, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'f_opt', def_fopt_ext, @(x) isnumeric(x));
    addParameter(par, 't_0', def_t0, @(x) isnumeric(x));
    addParameter(par, 'exit_at_tol', def_exit_at_tol, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'disp_j_step', def_disp_j_step, @(x) isnumeric(x) && (x>=0) && x==floor(x));
    addParameter(par, 'use_MFISTA', def_use_MFISTA, @(x) islogical(x) || x==1 || x==0);
    % Parse
    parse(par, f_eval, grad, R, varargin{:})
    % Rename
    f_eval = par.Results.f_eval;
    grad = par.Results.grad;
    R = sparse(par.Results.R);
    plus_op = par.Results.plus_op;
    x_0 = par.Results.x_0;
    k_min = par.Results.k_min;
    k_max = par.Results.k_max;
    j_max = par.Results.j_max;
    tol = par.Results.tol;
    restart = par.Results.restart;
    genHist = par.Results.genHist;
    verbose = par.Results.verbose;
    f_opt_ext = par.Results.f_opt;
    t0 = par.Results.t_0;
    exit_at_tol = par.Results.exit_at_tol;
    disp_j_step = par.Results.disp_j_step;
    use_MFISTA = par.Results.use_MFISTA;
    % Check arguments
    if isempty(x_0); x_0 = def_x0; end % Default x_0 if an empty array is provided
    if size(x_0, 2)>1; x_0 = x_0'; end % Make x_0 a column vector
%     if min(size(R))~=1 && ~isdiag(R) % Check if R is diagonal
%         warning('RFISTA:InputWarning',...
%                 'Matrix R must be a diagonal matrix or a vector. The diagonal elements will be used');
%     end
%     if min(size(R))~=1
%         R = diag(R); % Turn R into a vector if a matrix was given
%     end
%     if all(R < 0); error('RFISTA:InputError', 'Matrix R must be strictly positive definite'); end % R > 0
    if k_min >= k_max
        if verbose>=3
            warning('RFISTA:InputWarning', 'Value of k_min >= k_max. Algorithm will exit at k_max iterations allways');
        end
    end
    if tol==0; warning('RFISTA:InputWarning', 'Value of tol=0. Algorithm will never converge to this tolerance'); end
    if genHist>2
        %warning('RFISTA:InputWarning', 'Value of genHist greater that maximum of 2. Taking genHist = 2');
        genHist = 2;
    end
    if verbose>3
        %warning('RFISTA:InputWarning', 'Value of verbose greater that maximum of 3. Taking verbose = 3');
        verbose = 3;
    end
    % Get restart method
    if ischar(restart)
        if strcmp(restart, restart_chars{1});      restart = 0;
        elseif strcmp(restart, restart_chars{2});  restart = 1;
        elseif strcmp(restart, restart_chars{3});  restart = 2;
        elseif strcmp(restart, restart_chars{4});  restart = 3;
        elseif strcmp(restart, restart_chars{5});  restart = 4;
        elseif strcmp(restart, restart_chars{6});  restart = 5;
        elseif strcmp(restart, restart_chars{7});  restart = 6;
        elseif strcmp(restart, restart_chars{8});  restart = 7;
        else
            warning('RFISTA:InputWarning', 'Value of restart not recognized. Using default restart method: 1');
            restart = def_restart;
        end
    elseif isa(restart, 'function_handle')
        if genHist < 2
            warning('RFISTA:InputWarning', 'External restart method provided by a function handler is only available if genHist = 2. Setting genHist = 2');
            genHist = 2;
        end
    end
    if isempty(plus_op)
        plus_op = def_plus_op;
    end
    
    %% Initialization
    done = false; % Flag that indicates the end of the algorithm
    print_counter = 0;
    j = 0; % Counter for number of restarts
    k = 0; % Counter for the total number of FISTA iterations performed
    k_min_initial = k_min;
    r(:,1) = x_0; % Values of x_k at the restart points j
    %r_j_0 = x_0; % Value of r_j given to subFISTA algorithm as the initial condition of each j iteration
    n(1) = 0; % Number of iterations at each restart point j
    c(1) = 0; % Value of k at each restart j in which the restart condition would have exited had condition k>=k_min not existed
    fr(1) = f_eval(x_0);  % Value of the cost function f(r_j) at points r_j. This is stored independently of genHist value
    gr(1) = g_dual_norm_handler(x_0, plus_op, grad, R); % Value of ||g(r_j)||_* at points r_j. This is stored independently of genHist value
    k_res(1) = 0; % Value of k at each of the restart points j (helps in plotting fx alongside fx_k, for example)
    k_min_j(1) = k_min; % Value of k_min at the start of each iteration k
    doubling_j = NaN; % Indicates that a doubling happened at iteration j TODO: Add to Hist
    if genHist >= 1
        fx(1) = fr(1); % Value of the cost function at points x_k
        fy(1) = fr(1); % Value of the cost function at points y_k
        gx(1) = gr(1); % Value of ||g(x_k)||_* at points x_k
        gy(1) = gr(1); % Value of ||g(y_k)||_* at points y_k
    end
    if genHist >= 2
        x(:,1) = x_0; % Historic for x_k
        y(:,1) = x_0; % Historic for y_k
        t(1) = 0; % Historic for the values of t_k
    end
    
    % Verbose initial report to user
    if verbose >= 2
        fprintf('\tStarting RFISTA algorithm using restart scheme %s\n', restart_chars{restart+1});
    end
    if verbose >= 3
        fprintf('\tParamenets:\tk_min = %d\tk_max = %d\ttol = %g\n', k_min, k_max, tol);
        fprintf('\t\tNumber of decision variables: %d\n', size(R, 1));
        if restart == 5
            fprintf('\t\tf_opt = %g', f_opt_ext);
        end
        fprintf('\t\n');
        if exit_at_tol
            fprintf('\tAlgorithm will exit if ||g(x_k)||_* <= %g is satisfied (exit_at_tol = 1)\n', tol);
        else
            fprintf('\tAlgorithm will exit if ||g(r_j)||_* <= %g is satisfied (exit_at_tol = 0)\n', tol);
        end
        fprintf('\n');
        fprintf('\tIter_j Iter_k  n_j   c      g(r)          f(r)    f_flag k_min ');
        if restart == 1
            fprintf('Doubling');
        end
        fprintf('\n');
    end
    
    tic;
    %% Algorithm
    while ~done
        j = j+1;
        % Call FISTA algorithm
        if ~use_MFISTA
            [r(:,j+1), n(j+1), f_flag, H] = subFISTA(r(:,j), f_eval, grad, plus_op, R, k, k_min, k_max, restart, tol, genHist, verbose, f_opt_ext, t0, exit_at_tol);
        else
            [r(:,j+1), n(j+1), f_flag, H] = subMFISTA(r(:,j), f_eval, grad, plus_op, R, k, k_min, k_max, restart, tol, genHist, verbose, f_opt_ext, t0, exit_at_tol);
        end
        k = k + n(j+1);
        % Update historics
        c(j+1) = H.c;
        fr(j+1) = H.fx(end); % Value of the cost function f(r_j) at points r_j
        gr(j+1) = H.gx(end); % Value of ||g(r_j)||_* at points r_j
        k_res(j+1) = k; % Value of k at each of the restart points j
        k_min_j(j+1) = k_min; % Value of k_min at the start of each iteration k
        if genHist >= 1
            fx = [fx H.fx(2:end)]; % Value of the cost function at points x_k
            fy = [fy H.fy(2:end)]; % Value of the cost function at points y_k
            gx = [gx H.gx(2:end)]; % Value of ||g(x_k)||_* at points x_k
            gy = [gy H.gy(2:end)]; % Value of ||g(y_k)||_* at points y_k
        end
        if genHist >= 2
            x = [x H.x(:,2:end)]; % Historic for x_k
            y = [y H.y(:,2:end)]; % Historic for y_k
            t = [t H.t(2:end)]; % Historic for the values of t_k
        end
        f_m(j) = H.f_m;
        x_m(:,j) = H.x_m;
        k_m(j) = H.k_m;
        
        % Doubling step for LCR-FISTA: Compute the value of k_min
        k_min_prev = k_min;
        if restart == 1
            if j>1 && fr(j) - fr(j+1) > (fr(j-1) - fr(j))/exp(1) % (We only check this after the first iteration)
                doubling_j = 1;
                k_min = 2*k_min; % Double the value of k_min if the functional reduction condition is not satisfied
            else
                k_min = n(j+1); % Otherwise, take the new value of k_min as the number of iterations n(j+1)
                doubling_j = 0;
            end
        end
        
        % Exit conditions
        if f_flag > 0 && gr(j+1) <= tol % Exit if a value of r_j (or x_k - see argument exit_at_tol) is foud that satisfies the tolerance
            done = true;
            e_flag = 1;
        end
        if f_flag <= 0 % Exit because maximum number of iterations has been exceeded
            done = true;
            e_flag = -1;
        end
        if ~done && (j>=j_max)
            done = true;
            e_flag = -2;
        end
        
        % Verbose info to user
        if verbose >= 3
            if rem(j-1, disp_j_step)==0
                print_counter = print_counter + 1;
                if print_counter == 10
                    print_counter = 0;
                    fprintf('\tIter_j Iter_k  n_j   c      g(r)          f(r)    f_flag k_min ');
                    if restart == 1
                        fprintf('Doubling');
                    end
                    fprintf('\n');
                end
                fprintf('\t%6d %6d %4d %4d  %4.5e  %4.5e      %d  %4d', j, k, n(j+1), c(j+1), gr(j+1), fr(j+1), f_flag, k_min_prev);
                if restart == 1
                    fprintf('        %d', doubling_j);
                end
                fprintf('\n');
            end
        end
        
    end
    elapsed = toc;
    
    %% Construct and return results
    x_opt = r(:, end); % Optimal solution
    f_opt = fr(end); % Optimal cost function
    
    Hist.r = r; % Historic value of r_j
    Hist.n = n; % Historic value of n_j
    Hist.fr = fr; % Historic value of f(r_j)
    Hist.gr = gr; % Historic value of ||g(r_j)||_*
    if genHist >= 1
        Hist.fx = fx; % Historic value of f(x_k)
        Hist.fy = fy; % Historic value of f(y_k)
        Hist.gx = gx; % Historic value of ||g(x_k)||_*
        Hist.gy = gy; % Historic value of ||g(y_k)||_*
    else
        Hist.fx = []; Hist.fy = []; Hist.gx = []; Hist.gy = [];
    end
    if genHist >= 2
        Hist.x = x; % Historic value of x
        Hist.y = y; % Historic value of y
        Hist.t = t; % historic value of t
    else
        Hist. x = []; Hist.y = []; Hist.t = [];
    end
    Hist.f_m = f_m;
    Hist.x_m = x_m;
    Hist.k_m = k_m;
    Hist.c = c; % Value of c at each of the restart points
    Hist.k_res = k_res; % Value of k at each one of the restart points
    Hist.k_min = k_min_j; % Value of k_min at each of the outer iterations
    Hist.k = k; % Total number of inner iterations performed
    Hist.j = j; % Number of outer iterations performed
    Hist.x_opt = x_opt; % Value of the first return argument
    Hist.f_opt = f_opt; % Value of the second return argument
    Hist.e_flag = e_flag; % Value of the exit flag
    Hist.f_flag = f_flag; % Value of the exit flag of the last call to the inner FISTA algorithm
    Hist.time = elapsed; % Time taken by the algorithm
    Hist.param.k_max = k_max; % Value of k_max
    Hist.param.tol = tol; % Exit tolerance
    Hist.param.restart = restart; % Restart method used
    Hist.param.k_min = k_min_initial; % Initial value of k_min (at start of iteration 1)
    Hist.para.f_opt = f_opt_ext; % Optimal value of the objective function provided by the user (or the dafault value of 0)
    if restart == 1
        Hist.doubling_j = doubling_j;
    end
    
    % Verbose info to user
    if verbose >= 3
        fprintf('\n');
    end
    if e_flag <= 0
        if verbose >= 1
            if e_flag == -1
                fprintf('\tWARNING: Algorithm did not converge within the allowed number of iterations %d\n', k_max);
            elseif e_flag == -2
                fprintf('\tWARNING: Algorithm did not converge within the allowed number of restarts %d\n', j_max);
            end
            fprintf('\tRFISTA exit flag is %d. Exit flag of the last executed FISTA algorithm is %d\n', e_flag, f_flag);
            fprintf('\tElapsed time between start and end of algorithm: %fs\n', Hist.time);
        end
    else
        if verbose >=2
            fprintf('\tAlgorithm exited successfully with exit flag %d and the following results\n', e_flag);
            fprintf('\t\tNumber of outer loop iterations (j): %d\n', j);
            fprintf('\t\tNumber of inner loop iterations (k): %d (of a maximum of %d)\n', k, k_max);
            fprintf('\t\tf(x_opt) = %g\t\tg(x_opt) = %g\n', f_opt, gr(end));
            fprintf('\t\tElapsed time between start and end of algorithm: %fs\n', Hist.time);
        end
    end
    
    
end

%% Subfunctions

% FISTA algorithm
function [r, n, e_flag, H] = subFISTA(z, f_eval, grad, plus_op, R, k_total, k_min, k_max, restart, tol, genHist, verbose, f_opt, t0, exit_at_tol)
    % Initialization
    done = false; % Flag that indicates the end of the algorithm
    exit_cond = false; % Flag that is set to 1 if the exit condition is satisfied at iteration k
    k = 0; % Counter for number of iterations
    %yk = zeros(length(z), 1); % Value of y at iteration k. (No need to initialize it)
    yk1 = plus_op(z, grad, R); % Value of y at iteration k-1. This is the value of y_0
    xk = zeros(length(z), 1); % Value of x at iteration k
    xk1 = yk1; % Value of x at iteration k-1. This is the value of x_0
    %tk = 0; % Value of t at iteration k. (No need to initialize it)
    tk1 = t0; % Value of t at iteration k-1. We set the value of t_0 = 1
    % TODO: pre-allocate memory for fx, gy, fy and gx
    fx(1) = f_eval(xk1); % Value of the cost function at points x_k. This is stored independently of genHist value
    %gy(1) = g_norm_handler(yk1, plus_op, grad, R); % Value of ||g(y_k)||_* at points y_k. This is stored independently of genHist value
    c = NaN; % Value of k at which the exit condition of the restart method would have exited had k>=k_min not existed
    f_m = f_eval(z); % Minimum value of f found during the iterations
    x_m = z; % x for which the minimum value of f found during the iterations is attained
    k_m = 0; % Value of k where f_m is attained last
    
    if genHist >= 1
        fy(1) = fx(1); % Value of the cost function at points y_k. We have that y_0 = x_0
        gx(1) = g_dual_norm_handler(xk1, plus_op, grad, R); % Historic for the values of ||g(x_k)||_*
    end
    if genHist >= 2
        x(:,1) = xk1; % Historic of x
        y(:,1) = yk1; % Historic of y. We have that y_0 = x_0
        t(1) = tk1; % Historic of t
    end
    
    % Algorithm
    while ~done
        % Increase k and k_total counters
        k = k + 1; k_total = k_total + 1;
        
        % Main steps
        xk = plus_op(yk1, grad, R);
        tk = 0.5*(1 + sqrt(1+4*tk1^2));
        yk = xk + (tk1 - 1)*(xk - xk1)/tk;
        
        % Compute variables
        fx(k+1) = f_eval(xk);
        gy(k) = g_dual_norm(yk1, xk, R); % Compute || g(y_{k-1})||_*. It is computed here to take advantage of x_k = y_{k-1}+
        if genHist >= 1
            fy(k+1) = f_eval(yk); % Value of the cost function at points y_k
            gx(k+1) = g_dual_norm_handler(xk, plus_op, grad, R); % Historic for the values of ||g(x_k)||_*
        end
        if genHist >= 2
            x(:,k+1) = xk; % Historic of x
            y(:,k+1) = yk; % Historic of y
            t(:,k+1) = tk; % Historic of t
        end
        
        % Update f_min and x_min, is applicable
        if fx(k+1) <= f_m
            f_m = fx(k+1);
            x_m = xk;
            k_m = k;
        end
        
        % Compute exit condition
        if ~isa(restart, 'function_handle')
            % 0: No restart
            if restart == 0
                exit_cond = gy(k) <= tol;
                
            % 1: Linear Convergent scheme (T.Alamo et.al., 'Restart FISTA with Global Linear Convergence')
            elseif restart == 1
                m = floor(k/2)+1;
                exit_cond = (fx(m+1) - fx(k+1) <= (fx(1) - fx(m+1))/exp(1)) && (fx(k+1) <= fx(1)); 
                
            % 2: Functional scheme (B. O'Donoghue and E. Candes, 'Adaptive restart for accelerated gradient schemes')
            elseif restart == 2
                exit_cond = fx(k+1) > fx(k); 

            % 3: Gradient scheme (B. O'Donoghue and E. Candes, 'Adaptive restart for accelerated gradient schemes')
            elseif restart == 3
                exit_cond = comp_g(yk1, xk, R)'*(xk1 - xk) < 0;
                
            % 4: Optimal objective function (I. Necoara et.al., 'Linear convergence of first order methods for non-strongly convex optimization')
            elseif restart == 4
                exit_cond = fx(k+1) - f_opt <= (1/exp(1)^2)*(fx(1) - f_opt);
            
            % 5: Restart when the gradient mapping of y at iteration k g(yk) has decreased 1/e wrt the initial one (k=0)
            elseif restart == 5
                gyk = g_dual_norm_handler(yk, plus_op, grad, R);
                exit_cond = gyk <= gy(1)/exp(1);

            % 6: Restart when the gradient mapping of y at iteration k g(yk) has decreased k/(k+1) wrt the initial one (k=0)
            elseif restart == 6
                gyk = g_dual_norm_handler(yk, plus_op, grad, R);
                exit_cond = gyk <= k*gy(1)/(k+1);
                
            % 7: Delay. Similar to LCR. For the TAC article (Algorithm 3)
            elseif restart == 7
                m = floor(k/2);
                exit_cond = (fx(m+1) - fx(k+1) <= (fx(1) - fx(m+1))/3);
            
            end
            
        else % Use an external restart function
            exit_cond = restart(x, y, fx, fy, gx, gy, k, f_eval, grad, plus_op, @g_dual_norm_handler, R, k_min, k_max);
        end
        
        % Check if algorithm has finished
        % Algorithm exited by having satisfied the exit condition and k>=k_min
        if exit_cond
            if isnan(c) && k>=1; c = k; end % Save the first iteration in which the exit condition is satisfied. TODO: Do I need k>= ?
            if k >= k_min
                done = true;
                e_flag = 1;
            end
        end
        % Exited because a value of g(y_k) was found that was lower than the given tolerance (tol)
        if exit_at_tol && gy(k) <= tol% && ~done
            done = true;
            e_flag = 2;
        end
        % Algorithm exits due to exceeding the maximum allowed number of iteration
        if k_total >= k_max% && ~done
            done = true;
            e_flag = -1;
        end
        
        % Update variables
        xk1 = xk;
        yk1 = yk;
        tk1 = tk;
        
    end
    
    % Compute last value of gy
    gy(k+1) = g_dual_norm_handler(yk, plus_op, grad, R);
    
    % Create Historic and returned variables
    H.c = c;
    H.fx = fx;
    H.gy = gy;
    H.k_total = k_total;
    H.f_m = f_m;
    H.x_m = x_m;
    H.k_m = k_m;
    if genHist >= 1
        H.fy = fy;
        H.gx = gx;
    end
    if genHist >= 2
        H.x = x;
        H.y = y;
        H.t = t;
    end
    % TODO: Do I need to explicitly define these return variables? Can't i just use xk and k?
    r = xk; % Return variable
    n = k; % Return variable
    
end

% MFISTA algorithm
function [r, n, e_flag, H] = subMFISTA(x_0, f_eval, grad, plus_op, R, k_total, k_min, k_max, restart, tol, genHist, verbose, f_opt, t0, exit_at_tol)
    % Initialization
    done = false; % Flag that indicates the end of the algorithm
    exit_cond = false; % Flag that is set to 1 if the exit condition is satisfied at iteration k
    k = 0; % Counter for number of iterations
    yk1 = x_0; % Value of y at iteration k-1. This is the value of y_0
    %xk = zeros(length(z), 1); % Value of x at iteration k
    xk1 = yk1; % Value of x at iteration k-1. This is the value of x_0
    %tk = 0; % Value of t at iteration k. (No need to initialize it)
    tk1 = t0; % Value of t at iteration k-1. We set the value of t_0 = 1
    fx(1) = f_eval(xk1); % Value of the cost function at points x_k. This is stored independently of genHist value
    fz(1) = fx(1);
    %gy(1) = g_norm_handler(yk1, plus_op, grad, R); % Value of ||g(y_k)||_* at points y_k. This is stored independently of genHist value
    c = NaN; % Value of k at which the exit condition of the restart method would have exited had k>=k_min not existed
    %f_m = f_eval(z0); % Minimum value of f found during the iterations
    %x_m = z0; % x for which the minimum value of f found during the iterations is attained
    %k_m = 0; % Value of k where f_m is attained last
    if genHist >= 1
        fy(1) = fx(1); % Value of the cost function at points y_k. We have that y_0 = x_0
        gx(1) = g_dual_norm_handler(xk1, plus_op, grad, R); % Historic for the values of ||g(x_k)||_*
    end
    if genHist >= 2
        x(:,1) = xk1; % Historic of x
        y(:,1) = yk1; % Historic of y. We have that y_0 = x_0
        z(:,1) = xk1; % Historic of z
        t(1) = tk1; % Historic of t
    end
    
    % Algorithm
    while ~done
        % Increase k and k_total counters
        k = k + 1; k_total = k_total + 1;
    
        % Main steps
        zk = plus_op(yk1, grad, R);
        tk = 0.5*(1 + sqrt(1+4*tk1^2));
        fz(k+1) = f_eval(zk);
        if fz(k+1) <= fx(k)
            xk = zk;
        else
            xk = xk1;
        end
        yk = xk + tk1*(zk - xk)/tk + (tk1 - 1)*(xk - xk1)/tk;
        
        % Compute variables
        fx(k+1) = f_eval(xk);
        % TODO: Avoid computing another composite gradient mapping operator for || g(y_{k-1})||_*
        gy(k) = g_dual_norm_handler(yk1, plus_op, grad, R); % Compute || g(y_{k-1})||_*. It is computed here to take advantage of x_k = y_{k-1}+
        if genHist >= 1
            fy(k+1) = f_eval(yk); % Value of the cost function at points y_k
            gx(k+1) = g_dual_norm_handler(xk, plus_op, grad, R); % Historic for the values of ||g(x_k)||_*
        end
        if genHist >= 2
            x(:,k+1) = xk; % Historic of x
            y(:,k+1) = yk; % Historic of y
            z(:,k+1) = zk; % Historic of z
            t(:,k+1) = tk; % Historic of t
        end
    
        % Compute exit condition
        if ~isa(restart, 'function_handle')
            % 0: No restart
            if restart == 0
                exit_cond = gy(k) <= tol;
                
            % 1: Linear Convergent scheme (T.Alamo et.al., 'Restart FISTA with Global Linear Convergence')
            elseif restart == 1
                m = floor(k/2)+1;
                exit_cond = (fx(m+1) - fx(k+1) <= (fx(1) - fx(m+1))/exp(1)) && (fx(k+1) <= fx(1)); 
                
            % 2: Functional scheme (B. O'Donoghue and E. Candes, 'Adaptive restart for accelerated gradient schemes')
            elseif restart == 2
                exit_cond = fx(k+1) > fx(k);

            % 3: Gradient scheme (B. O'Donoghue and E. Candes, 'Adaptive restart for accelerated gradient schemes')
            elseif restart == 3
                exit_cond = comp_g(yk1, zk, R)'*(xk1 - xk) < 0;
                
            % 4: Optimal objective function (I. Necoara et.al., 'Linear convergence of first order methods for non-strongly convex optimization')
            elseif restart == 4
                exit_cond = fx(k+1) - f_opt <= (1/exp(1)^2)*(fx(1) - f_opt);

            % 5: Restart when the gradient mapping of y at iteration k g(yk) has decreased 1/e wrt the initial one (k=0)
            elseif restart == 5
                gyk = g_dual_norm_handler(yk, plus_op, grad, R);
                exit_cond = gyk <= gy(1)/exp(1);

            % 6: Restart when the gradient mapping of y at iteration k g(yk) has decreased k/(k+1) wrt the initial one (k=0)
            elseif restart == 6
                gyk = g_dual_norm_handler(yk, plus_op, grad, R);
                exit_cond = gyk <= k*gy(1)/(k+1);
                
            % 7: Delay. Similar to LCR. For the TAC article
            elseif restart == 7
                m = floor(k/2);
                exit_cond = (fx(m+1) - fx(k+1) <= (fx(1) - fx(m+1))/3);
            
            end
            
        else % Use an external restart function
            exit_cond = restart(x, y, fx, fy, gx, gy, k, f_eval, grad, plus_op, @g_dual_norm_handler, R, k_min, k_max);
        end
        
        % Check if algorithm has finished
        % Algorithm exited by having satisfied the exit condition and k>=k_min
        if exit_cond
            if isnan(c) && k>=1; c = k; end % Save the first iteration in which the exit condition is satisfied. TODO: Do I need k>= ?
            if k >= k_min
                done = true;
                e_flag = 1;
            end
        end
        % Exited because a value of g(y_k) was found that was lower than the given tolerance (tol)
        if exit_at_tol && gy(k) <= tol% && ~done
            done = true;
            e_flag = 2;
        end
        % Algorithm exits due to exceeding the maximum allowed number of iteration
        if k_total >= k_max% && ~done
            done = true;
            e_flag = -1;
        end
        
        % Update variables
        xk1 = xk;
        yk1 = yk;
        tk1 = tk;
        
    end
    
    % Compute last value of gy
    gy(k+1) = g_dual_norm_handler(yk, plus_op, grad, R);
    
    % Create Historic and returned variables
    H.c = c;
    H.fx = fx;
    H.gy = gy;
    H.k_total = k_total;
    H.f_m = fx(end);
    H.x_m = xk;
    H.k_m = k;
    if genHist >= 1
        H.fy = fy;
        H.gx = gx;
    end
    if genHist >= 2
        H.x = x;
        H.y = y;
        H.z = z;
        H.t = t;
    end
    % TODO: Do I need to explicitly define these return variables? Can't i just use xk and k?
    r = xk; % Return variable
    n = k; % Return variable
    
end
   
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
