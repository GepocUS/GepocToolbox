%% optim_RAFOM - Optimal restart algorithm for accelerated first order methods
% 
% [x_opt, f_opt, e_flag, Hist] = optim_RAFOM(alg, f_eval, x_0) computes the 
% minimizer x_opt of the optimization problem,
% 
%   x_opt = arg min f(x)
%              s.t. x in R^n
% 
% where f(x) is an extended real valued proper closed convex function, using
% the optimal restart algorithm for accelerated first order methods (AFOM)
% presented in (see Algorithm 2):
%
% "Restart of accelerated first order methods with linear convergence under
% a quadratic functional growth condition" T. Alamo, P. Krupa, D. Limon, in
% IEEE Transactions on Automatic Control, 2022.
% 
% OUTPUTS:
%  - x_opt: Solution of the optimization problem.
%  - f_opt: Value of F at the optimum, i.e. f_opt = F(x_opt)
%  - e_flag: Integer that indicates the exit condition of the algorithm.
%         1: Algorithm converged successfully to a solution within the given tolerance.
%         2: The inner algorithm found a solution within its own criteria. This exit
%            condition is only applied if the function input exit_at_tol is set to true.
%        -1: Algorithm did not converge before reaching the maximum allowed iterations.
%  - Hist: Structure containing additional data. See documentation for further details.
% 
% INPUTS:
%  - alg: Handler to a function that runs the optimization algorithm which computes
%         the iterates x_k. The function must have the following inputs and outputs:
%           x_n = alg(x_0, n_min, n_max),
%         where x_0 is the initial condition of the algorithm, n number of
%         iterations it will run for, and x_n is the value of x_k after n iterations.
%         This function also supports the following Gepoc standard:
%           [x_n, f_n, e_n, H_n] = alg(x_0, n),     
%         where f_n is the value of f(x_n), e_n is the exit flag of the algorithm (1 if it
%         has found an optimum, any other number if it has not), and H_n is a structure
%         that can contain any assortment of  historics and data. The contents of Hist
%         will be partly generated from H_n is they are available and if genHist is large
%         enough. In order for this to work input variable is_gepoc must be set to true.
%       
%  - f_eval: Handler to a function that returns the value of f at a point x. f(x) = f_eval(x)
%  - x0: Initial value of x (at iteration 0).
%        If no value is provided it defaults to a vector of zeros.
% 
% OPTIONAL INPUTS (must be passed using 'name' 'value' pairs):
%  - tol: Exit tolerance of the algorithm. Must be > 0. Default to 1e-9
%  - alpha_hat: Real number greater than 1. Design parameter of the algorithm. Default is 15.
%  - j_max: Maximum number of calls to alg. Integer > 0. Defaults to 100.
%  - genHist: Determines how much information is generated for output Hist. Can be 0, 1 or 2.
%             Higher values increase the amount of information returned in Hist. Defaults to 1.
%             A value of 2 can increase the computation time of the function significantly.
%  - verbose: Integer that can take value 0, 1, 2 or 3. It indicates how much information is
%             displayed to the user. No information is displayed if set to 0. Defaults to 1.
%  - exit_at_tol: If 1, the algorithm exits as soon as a value of x_k is found that satisfies the
%                 exit condition. If 0, the exit condition is only checked at points r_j.
%                 It defaults to 1.
%  - is_gepoc: Boolean variable that indicates if the alg function handler is Gepoc compatible.
%              Defaults to true.
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [x_opt, f_opt, e_flag, Hist] = optim_RAFOM (alg, f_eval, x_0, varargin)

    %% Default values
    def_tol = 1e-9; % Default value of the exit tolerance
    def_j_max = 100; % Default maximum number of outer iterations, i.e. calls to alg
    def_genHist = 0; % Default amount of data generated for Hist: Only stores fx, gy and basic data
    def_verbose = 1; % Default amount of information displayed: Displays end of algorithm info only 
    def_exit_at_tol = 1; % Default value of exit_at_rol. Determines if the exit flag of alg is taken as a valid successful exit condition
    def_only_alg_exit = 0; % Default value of variable only_alg_exit. If true, the algorithm only exits if alg returns an exit flag equal to 1
    def_k_max = 10000; % Default value for the liberated maximum number of iterations of alg
    
    %% Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'optim_RAFOM';
    % Required
    addRequired(par, 'alg', @(x) isa(x, 'function_handle'));
    addRequired(par, 'f_eval', @(x) isa(x, 'function_handle'));
    addRequired(par, 'x_0',  @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    % Name-value parameters
    addParameter(par, 'tol', def_tol, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'j_max', def_j_max, @(x) isnumeric(x) && (x>0) && x==floor(x));
    addParameter(par, 'genHist', def_genHist, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'verbose', def_verbose, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'exit_at_tol', def_exit_at_tol, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'only_alg_exit', def_only_alg_exit, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'k_max', def_k_max, @(x) isnumeric(x) && (x>0) && x==floor(x));
    
    % Parse
    parse(par, alg, f_eval, x_0, varargin{:})
    % Rename
    alg = par.Results.alg;
    f_eval = par.Results.f_eval;
    x_0 = par.Results.x_0;
    tol = par.Results.tol;
    j_max = par.Results.j_max;
    genHist = par.Results.genHist;
    verbose = par.Results.verbose;
    exit_at_tol = par.Results.exit_at_tol;
    only_alg_exit = par.Results.only_alg_exit;
    k_max = par.Results.k_max;
    
    % Check arguments
    if size(x_0, 2) > 1; x_0 = x_0'; end % Make x_0 a column vector
    if genHist > 2; genHist = 2; end
    if verbose > 3; verbose = 3; end
    if only_alg_exit; exit_at_tol = true; end

    %% Inicialization
    r = zeros(length(x_0), j_max+1);
    r(:,1) = x_0; % Values of x_k at the restart points j
    
    fr = zeros(1, j_max+1);
    fr(1) = f_eval(x_0); % Value of the cost function f(r_j) at points r_j. This is stored independently of genHist value
    
    m = 1; % m_0
    m_1 = 1; % m_{j-1}
    
    done = false; % Flag that indicates the end of the algorithm
    j = 0; % Counter for number of restarts
    k = 0; % Total number of alg iterations
    print_counter = 0; % Used for determining when to print the line showing what each column of the table represents
    
    % Data obtained from alg
    if genHist >= 1
        k_res(1) = 1;
        fx = []; % Historic of f(x_k)
        gx = []; % Historic of ||g(x_k)||_*
        sj = zeros(1, j_max); % Historic of s_j
        nj = zeros(1, j_max); % Historic of n_j
        mj = zeros(1, j_max); % Historic of m_j
        alg_e_flag = zeros(1, j_max+1); % Exit flags returned by alg at each iteration j
        alg_times = zeros(1, j_max+1); % Computation times of alg at each iteration j
    end
    if genHist >= 2
        x = []; % Historic of x_k
    end
    
    %% Initial verbose to user
    if verbose >= 2
        fprintf('\tStarting optim_RAFOM algorithm');
    end
    if verbose >= 3
        fprintf('\tParamenets:\tj_max = %d\ttol = %g\n', j_max, tol);
        if exit_at_tol
            fprintf('\toptim_RAFOM will also exit if exit condition of alg is satisfied.\n');
        end
        fprintf('\n');
        fprintf('\tIter_j  Iter_k    n_j    m_j    f(r_j)       res');
        fprintf('\n');
    end

    %% Algorithm
    start = tic;
    while ~done
        j = j + 1;
        
        % Compute s_j
        if j > 2
            s = sqrt( ( fr(j-1) - fr(j) )/( fr(j-2) - fr(j) ) );
        else
            s = 0;
        end
        
        % Compute n_j
        n = max([m, 4*s*m_1]);
        n_int = ceil(n);

        % Compute r_j
        [r(:,j+1), ~, e_n, H_n] = alg(r(:,j), n_int);
        m = H_n.k;
        k = k + H_n.k;
        if (e_n == -1) && (verbose>=1)
            warning('optim_RAFOM: alg did not converge within the maximum allowed iterations');
        end

        
        % Compute f(r_j)
        fr(j+1) = f_eval(r(:,j+1)); % Store new value of f(r_j)
            % Check that f(r_j) <= f(r_{j-1})
        if (fr(j+1) > fr(j)) && (verbose>=1)
            war_string = ['optim_RAFOM: The value of f(r_j) if greater than f(r_{j-1}) at iteration j=' num2str(j) '.'];
            warning(war_string);
        end
        
        % Generate historics
        if genHist >= 1
            k_res(j+1) = k + 1;
            fx = [fx H_n.fx(1:end-1)]; % Include all except the last
            gx = [gx H_n.gx(1:end-1)]; % Include all except the last
            sj(j) = s;
            nj(j) = n;
            mj(j) = m;
            alg_e_flag(j) = e_n;
            alg_times(j) = H_n.time;
        end
        if genHist >= 2
            x = [x H_n.x(:,1:end-1)]; % Include all except the last
        end
        
        % Exit conditions
            % Exit condition of optim_RAFOM (I only evaluate it once j>=2)
        res = fr(j) - fr(j+1);
        if j > 0
            if (res < tol) && ~only_alg_exit
                done = true;
                e_flag = 1;
            end
        end
            % Solution found by alg
        if e_n == 1 && ~done
            done = true;
            e_flag = 2;
        end
            % Maximum j attained
        if j >= j_max && ~done
            done = true;
            e_flag = -1;
        end
            % Maximum k attained
        if k >= k_max && ~done
            done = true;
            e_flag = -2;
        end
           
        % Print information (Iter_j, Iter_k, n_j, alpha_j)
        if verbose >= 3
            print_counter = print_counter + 1;
            if print_counter == 10
                print_counter = 0;
                fprintf('\tIter_j  Iter_k    n_j    m_j    f(r_j)       res');
                fprintf('\n');
            end
            fprintf('\t%6d  %6d %6d %6d  %2.2e   %2.2e', j, k, n_int, m, fr(j+1), res);
            fprintf('\n');
        end
        
    end
    elapsed = toc(start);
    
    %% Construct and return results
    x_opt = r(:, j+1); % Optimal solution
    f_opt = f_eval(r(:, j+1)); % Optimal cost function
    
    Hist.time = elapsed;
    Hist.r = r(:, 1:j+1);
    Hist.fr = fr(1:j+1);
    Hist.j = j;
    Hist.k = k;
    Hist.n = n;
    Hist.n_int = n_int;
    if genHist >= 1
        Hist.k_res = k_res - 1;
        Hist.fx = [fx H_n.fx(end)]; % Add the last one
        Hist.gx = [gx H_n.gx(end)]; % Add the last one
        Hist.gr = Hist.gx(k_res);
        Hist.sj = sj(1:j);
        Hist.nj = nj(1:j);
        Hist.mj = mj(1:j);
        Hist.alg_e_flag = alg_e_flag;
        Hist.alg_times = alg_times;
    end
    if genHist >= 2
        Hist.x = [x H_n.x(:,end)]; % Add the last one   
    end

    % Verbose info to user
    if verbose >= 3
        fprintf('\n');
    end
     if e_flag <= 0
        if verbose >= 1
            fprintf('\tWARNING: Algorithm did not converge within the allowed number of restarts %d\n', j_max);
            fprintf('\tElapsed time between start and end of algorithm: %fs\n\n', Hist.time);
        end
     elseif verbose>=2
        fprintf('\tAlgorithm exited successfully with exit flag %d and the following results\n', e_flag);
        fprintf('\t\tNumber of outer loop iterations (j): %d\n', j);
        fprintf('\t\tNumber of inner loop iterations (k): %d\n', k);
        fprintf('\t\tf(x_opt) = %g\n', f_opt);
        fprintf('\t\tElapsed time between start and end of algorithm: %fs\n\n', Hist.time);
     end
    
end
