%% FISTA_simpleQP - FISTA algorithm for solving simple QP problems of a particular form
% 
% This solver finds the solution of QP problems of the form
% 
%   min 1/2 z'*H*z + q'*z
%   s.t. A*z = b
%        LB <= z <= UB,
% 
% where H is a positive definite diagonal matrix, LB < UB and A is full row rank.
% 
% It uses the FISTA algorithm, where the QP problem is solved using its dual formulation.
% For a detailed explanation of the solver, see:
% "Implementation of model predictive control in programmable logic controllers", by
% P. Krupa, D. Limon and T. Alamo, in Transaction on Control Systems Technology, 2020.
% 
% [z_opt, f_opt, e_flag, Hist] = FISTA(varargin)
% 
% OUTPUTS:
%   - z_opt: Optimal decision variables of the QP problem
%   - f_opt: Optimal cost function of the QP problem
%   - e_flag: exit flag of thee algorithm. 1: optimal found. 0: maximum iterations attained.
%   - Hist: Structure containing additional data.
%
% INPUTS:
% There are two ways of calling the function.
% 1) FISTA_simpleQP(myQP): In this case, myQP is an instance of class QP, or extends it.
%                          Note thatthe constraints C*z<=d will be ignored if there are any.
% 2) FISTA_simpleQP(H, q, A, b, LB, UB)
%
% OPTIONAL INPUTS:
%
%

function [z_opt, f_opt, e_flag, Hist] = FISTA_simpleQP(varargin)
    %% Default values
    def_lambda0 = []; % Set default value of the dual variables to 0 by setting this to empty
    def_LSE_solver = []; % Default solver for the Linear System of Equations
    def_k_min = 0; % Default value for the minimum number of iterations
    def_k_max = 2000; % Default value for the maximum number of iterations
    def_tol = 1e-4; % Default value of the exit tolerance
    def_useSparse = false; % Default value of use_sparse. Determines if matrices are created as sparse
    def_genHist = 0; % Default amout of data generated and returned in Hist
    def_verbose = 1; % Default ampunt of information displayed to user on screen
    def_t0 = 1; % Default initial value of t_k
    def_initialize_Hist = true; % Default value og initialize_Hist
    def_print_counter_step = 10; % Number of iterations between verbose to user (only if verbose >= 3)
    
    %% Obtain QP variables from arguments
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'FISTA_simpleQP';
    % Required
    addRequired(par, 'H', @(x) isnumeric(x)); % Hessian H
    % Optional
    addOptional(par, 'q', [], @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    addOptional(par, 'A', [], @(x) isnumeric(x) || isempty(x));
    addOptional(par, 'b', [], @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    addOptional(par, 'LB', [], @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    addOptional(par, 'UB', [], @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    % Name-value parameters
    addParameter(par, 'lambda_0', def_lambda0, @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
    addParameter(par, 'LSE_solver', def_LSE_solver, @(x) isa(x, 'function_handle') || isempty(x));
    addParameter(par, 'use_sparse', def_useSparse, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'k_min', def_k_min, @(x) isnumeric(x) && (x>=0) && x==floor(x));
    addParameter(par, 'k_max', def_k_max, @(x) isnumeric(x) && (x>0) && x==floor(x));
    addParameter(par, 'tol', def_tol, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'genHist', def_genHist, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 'verbose', def_verbose, @(x) isnumeric(x) && (x>=0));
    addParameter(par, 't_0', def_t0, @(x) isnumeric(x));
    addParameter(par, 'initialize_Hist', def_initialize_Hist, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'print_counter_step', def_print_counter_step, @(x) isnumeric(x) && (x>=0));
    
    % Option 1: FISTA_simpleQP(myQP), where myQP is an instance of class QP
    if isa(varargin{1}, 'QP')
        H = varargin{1}.H;
        q = varargin{1}.q;
        A = varargin{1}.A;
        b = varargin{1}.b;
        LB = varargin{1}.LB;
        UB = varargin{1}.UB;
        % Parse
        parse(par, H, q, A, b, LB, UB, varargin{2:end});
        
    % Option 2: FISTA_simpleQP(H, q, A, b, LB, UB)
    else
        % Parse
        parse(varargin{:});
        % Rename main variables
        H = par.Results.H;
        q = par.Results.q;
        A = par.Results.A;
        b = par.Results.b;
        LB = par.Results.LB;
        UB = par.Results.UB;
    end
    
    % Rename the rest of variables
    lambda_0 = par.Results.lambda_0;
    LSE_solver = par.Results.LSE_solver;
    use_sparse = par.Results.use_sparse;
    k_min = par.Results.k_min;
    k_max = par.Results.k_max;
    tol = par.Results.tol;
    genHist = par.Results.genHist;
    verbose = par.Results.verbose;
    t_0 = par.Results.t_0;
    initialize_Hist = par.Results.initialize_Hist;
    print_counter_step = par.Results.print_counter_step;
    
    % Check arguments
    if isempty(lambda_0); lambda_0 = zeros(size(A, 1), 1); end
    if size(lambda_0, 2)>1; lambda_0 = lambda_0'; end % Make lambda_0 a column vector
    if k_min >= k_max
        if verbose>=3
            warning('FISTA_simpleQP:InputWarning', 'Value of k_min >= k_max. Algorithm will exit at k_max iterations allways');
        end
    end
    if tol==0; warning('FISTA_simpleQP:InputWarning', 'Value of tol=0. Algorithm will never converge to this tolerance'); end
    if genHist > 2; genHist = 2; end
    if verbose > 3; verbose = 3; end
    if t_0 <= 0; warning('FISTA_simpleQP:InputWarning', 'Value of t_0 must be > 0'); end
    
    % Get or compute W matrix
    Hinv = 1./diag(H);
    if isa(varargin{1}, 'QP')
        W = varargin{1}.Weq;
    else
        W = A*diag(Hinv)*A';
    end
    if isempty(LSE_solver)
        W = inv(W);
    %  Make  W sparse
    elseif use_sparse
        W = sparse(W);
    end 
    % Make A sparse
    if use_sparse
        A =  sparse(A);
    end
    
    %% Initialization
    done = false; % Flag that indicates the end of the algorithm
    dim = size(H, 1); % Number of decision variables
    k = 0; % Counter for the total number of iterations performed
    z =  zeros(dim, 1); % Decision variables of the QP problem
    lambda = lambda_0; % Dual variable
    mu_1 = zeros(size(A, 1), 1); % Value of the linearization point mu at the previous iteration
    tk1 = t_0; % Value of t at iteration k-1
    print_counter = 0; % Counter used to determine when to show info to user
    print_title = 0; % Conter used to determine when to reprint titles line
    
    
    %% Algorithm
    while ~done
        k = k + 1;
        
        %% Compute z_k
        
        % Compute A'*lambda
        zk_aux =  A'*lambda;
        
        % Compute zk
        for i = 1:dim
            z(i) = min(  max( (zk_aux(i) - q(i))*Hinv(i), LB(i)), UB(i));
        end
        
        %% Compute D_lambda
        
        %  Compute residual A*z - b
        res = A*z - b;
        
        % Solve linear system of equations
        if isempty(LSE_solver)
            D_lambda = -W*res;
        else
            D_lambda = LSE_solver(W, res);
        end
        
        %% Compute lambda
        mu = lambda + D_lambda;
        tk = 0.5*(1 + sqrt(1+4*tk1^2));
        lambda = mu + (tk1 - 1)*(mu - mu_1)/tk;
        
        % Compute exit condition
        norm_res = norm(res, 2);
        if (norm_res <= tol) && (k >= k_min)
            done = true;
            e_flag = 1;
        elseif k >= k_max
            done = true;
            e_flag = -1;
        end
        
        % Update variables
        tk1 = tk;
        mu_1 = mu;
        
    end
    
    %% Construct historic and return variables
    z_opt = z;
    f_opt = 0.5*z_opt'*H*z_opt + q'*z_opt;
    
    Hist.k = k;
    
end

%% TODOS:
%   - TODO: Add generation of Historic variables
%   - TODO: Add verbose
%   - TODO: Add function handler for computing A*z - b
