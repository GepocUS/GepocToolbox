%% compute_PIS_ellip - compute ellipsoidal PIS of the given system
% 
% Consider the system x+ = A*x + B*u + E*w, where w is a disturbance
% that is assumed to be contained in a compact set W.
% 
% This function computes matrices P and K such that the ellipsoid 
% {x: x'*P*x <= 1} is a positive invariant set (PIS) of the closed-loop
% system with the control law u = K*x. Additionally, the ellipsoid will 
% satisfy constraints Cx*x <= dx and Cu*u <= du for all x belonging to it.
% 
% Other conditions can also be imposed, such as a certain contraction
% factor, a reduction of the constraints on u or the imposition of the
% discrete algebraic Riccaty equation.
% Finally, it can be made robust for all disturbances in W.
% 
% The ellipsoid is computed by solving an LMI optimization problem, which
% is constructed using the YALMIP toolbox (https://yalmip.github.io/).
% 
% For information on the LMI problems and general procedure used here, we
% refer the user to Section 4 of:
% I. Alvarado, et.al. "Tractable robust MPC design based on nominal predictions",
% in Journal of Process Control, 2022.
%
% [P, K, found, status] = compute_PIS_ellip(sys, cons, [param], [opt])
% 
% INPUTS:
% - sys: A structure containing the system's information.
%        The fields of the structure must be matrices A, B and E.
% - cons: A structure containing the constraints that must be satisfied
%         by the ellipsoid, i.e., fields Cx, Cu, dx and du.
% - param (optional): Structure containing other information. In particular:
%       - Q, R: Matrices for the imposition of the Riccati equation.
%               Default to []. If [], the Riccati equation is not imposed.
%       - V: Matrix containing the vertices of W in each of its rows.
%            Defaults to []. If [], the ellipsoid does not consider disturbance.
%       - rho: Scalar that modifies the input constraints imposed on the ellipsoid.
%              In particular, we force Cu*K*x <= rho*du. Defaults to 1.
%       - mu: Scalar that imposes a certain contraction factor on the ellipsoid.
%             That is, we impose mu*P <= Ak'*P*Ak, where Ak = A + B*K.
%             Defaults to 0. If 0, the condition is not imposed.
% - opt(optional): A structure containing options. Its fields are:
%       - solver: String that determines the solver used to solve the LMI problem.
%                 Can be "SDPT3" or "sedumi", which correspond to the solvers:
%                 > SDPT3: https://github.com/SQLP/SDPT3
%                 > SeDuMi: https://github.com/SQLP/SeDuMi
%                 Defaults to "SDPT3".
%       - type: String. Selects the type of invariant set. The options are:
%               - 'none' (default): Returns any invariant set.
%               - 'min': Finds the smallest invariant set.
%               - 'max': Finds the largest invariant set.
%       - solver_verbose: Integer that sets the solver's verbose. Default: 0.
%       - verbose: Integer that sets this function's verbose. Default: 0.
%       - lambda_init: Scalar. Initial value of lambda. Default: 0.02.
%       - lambda_max: Maximum value of lambda. Default: 1.
%       - lambda_step: Step taken between values of lambda. Default: 0.02.
%         (We refer the user to the cited article for an explanation of lambda).
%       - optimize: Boolean. If true, then the function makes a search from
%                   lambda_init to lambda_max taking steps of size lambda_step
%                   and returns the best ellipsoid that it finds. If false,
%                   the function returns the LMI problem's solution for the
%                   value of lambda_init. Default: true.
%       - feas_gap: Scalar. Feasibility gap used to allow a certain degree of
%                   non-feasibility in the solver's solution. Default: 1e-14.
%       - normalize: Boolean. If true then the constraints are normalized so 
%                    that dx and du only contain 1s. This can help the solver
%                    by improving the numerical conditioning. Default: true.
% 
% OUTPUTS:
%   - P: Matrix P of the ellipsoid
%   - K: Gain matrix of the local control law
%   - found: Boolen. If true then the function found an ellipsoid satisfying
%            the selected conditions.
%   - status: Structure containing various additional fields. Only of use if
%             you know a lot about the inner workings of this function.
%
% This function requires a working installation of YALMIP (version 20200116
% onward) able to locate either the SDPT3 or the SeDuMi solver.
% 
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [P, K, found, status] = compute_PIS_ellip(sys, cons, param, opt)
    
    %% Default values
    def_opt.type = 'min';
    def_opt.solver = "SDPT3";
    def_opt.solver_verbose = 0;
    def_opt.lambda_init = 0.02;
    def_opt.lambda_step = 0.02;
    def_opt.lambda_max = 1.0;
    def_opt.optimize = true;
    def_opt.feas_gap = 1e-14;
    def_opt.verbose = 0;
    def_opt.normalize = true;
    def_param.Q = [];
    def_param.R = [];
    def_param.V = [];
    def_param.rho = 1;
    def_param.mu = 0;
    
    %% Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'compute_PIS_ellip';
    
    addRequired(par, 'sys');
    addRequired(par, 'cons');
    addOptional(par, 'param', def_param, @(x) isstruct(x));
    addOptional(par, 'opt', def_opt, @(x) isstruct(x));
    
    parse(par, sys, cons, param, opt);
    sys = par.Results.sys;
    cons = par.Results.cons;
    param = add_default_options_to_struct(par.Results.param, def_param);
    opt = add_default_options_to_struct(par.Results.opt, def_opt);

    % Scale constraints (make dx and du be ones)
    if opt.normalize
        Nx = diag(1./abs(cons.dx)); % Scaling for x
        Nu = diag(1./abs(cons.du)); % Scaling for u
        cons.Cx = Nx*cons.Cx;
        cons.Cu = Nu*cons.Cu;
        cons.dx = Nx*cons.dx;
        cons.du = Nu*cons.du;
    end

    % Check type
    if isempty(strfind(['min', 'max', 'none'], opt.type))
        error("opt.type must be 'min', 'max' or 'none'");
    end
    if strcmp(opt.type, 'none'); opt.optimize = false; end
    
    %% Solve problem
    done = false;
    found = false;
    lambda = opt.lambda_init;
    k = 0;
    if strcmp(opt.type, 'min')
        best_val = -inf;
    elseif strcmp(opt.type, 'max')
        best_val = inf;
    else
        best_val = NaN;
    end
    best_P = [];
    best_K = [];  
    best_status = [];
    best_lambda = [];
    
    if opt.verbose > 0
        fprintf("compute_PIS_ellip: Starting computation of ellipsoid. ");
        if strcmp(opt.type, 'min')
            fprintf("Computing minimal set.\n");
        elseif strcmp(opt.type, 'max')
            fprintf("Computing maximal set.\n");
        else
            fprintf("Computing invariant set (neither minimal nor maximal).\n");
        end
        
        fprintf("\tOptions: ");
        if opt.optimize
            fprintf("find best, ");
        else
            fprintf("find first, ");
        end
        if ~isempty(param.V)
            fprintf("robust, ");
        else
            fprintf("non-robust, ");
        end
        if param.mu > 0
            fprintf("contractive, ");
        else
            fprintf("non-contractive, ");
        end
        if ~isempty(param.Q) && ~isempty(param.R)
            fprintf("impose Riccati.\n");
        else
            fprintf("don't impose Riccati.\n");
        end
    end
    
    if opt.optimize
        num_lambda = length(opt.lambda_init:opt.lambda_step:opt.lambda_max) - 1;
        reverseStr = '';
    end

    while ~done
        k = k+1;
        
        if opt.verbose > 0 && opt.optimize
            percentDone = 100 * k / num_lambda;
            msg = sprintf('\tPercent done: %3.1f', percentDone); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        
        [P, K, e_flag, status] = local_compute_RPIS_ellip(sys, cons, param, opt, lambda);
        
        % Determine if we exit
        if ~opt.optimize
            
            if e_flag > 0
                found = true;
                status.k = k;
                status.lambda = lambda;
            else
                lambda = lambda + opt.lambda_step;
            end
            
            done = true;
            
        else
            
            if e_flag > 0
                
                val = trace(P);
                found = true;
                
                if strcmp(opt.type, 'min')
                    
                   if val >= best_val
                        best_P = P;
                        best_K = K;
                        best_status = status;
                        best_lambda = lambda;
                        best_val = val;
                   end
                   
                else
                    
                    if val <= best_val
                        best_P = P;
                        best_K = K;
                        best_status = status;
                        best_lambda = lambda;
                        best_val = val;
                    end

                end
                    
            end
            
            lambda = lambda + opt.lambda_step;
            
            if lambda >= opt.lambda_max
                P = best_P;
                K = best_K;
                status = best_status;
                status.lambda = best_lambda;
                done = true;
                status.k = k;
                if opt.verbose > 0
                    fprintf("\n");
                end
            end
            
        end

    end
    
    % Final verbose to user
    if opt.verbose > 0
        if found
            fprintf("\tFound a suitable ellipsoid for lambda = %2.4f\n\n", status.lambda);
        else
            fprintf("\tDid not find a suitable ellipsoid\n\n");
        end
    end

end

%% This function is the one that actually constructs and solves the LMI problem

function [P, K, e_flag, status] = local_compute_RPIS_ellip(sys, constraint, param, opt, lambda)
    
    %% Rename for convenience
    A = sys.A;
    B = sys.B;
    E = sys.E;
    n = size(A, 1); % Number of states
    m = size(B, 2); % Number of iputs
    Cx = constraint.Cx;
    dx = constraint.dx;
    Cu = constraint.Cu;
    du = constraint.du;
    Q = param.Q;
    R = param.R;
    V = param.V;
    
    %% Initialize YALMIP and declare variables
    yalmip('clear');

    S = sdpvar(n, n); % S = inv(P)
    Y = sdpvar(m, n); % Y = K*inv(P)
    if strcmp(opt.type, 'min'); gamma = sdpvar(1, 1); end
    
    %% Constraints
    cons = []; % List of constraints

    % Impose robustness w.r.t. W
    % Find pair P and K that is Schur stable and for which K is robust
    % A constraint is needed for each w \in vert(W), i.e. for each row of the argument V
    if ~isempty(V)
        for i = 1:size(V, 1)
            cons = [cons, [lambda*S, zeros(n,1), (A*S + B*Y)'; zeros(1,n), (1-lambda), V(i,:)*E'; (A*S + B*Y), E*V(i,:)', S] >=0];
        end
    end
    
    % Impose that P and K satisfy the relaxed discrete Riccati equation
    if ~isempty(Q) && ~isempty(R)
        Qsq = sqrtm(Q);
        Rsq = sqrtm(R);
        cons = [cons, [S, (A*S + B*Y)', (Qsq*S)', (Rsq*Y)';...
                       (A*S + B*Y), S, zeros(n ,n), zeros(n, m);...
                       (Qsq*S), zeros(n, n), eye(n), zeros(n, m);...
                       (Rsq*Y), zeros(m, n), zeros(m, n), eye(m)] >= 0];
    end

    % Impose state constraints
    for i = 1:size(Cx, 1)
        if strcmp(opt.type, 'min')
%             cons = [cons, [gamma*(dx(i))^2, Cx(i, :)*S;...
%                            S'*Cx(i, :)', S] >= 0];
            cons = [cons, Cx(i, :)*S*Cx(i, :)' <= gamma*(dx(i))^2];
        else
            cons = [cons, [(dx(i))^2, Cx(i, :)*S;...
                           S'*Cx(i, :)', S] >= 0];
        end
    end

    % Impose input constraints
    for i = 1:size(Cu, 1)
        cons = [cons, [(param.rho*du(i))^2, Cu(i, :)*Y;...
                       Y'*Cu(i, :)', S] >= 0];
    end

    % Make it mu-contractive
    if param.mu > 0
        cons = [cons, [param.mu*S (A*S + B*Y)'; (A*S + B*Y) S] >= 0];
    end
    
    %% Set solver settings
    if strcmp(opt.solver, 'SDPT3') % SDPT3 solver
        solver_options = sdpsettings('solver', 'SDPT3', 'cachesolver', 1, 'verbose', opt.solver_verbose); % Use SDPT3 solver
    
    elseif strcmp(opt.solver, 'sedumi') % SeDuMi solver
        solver_options = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-10, 'sedumi.cg.restol', 1e-6,...
            'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, 'sedumi.stepdif', 2,...
            'cachesolver', 1, 'verbose', opt.solver_verbose, 'debug', 1);
    end
    
    %% Solve optimization problem and extract values
    
    if strcmp(opt.type, 'min')
        status = optimize(cons, gamma, solver_options);
    elseif strcmp(opt.type, 'max')
        status = optimize(cons, -trace(S), solver_options);
    else
        status = optimize(cons, [], solver_options);
    end
    
    % Extract values
    Sk = value(S);
    Yk = value(Y);
    P = inv(Sk);
    K = (Sk\Yk')';
    
    % Diagnostics
    [primal_feas, dual_feas] = check(cons);
    if min(primal_feas) > 0 && min(dual_feas) > 0
        e_flag = 1; % Feasible solution found
    elseif abs(min(primal_feas)) < opt.feas_gap && abs(min(dual_feas)) < opt.feas_gap
        e_flag = 2; % Sub-feasible solution within the tolerance
    else
        e_flag = 0; % Non feasible solution found
    end
    
    % Return other status information
    if strcmp(opt.type, 'min')
        status.f_opt = value(gamma);
    else
        status.f_opt = -trace(Sk);
    end
    
    status.primal_feas = primal_feas;
    status.dual_feas = dual_feas;
    status.e_flag = e_flag;

end
