%% EqualityMPC
% Extends the ssMPC class
% This is a simple MPC with temirnal equality constraint and a single prediction horizon
%

% Author: Pablo Krupa (pkrupa@us.es)
% Creation: 2020/05/03
% Last update: 2020/05/03
% 
% Changelog: 
%   v0.1 (2020/05/03): Initial commit version
%

classdef EqualityMPC < ssMPC
    
    methods
    %% CONSTRUCTOR
    function self = EqualityMPC(model, Q, R, N, varargin) % (model, Q, R, N, x0, xr, ur)
        % Detect model class
        if isa(model, 'ssModel')
            is_ssModel = true;
        else
            is_ssModel = false;
        end
        % Extract from model
        if is_ssModel
            nx = model.n_x;
            nu = model.n_u;
            ny = model.n_yr;
        else
            nx = size(model.A, 1);
            nu = size(model.B, 2);
            ny = size(model.C, 1);
        end
        % Defaults
        def_x0 = zeros(nx,1);
        def_xr = zeros(nx,1);
        def_ur = zeros(nu,1);
        def_solver = 'quadprog';
        % Parser
        par = inputParser;
        par.CaseSensitive = false;
        par.FunctionName = 'EqualityMPC_Constructor';
        % Required
        addRequired(par, 'model');
        addRequired(par, 'Q', @(x) isnumeric(x) && (size(x,1)==size(x,2)));
        addRequired(par, 'R', @(x) isnumeric(x) && (size(x,1)==size(x,2)));
        addRequired(par, 'N', @(x) isnumeric(x) && isscalar(x) && (x>0));
        % Optional
        addOptional(par, 'x0', def_x0, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(par, 'xr', def_xr, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(par, 'ur', def_ur, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        % Name-value parameters
        addParameter(par, 'solver', def_solver);
        % Parse
        parse(par, model, Q, R, N, varargin{:})
        % Rename
        x0 = par.Results.x0;
        xr = par.Results.xr;
        ur = par.Results.ur;
        % Defauls
        if isempty(x0); x0 = def_x0; end
        if isempty(xr); xr = def_xr; end
        if isempty(ur); ur = def_ur; end
        % Make vectors column vectors
        if size(x0,2)>1; x0=x0'; end % Make x0 a column vector
        if size(xr,2)>1; xr=xr'; end % Make xr a column vector
        if size(ur,2)>1; ur=ur'; end % Make ur a column vector
        
        % Compute ingredients of QP problem
            % Hessian and q vector
        [H, q] = EqualityMPC.compute_H(par.Results.Q, par.Results.R, par.Results.N);
            % Equality constraints
        [A, b] = EqualityMPC.compute_eq(par.Results.model, par.Results.N);
            % Inequality constraints
        [C, d, LB, UB] = EqualityMPC.compute_ineq(par.Results.model, par.Results.N);
        
        % Construct MPC
        self@ssMPC(model, nx, nu, ny, H, q, A, b, C, d, LB, UB, par.Results.solver,...
            par.Results.Q, par.Results.R, N, N, x0, xr, ur);
        self = self.ref_update;
        self = self.x0_update;
        self.startup_MPC = false;
        
    end
    
    function u = extract_control(self, z)
            u = z(self.n+(1:self.m));
    end
    
    end
    
    %% PROTECTED METHODS
    methods (Access = protected)
       
        function self = ref_update(self)
            self.q = -kron(ones(self.N, 1), [self.Q*self.xr; self.R*self.ur]);
            self.b(self.N*self.n+1:end) = self.xr;
        end
        
        function self = x0_update(self)
            self.b(1:self.n) = self.x0;
        end
        
    end
    
    %% STATIC METHODS
    methods (Static)
        
        function [H, q] = compute_H(Q, R, N)
            if isobject(Q)
                R = Q.R;
                N = Q.N;
                Q = Q.Q;
            end
            % Calculate the Hessian H and q vectors for the nominal MPC formulation
            nx = size(Q, 1);
            nu = size(R, 1);
            H = kron(eye(N), [Q zeros(nx,nu); zeros(nu,nx) R]);
            q = zeros(N*(nx+nu),1);
        end
        
        function [Az, b] = compute_eq(model, N)
            if isobject(model) && ~isa(model, 'ssModel')
                    N = model.N;
                model = model.model;
            end
            A = model.A;
            if isa(model, 'ssModel')
                B = model.Bu;
            else
                B = model.B;
            end
            % Calculate the A and b equality matrix and vector, respectively, for the nominal MPC formulation
            nx = size(A, 1);
            nu = size(B, 2);
            Az = kron(eye(N), [A B]); % Diagonal de la matriz
            j = 0;
            for i=1:nx:nx*N-nx % Inserto las matrices -I en Az
                j = j+1;
                Az(i:i+nx-1,((j-1)*(nx+nu)+(nu+nx+1)):((j-1)*(nx+nu)+(nx+nu+1)+nx-1)) = -eye(nx);
            end
            Az = [eye(nx) zeros(nx, nx*N-nx+nu*N); Az]; % Initial condition
            b = zeros((N+1)*nx, 1);
        end
       
        function [C, d, LB, UB] = compute_ineq(model, N)
            if isobject(model) && ~isa(model, 'ssModel')
                N = model.N;
                model = model.model;
            end
            A = model.A;
            if isa(model, 'ssModel')
                B = model.Bu;
            else
                B = model.B;
            end
            % Calculate inequality matrix and vector C and d, as well as lower and upper bounds LB and UB for the nominal MPC formulation
            nx = size(A, 1);
            nu = size(B, 2);
            C = kron(eye(N), [eye(nx) zeros(nx, nu); zeros(nu, nx) eye(nu); -eye(nx) zeros(nx, nu); zeros(nu, nx) -eye(nu)]);
            d = kron(ones(N,1), [model.UBx; model.UBu; -model.LBx; -model.LBu]);
            LB = kron(ones(N,1), [model.LBx; model.LBu]); 
            UB = kron(ones(N,1), [model.UBx; model.UBu]);
        end 

    end
    
end

%% TODOS:
% TODO: Add documentation
% TODO: Add an example
% TODO: Test and debug
