%% EqualityMPC - MPC formulation with no terminal cost and an equality terminal constraint
% This class extends the ssMPC class (which itself extends the QP class)
%
% A detailed description of the MPCT formulation can be found in equation (8) of:
% P. Krupa, D. Limon, and T. Alamo, �Implementation of model predictive
% control in programmable logic controllers,� IEEE Transactions on
% Control Systems Technology, 2020.
% 
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
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
            u = z(1:self.m);
    end
    
    end
    
    %% PROTECTED METHODS
    methods (Access = protected)
       
        function self = ref_update(self)
            self.q = -[self.R*self.ur; kron(ones(self.N-1, 1), [self.Q*self.xr; self.R*self.ur])];
            self.b((self.N-1)*self.n+1:end) = self.xr;
        end
        
        function self = x0_update(self)
            self.b(1:self.n) = -self.model.A*self.x0;
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
            H = blkdiag(R, kron(eye(N-1), blkdiag(Q, R)));
            q = zeros(N*(nx+nu)-nx,1);
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
            Az = kron(eye(N-1), [A B]); % Diagonal of the matrix
            j = 0;
            for i=1:nx:nx*N-nx % Insert matrices -I in Az
                j = j+1;
                Az(i:i+nx-1,((j-1)*(nx+nu)+(nu+nx+1)):((j-1)*(nx+nu)+(nx+nu+1)+nx-1)) = -eye(nx);
            end
            Az = [B -eye(nx) zeros(nx, size(Az, 2) - nx); zeros(size(Az, 1), nu) Az]; % Initial condition
            Az = Az(:,1:end-nx);
            b = zeros(N*nx, 1);
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
            %C = kron(eye(N), [eye(nx) zeros(nx, nu); zeros(nu, nx) eye(nu); -eye(nx) zeros(nx, nu); zeros(nu, nx) -eye(nu)]);
            %d = kron(ones(N,1), [model.UBx; model.UBu; -model.LBx; -model.LBu]);
            C = [];
            d = [];
            LB = [model.LBu; kron(ones(N-1,1), [model.LBx; model.LBu])]; 
            UB = [model.UBu; kron(ones(N-1,1), [model.UBx; model.UBu])];
        end 

    end
    
end

