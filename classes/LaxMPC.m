%% Lax MPC - MPC formulation with no terminal constraint and a terminal cost ||x - xr||_P
% This class extends the ssMPC class (which itself extends the QP class)
%
% A detailed description of the MPCT formulation can be found in equation (9) of:
% P. Krupa, D. Limon, and T. Alamo, �Implementation of model predictive
% control in programmable logic controllers,� IEEE Transactions on
% Control Systems Technology, 2020.
%
% The decision variables of this controller are:
%   z = (u_0, x_1, u_1, ..., x_{N-1}, u_{N-1}, x_N}
% 
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

% Author: Pablo Krupa (pkrupa@us.es)
% Creation: 2020/05/07
% Last update: 2020/09/04
% 
% Changelog: 
%   v0.1 (2020/05/07): Initial commit version
%   v0.1 (2020/09/04): Added documentation
%   v0.2 (2020/09/14): Added constraints C z <= d for bounds on C x + D u
%

classdef LaxMPC < ssMPC
    
    %% PROPERTIES
    properties
        P % Terminal cost
    end
    
    properties (Hidden=true)
        C_x0 % For updating vector d with the new x0 (see method x0_update)
        d_x0 % For updating vector d with the new x0 (see method x0_update)
        Walpha % Alpha matrices of the Cholesky factorization of the Weq matrix
        Wbeta % Beta matrices of the Cholesky factorization of the Weq matrix (diagonal elemnts inverted)
    end
    
    methods
    %% CONSTRUCTOR
    function self = LaxMPC(model, Q, R, P, N, varargin) % (model, Q, R, P, N, x0, xr, ur)
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
        addRequired(par, 'P', @(x) isnumeric(x) && (size(x,1)==size(x,2)));
        addRequired(par, 'N', @(x) isnumeric(x) && isscalar(x) && (x>0));
        % Optional
        addOptional(par, 'x0', def_x0, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(par, 'xr', def_xr, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(par, 'ur', def_ur, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        % Name-value parameters
        addParameter(par, 'solver', def_solver);
        % Parse
        parse(par, model, Q, R, P, N, varargin{:})
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
        [H, q] = LaxMPC.compute_H(par.Results.Q, par.Results.R, par.Results.P, par.Results.N);
            % Equality constraints
        [A, b] = LaxMPC.compute_eq(par.Results.model, par.Results.N);
            % Inequality constraints
        [C, d, LB, UB, C_x0, d_x0] = LaxMPC.compute_ineq(par.Results.model, par.Results.N);
        
        % Construct MPC
        self@ssMPC(model, nx, nu, ny, H, q, A, b, C, d, LB, UB, par.Results.solver,...
            par.Results.Q, par.Results.R, N, N, x0, xr, ur, 'computeInverse', true, 'computeWeq', true);
        self.P = par.Results.P;
        self.C_x0 = C_x0;
        self.d_x0 = d_x0;
        self = self.ref_update;
        self = self.x0_update;
        self = self.Wc_update;
        self.startup_MPC = false;
 
    end
    
    function u = extract_control(self, z)
            u = z(1:self.m);
    end
    
    %% GETTERS and SETTERS
    
    % P
    function self = set.P(self, value)
        self.P = value;
        if ~self.startup_MPC
            H = compute_H(self);
            self.H = H;
            self = self.ref_update;
        end
    end
    
    end
    
    %% PROTECTED METHODS
    methods (Access = protected)
        
        function self = ref_update(self)
            self.q = -[self.R*self.ur; kron(ones(self.N-1, 1), [self.Q*self.xr; self.R*self.ur]); self.P*self.xr];
        end
        
        function self = x0_update(self)
            
            % Update vector b
            self.b(1:self.n) = -self.model.A*self.x0;
            
            % Update vector d
            if ~isempty(self.d_x0)
                self.d(1:length(self.d_x0)) = self.d_x0 + self.C_x0*self.x0;
            end
            
        end
        
        function self =  Wc_update(self)
            [self.Walpha,  self.Wbeta] = LaxMPC.sparse_Wc(self);
        end
        
    end
    
    %% STATIC METHODS
    methods (Static)
        
        function [H, q] = compute_H(Q, R, P, N)
            % Obtain variables from function inputs
            if isobject(Q) % This case assumes that the argument Q is in fact an instance of LaxMPC
                R = Q.R;
                P = Q.P;
                N = Q.N;
                Q = Q.Q;
            end
            % Compute the Hessian H and q vectors for the nominal MPC formulation
            nx = size(Q, 1);
            nu = size(R, 1);
            H = blkdiag(R, kron(eye(N-1), [Q zeros(nx,nu); zeros(nu,nx) R]), P);
            q = zeros(N*(nx+nu),1);
        end
        
        function [Az, b] = compute_eq(model, N)
            % Obtain variables from function inputs
            if isobject(model) && ~isa(model, 'ssModel') % This case assumes that the argument Q is in fact an instance of LaxMPC
                    N = model.N;
                model = model.model;
            end
            A = model.A;
            if isa(model, 'ssModel')
                B = model.Bu;
            else
                B = model.B;
            end
            % Compute the A and b equality matrix and vector, respectively.
            nx = size(A, 1);
            nu = size(B, 2);
            Az = kron(eye(N-1), [model.A model.Bu]); % Diagonal of the matrix
            j = 0;
            for i=1:nx:nx*N-nx % Insert matrices -I in Az
                j = j+1;
                Az(i:i+nx-1,((j-1)*(nx+nu)+(nu+nx+1)):((j-1)*(nx+nu)+(nx+nu+1)+nx-1)) = -eye(nx);
            end
            Az = [model.Bu -eye(nx) zeros(nx, size(Az, 2) - nx); zeros(size(Az, 1), nu) Az]; % Initial condition
            b = zeros(N*nx, 1);
        end
        
        function [C, d, LB, UB, C_x0, d_x0] = compute_ineq(model, N)
            % Obtain variables from function inputs
            if isobject(model) && ~isa(model, 'ssModel') % This case assumes that the argument Q is in fact an instance of LaxMPC
                N = model.N;
                model = model.model;
            end
            A = model.A;
            if isa(model, 'ssModel')
                B = model.Bu;
                Cc = model.Cc;
                Dc = model.Dcu;
                LBy_indx = find(~isinf(model.LBy));
                UBy_indx = find(~isinf(model.UBy));
                LBy = model.LBy(LBy_indx);
                UBy = model.UBy(UBy_indx);
            else
                B = model.B;
                if isfield(model, 'C')
                    Cc = model.C;
                else
                    Cc = [];
                end
                if isfield(model, 'D')
                    Dc = model.D;
                else
                    Dc = [];
                end
                if isfield(model, 'LBy')
                    LBy_indx = find(~isinf(model.LBy));
                    LBy = model.LBy(LBy_indx);
                else
                    LBy_indx = [];
                    LBy = [];
                end
                if isfield(model, 'UBy')
                    UBy_indx = find(~isinf(model.UBy));
                    UBy = model.UBy(UBy_indx);
                else
                    UBy_indx = [];
                    UBy = [];
                end
            end
            % Dimensions
            nx = size(A, 1);
            nu = size(B, 2);
            ny = size(Cc, 1);
            % Make Cc = 0 or Dc = 0 under certain conditions
            if ~isempty(Cc) && isempty(Dc)
                Dc = zeros(ny, nu);
            elseif isempty(Cc) && ~isempty(Dc)
                Cc = zeros(ny, nu);
            end
            
            % Compute the inequality constraint matrix and vector for C z <= d
            if (~isempty(Cc) || ~isempty(Dc))
                % Compute upper bounds
                if ~all(isinf(model.UBy))                   
                    Ccu = Cc(UBy_indx, :);
                    Dcu = Dc(UBy_indx, :);
                    Cu = kron(eye(N-1), [Ccu Dcu]); % Prediction horizon
                    Cu = [zeros(size(Cu, 1), nu), Cu]; % u_0
                    Cu = [Cu, zeros(size(Cu, 1), nx)]; % x_N
                    du = kron(ones(N-1, 1) , UBy);
                else
                    Ccu = [];
                    Dcu = [];
                    Cu = [];
                    du = [];
                end
                % Compute lower bounds
                if ~all(isinf(model.UBy))
                    Ccl = -Cc(LBy_indx, :);
                    Dcl = -Dc(LBy_indx, :);
                    Cl = kron(eye(N-1), [Ccl Dcl]); % Prediction horizon                    
                    Cl = [zeros(size(Cl, 1), nu), Cl]; % u_0             
                    Cl = [Cl, zeros(size(Cl, 1), nx)]; % x_N
                    dl = kron(ones(N-1, 1) , -LBy);
                else
                    Ccl = [];
                    Dcl = [];
                    Cl = [];
                    dl = [];
                end
                
                % Compute matrix C and vector d
                C = [Cu; Cl];
                d = [du; dl];
                if ~all(all(Dcu == 0))
                    C_aux = [Dcu zeros(nu, size(C, 2) - nu)];
                    d_x0 = UBy;
                    C_x0 = -Ccu;
                else
                    C_aux = [];
                    d_x0 = [];
                    C_x0 = [];
                end
                if ~all(all(Dcl == 0))
                    C_aux = [C_aux; Dcl zeros(nu, size(C, 2) - nu)];
                    d_x0 = [d_x0; -LBy];
                    C_x0 = [C_x0; -Ccl];
                end
                C = [C_aux; C];
                d = [d_x0; d];
                
            else
                C = [];
                d = [];
                C_x0 = [];
                d_x0 = [];
            end
            
            % Compute lower and upper bounds of the decision variables
            LB = [model.LBu; kron(ones(N-1,1), [model.LBx; model.LBu]); model.LBx]; 
            UB = [model.UBu; kron(ones(N-1,1), [model.UBx; model.UBu]); model.UBx];
        end
        
        function [alpha, beta] = sparse_Wc(self)
            
            % Initialize
            alpha = zeros(self.n, self.n, self.N-1);
            beta = zeros(self.n, self.n, self.N);
            Wc = chol(self.Weq);
            
            % Compute alpha matrices
            for i = 1:self.N-1
                alpha(:,:,i)  = Wc((i-1)*self.n+(1:self.n),i*self.n+(1:self.n));
            end
            
            % Compute beta matrices
            for i = 1:self.N
                beta(:,:,i) = Wc((i-1)*self.n+(1:self.n),(i-1)*self.n+(1:self.n));
                % Invert elements in the diagonal
                for j = 1:self.n 
                    beta(j,j,i) = 1/beta(j,j,i);
                end
            end
            
        end
        
    end
        
end

%% TODOS:
% TODO: Add update call to Walpha and Wbeta if any of its ingredients is changed (matrices, model or N). Extend parent methods.
