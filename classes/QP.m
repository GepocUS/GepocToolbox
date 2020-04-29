%% QP - Class for Quadratic Programming problems
%
% Defines a QP problem of the form
% 
%   min 1/2 z'*H*z + q'*z
%     s.t. A*z = b
%          C*z <= d
%          LB <= z <= UB
% 
% Properties:
% - H, q, A, b, C, d, LB, UB: Ingredients of the QP problem
% - dim: Number of decision variables (dimension of H)
% - n_eq: Numer of equality constraints
% - n_ineq: Number of inequality constraints
% - solver: String that indicates the solver used by method QP.solve
% 
% Useful hidden properties
% - options: variable (generally a structure) containing the options passed to the QP solver
% - Hinv: quick way to access the inverse of H. It is only computed the first time it is called
% - isDiag: Boolean indicating if the Hessian is diagonal
% - isPosDef: Boolean indicating if the Hessian is strictly positive definite
%
% Useful methods
% - QP.eval(z): evaluates the functional at point z
% - QP.gradient(z): Returnes the gradient of the functional at point z
% - QP.isConsistent: Returnes 1 if the dimensions of the ingredients are congruent. Returns 0 otherwise.
% 
% Constructor: myQP = QP(H, q, A, b, C, d, LB, UB)
%   - Required arguments: H, q
%   - Optional arguments: A, b, C, d, LB, UB. They default to []
%   - Parameter/value arguments:
%       - solver: string indicating the solver. Defaults to 'quadprog'
%       - computeInverse: Boolean. If 1, Hinv is computed and stored in the constructor method
% 
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

% Author: Pablo Krupa (pkrupa@us.es)
% Creation: 2020/04/29
% Last update: 2020/04/29
% 
% Changelog: 
%   v0.1 (2020/04/29): Initial commit version. Only basic attributes and methods included.
%

classdef QP
    %% PROPERTIES
    properties
        H % Hessian
        q % Independent vector
        A % Matrix for equality constraints
        b % Vector of equality constraints
        C % Matrix for inequality constraints
        d % Vector of inequality constraints
        LB % Lower bounds for decision variables
        UB % Upper bounds for decision variables
        dim % Number of decision variables (Dimension of the Hessian)
        n_eq % Number of equality constraints (rows of A)
        n_ineq % Number of inequality constraints (rows of C)
        solver % String containing the solver called by method 'solve'. Defaults to 'quadprog' % TODO: Add possibility to provide external function handler
    end
    properties (Dependent, Hidden=true)
        isDiag % Boolean that indicates if the Hessian is diagonal
        isPosDef % Boolean that indicates if the QP problem is strictly positive definite
    end
    properties (Hidden=true)
        computeInverse % Determines if the inverse of the Hessian is computed on construction. Defaults to 0.
        options % Options for solver
        debug % Debug mode (adds verbose). Defaults to 0
        Hinv % Inverve of the Hessian (if it is near singlular it is set to empty)
    end
    properties (Access = protected, Hidden=true)
        startup % Determines if the object if being created. It avoind updates while constructor runs (is set to 0 after the constructor finishes)
    end
   
   methods
   %% CONSTRUCTOR
        function self = QP(H, q, varargin) % myQP = QP(H, q, A, b, C, d, LB, UB, varargin)
            % Default values
            def_A = [];
            def_b = [];
            def_C = [];
            def_d = [];
            def_LB = [];
            def_UB = [];
            def_solver = 'quadprog';
            def_debug = false;
            def_computeInverse = false;
            % Parser
            p = inputParser;
            p.CaseSensitive = false;
            p.FunctionName = 'QP_Constructor';
            % Required
            addRequired(p, 'H', @(x) isnumeric(x) && (size(x,1)==size(x,2)));
            addRequired(p, 'q', @(x) isnumeric(x) && (min(size(x))==1));
            % Optional
            addOptional(p,'A', def_A, @(x) isnumeric(x) || isempty(x));
            addOptional(p,'b', def_b, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
            addOptional(p,'C', def_C,@(x) isnumeric(x) || isempty(x));
            addOptional(p,'d', def_d, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
            addOptional(p,'LB', def_LB, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
            addOptional(p,'UB', def_UB, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
            % Name-value parameters
            addParameter(p, 'solver', def_solver);
            addParameter(p, 'debug', def_debug, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'computeInverse', def_computeInverse, @(x) islogical(x) || x==1 || x==0);
            % Parse
            parse(p, H, q, varargin{:})
            % Make vectors column vectors
            if size(p.Results.q,2)>1; p.Results.q=p.Results.q'; end % Make f a column vector
            if size(p.Results.b,2)>1; p.Results.b=p.Results.b'; end % Make b a column vector
            if size(p.Results.d,2)>1; p.Results.d=p.Results.d'; end % Make b a column vector
            if size(p.Results.LB,2)>1; p.Results.LB=p.Results.LB'; end % Make LB a column vector
            if size(p.Results.UB,2)>1; p.Results.UB=p.Results.UB'; end % Make UB a column vector

            % Construct
            self.debug = p.Results.debug;
            self.startup = true; % Startup being initiated
            self.H = p.Results.H;
            self.q = p.Results.q;
            self.A = p.Results.A;
            self.b = p.Results.b;
            self.C = p.Results.C;
            self.d = p.Results.d;
            self.LB = p.Results.LB;
            self.UB = p.Results.UB;
            self.solver = p.Results.solver;
            self.computeInverse = p.Results.computeInverse;
            if strcmp(self.solver, 'quadprog')
                self.options = optimoptions('quadprog', 'Display', 'off');
            else
                self.options = {};
            end
            if self.computeInverse; self = self.updateHinv; end
            self = self.updateDim;
            self = self.updateNeq;
            self = self.updateNineq;
            if self.debug; isConsistent(self); end
        end
       
   %% METHODS
   
        function varargout = solve(self, varargin)
            % Obtain the solution z of the QP problem using the solver method set in property QP.solver
            % Returns the optimal z, optimal value f, exit_flag of the solver and a structure Hist with additional information
            % If a function handler to a solver is used, then the method returnes whatever the solver returnes.

            if isa(self.solver, 'function_handle')
                % Use external function handler. It must accept the arguments as shown here.
                varargout{:} = self.solver(self.H, self.q, self.C, self.d, self.A, self.b, self.LB, self.UB, varargin);
            else

                % quadprog
                if strcmp(self.solver, 'quadprog')
                    if isempty(varargin)
                        vars = {[], self.options};
                    else
                        vars = varargin;
                    end
                    [z, f, exit_flag, output, lambda] = quadprog(self.H, self.q, self.C, self.d, self.A, self.b, self.LB, self.UB, vars{:});
                    Hist = output; % Additional information returned by quadprog
                    Hist.lamdba = lambda ; % Optimal dual variables for equality and inequality constraints
                    Hist.k = output.iterations; % Number of iterations

                % Solver not available
                else
                    warning('QP:SolverNotRecognized', 'Solver not recognized. Returning empty arrays\n')
                    z = [];
                    k = [];
                    f = [];
                    exit_flag = [];
                    Hist = [];
                end

                % Function outputs
                varargout{1} = z; % Optimal solution
                varargout{2} = f; % Optimal cost
                varargout{3} = exit_flag; % Exit flag returned by the solver
                varargout{4} = Hist; % Additional information returned by the solver

            end
        end

        function f = eval(self, z)
            % Returns the evauation of the primal QP problem at point z
            f = zeros(1, size(z, 2));
            for i=1:size(z, 2)
                f(i) = 0.5*z(:,i)'*self.H*z(:,i) + self.q'*z(:,i);
            end
        end
       
        function g = gradient(self, z)
            % Returns the gradient of the primal QP problem at point z
            g = zeros(self.dim, size(z, 2));
            for i=1:size(z,2)
                g(:,i) = self.H*z(:,i) + self.q;
            end
        end

        function [r, message] = isConsistent(self)
            % Returns 1 if all dimensions are consistent. 0 if not. Errors detected are printed in the message string
            r = 1;
            message = '';
            if size(self.H, 1) ~= size(self.q, 1)
                r = 0;
                message = [message 'Dimensions of H and q not consistent\n'];
            end     
            if size(self.A, 1) ~= size(self.b, 1)
                r = 0;
                message = [message 'Dimensions of A and b not consistent\n'];
            end
            if size(self.C, 1) ~= size(self.d, 1)
                r = 0;
                message = [message 'Dimensions of C and d not consistent\n'];
            end
            if size(self.H, 1) ~= size(self.A, 2) && ~isempty(self.A)
                r = 0;
                message = [message 'Dimensions of H and A not consistent\n'];
            end
            if size(self.H, 1) ~= size(self.C, 2) && ~isempty(self.C)
                r = 0;
                message = [message 'Dimensions of H and C not consistent\n'];
            end
            if size(self.LB, 1) ~= size(self.UB, 1) && ~isempty(self.LB) && ~isempty(self.UB)
                r = 0;
                message = [message 'Dimensions of LB and UB not consistent\n'];
            end
            if size(self.H, 1) ~= size(self.LB, 1) && ~isempty(self.LB)
                r = 0;
                message = [message 'Dimensions of H and LB not consistent\n'];
            end
            if size(self.H, 1) ~= size(self.UB, 1) && ~isempty(self.UB)
                r = 0;
                message = [message 'Dimensions of H and UB not consistent\n'];
            end
            if self.debug && r == 0
                fprintf('Dimension mismatch in QP. Message errors:\n');
                fprintf(message);
            end
        end
       
   %% SETTERS and GETTERS
   
    function self = set.H(self, value)
        % Update the value of H
        self.H = value;
        if ~self.startup
            self = self.updateHinv;
            self = self.updateDim;
        end
    end
    
    function self = set.q(self, value)
        % Update the value of q
        if size(value,2)>1; value=p.value'; end % Make f a column vector
        self.q = value;
        if ~self.startup
            self = self.updateDim;
        end
    end

    function self = set.A(self, value)
        % Update the value of A
        self.A = value;
        if ~self.startup
            self.updateNeq;
        end
    end

    function self = set.b(self, value)
        % Update the value of b
        if size(value,2)>1; value=p.value'; end % Make f a column vector
        self.b = value;
        if ~self.startup
            self.updateNeq;
        end
    end

    function self = set.C(self, value)
        % Update the value of C
        self.C = value;
        if ~self.startup
            self.updateNineq;
        end
    end

    function self = set.d(self, value)
        % Update the value of d
        if size(value,2)>1; value=p.value'; end % Make f a column vector
        self.d = value;
        if ~self.startup
            self.updateNineq;
        end
    end

    function self = set.LB(self, value)
        % Update the value of LB
        if size(value,2)>1; value=p.value'; end % Make f a column vector
        self.LB = value;
    end

    function self = set.UB(self, value)
        % Update the value of UB
        if size(value,2)>1; value=p.value'; end % Make f a column vector
        self.UB = value;
    end

    function self = set.solver(self, solver_str)
        % Update the value of solver
        if ~ischar(solver_str) && ~isa(solver_str, 'function_handle')
            error('QP:InputError', 'QP.solver must be a string of characters or a function handler')
        end
        self.solver = solver_str;
        if ~strcmp(self.solver, solver_str) && ~self.startup
            self.options = {};
            if self.debug
                fprintf('Set solver: Erased options for solver');
            end
        end
    end
    
    function r = get.isPosDef(self)
        % Getter for dependant variable isPosDef
        r = min(eig(self.H)) > 0;
    end

    function r = get.isDiag(self)
        % Getter for dependant variable isDiag
        r = isdiag(self.H);
    end
    
    function r = get.Hinv(self)
        % Getter for Hinv. If it is empty, compute it first
        if isempty(self.Hinv)
            self = self.updateHinv; 
        end
        r = self.Hinv;
    end
    
   end

 %% PROTECTED METHODS
    methods (Access = protected)

        function self = updateHinv(self)
            % Updates the inverse of the Hessian
            lastwarn('');
            self.Hinv = inv(self.H);
            [~, warnId] = lastwarn;
            if strcmp(warnId, 'MATLAB:singularMatrix') % If we detect that H is singular we throw an error
                if self.debug
                    warning('QP:HsingularMatrix', 'Hessian H is singular to working precision. Cannot compute inverse')
                end
                self.Hinv = [];
            end
        end
       
       function self = updateDim(self)
           % Update the dim variable
           self.dim = size(self.H, 1);
       end
       
       function self = updateNeq(self)
           % Update the n_eq variable
           self.n_eq = size(self.A, 1);
       end
       
       function self = updateNineq(self)
           % Update the n_ineq variable
           self.n_ineq = size(self.C, 1);
        end
       
    end
end
