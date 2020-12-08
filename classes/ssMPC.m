%% ssMPC - Abstract class for state space MPC
% This is the abstract class for all the state space based MPC controller classes
% The non-abstract classes must define their own properties for additional cost function matrices
%
% This class extends the QP class. Therefore, it inherits its properties and methods.
%
% Properties
%   - model: an instance of ssModel class. The model of the system. Used as the prediction model.
%   - Q, R: Cost function matrices for the state and control input, respectively.
%   - Nc, Np: Control horizon and prediction horizon, respectively
%   - x0: Current value of the system state (in incremental units)
%   - x0Eng % Current value of the system state (in engineering units)
%   - xr % Current value of state reference (in incremental units)
%   - xrEng % Current value of state reference (in engineering units)
%   - ur % Current value of input reference (in incremental units)
%   - urEng % Current value of input reference (in engineering units)
%
% Methods
%   - ssMPC.solve(args): solves the MPC. It calls the solver defined by the QP property 'solver'
%       args can be any series of optional agruments accepted by the selected solver.
%       It returnes the following arguments: [u, z, f, exit_flag, Hist] = solve(args)
%           - u: Control action to be applied to the system (in incremental units)
%           - z: Optimal decision variables returned by the solver (if the solver converged)
%           - f: Optimal cost function returned by the solver
%           - exit_flag: exit flag returned by the solver
%           - Hist: Additional information returned by the solver or constructed by the QP class solve method
%   - ssMPC.setRef(self, xr, ur): Sets the reference (xr, ur) to the given values. In incremental units
%   - ssMPC.setRefEng(self, xrEng, urEng): Sets the reference (xr, ur) to the given values. In engineering units
%
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

% Author: Pablo Krupa (pkrupa@us.es)
% Creation: 2020/05/03
% Last update: 2020/09/04
% 
% Changelog: 
%   v0.1 (2020/05/03): Initial commit version
%   v0.1 (2020/09/04): Added documentation
%

classdef (Abstract) ssMPC < QP
    %% PROPERTIES
    properties
        model % Matrices and vectors of the system model. It must be a state space model
        Q % Cost for the system states
        R % Cost for the system inputs
        Nc % Control horizon
        Np % Prediction horizon
        x0 % Current value of system state (in incremental units)
        x0Eng % Current value of the system state (in engineering units)
        xr % Current value of state reference (in incremental units)
        xrEng % Current value of state reference (in engineering units)
        ur % Current value of input reference (in incremental units)
        urEng % Current value of input reference (in engineering units)
    end
    properties (SetAccess = private)
        n % Number of system states
        m % Number of system inputs
        p % Number of tracked system outputs
        N % Short access for the control horizon (used in MPC controllers that only have one horizon)
    end
    properties (SetAccess = protected, Hidden = true)
        startup_MPC % Determines if the object if being created. It avoind updates while constructor runs (is set to 0 after the constructor finishes)
        lockSetter % Boolean used to avoid reciprocal setter loops
    end
    
    methods       
    %% CONSTRUCTOR
    function self = ssMPC(model, nx, nu, ny, H, q, A, b, C, d, LB, UB, solver, Q, R, Nc, Np, x0, xr, ur, varargin)
        self@QP(H, q, A, b, C, d, LB, UB, 'solver', solver, varargin{:});
        self.startup_MPC = true;
        self.lockSetter = false;
        self.model = model;
        self.Q = Q;
        self.R = R;
        self.N = Nc;
        self.Nc = Nc;
        self.Np = Np;
        self.n = nx;
        self.m = nu;
        self.p = ny;
        self.x0 = x0;
        self.xr = xr;
        self.ur = ur;       
    end
    
    %% METHODS
    function [u, z, f, exit_flag, Hist] = solve(self, varargin)
        % Calls the QP solver
        [z, f, exit_flag, Hist] =solve@QP(self, varargin{:});
        u = self.extract_control(z);
    end
    
    function self = setRef(self, xr, ur)
        % Sets the state reference xr and input reference ur to the given values in incremental units. Updates the QP problem
        self.xr = xr;
        self.ur = ur;
    end
    
    function self = setRefEng(self, xrEng, urEng)
        % Sets the state reference xr and input reference ur to the given values in engineeging units. Updates the QP problem
        self.xrEng = xrEng;
        self.urEng = urEng;
    end 
    
    %% SETTERS and GETTERS
    
    % Q
    function self = set.Q(self, value)
        self.Q = value;
        if ~self.startup_MPC
            H = compute_H(self);
            self.H = H;
            self = self.ref_update;
        end
    end
    
    % R
    function self = set.R(self, value)
        self.R = value;
        if ~self.startup_MPC
            H = compute_H(self);
            self.H = H;
            self = self.ref_update;
        end
    end
    
    % N
    function self = set.N(self, value)
        if ~self.lockSetter
            self.lockSetter = true;
            if ~self.startup_MPC
                self.Nc = value;
                self.N = value;
                [H, q, A, b, C, d, LB, UB] = compute_QP(self);
                self.H = H;
                self.q = q;
                self.A = A;
                self.b = b;
                self.C = C;
                self.d = d;
                self.LB = LB;
                self.UB = UB;
            else
                self.N = value;
            end
        end
        self.lockSetter = false;
    end
    
    % Nc
    function self = set.Nc(self, value)
        if ~self.lockSetter
            self.lockSetter = true;
            if ~self.startup_MPC
                self.Nc = value;
                self.N = value;
                [H, q, A, b, C, d, LB, UB] = compute_QP(self);
                self.H = H;
                self.q = q;
                self.A = A;
                self.b = b;
                self.C = C;
                self.d = d;
                self.LB = LB;
                self.UB = UB;
            else
                self.Nc = value;
            end
        end
        self.lockSetter = false;
    end
    
    % Np
    function self = set.Np(self, value)
        if ~self.startup_MPC
            self.Np = value;
            [H, q, A, b, C, d, LB, UB] = compute_QP(self);
            self.H = H;
            self.q = q;
            self.A = A;
            self.b = b;
            self.C = C;
            self.d = d;
            self.LB = LB;
            self.UB = UB;
        else
            self.Np = value;
        end
    end
    
    % x0
    function self = set.x0(self, value)
        if size(value,2)>1; value=value'; end % Make column vector
        self.x0 = value;
        if ~self.lockSetter
            self.lockSetter = true;
            if isa(self.model, 'ssModel')
                self.x0Eng = x2Eng(self.model, value);
            else
                self.x0Eng = value;
            end
            if ~self.startup_MPC
                self = self.x0_update; 
            end
            self.lockSetter = false;
        end
    end
    function self = set.x0Eng(self, value)
        if size(value,2)>1; value=value'; end % Make column vector
        self.x0Eng = value;
        if ~self.lockSetter
            self.lockSetter = true;
            if isa(self.model, 'ssModel')
                self.x0 = x2Inc(self.model, value);
            else
                self.x0 = value;
            end
            if ~self.startup_MPC
                self = self.x0_update;
            end
            self.lockSetter = false;
        end
    end
    
    % xr
    function self = set.xr(self, value)
        if size(value,2)>1; value=value'; end % Make column vector
        self.xr = value;
        if ~self.lockSetter
            self.lockSetter = true;
            if isa(self.model, 'ssModel')
                self.xrEng = x2Eng(self.model, value);
            else
                self.xrEng = value;
            end
            if ~self.startup_MPC
                self = self.ref_update;
            end
            self.lockSetter = false;
        end
    end
    function self = set.xrEng(self, value)
        if size(value,2)>1; value=value'; end % Make column vector
        self.xrEng = value;
        if ~self.lockSetter
            self.lockSetter = true;
            if isa(self.model, 'ssModel')
                self.xr = x2Inc(self.model, value);
            else
                self.xr = value;
            end
            if ~self.startup_MPC
                self = self.ref_update;
            end
            self.lockSetter = false;
        end
    end
    
    % ur
    function self = set.ur(self, value) % Setter. Updates the QP problem
        if size(value,2)>1; value=value'; end % Make column vector
        self.ur = value;
        if ~self.lockSetter
            self.lockSetter = true;
            if isa(self.model, 'ssModel')
                self.urEng = u2Eng(self.model, value);
            else
                self.urEng = value;
            end
            if ~self.startup_MPC
                self = self.ref_update;
            end
            self.lockSetter = false;
        end
    end
    function self = set.urEng(self, value)
        if size(value,2)>1; value=value'; end % Make column vector
        self.urEng = value;
        if ~self.lockSetter
            self.lockSetter = true;
            if isa(self.model, 'ssModel')
                self.ur = u2Inc(self.model, value);
            else
                self.ur = value;
            end
            if ~self.startup_MPC
                self = self.ref_update;
            end
            self.lockSetter = false;
        end
    end
    
    end % End of 'methods'
    
    %% PROTECTED METHODS
    methods (Access = protected)
        
        function [H, q, A, b, C, d, LB, UB] = compute_QP(self)
            % Computes the ingredients of the QP problem
            [H, q] =  compute_H(self);
            [A, b] = compute_eq(self);
            [C, d, LB, UB] = compute_ineq(self);
        end
        
    end

    %% ABSTRACT METHODS
    methods (Abstract, Access = protected)
        
        ref_update(self) % Updates the QP ingredients that depend on the reference
        x0_update(self) % Updates the QP ingredients that depend on the current state
        
    end
    
    methods (Abstract)
        
        extract_control(varargin) % Extracts the control action from the given vector of decision variables 
        
    end
    
    methods (Abstract, Static)
       
        compute_H(varargin) % Compute the H and q ingredients of the QP
        compute_eq(varargin) % Compute the A and b ingredients of the QP (equality constraints)
        compute_ineq(varargin) % Compute the C and d ingredients of the QP (inequality constraints)

    end

end

