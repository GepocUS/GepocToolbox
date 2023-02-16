%% Oracle - A class for generating function oracles
%
% This class implements an oracle for a given function
% f: R^n -> R^m.
% The class provides methods to evaluate the function
% and its derivative.
% The objective of this class is to provide a common
% interface for working with functions in the GepocToolbox.
% In particular for its use in optimization algorithms.
%
% orac = Oracle(param, n, 'f_eval', f, 'g_eval', g)
%
% Constructor arguments:
%   - param (required): Structure containing the parameters
%                       that define the function.
%   - n (required): Dimension of the input space of the function.
%                   It can be a scalar (function f: R^n -> R^m) or
%                   An array n = [n_1, n_2, ... n_p], in which case
%                   f: R^n_1 x R^n_2 x ... x R^n_p -> R^m.
%   - m (optional): Dimension of the output space of the function.
%                   If not provided, it will be computed internally.
%   - 'f_eval': Function handler to a function that evaluates the
%               oracle function value. The function handler must
%               have the following format:
%               f = @(x, p) f_eval(x, p)
%               where x is the point where the function will be 
%               evaluated and p is the structure of parameters
%               that define the function (constructor's param input).
%   - 'g_eval': Function handler to a function that evaluates the
%               oracle function gradient. The function handler must
%               have the following format:
%               g = @(x, p) g_eval(x, p)
%               where x is the point where the gradient will be 
%               evaluated and p is the structure of parameters
%               that define the function (constructor's param input).
%   - 'L': Lipschitz constant of the function.
%   - 'mu': Strong convexity parameter of the function.
%
% Methods:
%   - f(x_1,.. x_p): Evaluates the function value at (x_1,.. x_p)
%   - g(x_1,.. x_p): Evaluates the gradient at (x_1,.. x_p)
%   - set_f(f): Resets the f_eval function handler to the new handler f.
%   - set_g(g): Resets the g_eval gradient handler to the new handler g.
%   - get_f_hand(): Returns the function evaluation function handler.
%   - get_g_hand(): Returns the gradient evaluation function handler.
%
% Simple example: Oracle for function f(x) = 0.5*x'*Q*x + c'*x
%   n = 4; % Declare dimension
%   param = struct('Q', eye(n), 'c', ones(n,1)); % Structure with parameters
%   f = @(x, p) 0.5*x'*p.Q*x + p.c'*x; % Handler for function value evaluation
%   g = @(x, p) p.Q*x + c.p; % Handler for gradient evalation
%   L = max(eig(param.Q)); % Compute Lipschitz constant
%   qpOracle = Oracle(param, n, 'f_eval', f, 'g_eval', g, 'L', L);
%   qpOracle.f(rand(n, 1)) % Print f(x) for a random x
% 
% For additional help on a method call: help Oracle.method_name
% 
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

classdef Oracle < handle

    %% PROPERTIES
    properties
        param % Structure containing the parameters of the function
        n % Dimension of each input variable
        m % Output dimension of the function
        dim % Number of input variables of the function
        L % Lipschitz smoothness constant
        mu % Strong convexity constant
    end
    properties (Hidden = true, SetAccess=protected, GetAccess=protected)
        f_eval_hand_ % Function handler to evaluate the function value. See Oracle.f()
        g_eval_hand_ % Function hendler to evaluate the gradient value. See Oracle.g()
    end

    methods

    %% CONSTRUCTOR
    
    function self = Oracle(varargin)

        % Default values
        def_n = [];
        def_m = [];
        def_f_eval = [];
        def_g_eval = [];
        def_L = [];
        def_mu = [];

        % Parser
        par = inputParser;
        par.CaseSensitive = false;
        par.FunctionName = 'Oracle.constructor()';

        % Required
        addRequired(par, 'param', @(x) isa(x, 'struct') || isempty(x));
        addRequired(par, 'n', @(x) isnumeric(x) || isempty(x));
        % Optional
        addOptional(par, 'm', def_m, @(x) isnumeric(x) || isempty(x));
        % Name-value parameters
        addParameter(par, 'f_eval', def_f_eval, @(x) isa(x, 'function_handle') || isempty(x));
        addParameter(par, 'g_eval', def_g_eval, @(x) isa(x, 'function_handle') || isempty(x));
        addParameter(par, 'L', def_L, @(x) (isnumeric(x) && (x >= 0)) || isempty(x));
        addParameter(par, 'mu', def_mu, @(x) (isnumeric(x) && (x >= 0)) || isempty(x));

        % Parse
        parse(par, varargin{:})

        % Set properties
        self.param = par.Results.param;
        self.n = par.Results.n;
        self.dim = length(self.n);
        self.m = par.Results.m;
        self.L = par.Results.L;
        self.mu = par.Results.mu;
        self.f_eval_hand_ = par.Results.f_eval;
        self.g_eval_hand_ = par.Results.g_eval;

        % Try to set other properties
        if isempty(self.m)
            self.update_m_(); % Automatically set m
        end

    end

    %% GETTERS and SETTERS

    function set_f(self, value)
        % Oracle.set_f() - Sets the function handler for function evaluation
        % See also: Oracle.f()
        if isa(value, 'function_handle')
            self.f_eval_hand_ = value;
            self.update_m_();
        else
            throw(MException("Oracle:NotaHandler","Oracle.set_f() must receive a function handler"));
        end
    end

    function set_g(self, value)
        % Oracle.set_g() - Sets the function handler for gradient evaluation
        % See also: Oracle.g()
        if isa(value, 'function_handle')
            self.g_eval_hand_ = value;
        else
            throw(MException("Oracle:NotaHandler","Oracle.set_g() must receive a function handler"));
        end
    end

    function f_hand = get_f_hand(self)
        % Oracle.f_hand() - Returns function handler for the function evaluation
        f_hand = @(varargin) self.f_eval_hand_(varargin{:}, self.param);
    end

    function g_hand = get_g_hand(self)
        % Oracle.f_hand() - Returns function handler for the gradient evaluation
        g_hand = @(varargin) self.g_eval_hand_(varargin{:}, self.param);
    end

    %% PUBLIC METHODS

    function f_val = f(self, varargin)
        % Oracle.f() - Evaluate function value
        if ~isempty(self.f_eval_hand_)
            f_val = self.f_eval_hand_(varargin{:}, self.param);
        else
            f_val = [];
        end
    end

    function g_val = g(self, varargin)
        % Oracle.f() - Evaluate gradient value
        if ~isempty(self.g_eval_hand_)
            g_val = self.g_eval_hand_(varargin{:}, self.param);
        else
            g_val = [];
        end
    end

    end

    %% PROTECTED METHODS

    methods (Access = protected)

    function [ varargout ] = rand_point_(self)
        for i = 1:self.dim
            varargout{i} = rand(self.n(i), 1);
        end
    end

    function update_m_(self)
        if ~isempty(self.f_eval_hand_)
            point = cell(1, self.dim);
            [point{:}] = self.rand_point_();
            self.m = length(self.f(point{:}));
        end
    end

    end

end

