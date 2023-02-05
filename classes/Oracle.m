%% Oracle - A class for generating function oracles
%
% 

classdef Oracle < handle

    %% PROPERTIES
    properties
        param
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
        self.f_eval_hand_ = par.Results.f_eval;
        self.g_eval_hand_ = par.Results.g_eval;

        % Try to set other properties
        if isempty(self.m)
            self.update_m_();
            % if ~isempty(self.f_eval_hand_)
            %     point = cell(1, self.dim);
            %     [point{:}] = self.rand_point_();
            %     self.m = length(self.f(point{:}));
            % end
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

