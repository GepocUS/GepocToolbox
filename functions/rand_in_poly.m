%% rand_in_poly - random point generator inside a polyhedron
%
% x = rand_in_poly(P) returns a random point x inside the polyhedron P,
% which is an instance of the MPT3 toolbox class Polyhedron.
% 
% x = rand_in_poly(P, 'opt_name', opt_value, ...) is used to set options.
%
% A random point is generated using a simple hit-and-run sampler starting
% from an initial point x0 and iterating for num_iters iterations.
% The hit-and-run sampler converges to a uniform distribution sampler in
% the polyhedron as num_iters increases.
%
% This class requires the use of the MPT3 toolbox: https://www.mpt3.org/
%
% INPUTS:
%   - P: Instance of the Polyhedron class.
% 
% OPTIONS (as name-value arguments):
%   - 'x0': Initial point of the hit-and-run algorithm. (def: zeros)
%   - 'on_bound': If true, x lies on the bound of P. (def: false)
%   - 'num_iters': Number of iterations of the hit-and-run algorithm.
%   - 'max_exit': If true, an additional exit criterion that looks at the
%                 norm of the current value of x is used. Useful if P is
%                 unbounded, to avoid x becoming too large. (def: false)
%   - 'seed': Reset Matlab's rng to the given seed number.
%   - 'max_norm_penalty' and 'max_norm': Used to determine when the exit
%      condition used if max_exit=true is satisfied. Only modify if you
%      know what you are doing.
%
% OUTPUT:
%   - x: random point in the polyhedron P.
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function x = rand_in_poly(P, varargin)
    
    % Dimensions of the A matrix
    [n, m] = size(P.A);
    
    % Default values
    def_x0 = [];
    def_on_bound = false;
    def_num_iters = (n+m)*10;
    def_seed = [];
    def_max_exit = false;
    def_max_norm_penalty = 1e2;
    def_max_norm = [];
    
    % Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'rand_in_poly';
    
    % Name-value parameters
    addParameter(par, 'x0', def_x0, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
    addParameter(par, 'on_bound', def_on_bound, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'num_iters', def_num_iters, @(x) isnumeric(x) && (x>0) && x==floor(x));
    addParameter(par, 'seed', def_seed, @(x) (isnumeric(x) && (x>0)) || isempty(x));
    addParameter(par, 'max_exit', def_max_exit, @(x) islogical(x) || x==1 || x==0);
    addParameter(par, 'max_norm_penalty', def_max_norm_penalty, @(x) (isnumeric(x) && (x>0)));
    addParameter(par, 'max_norm', def_max_norm, @(x) (isnumeric(x) && (x>0)) || isempty(x));
    
    % Parse
    parse(par, varargin{:})
    
    % Get initial point in the polyhedron
    if ~isempty(par.Results.seed)
        rng(par.Results.seed);
    end
    if isempty(par.Results.x0)
        cheby = P.chebyCenter(); % Compute Chebyshev center
        if cheby.r == Inf
            error("Polyhedron must be bounded");
        end
        x = cheby.x; 
    else
        x = par.Results.x0;
    end
    if par.Results.max_exit
        if isempty(par.Results.max_norm)
            max_norm = par.Results.max_norm_penalty*norm(P.A*ones(m, 1) - P.b, 2);
        else
            max_norm = par.Results.max_norm;
        end
    end
    
    % TODO: Delete (for debugging a two-dimensional example)
%     figure(1); clf(1);
%     P.plot('wire', true, 'linewidth', 2);
%     hold on;
%     plot(x(1), x(2), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
    
    % Hit-and-run-like algorithm
    for i=1:par.Results.num_iters
        
        % Pick random direction
        d = randn(m, 1);
        d = d/norm(d, 2);
        
        % Find maximum and minimum movement in the direction
        z = P.A*d;
        c = (P.b - P.A*x)./z;
        z_neg = z<0;
        tmin = max(c(z_neg));
        tmax = min(c(~z_neg));
        
        % Check for recesion direction
        if isempty(tmax) && isempty(tmin)
            error("Unbounded movement in all directions");
        end
        
        % If only one of the direction is of recession we try to stay close to the bounds
        if isempty(tmax)
            tmax = -tmin;
        elseif isempty(tmin)
            tmin = -tmax;
        end
        
        % Move to a random point in the segment
        step = rand(1);
        x = x + ( tmin + (tmax-tmin) * step )*d;
        
        if par.Results.max_exit && norm(x, 2) > max_norm
            break;
        end
        
        % TODO: Delete (for debugging a two-dimensional example)
%         plot(x(1), x(2), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
        
    end
    
    % Take x on the bound of the set
    if par.Results.on_bound
       x = x - ( tmin + (tmax-tmin) * step )*d; % Reverse last step
       x = x + tmax*d; % Take the step to reach a point on the bound
    end
    
    % TODO: Delete (for debugging a two-dimensional example)
%     plot(x(1), x(2), 'mo', 'MarkerSize', 3, 'MarkerFaceColor', 'm');
    
end
