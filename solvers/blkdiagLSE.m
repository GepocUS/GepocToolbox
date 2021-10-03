%% blkdiagLSE - Linear system of equations solver for block diagonal matrices
% 
% This is a sparse solver is for systems of equations W*z = c, where z \in \R^{N  n},
% c \in \R^{N n}, and  is a block diagonal matrix whose Cholesky factorization
% has the following structure:
% 
%   Wc = (beta_1, alpha_1, 0,       0,   ...         0;
%         0,      beta_2,  alpha_2, 0,   ...         0;
%         0,      ...      ...      ...  ...         0;
%                                        beta_{N-1}, alpha_{N-1};
%         0,     ...       ...      ...  ...         beta_N},
%
%   where the matrices beta_i \in \R^{n x n} are upper triangular and
%   the matrices alpha \in \R^{n x n}.
% 
% OUTPUTS:
%   - z: vector that solver the system of equations.
% 
% INPUTS (there are two ways of calling the function):
% 
% Option 1: z = blkdiagLSE(Beta,  Alpha,  c, is_inverted)
%   - Beta: Matrix of dimensions (n, n, N). Each (:, :, i) matrix corresponds to beta_i.
%   - Alpha: Matrix of dimensions (n, n, N-1). Each (:, :, i) matrix corresponds to alpha_i.
%   - c: vector c of the system of equations to be solved.
%   - is_inverted (optional, defaults to false): If this option is set to true, then
%       the function assumes that the matrices contained in Beta have had their diagonal 
%       elements inverted.
% 
% Option 2: z = blkdiagLSE(W,  c,  n)
%   - W, c: Matrix W and vector c of the system of equations.
%   - n: dimension of the matrrices alpha and beta
%
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function z = blkdiagLSE(varargin)
    %% Default values
    def_is_inverted = false;
    
    %% Parser
    par = inputParser;
    par.CaseSensitive = false;
    par.FunctionName = 'blkdiagLSE';
    
    % Option 1: z = blkdiagLSE(Alpha,  Beta,  c, is_inverted)
    if size(varargin{1}, 3) > 1
        
        % Required
        addRequired(par, 'Beta', @(x) isnumeric(x));
        addRequired(par, 'Alpha', @(x) isnumeric(x));
        addRequired(par, 'c', @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
        % Optional
        addOptional(par, 'is_inverted', def_is_inverted, @(x) islogical(x) || x==1 || x==0);
        
        % Parse
        parse(par,  varargin{:});
        
        % Rename variables
        Alpha = par.Results.Alpha;
        Beta = par.Results.Beta;
        c = par.Results.c;
        is_inverted = par.Results.is_inverted;
        
        % Check arguments
        if ~(size(Alpha, 3) == size(Beta, 3)-1)
            error('There should be one less matrix alpha than beta');
        end
        
        % Get n and N
        n = size(Alpha, 1);
        N = size(Beta, 3);
        
    else
        
        % Required
        addRequired(par, 'W', @(x) isnumeric(x));
        addRequired(par, 'c', @(x) isnumeric(x) && (min(size(x))==1) || isempty(x));
        addRequired(par, 'n', def_verbose, @(x) isnumeric(x) && (x>=0) && x==floor(x) && isfinite(x));
        
        % Parse
        parse(par, varargin{:});
        
        % Rename variables
        W = par.Results.W;
        c =  par.Results.c;
        n = par.Results.n;
        is_inverted = true;
        
        % Compute Cholesky factorization
        Wc = chol(W);
        
        % Compute N
        N = length(c)/n;
        if N ~= floor(N)
            error('length(c)/n is not an integer');
        end
        
        % Compute matrices Alpha and Beta
        Alpha = zeros(self.n, self.n, self.N-1);
        Beta = zeros(self.n, self.n, self.N);

        % Compute Alpha matrices
        for i = 1:N-1
            Alpha(:,:,i)  = Wc((i-1)*n+(1:n),i*n+(1:n));
        end

        % Compute Beta matrices
        for i = 1:N
            Beta(:,:,i) = Wc((i-1)*n+(1:n),(i-1)*n+(1:n));
            % Invert elements in the diagonal
            for j = 1:n
                Beta(j,j,i) = 1/Beta(j,j,i);
            end
        end
        
    end
    
    %% Initialize
    if size(c, 2) == 1
        z = c;
    else
        z = c';
    end
    
    %% Forward substitution
    
    % Compute first n elements
    for j = 1:n 
        for i = 1:j-1
            z(j) = z(j) - Beta(i, j, 1)*z(i);
        end
        if is_inverted
            z(j) = z(j)*Beta(j, j, 1);
        else
            z(j) = z(j)/Beta(j, j, 1);
        end
    end
    
    % Compute the rest of the elements
    for k = 1:N-1
        for j = 1:n
            for i = 1:n
                z(j + n*k) = z(j + n*k) - Alpha(i, j, k)*z(i + n*(k-1));
            end
            for  i = 1:j-1
                z(j + n*k) = z(j + n*k) - Beta(i, j, k+1)*z(i + n*k);
            end
            if is_inverted
                z(j + n*k) = Beta(j, j, k+1)*z(j + n*k);
            else
                z(j + n*k) = Beta(j, j, k+1)/z(j + n*k);
            end
        end
    end
    
    %% Backwards substitution
    
    % Compute last n elements
    for j = n:-1:1
        for i = n:-1:j+1
            z(j + (N-1)*n) = z(j + (N-1)*n) - Beta(j, i, N)*z(i + (N-1)*n);
        end
        if is_inverted
            z(j + (N-1)*n) = z(j + (N-1)*n)*Beta(j, j, N);
        else
            z(j + (N-1)*n) = z(j + (N-1)*n)/Beta(j, j, N);
        end
    end
    
    % Compute the rest of the elements
    for k = N-2:-1:0
        for j = n:-1:1
            for i = n:-1:1
                z(j + n*k) = z(j + n*k) - Alpha(j, i, k+1)*z(i + n*(k+1));
            end
            for i = n:-1:j+1
                z(j + n*k) = z(j + n*k) - Beta(j, i, k+1)*z(i + n*k);
            end
            if is_inverted
                z(j + n*k) = Beta(j, j, k+1)*z(j + n*k);
            else
                z(j + n*k) = Beta(j, j, k+1)/z(j + n*k);
            end
        end
    end
    
end          

