%% pontryagin_poly_ellip - Pontryagin difference between polyhedron and ellipsoid
%
% Computes the pontryagin difference between X-K*E, where X = {x: A*x <= b},
% E = {x: (x - c)'*P*(x - c) <= r^2}, and K is a matrix.
% The function returns the vector bhat such that Xhat = {x: A*x <= bhat} is
% the pontryagin difference.
%
% INPUTS:
%   - A: m-by-n matrix defining the polyhedron X.
%   - b: m-by-1 vector defining the polyhedron X.
%   - P: n_e-by-n_e matrix defining the geometry of the ellipsoid E.
%   - c: n_e-by-1 vector defining the center of the ellipsoid E.
%   - r: positive scalar defining the size of the ellipsoid E.
%   - K: (optional): n-by-n_e matrix. Defaults to the identity matrix.
% 
% OUTPUTS:
%   - bhat: Vector defining the polyhedron Xhat.
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function [bhat] = pontryagin_poly_ellip(A, b, P, c, r, K)
    
    % Get dimensions and check arguments
    [m, n] = size(A);
    n_e = size(P,1);
    if diff(size(P))
       error("P must be a suqre matrix");
    end
      
    if nargin == 5
        K = eye(n, n_e);
        if n ~= n_e
            error("Columns of A must be equal to the dimension of P");
        end
    else
        [mk, nk] = size(K);
        if mk~=n || nk~=n_e
            error("Dimension of K is not correct")
        end
    end
    Ak = A*K;
    
    bhat = zeros(m, 1);
    Pi = inv(P);
    
    for i = 1:m
        bhat(i) = b(i) - r*sqrt(Ak(i,:)*Pi*Ak(i,:)') - Ak(i,:)*c;      
    end

end
