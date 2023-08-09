%% compute_ellip_inner_poly
% 
% [Poly, H] = compute_ellip_inner_poly(P)
% 
% Computes an inner polyhedral approximation of the
% ellipsoid {x : x'*P*x <= 1}
% 
% Returns the polyhedral in two forms:
%   - Poly: Polyhedral object from the MPT3 toolbox
%   - H: Matrix which describes the polyhedron as
%        || H*x ||_inf <= 1
% 
% Procedure obtained from Lemma 2 in:
% T. Alama, et.al., "Estimation of the domain of attraction for saturated
% discrete-time systems", in International Journal of Systems Science, 2006.
%
% This function requires a working installation of the MPT3 toolbox.
%
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
%

function [Poly, H] = compute_ellip_inner_poly(P)

    % Compute eigenvalues and eigenvectors of P
    [V, D] = eig(P);
    n = size(P, 1);
    
    % Compute the matrix of the polyhedron
    H = zeros(n, n);
    for i = 1:n
        H(i, :) = sqrt(n*D(i, i))*V(:, i)';
    end
    
     % Instantiate Polyhedron object
    Poly = Polyhedron([H; -H], [ones(n, 1); ones(n, 1)]);

end
