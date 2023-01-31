%% OscMass_gen_ss - Generate continuous-time model of the oscillating masse system
% 
% This function generates the continuous-time state space model
% of p consecutive masses connected by springs. The state of the
% system is given by x = (x1, x2, .. xp, v1, v2, ... vp), where
% xi is the relative position of mass i, and vi is its velocity.
% The input is given by u = (f1, f2, ... fm), where each fi is
% a force applied to one of the masses.
% The first and last mass are connected via springs to a solid
% an immovable wall.
% 
% INPUTS:
%   - M: Vector of dimension p. Each component of the vector is
%        the mass of each mass, from left to right.
%   - K: Vector of dimension p+1. Each component is the srping
%        constant of each spring, from left to right.
%   - F: Vector of booleans of dimension p. Each component 
%        indicates if the associated mass has an external force
%        acting on it.
%
% OUTPUTS:
%   - sys: Instance of the ss class. Continuous-time state space model.
% 
% For an example of a system of this kind see:
%   M. Kögel, R. Findeisen, "A fast gradient method for embedded
%   linear predictive control", in IFAC World Congress, 2011.
% 
% See also: OscMass_benchmark
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function sys = OscMass_gen_ss(M, K, F)
    
    %% Check arguments
    p = length(M);
    
    %% Generate matrix A
    
    % We start by generating the part corresponding to the velocities
    Av = zeros(p);
    Av(1, 1:2) = [-(K(1) + K(2)), K(2)]; % First row
    Av(p, p-1:p) = [K(p), -(K(p) + K(p+1))]; % Last row
    % Fill al other rows
    for i = 2:p-1
        Av(i, i-1:i+1) = [K(i), -(K(i) + K(i+1)), K(i+1)];
    end
    
    % Divide by the masses
    for i = 1:p
        Av(i, :) = Av(i, :)/M(i);
    end
    
    % Generate A
    A = [zeros(p), eye(p); Av, zeros(p)];
    
    %% Generate matrix B
    B = [zeros(p); diag(1./M)];
    B = B(:, F == true); % Only keep the columns with external forces
    
    %% Generate system
    sys = ss(A, B, eye(2*p), zeros(2*p, length(find(F))));

end
