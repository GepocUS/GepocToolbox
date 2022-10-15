%% BaP_gen_ss - Generate continuous-time model of the ball and plate system
%
% This function generates the continuous-time state space model
% of the ball and plate system described in:
% 
% "Harmonic based model predictive control for set-point tracking", P. Krupa,
% D. Limon and T. Alamo, in Transactions on Automatic control, 2021.
%
% INPUTS:
%   - param: Structure containing the parameters of the system.
%       - .m: Mass of the ball.
%       - .r: Radius of the ball.
%       - .g: Gravitational constant.
%       - .I: Mass momentum of inertia of the ball.
% 
% OUTPUTS:
%   - sys: Instance of the ss class. Continuous-time state space model.
%
% See also: BaP_ode, BaP_bechmark
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function sys = BaP_gen_ss(param)

    % Compute the variable that relates the parameters of te system with the dynamics
    M = param.m/(param.m + param.I/param.r^2);

    % Compute A matrix
    A = zeros(8, 8);
        % x(1): Position of ball in axis 1
    A(1, 2) = 1;
        % x(2): Velocity of ball in axis 1
    A(2, 3) = M*param.g;
        % x(3): Angle of axis 1
    A(3, 4) = 1;
        % x(5): Position of ball in axis 2
    A(5, 6) = 1;
        % x(6): Velocity of ball in axis 2
    A(6, 7) = M*param.g;
        % x(7): Angle of axis 2
    A(7, 8) = 1;

    % Compute B matrix
    B = zeros(8, 2);
        % x(4): Velocity of angle of axis 1
    B(4, 1) = 1;
        % x(8): Velocity of angle of axis 2
    B(8, 2) = 1;

    % Compute C matrix
    C = eye(8);

    % Compute D matrix
    D = zeros(8, 2);

    % Construct continuos-time state space model
	sys = ss(A, B, C, D);
    
end
