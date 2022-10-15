%% BaP_ode - evaluation of the ODE of the ball and plate system
%
% This function evaluates the ODE equations of the ball and plate
% system described in:
% 
% "Harmonic based model predictive control for set-point tracking", P. Krupa,
% D. Limon and T. Alamo, in Transactions on Automatic control, 2021.
%
% INPUTS:
%   - t: Current time. It is not used in the ode, but is included
%        to be compliant with the requirements of Matlab's ode functions.
%   - x: Current state of the system.
%        x = (p_1, dp_1, teta_1, dteta_1, p_2, dp_2, teta_2, dteta_2)
%   - u: Current control input.
%        u = (ddp_1, ddp_2)
%   - param: Structure containing the parameters of the system.
%       - .m: Mass of the ball.
%       - .r: Radius of the ball.
%       - .g: Gravitational constant.
%       - .I: Mass momentum of inertia of the ball.
%
% OUTPUTS:
%   - dx: Time derivative of the system state.
%
% See also: BaP_gen_ss, BaP_bechmark
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function dx = BaP_ode(t, x, u, param)

    % Compute the variable that relates the parameters of the system with the dynamics
    M = param.m/(param.m + param.I/param.r^2);
    
    % Compute the derivative of the system state
    dx = zeros(8, 1);

    dx(1) = x(2);
    dx(2) = M*(x(1)*x(4)^2 + x(5)*x(4)*x(8) + param.g*sin(x(3)));
    dx(3) = x(4);
    dx(4) = u(1);
    dx(5) = x(6);
    dx(6) = M*(x(4)*x(8)^2 + x(1)*x(4)*x(8) + param.g*sin(x(7)));
    dx(7) = x(8);
    dx(8) = u(2);

end
