%% rk4 - Runge-Kutta 4 integrator
% 
% Performs a Runge-Kutta 4 integration step.
% xk1 = rk4_(f, xk, u, t, h) returns the state at time t+h of the system x_dot = f(x, u, t).
% 
% INPUTS:
%   - f: Function handler that returns the derivative of x, i.e. x_dot = f(x, u, t)
%   - xk: Current state
%   - u: Input to be applied during the integration step
%   - t: Current time
%   - h: Integration step
% 
% OUTPUTS:
%   - xk1: Value of the state at time t+h
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function xk1 = rk4(f, xk, u, t, h)

    k1 = f(xk,      u, t    )*h;
    k2 = f(xk+k1/2, u, t+h/2)*h;
    k3 = f(xk+k2/2, u, t+h/2)*h;
    k4 = f(xk+k3,   u, t+h  )*h;

    xk1 = xk + (k1 + 2*k2 + 2*k3 + k4)/6;

end
