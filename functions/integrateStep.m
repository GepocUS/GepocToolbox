%% integrateStep - Integrates a system from time t_0 to t_f
%
% x_f = integrateStep(x_0, u, t_0, t_f, h, integrator)
% 
% This function returns the state x_f of a system at time t_f
% by numerical integration from the given state x_0 at time t_0.
% 
% INPUTS:
%   - x_0: State at time t_0
%   - u: System input to be applied during the time interval
%   - t_0: Initial time
%   - t_f: Time at the end of the integration
%   - h: Step time of the integrator
%   - integrator: Function handler to an integrator function: x_{t+h} = integrator(x_t, u, t, h)
% 
% OUTPUTS:
%   - x_f: State at time t_f
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function x_f = integrateStep(x_0, u, t_0, t_f, h, integrator)

    % Initial steps
    t = t_0;
    N = (t_f - t_0)/h;

    % Compute x_f
    x_f = x_0;
    for k = 1:N
        x_f = integrator(x_f, u, t, h);
        t = t + h;
    end

end
