%% Reactors_ode - evaluation of the ODE of the reactors system
%
% This function evaluates the ODE equations of the reactors
% system described in:
% 
% "Implementation of model predictive control in programmable logic controllers",
% P. Krupa, D. Limon and T. Alamo, in Transactions on Control Systems Technology, 2021.
%
% The continuous-time model of the system is obtained using Reactors_gen_ss()
% The operating point of the system is obtained using Reactors_ode()
% 
% INPUTS:
%   - t: Current time. It is not used in the ode, but is included
%        to be compliant with the requirements of Matlab's ode functions.
%   - x: Current state of the system.
%        x = (h_1, c_A1, c_B1, T_1, h_2, c_A2, c_B2, T_2, h_3, c_A3, c_B3, T_3)
%   - u: Current control input.
%        u = (Q_1, Q_2, Q_3, F_f1, F_f2, F_R)
%   - param: Structure containing the parameters of the system.
%            The fields of the strucure are: xA0, xB0, xC0, T0, rho, Cp, alphaA,
%            alphaB, alphaC, alphaD, kA, kB, EAR, EBR, DeltaHA, DeltaHB, A1, A2,
%            A3, kv1, kv2, kv3
%            The meaning of each field is self explanatory from the description 
%            of the system parameters provided in the above reference.
%
% OUTPUTS:
%   - dx: Time derivative of the system state.
%
% See also: Reactors_gen_ss, Reactors_benchmark
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function dx = ode_Reactors(t, x, u, param)

    %% Estract state and input variables into readible names

    % States
    H1 = x(1);
    xA1 = x(2);
    xB1 = x(3);
    T1 = x(4);
    H2 = x(5);
    xA2 = x(6);
    xB2 = x(7);
    T2 = x(8);
    H3 = x(9);
    xA3 = x(10);
    xB3 = x(11);
    T3 = x(12);

    % Inputs
    Q1 = u(1);
    Q2 = u(2);
    Q3 = u(3);
    Ff1 = u(4);
    Ff2 = u(5);
    FR = u(6);
    
    %% Get variable values
    F1 = param.kv1*H1;
    F2 = param.kv2*H2;
    F3 = param.kv3*H3;
    kA1 = param.kA*exp(-param.EAR/T1);
    kA2 = param.kA*exp(-param.EAR/T2);
    kB1 = param.kB*exp(-param.EBR/T1);
    kB2 = param.kB*exp(-param.EBR/T2);
    xC3 = 1 - xA3 - xB3;
    x3hat = param.alphaA*xA3 + param.alphaB*xB3 + param.alphaC*xC3;
    xAR = (param.alphaA*xA3)/x3hat;
    xBR = (param.alphaB*xB3)/x3hat;
    FD = param.alphaD*FR;
    TR = T3;
    
    %% First order differential equations
    H1d = (1/(param.rho*param.A1))*(Ff1 + FR - F1);

    xA1d = (1/(param.rho*param.A1*H1))*(Ff1*(param.xA0 - xA1) + FR*(xAR - xA1)) - kA1*xA1;

    xB1d = (1/(param.rho*param.A1*H1))*(Ff1*(param.xB0 - xB1) + FR*(xBR - xB1)) + kA1*xA1 - kB1*xB1;

    T1d = (1/(param.rho*param.A1*H1))*(Ff1*(param.T0 - T1) + FR*(TR - T1)) - (1/param.Cp)*(kA1*xA1*param.DeltaHA +...
          kB1*xB1*param.DeltaHB) + Q1/(param.rho*param.A1*param.Cp*H1);

    H2d = (1/(param.rho*param.A2))*(Ff2 + F1 - F2);

    xA2d = (1/(param.rho*param.A2*H2))*(Ff2*(param.xA0 - xA2) + F1*(xA1 - xA2)) - kA2*xA2;

    xB2d = (1/(param.rho*param.A2*H2))*(Ff2*(param.xB0 - xB2) + F1*(xB1 - xB2)) + kA2*xA2 - kB2*xB2;

    T2d = (1/(param.rho*param.A2*H2))*(Ff2*(param.T0 - T2) + F1*(T1 - T2)) - (1/param.Cp)*(kA2*xA2*param.DeltaHA +...
          kB2*xB2*param.DeltaHB) + Q2/(param.rho*param.A2*param.Cp*H2);

    H3d = (1/(param.rho*param.A3))*(F2 - FD - FR - F3);

    xA3d = (1/(param.rho*param.A3*H3))*(F2*(xA2 - xA3) - (FD + FR)*(xAR - xA3));

    xB3d = (1/(param.rho*param.A3*H3))*(F2*(xB2 - xB3) - (FD + FR)*(xBR - xB3));

    T3d = (1/(param.rho*param.A3*H3))*(F2*(T2 - T3) - (FD + FR)*(TR - T3)) + Q3/(param.rho*param.A3*param.Cp*H3);

    %% Construct vector dx
    dx = [H1d; xA1d; xB1d; T1d; H2d; xA2d; xB2d; T2d; H3d; xA3d; xB3d; T3d];

end
