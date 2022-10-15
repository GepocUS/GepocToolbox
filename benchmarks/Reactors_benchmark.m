%% Reactors_benchmark - Generate benchmark of the reactors system 
% 
% Returns a specific instance of the reactors system.
% In particular, it returns the reactors system described in:
%
% "Implementation of model predictive control in programmable logic controllers",
% P. Krupa, D. Limon and T. Alamo, in Transactions on Control Systems Technology, 2021.
%
% The continuous-time model of the system is obtained using Reactors_gen_ss()
%
% The state of the system is given by:
%   x = (h_1, c_A1, c_B1, T_1, h_2, c_A2, c_B2, T_2, h_3, c_A3, c_B3, T_3)
% The input of the system is given by:
%   u = (Q_1, Q_2, Q_3, F_f1, F_f2, F_R)
%
% OUTPUTS:
%   - sysReac: Instance of ssModel. Discrete-time model of the system
%
% See also: Reactors_gen_ss, Reactors_ode
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function sysReac = Reactors_benchmark()

    %% Define system

    % Parameters
        % External conditions
    xA0 = 1; % Concentration of A in the input flows
    xB0 = 0; % Concentration of A in the input flows
    xC0 = 0; % Concentration of A in the input flows
    T0 = 313; % Temperature of the input flows
        % Tanks
    A1 = 1; % Area of reactor 1
    A2 = 1; % Area of reactor 2
    A3 = 1; % Area of separator
    kv1 = 50; % Emptying speed of reactor 1
    kv2 = 50; % Emptying speed of reactor 2
    kv3 = 30; % Emptying speed of separator
        % Material (we assume it is the same regardless of the relative concentrations of the agents)
    rho = 1100; % Density
    Cp = 4; % Specific heat
    alphaA = 3.5;
    alphaB = 1.1;
    alphaC = 0.5;
    alphaD = 1e-3;
        % Reactions (speeds, etc)
    kA = 1e-5;
    kB = 0.5e-5;
    EAR = -2840;
    EBR = -2077;
    DeltaHA = -100;
    DeltaHB = -39;

    param_sys = struct('xA0', xA0, 'xB0', xB0, 'xC0', xB0, 'T0', T0, 'rho', rho, 'Cp', Cp, 'alphaA', alphaA, 'alphaB', alphaB,...
                       'alphaC', alphaC, 'alphaD', alphaD, 'kA', kA, 'kB', kB, 'EAR', EAR, 'EBR', EBR, 'DeltaHA', DeltaHA,...
                       'DeltaHB', DeltaHB, 'A1', A1, 'A2', A2, 'A3', A3, 'kv1', kv1, 'kv2', kv2, 'kv3', kv3);

    % Constraints (in engineering units)
        % State constraints
    LBxEng = kron(ones(3,1), [0; 0; 0; 320]); % Lower bound for the states
    UBxEng = kron(ones(3,1), [2; 1; 1; 348]); % Upper bound for the states
    UBxEng(12) = 338; % Make upper bound of T3 smaller (75ºC)
        % Input constraints
    LBuEng = -[5000; 5000; 5000; 0; 0; 0]; % Lower bound for the control inputs
    UBuEng = [5000; 5000; 5000; 50; 50; 50]; % Upper bound for the control inputs

    % Operating point
        % Inputs
    Q1_0 = 0;
    Q2_0 = 0;
    Q3_0 = 0;
    Ff1_0 = 30;
    Ff2_0 = 10;
    FR_0 = 5;
    uOpPoint = [Q1_0; Q2_0; Q3_0; Ff1_0; Ff2_0; FR_0];
        % State
    fprintf('Computing operating point of the reactors system\n')
    options_fsolve = optimoptions('fsolve', 'StepTolerance', 1e-9, 'OptimalityTolerance', 1e-9, 'FunctionTolerance', 1e-9, 'MaxIterations', 1000);
    xinit_fsolve = kron( ones(3,1), [10; 0.5; 0.5; 300]);
    xOpPoint = fsolve(@(xx) Reactors_ode(0, xx, uOpPoint, param_sys), xinit_fsolve, options_fsolve);

    %% Generate the continuous-time model
    sysC = Reactors_gen_ss(param_sys, xOpPoint, uOpPoint);

    %% Generate ssModel
    Ts = 3; % Sampling time

    % Scaling matrices
    Nx = 1./kron(ones(3,1), [1; 1; 1; 10]);
    Nu = 1./[1000*ones(3, 1); 10*ones(3,1)];
    Ny = Nx;

    % Create ssModel
    sysD = c2d(sysC, Ts);
    sysReac = ssModel(sysD.A, sysD.B, sysD.C, sysD.D, Ts);
    sysReac.X0 = xOpPoint;
    sysReac.IN0 = uOpPoint;
    sysReac.Nx = Nx;
    sysReac.Nu = Nu;
    sysReac.Ny = Ny;
    sysReac.LBxEng = LBxEng;
    sysReac.UBxEng = UBxEng;
    sysReac.LBuEng = LBuEng;
    sysReac.UBuEng = UBuEng;
    sysReac.LByEng = sysReac.Cc*sysReac.LBxEng;
    sysReac.UByEng = sysReac.Cc*sysReac.UBxEng;

end
