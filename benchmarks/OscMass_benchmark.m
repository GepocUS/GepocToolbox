%% OscMass_benchmark - Generate benchmark of the ball and plate system 
% 
% Returns a specific instance of the oscillating masses system.
% In particular, it returns the oscillating masses system described in:
% 
% "A sparse ADMM-based solver for linear MPC subject to terminal quadratic constraint",
% P. Krupa, R. Jaouani, D. Limon and T. Alamo, arXiv preprint, 2021.
% 
% The continuous-time model of the system is obtained using OscMass_gen_ss()
% 
% The state of the system is given by:
%   x = (p_1, p_2, p_3, v_1, v_2, v_3)
% The input of the system is given by:
%   u = (F_f, F_l)
%
% OUTPUTS:
%   - sys: Instance of ssModel. Discrete-time model of the system
%
% See also: OscMass_gen_ss
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function sysOscMass = OscMass_benchmark()

    %% Define system

    % Parameters
    p = 3; % Number of objects
    M = [1; 0.5*ones(p-2, 1); 1]; % Mass of each object
    K = 2*ones(p+1, 1); % Spring constants
    F = [1; zeros(p-2, 1); 1]; % Objects on which external forces are applied

    % Constraints (in engineering units)
        % State constraints
    LBxEng = -[ones(p, 1); 1000*ones(p, 1)]; % Lower bound for the states
    UBxEng = [0.3; 0.3; 0.3; 1000*ones(p, 1)]; % Upper bound for the states
        % Input constraints
    LBuEng = -[0.8; 0.8]; % Lower bound for the control inputs
    UBuEng = [0.8; 0.8]; % Upper bound for the control inputs

    % Operating point
    xOpPoint = zeros(2*p, 1);
    uOpPoint = zeros(sum(F), 1);

    %% Generate continuous-time model
    sysC = OscMass_gen_ss(M, K, F);

    %% Generate ssModel
    Ts = 0.2; % Sampling time

    % Scaling matrices
    Nx = [10*ones(p, 1); ones(p, 1)];
    Nu = ones(sum(F), 1);
    Ny = Nx;

    % Create ssModel
    sysD = c2d(sysC, Ts);
    sysOscMass = ssModel(sysD.A, sysD.B, sysD.C, sysD.D, Ts);
    sysOscMass.X0 = xOpPoint;
    sysOscMass.IN0 = uOpPoint;
    sysOscMass.Nx = Nx;
    sysOscMass.Nu = Nu;
    sysOscMass.Ny = Ny;
    sysOscMass.LBxEng = LBxEng;
    sysOscMass.UBxEng = UBxEng;
    sysOscMass.LBuEng = LBuEng;
    sysOscMass.UBuEng = UBuEng;
    sysOscMass.yc = [1 2 3];
    sysOscMass.LByEng = sysOscMass.Cc*sysOscMass.LBxEng;
    sysOscMass.UByEng = sysOscMass.Cc*sysOscMass.UBxEng;

end
