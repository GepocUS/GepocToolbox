%% BaP_benchmark - Generate benchmark of the ball and plate system 
%
% Returns a specific instance of the ball and plate system.
% In particular, it returns the ball and plate system described in:
%
% "Harmonic based model predictive control for set-point tracking",
% P. Krupa, D. Limon and T. Alamo, in Transactions on Automatic control, 2021.
%
% The continuous-time model of the system is obtained using gen_sys_BaP.m
%
% The state of the system is given by:
%   x = (p_1, dp_1, teta_1, dteta_1, p_2, dp_2, teta_2, dteta_2)
% The input of the system is given by:
%   u = (ddp_1, ddp_2)
%
% OUTPUTS:
%   - sysBaP: Instance of ssModel. Discrete-time model of the system
%
% See also: BaP_gen_ss, BaP_ode
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function sysBaP = BaP_benchmark()

    %% Define system

    % Parameters
    param_sys.m = 0.05; % Mass of the ball [Kg]
    param_sys.r = 0.01; % Radius of the ball [m]
    param_sys.g = 9.81; % Gravitational constant [m/s^2]
    param_sys.I = 2*param_sys.m*param_sys.r^2/5; % Mass momentum of inertia of the ball [Kg m^2]. This one is for a solid ball
    %param.I = 2*param.m*param.r^2/3; % Mass momentum of inertia of the ball [Kg m^2]. This one is for a hollow ball

    % Constraints (in engineering units)
    accel_limit = 0.4; % Bound on the acceleration of the ball
    speed_limit = 0.5; % Bound on the speed of the ball
    angle_limit = pi/4; % Bound on the angle of the plate

        % State constraints
    LBxEng = -[10000; speed_limit; angle_limit; 10000;  10000; speed_limit; angle_limit; 10000]; % Lower bound for the states
    UBxEng = [10000; speed_limit; angle_limit; 10000;  10000; speed_limit; angle_limit; 10000]; % Upper bound for the states
        % Input constraints
    LBuEng = -accel_limit*ones(2,1); % Lower bound for the control inputs
    UBuEng = accel_limit*ones(2, 1); % Upper bound for the control inputs

    % Operating point
    xOpPoint = zeros(8, 1);
    uOpPoint = zeros(2, 1);

    %% Generate continuous-time model
    sysC = gen_sys_BaP(param_sys);

    %% Generate ssModel
    Ts = 0.2; % Sampling time

    % Scaling matrices
    Nx = [0.1; 1; 1; 1; 0.1; 1; 1; 1]; % Positions are scaled by 0.1 (units in dm)
    Nu = 1*ones(2, 1);
    Ny = Nx;

    % Create ssModel
    sysD = c2d(sysC, Ts);
    sysBaP = ssModel(sysD.A, sysD.B, sysD.C, sysD.D, Ts);
    sysBaP.X0 = xOpPoint;
    sysBaP.IN0 = uOpPoint;
    sysBaP.Nx = Nx;
    sysBaP.Nu = Nu;
    sysBaP.Ny = Ny;
    sysBaP.LBxEng = LBxEng;
    sysBaP.UBxEng = UBxEng;
    sysBaP.LBuEng = LBuEng;
    sysBaP.UBuEng = UBuEng;
    sysBaP.yc = [2 3 6 7];
    sysBaP.LByEng = sysBaP.Cc*sysBaP.LBxEng;
    sysBaP.UByEng = sysBaP.Cc*sysBaP.UBxEng;
    sysBaP.yr = [1 5]; % The controlled outputs are the position of the ball in each one of the axis

end
