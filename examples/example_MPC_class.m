%% This script is an example of use of various MPC controller classes
% Example script for the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 
% We define a simple system and control it using two MPC controllers, the
% MPC for tracking (TrackingMPC.m) and the Lax MPC (LaxMPC.m).
% We show how to construct the classes, set the reference and current system
% state and solve the optimization problem.
% We use the ssModel class to define the system
%
clear; clc;

%% Define the system
% We define a simple system diven by the state space model x_{k+1} = A x_k + B u_k
A = [0.9 0.8; 0 0.8];
B = [0; 0.95];
C = eye(size(A, 1));
D = 0;
Ts = 1;

% Subject the system to box constraints
LBx = -[10; 10]; % Lower bound constraints for the state
UBx = -LBx; % Uppper bound constraints for the state
LBu = -1; % Lower bound constraints for the input
UBu = 1; % Upper bound constraints for the input

% Construct the system model as an ssModel object
sys = ssModel(A, B, C, D, Ts);

% Set the box constraints
sys.LBx = LBx;
sys.UBx = UBx;
sys.LBu = LBu;
sys.UBu = UBu;

%% Define the initial system state and reference
x0 = rand(sys.n_x, 1); % Initial system state
xr = rand(sys.n_x, 1); % Reference for the system state
ur = rand(sys.n_u, 1); % Reference for the system input

%% Define the ingredients of the MPC controllers

% Both controllers
N = 10;
Q = 5*eye(sys.n_x);
R = eye(sys.n_u);

% MPC for tracking
T = 2*N*Q;
S = R;

% Lax MPC
P = 2*N*Q;

%% Construct the MPC objects

% Tracking MPC
MPCT = TrackingMPC(sys, Q, R, T, S, N);

% Lax MPC
LMPC = LaxMPC(sys, Q, R, P, N);

%% Set the current state and reference

% Tracking MPC
MPCT.x0 = x0;
MPCT.xr = xr;
MPCT.ur = ur;

% Lax MPC
LMPC.x0 = x0;
LMPX.xr = xr;
LMPC.ur = ur;

%% Solve the MPC problems

% Tracking MPC
[ut, zt, ft, et, Ht] = MPCT.solve;

% Lax MPC
[ul, zl, fl, el, Hl] = LMPC.solve;
