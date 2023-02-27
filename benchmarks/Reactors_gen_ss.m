%% Reactors_gen_ss - Generate continuous-time model of the reactors system
%
% This function generates the continuous-time state space model
% of the reactors system described in:
% 
% "Implementation of model predictive control in programmable logic controllers",
% P. Krupa, D. Limon and T. Alamo, in Transactions on Control Systems Technology, 2021.
%
% INPUTS:
%   - p: Structure containing the parameters of the system.
%        The fields of the strucure are: xA0, xB0, xC0, T0, rho, Cp, alphaA,
%        alphaB, alphaC, alphaD, kA, kB, EAR, EBR, DeltaHA, DeltaHB, A1, A2,
%        A3, kv1, kv2, kv3
%        The meaning of each field is self explanatory from the description 
%        of the system parameters provided in the above reference.
%   - x0: Operating point of the system state.
%   - u0: Operating point of the system input.
%
% OUTPUTS:
%   - sys: Instance of the ss class. Continuous-time state space model.
%
% See also: Reactors_ode, Reactors_benchmark
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function sys = Reactors_gen_ss(p, x0, u0)

    %% Rename variables
    
    % State variables
    H1_0 = x0(1); xA1_0 = x0(2); xB1_0 = x0(3); T1_0 = x0(4);
    H2_0 = x0(5); xA2_0 = x0(6); xB2_0 = x0(7); T2_0 = x0(8);
    H3_0 = x0(9); xA3_0 = x0(10); xB3_0 = x0(11); T3_0 = x0(12);
    
    % Input variables
    Q1_0 = u0(1); Q2_0 = u0(2); Q3_0 = u0(3);
    Ff1_0 = u0(4); Ff2_0 = u0(5); FR_0 = u0(6);  
    
    %% Calculate OpPoint value of inner variables
    den_AR_0 = (p.alphaA - p.alphaC)*xA3_0 + (p.alphaB - p.alphaC)*xB3_0 + p.alphaC;
    xAR_0 = (p.alphaA*xA3_0)/den_AR_0;
    xBR_0 = (p.alphaB*xB3_0)/den_AR_0;
    d_xAR_xA3_0 = -(p.alphaA - p.alphaC)*p.alphaA*xA3_0/(den_AR_0^2) + p.alphaA/den_AR_0;
    d_xAR_xB3_0 = -(p.alphaB - p.alphaC)*p.alphaA*xA3_0/(den_AR_0^2);
    d_xBR_xA3_0 = -(p.alphaA - p.alphaC)*p.alphaB*xB3_0/(den_AR_0^2);
    d_xBR_xB3_0 = -(p.alphaB - p.alphaC)*p.alphaB*xB3_0/(den_AR_0^2) + p.alphaB/den_AR_0;
    kA1_0 = p.kA*exp(-p.EAR/T1_0);
    kB1_0 = p.kB*exp(-p.EBR/T1_0);
    kA2_0 = p.kA*exp(-p.EAR/T2_0);
    kB2_0 = p.kB*exp(-p.EBR/T2_0);
    
    %% Calculate the elements of matrix A
    % The order of the states is: H1, xA1, xB1, T1, H2, xA2, xB2, T2, H3, xA3, xB3, T3
    A = zeros(8);
    % H1
    A(1,1)  = (-p.kv1)/(p.rho*p.A1); % H1 <- H1
    % xA1
    A(2,1)  = -p.rho*p.A1*(Ff1_0*(p.xA0 - xA1_0) + FR_0*(xAR_0 - xA1_0))/(p.rho*p.A1*H1_0)^2; % xA1 <- H1
    A(2,2)  = (-Ff1_0 - FR_0)/(p.rho*p.A1*H1_0) - kA1_0; % xA1 <- xA1
    A(2,4)  = -xA1_0*p.kA*exp(-p.EAR/T1_0)*p.EAR/(T1_0^2); % xA1 <- T1
    A(2,10) = FR_0*d_xAR_xA3_0/(p.rho*p.A1*H1_0); % xA1 <- xA3
    A(2,11) = FR_0*d_xAR_xB3_0/(p.rho*p.A1*H1_0); % xA1 <- xB3
    % xB1
    A(3,1)  = -p.rho*p.A1*(Ff1_0*(p.xB0 - xB1_0) + FR_0*(xBR_0 - xB1_0))/(p.rho*p.A1*H1_0)^2; % xB1 <- H1
    A(3,2)  = kA1_0; % xB1 <- xA1
    A(3,3)  = (-Ff1_0 - FR_0)/(p.rho*p.A1*H1_0) - kB1_0;
    A(3,4)  = -xB1_0*p.kB*exp(-p.EBR/T1_0)*p.EBR/(T1_0^2) + xA1_0*p.kA*exp(-p.EAR/T1_0)*p.EAR/(T1_0^2); % xB1 <- T1
    A(3,10) = FR_0*d_xBR_xA3_0/(p.rho*p.A1*H1_0); % xB1 <- xA3
    A(3,11) = FR_0*d_xBR_xB3_0/(p.rho*p.A1*H1_0); % xB1 <- xB3
    % T1
    A(4,1)  = -p.rho*p.A1*(Ff1_0*(p.T0 - T1_0) + FR_0*(T3_0 - T1_0))/(p.rho*p.A1*H1_0)^2 - p.rho*p.A1*p.Cp*Q1_0/(p.rho*p.A1*H1_0*p.Cp)^2; % T1 <- H1
    A(4,2)  = -kA1_0*p.DeltaHA/p.Cp; % T1 <- xA1
    A(4,3)  = -kB1_0*p.DeltaHB/p.Cp; % T1 <- xB1
    A(4,4)  = (-Ff1_0 - FR_0)/(p.rho*p.A1*H1_0) - (1/p.Cp)*(xA1_0*p.DeltaHA*p.kA*exp(-p.EAR/T1_0)*p.EAR/(T1_0^2) + xB1_0*p.DeltaHB*p.kB*exp(-p.EBR/T1_0)*p.EBR/(T1_0^2)); % T1 <- T1
    A(4,12) = FR_0/(p.rho*p.A1*H1_0); % T1 <- T3
    % H2
    A(5,1)  = p.kv1/(p.rho*p.A2); % H2 <- H1
    A(5,5)  = -p.kv2/(p.rho*p.A2);
    % xA2
    A(6,1)  = p.kv1*(xA1_0 - xA2_0)/(p.rho*p.A2*H2_0); % xA2 <- H1
    A(6,2)  = p.kv1*H1_0/(p.rho*p.A2*H2_0); % xA2 <- xA1
    A(6,5)  = -p.rho*p.A2*(Ff2_0*(p.xA0 - xA2_0) + p.kv1*H1_0*(xA1_0 - xA2_0))/(p.rho*p.A2*H2_0)^2; % xA2 <- H2
    A(6,6)  = (-Ff2_0 - p.kv1*H1_0)/(p.rho*p.A2*H2_0) - kA2_0; % xA2 <- xA2
    A(6,8)  = -xA2_0*p.kA*exp(-p.EAR/T2_0)*p.EAR/(T2_0^2);  % xA2 <- T2
    % xB2
    A(7,1)  = p.kv1*(xB1_0 - xB2_0)/(p.rho*p.A2*H2_0); % xB2 <- H1
    A(7,3)  = p.kv1*H1_0/(p.rho*p.A2*H2_0); % xB2 <- xB1
    A(7,5)  = -p.rho*p.A2*(Ff2_0*(p.xB0 - xB2_0) + p.kv1*H1_0*(xB1_0 - xB2_0))/(p.rho*p.A2*H2_0)^2; % xB2 <- H2
    A(7,6)  = kA2_0; % xB2 <- xA2
    A(7,7)  = (-Ff2_0 - p.kv1*H1_0)/(p.rho*p.A2*H2_0) - kB2_0; % xB2 <- xB2
    A(7,8)  = -xB2_0*p.kB*exp(-p.EBR/T2_0)*p.EBR/(T2_0^2) + xA2_0*p.kA*exp(-p.EAR/T2_0)*p.EAR/(T2_0^2); % xB2 <- T2
    % T2
    A(8,1)  = p.kv1*(T1_0 - T2_0)/(p.rho*p.A2*H2_0); % T2 <- H1
    A(8,4)  = p.kv1*H1_0/(p.rho*p.A2*H2_0); % T2 <- T1
    A(8,5)  = -p.rho*p.A2*(Ff2_0*(p.T0 - T2_0) + p.kv1*H1_0*(T1_0 - T2_0))/(p.rho*p.A2*H2_0)^2 - p.rho*p.A2*p.Cp*Q2_0/(p.rho*p.A2*H2_0*p.Cp)^2; % T2 <- H2
    A(8,6)  = -kA2_0*p.DeltaHA/p.Cp; % T2 <- xA2
    A(8,7)  = -kB2_0*p.DeltaHB/p.Cp; % T2 <- xB2
    A(8,8)  = (-Ff2_0 - p.kv1*H1_0)/(p.rho*p.A2*H2_0) - (1/p.Cp)*(xA2_0*p.DeltaHA*p.kA*exp(-p.EAR/T2_0)*p.EAR/(T2_0^2) + xB2_0*p.DeltaHB*p.kB*exp(-p.EBR/T2_0)*p.EBR/(T2_0^2)); % T2 <- T2
    % H3
    A(9,5)  = p.kv2/(p.rho*p.A3); % H3 <- H2
    A(9,9)  = -p.kv3/(p.rho*p.A3); % H3 <- H3
    % xA3
    A(10,5) = p.kv2*(xA2_0 - xA3_0)/(p.rho*p.A3*H3_0); % xA3 <- H2
    A(10,6) = p.kv2*H2_0/(p.rho*p.A3*H3_0); % xA3 <- xA2
    A(10,9) = -p.rho*p.A3*(p.kv2*H2_0*(xA2_0 - xA3_0) - (1+p.alphaD)*FR_0*(xAR_0 - xA3_0))/(p.rho*p.A3*H3_0)^2; % xA3 <- H3
    A(10,10)= (-p.kv2*H2_0 + (1+p.alphaD)*FR_0 - (1+p.alphaD)*FR_0*d_xAR_xA3_0)/(p.rho*p.A3*H3_0); % xA3 <- xA3
    A(10,11)= -(1+p.alphaD)*FR_0*d_xAR_xB3_0/(p.rho*p.A3*H3_0); % xA3 <- xB3
    % xB3
    A(11,5) = p.kv2*(xB2_0 - xB3_0)/(p.rho*p.A3*H3_0); % xB3 <- H2
    A(11,7) = p.kv2*H2_0/(p.rho*p.A3*H3_0); % xB3 <- xB2
    A(11,9) = -p.rho*p.A3*(p.kv2*H2_0*(xB2_0 - xB3_0) - (1+p.alphaD)*FR_0*(xBR_0 - xB3_0))/(p.rho*p.A3*H3_0)^2; % xB3 <- H3
    A(11,10)= -(1+p.alphaD)*FR_0*d_xBR_xA3_0/(p.rho*p.A3*H3_0); % xB3 <- xA3
    A(11,11)= (-p.kv2*H2_0 + (1+p.alphaD)*FR_0 - (1+p.alphaD)*FR_0*d_xBR_xB3_0)/(p.rho*p.A3*H3_0); % xB3 <- xB3
    % T3
    A(12,5) = p.kv2*(T2_0 - T3_0)/(p.rho*p.A3*H3_0); % T3 <- H2
    A(12,8) = p.kv2*H2_0/(p.rho*p.A3*H3_0); % T3 <- T2
    A(12,9) = -p.rho*p.A3*p.kv2*H2_0*(T2_0 - T3_0)/(p.rho*p.A3*H3_0)^2 - p.rho*p.A3*p.Cp*Q3_0/(p.rho*p.A3*H3_0*p.Cp)^2; % T3 <- H3
    A(12,12)= -p.kv2*H2_0/(p.rho*p.A3*H3_0); % T3 <- T3
    
    %% Calculate the elements of matrix B
    % The order of the outputs is: (Q_1, Q_2, Q_3, F_f1, F_f2, F_R)
    B = zeros(8, 6);
    % H1
    B(1,4) = 1/(p.rho*p.A1); % H1 <- Ff1
    B(1,6) = 1/(p.rho*p.A1); % H1 <- FR
    % xA1
    B(2,4) = (p.xA0 - xA1_0)/(p.rho*p.A1*H1_0); % xA1 <- Ff1
    B(2,6) = (xAR_0 - xA1_0)/(p.rho*p.A1*H1_0); % xA1 <- FR
    % xB1
    B(3,4) = (p.xB0 - xB1_0)/(p.rho*p.A1*H1_0); % xB1 <- Ff1
    B(3,6) = (xBR_0 - xB1_0)/(p.rho*p.A1*H1_0); % xB1 <- FR
    % T1
    B(4,4) = (p.T0 - T1_0)/(p.rho*p.A1*H1_0); % T1 <- Ff1
    B(4,1) = 1/(p.rho*p.A1*H1_0*p.Cp); % T1 <- Q1
    B(4,6) = (T3_0 - T1_0)/(p.rho*p.A1*H1_0); % T1 <- FR
    % H2
    B(5,5) = 1/(p.rho*p.A2); % H2 <- Ff2
    % xA2
    B(6,5) = (p.xA0 - xA2_0)/(p.rho*p.A2*H2_0); % xA2 <- Ff2
    % xB2
    B(7,5) = (p.xB0 - xB2_0)/(p.rho*p.A2*H2_0); % xB2 <- Ff2
    % T2
    B(8,5) = (p.T0 - T2_0)/(p.rho*p.A2*H2_0); % T2 <- Ff2
    B(8,2) = 1/(p.rho*p.A2*H2_0*p.Cp); % T2 <- Q2
    % H3
    B(9,6) = -(1+p.alphaD)/(p.rho*p.A3); % H3 <- FR
    % xA3
    B(10,6)= -(1+p.alphaD)*(xAR_0 - xA3_0)/(p.rho*p.A3*H3_0); % xA3 <- FR
    % xB3
    B(11,6)= -(1+p.alphaD)*(xBR_0 - xB3_0)/(p.rho*p.A3*H3_0); % xB3 <- FR
    % T3
    B(12,3)= 1/(p.rho*p.A3*H3_0*p.Cp); % T3 <- Q3

    %% Calculate the elements of matrix C
    C = eye(12);

    %% Calculate the elements of matrix D
    D = zeros(12, 6);

    %% Construct continuos-time state space model
    sys = ss(A, B, C, D);
    
end
