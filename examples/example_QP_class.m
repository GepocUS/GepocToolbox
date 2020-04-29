%% This script is an example of the use of the class QP
%
% We construct and solve a simple QP problem using quadprog. 
% We also show how to change the ingredients of the QP problem after construction
% and we show how to use a solver that is not listed among the default solvers of the class.
% 
clear; clc;

%% Initialize ingredients of the QP problem
% We are going to be constructing and solving QP problems of the form
% 
%   min 1/2 z'*H*z + q'*z
%     s.t. A*z = b
%          C*z <= d
%          LB <= z <= UB
%
% We first need to initialize the ingredients of our example
% For now, we will only define a QP problem with box constraints

n = 10; % Number of decision variables

% Obtain a positive definite matrix H and a vector q
H = 1*rand(n, n);
H = 0.5*(H'*H);
q = -1*rand(n, 1);

% Obtain box constraints
UB = rand(n, 1);
LB = -UB;

%% Construct instance of QP class
% Fist we construct the QP problem using the constructor method myQP = QP(H, q, A, b, C, d, LB, UB)
myQP = QP(H, q, [], [], [], [], LB, UB);

% We can also obtain a QP problem with no box constraints by using the constructos myQP_unconstrained = QP(H, q)
myQP_unconstrained = QP(H, q);
% Note that only the H and q arguments are required in the constructor. All other arguments are optional

%% Solving the QP
% Solving the QP is as simple as calling its QP.solve method as follows
[z_opt, f_opt, e_flag] = myQP.solve;
% z_opt is the optimal solution, f_opt the value of the optimal functional and e_flag is a flag returned by the solver
% In this case, since we are using rhe default solver 'quadprog', e_flag = 1 indicates that the problem was successfully solved

% We can do the same for the unconstrained QP
[z_opt_unconstrained, f_opt_unconstrained, e_flag_unconstrained] = myQP_unconstrained.solve;

%% Calling other methods of the QP solver
% We show how calls to the methods of the QP class are done
% As an example, let us call the method QP.eval, which returned the value of the cost function at the given point z
f = myQP.eval(z_opt);
% Since we are evaluating the functional at point z_opt, its value should be exaclty equal to the f_opt obtained above

%% Changing parameters of myQP
% The properties of the QP instance can be changed after its construction
% For example, you could change the value of q as follows:
% myQP.q = -10*rand(n, 1);
% If you did so, the variable myQP in your workspace would update its q to the given value
% You could also add an equality constraint as follows
% myQP.A = rand(1, n);
% myQP.b = 0;

%% Using non-supported solvers
% The class QP supports a list of solvers, which can be selected with the property 'solver'
% However, it also supports the user providing his own function handler to any other solver.
% The only requirement is that the solver must accept the following arguments, in this order
% z_opt = my_QP_solver(H, q, c, d, A, b, LB, UB, varargin);
% where varargin can be any amount of additional arguments

% Lets get a function handler to our own solver. This solver is for an unconstrained QP
solver_unconstrained_QP = @(H, q, c, d, A, b, LB, UB, varargin) -inv(H)*q;

% We then set the solver property of my instance myQP_unconstrained to our function handler
myQP_unconstrained.solver = solver_unconstrained_QP;

% We can now solve the problem as usual
z_opt_mysolver = myQP_unconstrained.solve();
% z_opt_mysolver should be equal to the value returned when we used the quadprog solver, i.e. z_opt_mysolver = z_opt_unconstrained
