function J = Escape_Cost_Function(s, x0, t0, t1)
mu_e = 3.986004418*1e14; % in m^3/s^-2

% Solve ODEs using the guessed parameters
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
%[~, y] = ode45(@(t,y) Escape_ODE(t, y), [t0, t1], initial_conditions(s, x0), options); 

[~, y] = ode45(@(t,y) VariableThrustODE(t, y), [t0, t1], initial_conditions(s, x0), options);


% Define the cost as the sum of squared errors in terminal conditions
J = 0;