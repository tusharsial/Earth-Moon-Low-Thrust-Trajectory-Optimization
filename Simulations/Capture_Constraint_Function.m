%% Constraint Function (Ensures Boundary Conditions are Satisfied)
function [c, ceq] = Capture_Constraint_Function(s, x0, t0, t1)
mu_m = 4.9048695*1e12; % in m^3/s^-2

% Solve ODEs using the guessed parameters
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[~, y] = ode45(@(t,y) Capture_ODE(t, y), [t0, t1], initial_conditions(s, x0), options); 
%[~, y] = ode45(@(t,y) VariableThrustODE(t, y), [t0, t1], initial_conditions(s, x0), options);

% Extract final values at t1
y_final = y(end, :);

% No inequality constraints (c = [])
c = [];

% Equality constraints (terminal conditions)
ceq = [y_final(5) + mu_m/y_final(1)^2; y_final(6) + y_final(2); y_final(7) + y_final(3)]; % for moon
