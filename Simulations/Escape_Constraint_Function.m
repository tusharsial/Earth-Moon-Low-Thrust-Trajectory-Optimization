%% Constraint Function (Ensures Boundary Conditions are Satisfied)
function [cineq, ceq] = Escape_Constraint_Function(s, x0, t0, t1)
mu_e = 3.986004418*1e14; % in m^3/s^-2

% Solve ODEs using the guessed parameters
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
%[~, y] = ode45(@(t,y) Escape_ODE(t, y), [t0, t1], initial_conditions(s, x0), options); 
[~, y] = ode45(@(t,y) VariableThrustODE(t, y), [t0, t1], initial_conditions(s, x0), options);

% Extract final values at t1
y_final = y(end, :);

% No inequality constraints (c = [])
cineq = [];

% Equality constraints (terminal conditions)
%ceq = [y_final(5) + mu_e/(y_final(1)^2); y_final(6) + y_final(2); y_final(7) + y_final(3)]; % for earth

%ceq = [y_final(6) + mu_e/(y_final(1)^2); y_final(7) + y_final(2); y_final(8) + y_final(3); y_final(9) + (T_final/y_final(5)^2)*(y_final(7)*sin(alph_final) + y_final(8)*cos(alph_final))]; % for earth with variable thrust
ceq = [y_final(6) + mu_e/(y_final(1)^2); y_final(7) + y_final(2); y_final(8) + y_final(3); y_final(9)]; % for earth with variable thrust