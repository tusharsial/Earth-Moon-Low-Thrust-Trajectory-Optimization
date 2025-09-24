function J = Capture_Cost_Function(s, x0, t0, t1)

mu_m = 4.9048695*1e12; % in m^3/s^-2

%g = 9.80665;
%Isp = 10046.5287022250; % Constant Specific Impulse value 
%c = Isp*g; % Exhaust velocity (in m/s)
%T_max = 2942;

% Solve ODEs using the guessed parameters
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[~, y] = ode45(@(t,y) Capture_ODE(t, y), [t0, t1], initial_conditions(s, x0), options); 

%[~, y] = ode45(@(t,y) VariableThrustODE(t, y), [t0, t1], initial_conditions(s, x0), options);


% Define the cost as the sum of squared errors in terminal conditions
J = 0;

%% Define the cost as the sum of squared errors in terminal conditions (Variable Thrust)
% alph_final = atan2(-y_final(7),-y_final(8)); % Optimal Thrust angle
% %m_inv_final = (y_final(9)/c)*(1/(y_final(7)*sin(alph_final) + y_final(8)*cos(alph_final)));
% X_final = (y_final(7)*sin(alph_final) + y_final(8)*cos(alph_final))/y_final(5) - y_final(9)/c;
% 
% if X_final < 0
%     T_final = T_max;
% elseif X_final > 0
%     T_final = 0;
% else
%     T_final = (c^2/(y_final(9)*(c-1)))*((-y_final(7)+y_final(8)*y_final(3)/y_final(1))*sin(alph_final) + (-2*y_final(7)*y_final(3)/y_final(1) + y_final(8)*y_final(2)/y_final(1))*cos(alph_final));
% end
% 
% J = (y_final(6) + mu_e/y_final(1)^2)^2 + (y_final(9) + (T_final/y_final(5)^2)*(y_final(7)*sin(alph_final) + y_final(8)*cos(alph_final)))^2;