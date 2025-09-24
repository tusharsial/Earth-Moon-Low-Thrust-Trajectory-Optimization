function dydt = Escape_ODE(t, y)

mu_e = 3.986004418*1e14; % in m^3/s^-2

T = 2942; % Constant Thrust magnitude (in N)
m_leo = 1e5; % initial mass at LEO (in kg) 
mdot = 107.5/3600; % time rate of change of mass due to thrust (in kg/s)

m = m_leo - mdot*t;
a_T = T/m; 

% Define your system of 7 differential equations
dydt = zeros(7,1);

u = atan2(-y(6), -y(7)); % Optimal control 

%For escape
dydt(1) = y(2);
dydt(2) = y(3)^2/y(1) - mu_e/y(1)^2 + a_T*sin(u);
dydt(3) = -y(2)*y(3)/y(1) + a_T*cos(u);
dydt(4) = y(3)/y(1);
dydt(5) = y(6)*(y(3)^2/y(1)^2 - 2*mu_e/y(1)^3) - y(7)*y(2)*y(3)/(y(1)^2);
dydt(6) = -y(5) + y(7)*y(3)/y(1);
dydt(7) = -2*y(6)*y(3)/y(1) + y(7)*y(2)/y(1);