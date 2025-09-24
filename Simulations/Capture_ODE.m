function dydt = Capture_ODE(t, y)
global m_c;

mu_m = 4.9048695*1e12; % in m^3/s^-2

T = 2942; % Constant Thrust magnitude (in N)
mdot = -107.5/3600; % time rate of change of mass due to thrust (in kg/s)

m = m_c - mdot*t;
a_T = T/m; 

% Define your system of 7 differential equations
dydt = zeros(7,1);

u = atan2(-y(6), -y(7))+pi; % Optimal control 

%for capture
dydt(1) = -y(2);
dydt(2) = -(y(3)^2/y(1) - mu_m/y(1)^2 + a_T*sin(u));
dydt(3) = -(-y(2)*y(3)/y(1) + a_T*cos(u));
dydt(4) = -y(3)/y(1);
dydt(5) = -(y(6)*(y(3)^2/y(1)^2 - 2*mu_m/y(1)^3) - y(7)*y(2)*y(3)/(y(1)^2));
dydt(6) = -(-y(5) + y(7)*y(3)/y(1));
dydt(7) = -(-2*y(6)*y(3)/y(1) + y(7)*y(2)/y(1));