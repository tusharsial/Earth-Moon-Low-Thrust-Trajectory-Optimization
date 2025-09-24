function dydt = Hybrid_ODE_2(t, y, t_capt, m_llo)

% constant values for Earth-Moon system
mu_e = 3.986004418*1e14; % in m^3/s^-2
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu = mu_e + mu_m; % in m^3/s^-2
D = 3.84399*1e8; % Distance between Earth and Moon (in m)
omg = sqrt(mu/D^3); % Constant angular rate of Earth-Moon system

% Thrusting only upto escape phase
if t<= t_capt
    T = 2942; % Constant Thrust magnitude (in N)
    u = atan2(-y(6), -y(7))+pi; % Optimal control for escape phase
else
    T = 0; % coast phase
    u = 0; % coast phase
end

mdot = -107.5/3600; % time rate of change of mass due to thrust (in kg/s)
m = m_llo - mdot*t;

a_T = T/m; 

% CR3BP EoM in Earth-centered polar rotating frame 
dydt = zeros(8,1);

r_esc = (y(1)^2 - 2*D*y(1)*cos(y(4)) + D^2)^(1/2); % r_earth-s/c

% EoM
dydt(1) = -y(2);
dydt(2) = -(y(3)^2/y(1) - mu_m/y(1)^2 - mu_e*(y(1) - D*cos(y(4)))/r_esc^3 - mu_e*cos(y(4))/D^2 + a_T*sin(u) + 2*omg*y(3) + omg^2*y(1));
dydt(3) = -(-mu_e*D*sin(y(4))/r_esc^3 + mu_e*sin(y(4))/D^2 + a_T*cos(u) - 2*omg*y(2) - y(2)*y(3)/y(1));
dydt(4) = -y(3)/y(1);

% Costate 
dydt(5) = -(y(6)*(y(3)^2/y(1)^2 + mu_e/r_esc^3 - 3*mu_e*(y(1) - D*cos(y(4)))^2/r_esc^5 - 2*mu_m/y(1)^3 - omg^2) + y(7)*(-y(2)*y(3)/y(1)^2 + 3*mu_e*(y(1) - D*cos(y(4)))*D*sin(y(4))/r_esc^5) + y(8)*y(3)/y(1)^2);
dydt(6) = -(-y(5) + y(7)*(y(3)/y(1) + 2*omg)); 
dydt(7) = -(y(6)*(-2*omg - 2*y(3)/y(1)) + y(7)*y(2)/y(1) - y(8)/y(1));
dydt(8) = -(y(6)*(mu_e*D*sin(y(4))/r_esc^3 - 3*mu_e*(y(1) - D*cos(y(4)))*D*y(1)*sin(y(4))/r_esc^5 - mu_e*sin(y(4))/D^2) + y(7)*(mu_e*D*cos(y(4))/r_esc^3 - 3*mu_e*D^2*y(1)*(sin(y(4)))^2/r_esc^5 - mu_e*cos(y(4))/D^2));
