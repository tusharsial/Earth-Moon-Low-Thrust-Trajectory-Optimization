function dydt = VariableThrustODE(t, y)

mu_e = 3.986004418*1e14; % in m^3/s^-2
g = 9.80665;
%Isp = 10046.5287022250; % Constant Specific Impulse value 
Isp = 700;
%T_max = 2942;
T_max = 2942;
c = Isp*g; % Exhaust velocity (in m/s)
%m_leo = 1e5; % initial mass at LEO (in kg) 

alph = atan2(-y(7),-y(8)); % Optimal Thrust angle
X = (y(7)*sin(alph) + y(8)*cos(alph))/y(5) - y(9)/c;

if X < 0
    T = T_max;
elseif X > 0
    T = 0;
%else
%    T = (c^2/(y(9)*(c-1)))*((-y(7)+y(8)*y(3)/y(1))*sin(alph) + (-2*y(7)*y(3)/y(1) + y(8)*y(2)/y(1))*cos(alph));
end

disp(X);

a_T = T/y(5);

% Define your system of 7 differential equations
dydt = zeros(9,1);

dydt(1) = y(2);
dydt(2) = y(3)^2/y(1) - mu_e/y(1)^2 + a_T*sin(alph);
dydt(3) = -y(2)*y(3)/y(1) + a_T*cos(alph);
dydt(4) = y(3)/y(1);
dydt(5) = -T/c;
dydt(6) = y(7)*(y(3)^2/y(1)^2 - 2*mu_e/y(1)^3) - y(8)*(y(2)*y(3)/y(1)^2);
dydt(7) = -y(6) + y(8)*y(3)/y(1);
dydt(8) = -2*y(7)*y(3)/y(1) + y(8)*y(2)/y(1);
dydt(9) = (a_T/y(5))*(y(7)*sin(alph) + y(8)*cos(alph));