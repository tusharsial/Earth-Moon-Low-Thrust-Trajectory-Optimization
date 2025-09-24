function dydt = Coasting_ODE(t, y)

% constant values for Earth-Moon system
mu_e = 3.986004418*1e14; % in m^3/s^-2
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu = mu_e + mu_m; % in m^3/s^-2
D = 3.84399*1e8; % Distance between Earth and Moon (in m)
omg = sqrt(mu/D^3); % Constant angular rate of Earth-Moon system

%% Set Thrust


% CR3BP EoM in Earth-centered polar rotating frame
dydt = zeros(4,1);

dydt(1) = y(2);
dydt(2) = y(3)^2/y(1) - mu_e/y(1)^2 - mu_m*(y(1) + D*cos(y(4)))/(y(1)^2 + 2*D*y(1)*cos(y(4)) + D^2)^(3/2) + mu_m*cos(y(4))/D^2 + 2*omg*y(3) + omg^2*y(1);
dydt(3) = mu_m*D*sin(y(4))/(y(1)^2 + 2*D*y(1)*cos(y(4)) + D^2)^(3/2) - mu_m*sin(y(4))/D^2 - 2*omg*y(2) - y(2)*y(3)/y(1);
dydt(4) = y(3)/y(1);

