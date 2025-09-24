function [value, isterminal, direction] = soi_event(t, x)
    soi_m = 66193923.1913185; % Moon's SOI radius (km)
    r = x(1);  % Radial distance from Moon
    value = r - soi_m;  % Stop when r >= soi_moon
    isterminal = 1;        % Terminate integration
    direction = 1;         % Detect positive crossing (exiting SOI)
end