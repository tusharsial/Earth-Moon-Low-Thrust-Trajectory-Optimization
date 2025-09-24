% Clear workspace and close figures
clc; clear; close all;

% Define Earth radius (for scaling)
R_earth = 6371; % km

% Generate sample trajectory data
% Thrust arc (e.g., powered flight)
theta_thrust = linspace(0, pi/2, 100); % Angular range for thrust arc
x_thrust = -70 + 60 * cos(theta_thrust); % X-coordinates (Earth radii)
y_thrust = -15 + 30 * sin(theta_thrust); % Y-coordinates (Earth radii)

% Coast arc (e.g., ballistic flight)
theta_coast = linspace(pi/2, pi, 100); % Angular range for coast arc
x_coast = -10 + 10 * cos(theta_coast); % X-coordinates (Earth radii)
y_coast = 15 + 5 * sin(theta_coast);   % Y-coordinates (Earth radii)

% Combine trajectories
x_total = [x_thrust, x_coast];
y_total = [y_thrust, y_coast];

% Plot the trajectory
figure;
plot(x_thrust, y_thrust, 'r-', 'LineWidth', 2); % Thrust arc (red)
hold on;
plot(x_coast, y_coast, 'b--', 'LineWidth', 2);  % Coast arc (blue dashed)
grid on;
axis equal;
xlabel('X, Earth radii');
ylabel('Y, Earth radii');
title('Spacecraft Trajectory: Thrust Arc and Coast Arc');
legend('Thrust Arc', 'Coast Arc', 'Location', 'best');

% Add Earth and Moon markers (approximate positions)
plot(0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Earth at origin
plot(-60, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); % Moon (scaled position)

% Adjust axes to match the reference image
xlim([-70, 10]);
ylim([-15, 15]);