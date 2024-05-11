% Define the system of equations
function dydt = sys_ode(t, y)
    dydt = zeros(2,1);    % Initialize the output derivative vector
    dydt(1) = 3*y(1) + 4*y(2);  % Equation for dy/dt
    dydt(2) = -4*y(1) + 3*y(2); % Equation for dz/dt
end

% Initial conditions
y0 = [1; 0];  % y(0) = 1, z(0) = 0

% Time span
tspan = linspace(0, 10, 100);  % From t = 0 to t = 10 with 100 points

% Solve the ODE using lsode
[t, y] = lsode(@sys_ode, y0, tspan);

% Plotting the results
figure;
plot(t, y(:, 1), 'r-', t, y(:, 2), 'b--');
legend('y(t)', 'z(t)');
xlabel('Time t');
ylabel('Solutions y and z');
title('Solutions of the system of ODEs');
grid on;