function [u, x, t] = solve_wave_equation(c, L, T, Nx, Nt, u0, v0)
    % c: Wave speed
    % L: Length of the domain
    % T: Total time of simulation
    % Nx: Number of spatial points
    % Nt: Number of time steps
    % u0: Initial displacement function handle
    % v0: Initial velocity function handle

    % Spatial and temporal grid sizes
    dx = L / (Nx - 1);
    dt = T / Nt;

    % Check the CFL condition for stability
    if (c * dt / dx) > 1
        error('CFL condition not met. Increase Nx or decrease dt.');
    end
    
    % Initialize grids
    x = linspace(0, L, Nx);
    t = linspace(0, T, Nt);
    u = zeros(Nx, Nt);

    % Set initial conditions
    u(:, 1) = u0(x);   % Initial displacement
    % Initial velocity converted to first time step using central difference
    u(:, 2) = u(:, 1) + v0(x) * dt + 0.5 * (c^2 * dt^2 / dx^2) * ...
              ([u(2, 1); u(3:Nx, 1) - 2*u(2:Nx-1, 1) + u(1:Nx-2, 1); u(Nx-1, 1)] - 2 * u(:, 1) + [u(2, 1); u(3:Nx, 1); 0]);

    % Finite difference solution
    for n = 2:Nt-1
        u(2:Nx-1, n+1) = 2 * u(2:Nx-1, n) - u(2:Nx-1, n-1) + ...
                         (c^2 * dt^2 / dx^2) * (u(3:Nx, n) - 2 * u(2:Nx-1, n) + u(1:Nx-2, n));
    end
end

% Example usage
L = 10;      % Length of the spatial domain
T = 10;      % Total time
Nx = 100;    % Number of spatial points
Nt = 1000;   % Number of time steps
c = 1;       % Speed of wave propagation

% Initial conditions
u0 = @(x) sin(pi*x/L);  % Initial displacement
v0 = @(x) 0;            % Initial velocity

[u, x, t] = solve_wave_equation(c, L, T, Nx, Nt, u0, v0);

% Plot the results
mesh(t, x, u);
xlabel('Time');
ylabel('Position');
zlabel('Displacement');
title('Wave Equation Solution');