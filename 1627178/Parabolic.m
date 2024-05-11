function u = heat_equation_explicit(u0, alpha, L, T, Nx, Nt)
    % u0: Initial temperature distribution, a vector of length Nx
    % alpha: Thermal diffusivity
    % L: Length of the rod
    % T: Total time for simulation
    % Nx: Number of spatial points
    % Nt: Number of time steps

    % Spatial and temporal grid sizes
    dx = L / (Nx - 1);
    dt = T / Nt;

    % Stability check
    s = alpha * dt / dx^2;
    if s > 0.5
        error('The scheme is unstable! Reduce dt or increase dx.');
    end

    % Initialize temperature matrix
    u = zeros(Nx, Nt);
    u(:, 1) = u0;  % Set initial condition

    % Time-stepping loop
    for n = 1:Nt-1
        for i = 2:Nx-1
            u(i, n+1) = u(i, n) + s * (u(i+1, n) - 2*u(i, n) + u(i-1, n));
        end
        % Boundary conditions (fixed in this example)
        u(1, n+1) = u(1, n);  % Left boundary
        u(Nx, n+1) = u(Nx, n);  % Right boundary
    end
end

% Example usage
Nx = 50;
Nt = 1000;
L = 10;
T = 1;
alpha = 0.01;
u0 = sin(pi * linspace(0, L, Nx));  % Initial condition

u = heat_equation_explicit(u0, alpha, L, T, Nx, Nt);

% Plot the final temperature distribution
plot(linspace(0, L, Nx), u(:, end));
xlabel('Position');
ylabel('Temperature');
title('Temperature distribution at final time');