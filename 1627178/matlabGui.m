function heat_eq_gui()
    % Create the main figure
    hFig = figure('Name', 'Heat Equation Solver', 'NumberTitle', 'off', ...
                  'Position', [100, 100, 500, 300], 'Resize', 'off');

    % Add input fields
    uicontrol('Style', 'text', 'Position', [10, 260, 180, 20], 'String', 'Diffusion Coefficient (c):');
    hEditC = uicontrol('Style', 'edit', 'Position', [200, 260, 100, 20], 'String', '1');

    uicontrol('Style', 'text', 'Position', [10, 230, 180, 20], 'String', 'Total Time (T):');
    hEditT = uicontrol('Style', 'edit', 'Position', [200, 230, 100, 20], 'String', '1');

    uicontrol('Style', 'text', 'Position', [10, 200, 180, 20], 'String', 'Spatial Points (Nx):');
    hEditNx = uicontrol('Style', 'edit', 'Position', [200, 200, 100, 20], 'String', '50');

    uicontrol('Style', 'text', 'Position', [10, 170, 180, 20], 'String', 'Time Steps (Nt):');
    hEditNt = uicontrol('Style', 'edit', 'Position', [200, 170, 100, 20], 'String', '1000');

    % Add buttons for methods
    hButtonMethod1 = uicontrol('Style', 'pushbutton', 'Position', [10, 100, 150, 30], ...
                               'String', 'Explicit Method', 'Callback', {@solve_method, 1});
    hButtonMethod2 = uicontrol('Style', 'pushbutton', 'Position', [170, 100, 150, 30], ...
                               'String', 'Implicit Method', 'Callback', {@solve_method, 2});
    hButtonMethod3 = uicontrol('Style', 'pushbutton', 'Position', [330, 100, 150, 30], ...
                               'String', 'Crank-Nicolson', 'Callback', {@solve_method, 3});

    % Nested function for button callbacks
    function solve_method(hObject, eventData, methodID)
        c = str2double(get(hEditC, 'String'));
        T = str2double(get(hEditT, 'String'));
        Nx = str2double(get(hEditNx, 'String'));
        Nt = str2double(get(hEditNt, 'String'));
        dx = 1 / (Nx - 1);
        dt = T / Nt;
        x = linspace(0, 1, Nx);
        t = linspace(0, T, Nt);
        u = sin(pi * x); % Initial condition
        
        switch methodID
            case 1
                result = explicitMethod(c, dx, dt, Nx, Nt, u);
            case 2
                result = implicitMethod(c, dx, dt, Nx, Nt, u);
            case 3
                result = crankNicolsonMethod(c, dx, dt, Nx, Nt, u);
        end
        
        % Plot results
        figure;
        plot(x, result);
        xlabel('x');
        ylabel('Temperature');
        title(sprintf('Solution at T=%.2f using Method %d', T, methodID));
    end
end

function u_new = explicitMethod(c, dx, dt, Nx, Nt, u)
    r = c * dt / dx^2;
    for n = 1:Nt
        u_new = u;
        for i = 2:Nx-1
            u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
        end
        u = u_new;  % Update the old temperature profile for the next time step
    end
end

function u_new = implicitMethod(c, dx, dt, Nx, Nt, u)
    r = c * dt / dx^2;
    A = (1 + 2*r) * eye(Nx) + diag(-r * ones(1, Nx-1), 1) + diag(-r * ones(1, Nx-1), -1);
    A(1,1) = 1; A(1,2) = 0; A(Nx,Nx) = 1; A(Nx,Nx-1) = 0;  % Boundary conditions
    for n = 1:Nt
        u = A \ u;  % Solve the linear system
    end
    u_new = u;
end

function u_new = crankNicolsonMethod(c, dx, dt, Nx, Nt, u)
    r = c * dt / (2 * dx^2);
    A = (1 + r) * eye(Nx) + diag(-0.5*r * ones(1, Nx-1), 1) + diag(-0.5*r * ones(1, Nx-1), -1);
    B = (1 - r) * eye(Nx) + diag(0.5*r * ones(1, Nx-1), 1) + diag(0.5*r * ones(1, Nx-1), -1);
    A(1,1) = 1; A(1,2) = 0; A(Nx,Nx) = 1; A(Nx,Nx-1) = 0;  % Boundary conditions
    B(1,1) = 1; B(1,2) = 0; B(Nx,Nx) = 1; B(Nx,Nx-1) = 0;  % Boundary conditions
    for n = 1:Nt
        u = A \ (B * u);  % Solve the linear system
    end
    u_new = u;
end