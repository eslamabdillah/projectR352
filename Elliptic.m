function phi = laplace_jacobi(N, maxIter, tol)
    % N Number of points along each dimension (including boundaries)
    % maxIter Maximum number of iterations
    % tol Tolerance for stopping criterion based on max change in the grid

    % Initialize the potential array
    phi = zeros(N, N);

    % Set boundary conditions (example phi = 0 on the boundaries)
    phi(1, ) = 1;   % Top boundary
    phi(N, ) = 1;   % Bottom boundary
    phi(, 1) = 1;   % Left boundary
    phi(, N) = 1;   % Right boundary

    % Iterative solver
    for iter = 1maxIter
        phi_old = phi;  % Copy of the old phi for convergence check

        % Update the potential at each interior point
        for i = 2N-1
            for j = 2N-1
                phi(i, j) = 0.25  (phi_old(i+1, j) + phi_old(i-1, j) + ...
                                     phi_old(i, j+1) + phi_old(i, j-1));
            end
        end

        % Check for convergence
        if max(max(abs(phi - phi_old)))  tol
            fprintf('Convergence achieved after %d iterations.n', iter);
            break;
        end
    end

    if iter == maxIter
        fprintf('Max iterations reached without convergence.n');
    end
end

% Example usage
N = 50;         % Grid size
maxIter = 1000; % Maximum iterations
tol = 1e-5;     % Tolerance
phi_result = laplace_jacobi(N, maxIter, tol);

% Visualize the result
surf(phi_result);
title('Potential distribution (phi)');
xlabel('x');
ylabel('y');
zlabel('phi');