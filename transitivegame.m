function transitivegame()
    % Parameters
a=1.0;
b=1.5;
c=1.2;
sig1=1.5;
sig2=1.7;
sig3=1.25;
alpha=0.15;
beta=0.4;
xi=0.25;
p=0.0;
  %  p = 0.5;
    
    dt = 0.01; % Time step
    tspan = 0:dt:50000; % Time span
    %fprintf('%d\n', tspan); 
    % Initial conditions
    
    x0 = 0.3;
    y0 = 0.3;
    z0 = 0.3;

    % Numerical integration using the Runge-Kutta45
    [~, xyz] = ode45(@(t, xyz) transitivegameequations(xyz, a, b, c, sig1, sig2, sig3, alpha, beta, xi, p), tspan, [x0, y0, z0]);

    x = xyz(:, 1);
    y = xyz(:, 2);
    z = xyz(:, 3);

    % Plotting x,y,z without transcient
    lastSteps = 10000;
    if length(tspan) > lastSteps
        timeIdx = length(tspan) - lastSteps + 1 : length(tspan);
    else
        timeIdx = 1 : length(tspan);
    end

    figure;
    plot(tspan(timeIdx), x(timeIdx), 'b', 'DisplayName', 'x');
    hold on;
    plot(tspan(timeIdx), y(timeIdx), 'r', 'DisplayName', 'y');
    plot(tspan(timeIdx), z(timeIdx), 'g', 'DisplayName', 'z');
    xlabel('Time');
    ylabel('x,y,z');
    legend show;
  % figure;
   %  plot3(x(timeIdx), y(timeIdx), z(timeTdx));
%     figure;
%     plot(y(timeIdx), z(timeIdx));
%     figure;
%     plot(x(timeIdx),z(timeIdx));
    %title('Transitive Game Dynamics - Last 2000 Time Steps');
   % grid on;
end

function dxdt = transitivegameequations(xyz, a, b, c, sig1, sig2, sig3, alpha, beta, xi, p)
    % Extract variables
    x = xyz(1);
    y = xyz(2);
    z = xyz(3);

    % Eco-evolutionary equations
    dxdt = [x * (-(sig1 * x) + (a - sig1) * y + (b - 2.0 * p * b - sig1) * z + (sig1 - alpha));
            y * (-(a + sig2) * x - (sig2 * y) + (c - sig2) * z + (sig2 - beta));
            z * ((2.0 * p * b - b - sig3) * x - (c + sig3) * y - (sig3 * z) + (sig3 - xi))];
end