function transitivegame()
    % Parameters
    a = 0.60;
    b=0.10;
    c = 0.40;
    %sig1 = 0.75;
    sig2 = 0.8;
    sig3 = 0.90;
    alpha = 0.35;
    beta = 0.20;
    xi = 0.08;
    p = 0.5 ;
    
    dt = 0.01; % Time step
    tspan = 0:dt:70000; % Time span

    % Initial conditions
    x0 = 0.3;
    y0 = 0.3;
    z0 = 0.3;

    fp1 = fopen('bifurx1.dat', 'w');
    fp2 = fopen('bifury1.dat', 'w');
    fp3 = fopen('bifurz1.dat', 'w');
    fp4 = fopen('bifurx11.dat', 'w');
    fp5 = fopen('bifury11.dat', 'w');
    fp6 = fopen('bifurz11.dat', 'w');


    for j = 0:100
        sig1 = 0.0 + (2.0 / 100.0) * j;
        fprintf('%d\t%f\n', j, sig1);

        % Numerical integration using the Runge-Kutta45
        [~, xyz] = ode45(@(t, xyz) transitivegameequations(t, xyz, a, b, c, sig1, sig2, sig3, alpha, beta, xi, p), tspan, [x0, y0, z0]);

        x = xyz(:, 1);
        y = xyz(:, 2);
        z = xyz(:, 3);

        % Plotting x,y,z without transient
        lastSteps = 10000;
        if length(x) > lastSteps
            xx1 = x(end-lastSteps:end);
            yy1 = y(end-lastSteps:end);
            zz1 = z(end-lastSteps:end);

            if max(xx1) >= xx1(round(end/2)) && max(xx1) >= xx1(end)
                fprintf(fp1, '%f\t%f\n', sig1, max(xx1));
            end
            if max(yy1) >= yy1(round(end/2)) && max(yy1) >= yy1(end)
                fprintf(fp2, '%f\t%f\n', sig1, max(yy1));
            end
            if max(zz1) >= zz1(round(end/2)) && max(zz1) >= zz1(end)
                fprintf(fp3, '%f\t%f\n', sig1, max(zz1));
            end
        end
        if length(x) > lastSteps
            xx1 = x(end-lastSteps:end);
            yy1 = y(end-lastSteps:end);
            zz1 = z(end-lastSteps:end);

            if min(xx1) <= xx1(round(end/2)) && min(xx1) <= xx1(end)
                fprintf(fp4, '%f\t%f\n', sig1, min(xx1));
            end
            if min(yy1) <= yy1(round(end/2)) && min(yy1) <= yy1(end)
                fprintf(fp5, '%f\t%f\n', sig1, min(yy1));
            end
            if min(zz1) <= zz1(round(end/2)) && min(zz1) <= zz1(end)
                fprintf(fp6, '%f\t%f\n', sig1, min(zz1));
            end
        end
        
    end

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    fclose(fp6);
    
p1=load('bifurx1.dat');
p2=load('bifury1.dat');
p3=load('bifurz1.dat');
p4=load('bifurx11.dat');
p5=load('bifury11.dat');
p6=load('bifurz11.dat');

p11=unique(p1,'rows');
p22=unique(p2,'rows');
p33=unique(p3,'rows');
p44=unique(p4,'rows');
p55=unique(p5,'rows');
p66=unique(p6,'rows');

figure;
for i = 1:length(p11)-1
    if(p11(i,1)==p11(i+1,1))
        p11(i,2) = max(p11(i,2),p11(i+1,2));
        p11(i+1,2) = max(p11(i,2),p11(i+1,2));
    end
end
p11 = unique(p11,'rows');
plot(p11(:,1),p11(:,2))
hold on

for i = 1:length(p22)-1
    if(p22(i,1)==p22(i+1,1))
        p22(i,2) = max(p22(i,2),p22(i+1,2));
        p22(i+1,2) = max(p22(i,2),p22(i+1,2));
    end
end
p22 = unique(p22,'rows');
plot(p22(:,1),p22(:,2))

for i = 1:length(p33)-1
    if(p33(i,1)==p33(i+1,1))
        p33(i,2) = max(p33(i,2),p33(i+1,2));
        p33(i+1,2) = max(p33(i,2),p33(i+1,2));
    end
end
p33 = unique(p33,'rows');
plot(p33(:,1),p33(:,2))

for i = 1:length(p44)-1
    if(p44(i,1)==p44(i+1,1))
        p44(i,2) = min(p44(i,2),p44(i+1,2));
        p44(i+1,2) = min(p44(i,2),p44(i+1,2));
    end
end
p44 = unique(p44,'rows');
plot(p44(:,1),p44(:,2))

for i = 1:length(p55)-1
    if(p55(i,1)==p55(i+1,1))
        p55(i,2) = min(p55(i,2),p55(i+1,2));
        p55(i+1,2) = min(p55(i,2),p55(i+1,2));
    end
end
p55 = unique(p55,'rows');
plot(p55(:,1),p55(:,2))

for i = 1:length(p66)-1
    if(p66(i,1)==p66(i+1,1))
         p66(i,2) = min(p66(i,2),p66(i+1,2));
        p66(i+1,2) = min(p66(i,2),p66(i+1,2));
    end
end
p66 = unique(p66,'rows');
plot(p66(:,1),p66(:,2))

end

function plotUniqueData(data)
    for i = 1:size(data, 1)-1
        if data(i, 1) == data(i+1, 1)
            data(i, 2) = max(data(i, 2), data(i+1, 2));
            data(i+1, 2) = data(i, 2);
        end
    end
    data = unique(data, 'rows');
    plot(data(:, 1), data(:, 2));
end

function dxdt = transitivegameequations(~, xyz, a, b, c, sig1, sig2, sig3, alpha, beta, xi, p)
    % Extract variables
    x = xyz(1);
    y = xyz(2);
    z = xyz(3);

    % Eco-evolutionary equations
    dxdt = [x * (-(sig1 * x) + (a - sig1) * y + (b - 2.0 * p * b - sig1) * z + (sig1 - alpha));
            y * (-(a + sig2) * x - (sig2 * y) + (c - sig2) * z + (sig2 - beta));
            z * ((2.0 * p * b - b - sig3) * x - (c + sig3) * y - (sig3 * z) + (sig3 - xi))];
end
