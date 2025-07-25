%% implementation of Report's var step size method 

x0 = [0.1; 0.1; 0.1; 0.1; 0.1];

t_start = 0;
t_end = 20;

atol = 1e-1;
rtol = 1e-2;

% default h0 choice based on user defined atol and rtol
% h0 = 0.1 * (atol+ rtol*max(x0));

% user defined constant 
h0 = 0.95;

t = [];

X = [];
X(:,1) = x0;

t = [];
h_vals = [];
t(1) = t_start;

ti = t_start;
i = 1;
tic      

while ti< t_end
    [X(:,i+1), h] = runge34(h0, X(:,i), atol, rtol);
    
    ti = ti+h;
    t(i+1) = ti;
    h_vals(i) = h;

    i = i+1;
end
elapsedTime = toc;  
t_4 = t; 
x_4 = X;
fprintf('Elapsed time: %.4f seconds\n', elapsedTime);

%% Plots 
figure;
plot(t, X(1, :));
xlabel('Time');
ylabel('n(t)');
title('Runge-Kutta 4th Order for Various Step Sizes');
grid on;

epsilon = 0.1*h0;
figure;
plot(h_vals);
xlabel('Simulation Time Step');
ylabel('Final Simulation Step Size');
title('Step Size values over Simulation Time Steps');
ylim([-epsilon, h0 + epsilon]);
grid on;
%% functions 
function xip1 = runge4(hi, xi)
    k1 = f(xi);
    k2 = f(xi + 0.5 * hi * k1);
    k3 = f(xi + 0.5 * hi * k2);
    k4 = f(xi + hi * k3);

    phi = (1/6)*k1+(2/6)*k2+(2/6)*k3+(1/6)*k4;
    xip1 = xi+hi*phi;    
end 

function [xf, hf] = runge34(h0, xi, atol, rtol)
    E = 2;
    hi = h0;

    while E>1
        % p+1 method: RK4
            k1 = f(xi);
            k2 = f(xi + 0.5 * hi * k1);
            k3 = f(xi + 0.5 * hi * k2);
            k4 = f(xi + hi * k3);
        
            phi = (1/6)*k1+(2/6)*k2+(2/6)*k3+(1/6)*k4;
            x_4 = xi+hi*phi;  

            k1 = f(xi);
            k2 = f(xi + 1/3 * hi * k1);
            k3 = f(xi + 2/3 * hi * k2);

            phi = (1/4)*k1+(3/4)*k3;
            x_3 = xi+hi*phi;  

            err = abs(x_4-x_3);
            tol = atol+rtol*abs(xi);
            
            E = err ./ tol;

            E = norm(E) / sqrt(length(E));

            if E<=1
                xf = x_4;
                hf = hi;
                break
            end

            hi = hi*min(5, max(0.1, (1/E)^(1/4)));
    end 
end