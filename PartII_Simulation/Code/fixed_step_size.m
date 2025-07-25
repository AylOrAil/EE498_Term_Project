%% implementation of Report's var step size method 

x0 = [0.1; 0.1; 0.1; 0.1; 0.1];

t_start = 0;
t_end = 20;

hi = 1e-1;
ti = t_start:hi:t_end;

X = zeros(5, length(ti));
X(:,1) = x0;

tic 
for i = 1:length(ti)-1
    X(:,i+1) = runge4(hi, X(:,i));
end

elapsedTime = toc;  
fprintf('Elapsed time: %.4f seconds\n', elapsedTime);

t_2 = ti;
X_2 = X;

%% Plots 
figure;
plot(ti, X(1,:));
xlabel('Time');
ylabel('n(t)');
title(['Runge-Kutta 4th Order for Step Size ', num2str(hi)]);
grid on;

%% functions 
function xip1 = runge4(hi, xi)
    k1 = (xi);
    k2 = f(xi + 0.5 * hi * k1);
    k3 = f(xi + 0.5 * hi * k2);
    k4 = f(xi + hi * k3);

    phi = (1/6)*k1+(2/6)*k2+(2/6)*k3+(1/6)*k4;
    xip1 = xi+hi*phi;    
end