% used for storing params simulink results 

q = 5;
pi0 = 0.1; 
delta = 1;
alpha0 = 0.5;
beta1 = 1;
theta1 = 1;
alpha1 = 0.5;
beta10 = 2;
beta2 = 1;
theta2 = 1;
alpha2 = 0.5;
beta20 = 2;
pi1 = 0.9;
pi2 = 0.9;
qc = 10;
alpha3 = 2;
lambda11 = 0.25;
delta1 = 10;

save('sim_data.mat', 'q', 'pi0', 'delta', 'alpha0', ...
    'beta1', ...
    'theta1', 'alpha1', 'beta10', 'beta2', 'theta2', ...
    'alpha2', ...
    'beta20', 'pi1', 'pi2', 'qc', 'alpha3', ...
    'lambda11', 'delta1');