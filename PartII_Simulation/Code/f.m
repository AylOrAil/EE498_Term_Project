function dx = system_rhs(x)
    % State vector: x = [n; a; Z; S; C]

persistent q pi0 delta alpha0 beta1 theta1 alpha1 beta10 beta2 theta2 alpha2 beta20 pi1 pi2 qc alpha3 lambda11 delta1

if isempty(q)
s = coder.load('sim_data.mat',  'q', 'pi0', 'delta', 'alpha0', 'beta1', ...
    'theta1', 'alpha1', 'beta10', 'beta2', 'theta2', 'alpha2', ...
    'beta20', 'pi1', 'pi2', 'qc', 'alpha3', 'lambda11', 'delta1');

q = s.q; pi0 = s.pi0; delta = s.delta; alpha0 = s.alpha0; beta1 = s.beta1;
theta1 = s.theta1; alpha1 = s.alpha1; beta10 = s.beta10; beta2 = s.beta2; theta2 = s.theta2; alpha2 = s.alpha2;
beta20 = s.beta20; pi1 = s.pi1; pi2 = s.pi2; qc = s.qc; alpha3 = s.alpha3; lambda11 = s.lambda11; delta1 = s.delta1;
end

    n = x(1);
    a = x(2);
    Z = x(3);
    S = x(4);
    C = x(5);

    dn = q+ pi0*delta*S-alpha0*n-beta1*n*a;
    da = theta1*beta1*n*a - alpha1*a-beta10*a^2-beta2*a*Z;
    dZ = theta2*beta2*a*Z - alpha2*Z-beta20*Z^2;
    dS = pi1*alpha1*a+pi2*alpha2*Z-delta*S;
    dC = qc - alpha3*C-lambda11*a-delta1*S;
    dx = [dn; da; dZ; dS; dC];
end