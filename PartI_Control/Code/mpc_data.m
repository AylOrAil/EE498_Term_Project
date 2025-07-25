% system state space parameters

T = 0.2;
th = 1.5;
tau = 0.5;
d0 = 5;

% x:= [d v_acc v_rel a j]
% u:= torque 
% y:= [delta_x v_rel a j]  

% x(k+1) = Ax(k) + Bu(k) 
% y(k) = Cx(k) + F


A = [1 0 T -0.5*T^2 0;
     0 1 0 T 0;
     0 0 1 -T 0; 
     0 0 0 1-T/tau 0;
     0 0 0 -1/tau 0];

B = [0; 0; 0; T/tau; 1/tau];
C = [1 -th 0 0  0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];


E = [0.5*T^2; 
    0;
    T;
    0;
    0];

% constraints to be considered 

% state constraints: 
% d0 <= d 
% v_min <= v <= v_max 
% a_min <= a <= a_max
% j_min <= j <= j_max

% input constraints:
% u_min <= u <= u_max
v_min = 0;
v_max = 30;
a_min = -5;
a_max = 2;
j_min = -5;
j_max = 2;
u_min = -5.5;
u_max = 2.5;


% defining N= Np= Nc, G and H

N = 10;
z = zeros(size(B));
G = [z z z z z z z z z z ; 
    B z z z z z z z z z;  
    A*B B z z z z z z z z;  
    A^2*B A*B B  z z z z z z z;  
    A^3*B A^2*B A*B B z z z z z z;  
    A^4*B A^3*B A^2*B A*B B  z z z z z; 
    A^5*B A^4*B A^3*B A^2*B A*B B  z z z z;  
    A^6*B A^5*B A^4*B A^3*B A^2*B A*B B z z z;  
    A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B  z z;  
    A^8*B A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B  z;  
    A^9*B A^8*B A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B]; % (11x5) x 10
H = [eye(size(A));A; A^2; A^3; A^4; A^5; A^6; A^7; A^8; A^9; A^10]; % (11x5) x 5


c1 = [1 -th 0 0  0];
C1 = kron(eye(N+1),c1);

s2 = [0 1 0 0 0];
S2 = kron(eye(N+1), s2);

s4 = [0 0 0 1 0];
S4 = kron(eye(N+1), s4);

s5 = [0 0 0 0 1];
S5 = kron(eye(N+1), s5); % (11) x (11 x 5)


% defining input constraints F1 U <= V1

F1 = [eye(N); -eye(N)]; % (2x10) x 10
V1 =  [u_max * ones(N, 1); -u_min * ones(N,1)]; % (2 x 10) x 1

% state constraints will have the form S*G U <= V' - SHx_{0|k} = V -> V online 
F2 = -C1*G;
V2_p = -(-4)*ones(N+1, 1);

F3 = [S2; -S2]*G;
V3_p = [v_max*ones(N+1,1); -v_min*ones(N+1,1)];

F4 = [S4; -S4]*G;
V4_p = [a_max*ones(N+1,1); -a_min*ones(N+1,1)];

F5 = [S5; -S5]*G;
V5_p = [j_max*ones(N+1,1); -j_min*ones(N+1,1)];


% tuning Q and R

%Q = blkdiag(0, 1000, 1000, 0, 0);
Q =blkdiag(5,10,1,1);
R = 0.001;
Qf = Q; 

Qbar =  blkdiag(Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Qf);
Rbar = blkdiag(R, R, R, R, R, R, R, R, R, R);
Cbar = blkdiag(C, C, C, C, C, C, C, C, C, C, C);
    
%M = G'*Qbar*G + Rbar;

x0 = [100; 30; -10; 0; 0];

x_k = [];
u_k = [];
y_k = [];
i_k = [];

x_0_k = x0;
t_start = 0;
t_end = 40; 
Ts = T;


k_max = (t_end-t_start)/Ts-1;

%% saving weights and constant parameters 
save('mpc_data.mat', 'G', 'H', 'C1', 'Qbar', 'Rbar', ...
     'F1', 'F2', 'F3', 'F4', 'F5', ...
     'Cbar', 'S2', 'S4', 'S5', ...
     'V1', 'V2_p', 'V3_p', 'V4_p', 'V5_p', 'd0');