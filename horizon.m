function x_mpc = horizon(model, x, contact_schedule, k, p1x, p1z, p2x, p2z)
% Horizon of 20, discretization of .02s = 10 * .002
% x is 7 long vector, x,y,theta, zdot,ydot,thetadot, g

m = 9;
I = .0638;
Fmin = 10;
Fmax = 250;
alphas = contact_schedule(k:k+200,:);
alphas = alphas(1:10:end,:);

Ac = [0 0 0 1 0 0 0;
    0 0 0 0 1 0 0;
    0 0 0 0 0 1 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 1;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 ];
Ad = Ac + .02*eye(7);

B_horizon = cell(20,1);
for i = 1:length(alphas) % hopefully 1:20
    Bc = [0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        alphas(i, 1)/m, 0, alphas(i, 2)/m, 0;
        0, alphas(i, 1)/m, 0, alphas(i, 2)/m;
        -p1z*alphas(i,1)/I, p1x*alphas(i,1)/I, -p2z*alphas(i,2)/I, p2x*alphas(i,2)/I;
        0, 0, 0, 0];
    Bd = Bc*.02;
    B_horizon{i} = Bd;
end

%% Calculate Torso Trajectories Now
vxd = .5;
x_traj = x(1)*ones(1,20)+linspace(.2,4,20)*vxd;
y_traj = .55*ones(1,20);
theta_traj = 0*ones(1,20);
xd_traj = .5*ones(1,20);
yd_traj = 0*ones(1,20);
thetad_traj = 0*ones(1,20); 
g_traj = 9.81*ones(1,20);
x_traj = [x_traj; y_traj; theta_traj; xd_traj; yd_traj; thetad_traj; g_traj];

%% h and f' Now
Q = diag([1, 1, 1, 1, 1, 1, 0]);
R =  diag([1, 1, 1, 1]);
h = 2*blkdiag(Q, Q, Q, Q, Q, Q, Q, Q, Q ,Q ,Q ,Q ,Q, Q, Q, Q, Q, Q, Q, Q, ...
    R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R,R);
ft = [];
for i = 1:length(alphas)
    ft = [ft; x_traj(:,i)'*Q];
end
ft = -2*ft;
f = ft';

%% Aeq and beq Now

I20 = eye(20);              
E20 = diag(ones(19,1), -1); 
Aeq_left = kron(I20, eye(7)) + kron(E20, -Ad);
C = cellfun(@(X) -X, B_horizon, 'UniformOutput', false);
Aeq_right = blkdiag(C{:});
Aeq = [Aeq_left Aeq_right];
beq = [Ad*x;zeros(133,1)];

%% Inequalities Aqp bqp now
Aqp = [];
bqp = [];
for i = 1:length(alphas) % hopefully 1:20
    if alphas(i,1) == 0 && alphas(i, 2) == 1
        Ai_right = [1 0 0 0;
            -1 0 0 0;
            0 1 0 0;
            0 -1 0 0;
            0 0 0 1;
            0 0 0 -1;
            0 0 1 -mu;
            0 0 -1 -mu];
        
    elseif alphas(i,1) == 1 && alphas(i,2) == 0
        Ai_right = [0 0 1 0;
            0 0 -1 0;
            0 0 0 1;
            0 0 0 -1;
            0 1 0 0;
            0 -1 0 0;
            1 -mu 0 0;
            -1 -mu 0 0;];
    end
    Bi = [0 0 0 0 Fmax -Fmin 0 0];
    Aqp = [Aqp; zeros(8,7) Ai_right];
    bqp = [bqp; Bi];
end

%% Finally call the QP function
x_mpc = quadprog(h, f, Aqp, bqp, Aeq, beq);
