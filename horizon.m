function x_mpc = horizon(model, x, contact_schedule, k, p1x, p1z, p2x, p2z)
% Horizon of 20, discretization of .02s = 10 * .002
% x is 7 long vector, x,y,theta, zdot,ydot,thetadot, g

m = 9;
I = .0638;
Fmin = 10;
Fmax = 300;
mu = .7;
alphas = contact_schedule(k:10:k+190,:); % Extract the 20 horizons out of 200 contact_schedule
Ac = [0 0 0 1 0 0 0; % x
    0 0 0 0 1 0 0; % y
    0 0 0 0 0 1 0; % theta
    0 0 0 0 0 0 0; % x
    0 0 0 0 0 0 -1; % y
    0 0 0 0 0 0 0; % theta
    0 0 0 0 0 0 0 ]; % g
Ad = Ac*.02 + eye(7);

B_horizon = cell(20,1);
% for i = 1:length(alphas) % hopefully 1:20
%     Bc = [0 0 0 0;
%         0 0 0 0;
%         0 0 0 0;
%         alphas(i, 1)/m, 0, alphas(i, 2)/m, 0;
%         0, alphas(i, 1)/m, 0, alphas(i, 2)/m;
%         -p1z*alphas(i,1)/I, p1x*alphas(i,1)/I, -p2z*alphas(i,2)/I, p2x*alphas(i,2)/I;
%         0, 0, 0, 0];
%     Bd = Bc*.02;
%     B_horizon{i} = Bd;
% end

for i = 1:length(alphas) % hopefully 1:20
    Bc = [0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        1/m, 0, 1/m, 0;
        0, 1/m, 0, 1/m;
        p1z/I, -p1x/I, p2z/I, -p2x/I;
        0, 0, 0, 0];
    Bd = Bc*.02;
    %Bd(3,:) = Bd(3,:) - [.003627, 0, .003627, 0]; % c2d gives THIS for the discrete B matrix...
    B_horizon{i} = Bd;
end

% Build Ap and Bp now (these are NOT the inequality matrices)

N = 20;
n = size(Ad,1);
m = size(B_horizon{1},2);

M = zeros(N*n, N*m);

for j = 1:N
    Apow = eye(n);                 % Apow = Ad^(i-j)
    for i = j:N
        rows = (i-1)*n + (1:n);
        cols = (j-1)*m + (1:m);
        M(rows, cols) = Apow * B_horizon{j};
        Apow = Ad * Apow;          % next power
    end
end
Bp = M;
Ap = cell(20,1);
for i = 1:20
    Ap{i} = Ad^i;
end
Ap = cell2mat(Ap);

%% Calculate Torso Trajectories Now
vxd = 0;
x_traj = zeros(1,20);
y_traj = .45*ones(1,20);
theta_traj = ones(1,20)*10/180*pi;
xd_traj = zeros(1,20);
yd_traj = zeros(1,20);
thetad_traj = zeros(1,20); 
g_traj = 9.81*ones(1,20);
x_traj = [x_traj; y_traj; theta_traj; xd_traj; yd_traj; thetad_traj; g_traj];
x_ref = x_traj(:);

% vxd = .5;
% x_traj = x(1)*ones(1,20)+linspace(.2,4,20)*x(4);
% y_traj = .55*ones(1,20);
% theta_traj = zeros(1,20);
% xd_traj = linspace(x(4),vxd,20);
% yd_traj = zeros(1,20);
% thetad_traj = zeros(1,20); 
% g_traj = 9.81*ones(1,20);
% x_traj = [x_traj; y_traj; theta_traj; xd_traj; yd_traj; thetad_traj; g_traj];
% x_ref = x_traj(:);

%% h and f' Now
Q = diag([1, 10, 3, 1, 1, 1, 0]);
R =  diag([.0001, .0001, .0001, .0001]);
L = blkdiag(Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q);
K = blkdiag(R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R);
h = 2*Bp'*L*Bp+K;
ft = 2*x'*Ap'*L*Bp - 2*x_ref'*L*Bp;
f = ft';

%% Inequalities Constraints Aqp bqp. These are essentially friction pyramid and Fmin < Fy < Fmax constraints
% Aqp = [];
% bqp = [];
% Aqp_k = cell(1,20);
% for i = 1:length(alphas) % hopefully 1:20
%     % if alphas(i,1) == 0 && alphas(i, 2) == 1
%     %     Ai_right = [1 0 0 0;
%     %         -1 0 0 0;
%     %         0 1 0 0;
%     %         0 -1 0 0;
%     %         0 0 0 1;
%     %         0 0 0 -1;
%     %         0 0 1 -mu;
%     %         0 0 -1 -mu];
%     %     bi = [0; 0; 0; 0; Fmax; -Fmin; 0; 0]; % First four are saying F1 = 0 or F2 = 0. Can relax those a little bit. Next 2 is saying Fmin < Fy < Fmax. Last two are friction constraints.
%     %     bqp = [bqp; bi];
%     % 
%     % elseif alphas(i,1) == 1 && alphas(i,2) == 0
%     %     Ai_right = [0 0 1 0;
%     %         0 0 -1 0;
%     %         0 0 0 1;
%     %         0 0 0 -1;
%     %         0 1 0 0;
%     %         0 -1 0 0;
%     %         1 -mu 0 0;
%     %         -1 -mu 0 0;];
%     %     bi = [0; 0; 0; 0; Fmax; -Fmin; 0; 0]; % First four are saying F1 = 0 or F2 = 0. Can relax those a little bit. Next 2 is saying Fmin < Fy < Fmax. Last two are friction constraints.
%     %     bqp = [bqp; bi];
%     % 
%     % elseif alphas(i,1) == 1 && alphas(i,2) == 1
%     Ai_right = [0 -1 0 0;
%         0 0 0 -1;
%         0 1 0 0;
%         0 0 0 1;
%         1 -mu 0 0;
%         -1 -mu 0 0;
%         0 0 1 -mu;
%         0 0 -1 -mu];
%     bi = [-Fmin; -Fmin; Fmax; Fmax; 0; 0; 0; 0];
%     bqp = [bqp; bi];
%     Aqp_k{i}= [Ai_right];
% end
% Aqp = blkdiag(Aqp_k{1}, Aqp_k{2}, Aqp_k{3}, Aqp_k{4}, Aqp_k{5}, Aqp_k{6}, Aqp_k{7}, Aqp_k{8}, Aqp_k{9}, Aqp_k{10}, ...
%     Aqp_k{11}, Aqp_k{12}, Aqp_k{13}, Aqp_k{14}, Aqp_k{15}, Aqp_k{16}, Aqp_k{17}, Aqp_k{18}, Aqp_k{19}, Aqp_k{20});

%% Equality and inequality constraints. This is telling the MPC to not make any contact forces for the swing foot, and to respect friction in x and max/min in y.

Aeq = [];
beq = [];
Aqp = [];
bqp = [];

for i = 1:length(alphas) % hopefully 1:20
    if alphas(i,1) == 0 && alphas(i, 2) == 1
        % Define equality constraints for this time step for the left foot
        tok = zeros(2,80);
        tok(:, 4*i-3: 4*i) = [1 0 0 0;
            0 1 0 0];
        Aeq = [Aeq; tok];
        bi = [0; 0]; % Set left foot forces equal to 0
        beq = [beq; bi];

        % Define inequality constraints for this time step for the right
        % foot
        tok = zeros(4,80);
        tok(:, 4*i-3: 4*i) = [0 0 0 1;
             0 0 0 -1;
             0 0 1 -mu;
             0 0 -1 -mu];
        Aqp = [Aqp; tok];
        bi = [Fmax; -Fmin; 0; 0];
        bqp = [bqp; bi];

    elseif alphas(i,1) == 1 && alphas(i,2) == 0
        % Define equality constraints for this time step for the right foot
        tok = zeros(2,80);
        tok(:, 4*i-3: 4*i) = [0 0 1 0;
            0 0 0 1];
        Aeq = [Aeq; tok];
        bi = [0; 0]; % Set left foot forces equal to 0
        beq = [beq; bi];

        % Define inequality constraints for this time step for the left
        % foot
        tok = zeros(4,80);
        tok(:, 4*i-3: 4*i) = [0 1 0 0;
             0 -1 0 0;
             1 -mu 0 0;
            -1 -mu 0 0;];
        Aqp = [Aqp; tok];
        bi = [Fmax; -Fmin; 0; 0];
        bqp = [bqp; bi];

    elseif alphas(i,1) == 1 && alphas(i,2) == 1
        % Define inequality constraints for this time step for the right
        % foot and left foot
        tok = zeros(8,80);
        tok(:, 4*i-3: 4*i) = [0 -1 0 0;
            0 0 0 -1;
            0 1 0 0;
            0 0 0 1;
            1 -mu 0 0;
            -1 -mu 0 0;
            0 0 1 -mu;
            0 0 -1 -mu];
        Aqp = [Aqp; tok];
        bi = [-Fmin; -Fmin; Fmax; Fmax; 0; 0; 0; 0];
        bqp = [bqp; bi];
    end
end

%% Finally call the QP function
%error('wee')
x_mpc = quadprog(h, f, Aqp, bqp, Aeq, beq);
% [ineqs,eqs,ncalls] = deletionfilter(Aqp,bqp,Aeq,beq,[],[]);
% Aeq(find(eqs), :) = 0;
% beq(find(eqs), :) = 0;
% Aqp(find(ineqs), :) = 0;
% bqp(find(ineqs), :) = 0;
%error('wee')
%x_mpc = quadprog(h, f, Aqp, bqp, Aeq, beq);
end