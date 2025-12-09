% mj_sim_runner.m
clear
clc
close all

function tau = footForcesToTorques(model, data, F_des)
% footForcesToTorques  Map desired foot forces (2D) to joint torques. We
% use QP or MPC algorithms to calculate the desired foot contact forces
% in order to achieve specific dynamics for the robot. Then, we use this
% force-torque-mapping function to calculate which torques need to be
% applied to which joints in order to instantenously achieve those desired
% feet (end effector) forces.
%
%   tau = footForcesToTorques(model, data, F_des)
%
%   Inputs:
%     model  - py.mujoco.MjModel object
%     data   - py.mujoco.MjData  object (at current state)
%     F_des  - 4x1 vector of desired foot forces in world X/Z:
%                F_des = [Fx_L; Fz_L; Fx_R; Fz_R]
%
%   Output:
%     tau    - 4x1 vector of joint torques:
%                tau = [tau_LH; tau_LK; tau_RH; tau_RK]
%
%   Assumptions:
%     - Model is your pseudo_2D_humanoid.
%     - Joint order in qvel is:
%         [slide_x, slide_z, torso_pitch, left_hip, left_knee, right_hip, right_knee]
%     - You only actuate the last 4 joints (hips & knees).
%     - left_foot and right_foot bodies are your end-effectors, with body IDs:
%         world(0), torso(1), left_thigh(2), left_shin(3), left_foot(4),
%         right_thigh(5), right_shin(6), right_foot(7)
%       (MuJoCo 0-based indexing).
%
%   Adjust LEFT_FOOT_ID / RIGHT_FOOT_ID below if your IDs differ.

    arguments
        model
        data
        F_des (4,1) double
    end

    %--- Hard-coded body IDs (MuJoCo 0-based) ---
    LEFT_FOOT_ID  = int32(4);   % left_foot
    RIGHT_FOOT_ID = int32(7);   % right_foot

    %--- Import Python modules ---
    mujoco = py.importlib.import_module('mujoco');
    np     = py.importlib.import_module('numpy');

    % Number of velocity DOFs
    nv = int32(model.nv);

    %--- Helper to get 4x4 Jacobian (X,Z forces and Y moment) for a given foot body ---
    function J_foot = getFootJacobian2D(body_id)
        % Allocate numpy arrays (3 x nv) for translational and rotational Jacobians
        jacp = np.zeros(py.tuple({int32(3), nv}), pyargs('dtype', np.float64));
        jacr = np.zeros(py.tuple({int32(3), nv}), pyargs('dtype', np.float64));

        % Fill Jacobian at body origin (world frame)
        mujoco.mj_jacBody(model, data, jacp, jacr, body_id);

        % Convert to MATLAB double
        Jp = double(jacp);    % 3 x nv translational J
        Jr = double(jacr);    % 3 x nv rotational J

        % Extract planar rows: X (row 1) and Z (row 3) for forces
        % Extract Y rotation (row 2) for moment
        % Extract columns for the 4 actuated joints: LH, LK, RH, RK
        % Those are qvel indices 4:7 in MATLAB (1-based).
        % For 2D motion in XZ plane: forces in X,Z and moment about Y axis
        J_foot = [Jp(1,4:7);      % Fx row (translational X)
                  Jp(3,4:7);      % Fz row (translational Z)
                  Jr(2,4:7)];     % My row (rotational Y)    --> 3 x 4
    end

    %--- Compute Jacobians for left and right feet (3 x 4 each) ---
    J_L = getFootJacobian2D(LEFT_FOOT_ID);   % left foot
    J_R = getFootJacobian2D(RIGHT_FOOT_ID);  % right foot

    %--- Stack into full task Jacobian (6 x 4) ---
    % task vector is [Fx_L; Fz_L; My_L; Fx_R; Fz_R; My_R]
    % But we only have 4 DOF in actuated joints, so this is overdetermined
    % Use pseudo-inverse to find best-fit solution
    J = [J_L;
         J_R];   % 6 x 4

    %--- Map desired foot forces to joint torques ---
    % tau = J^+ * F  (pseudo-inverse since system is overdetermined)
    % tau = pinv(J) * F
    F_full = [F_des(1); F_des(2); 0; F_des(3); F_des(4); 0];  % [Fx_L; Fz_L; 0; Fx_R; Fz_R; 0]
    tau = pinv(J) * F_full;  % 4x1

    % Optional: handle singularities with a tiny damping (uncomment if needed)
    % lambda = 1e-6;
    % tau = J.' * ((J*J.' + lambda*eye(4)) \ F);
end


% --- ensure Python can see mj_sim.py ---
scriptFolder = fileparts(mfilename('fullpath'));
if count(py.sys.path, scriptFolder) == 0
    insert(py.sys.path, int32(0), scriptFolder);
end

% Import module into MATLAB
py.mj_sim.init("C:\Users\Zargai\Documents\AME_556\.venv\AME-556-Project\hector.xml");
module = py.importlib.import_module('mj_sim');
np = py.importlib.import_module("numpy");
% Load model & data from your XML file
model = module.model;
data  = module.data;

% Hardcode useful body IDs. Add 1 because MATLAB indexing. DOUBLE CHECK THESE IDS
% IF NEW BODIES ARE CREATED
% NOTE: NEXT SECTION  CREATES DICTIONARIES FOR GEOM/BODY IDS. USE THAT
% INSTEAD
torso_id     = 2 + 1;
leftfoot_id  = 4 + 1;
rightfoot_id = 7 + 1;

% Create dictionaries for mapping geom and body names to their IDs
maps = module.get_mujoco_name_maps(model);
body_map = maps{'body_name_to_id'};
geom_map = maps{'geom_name_to_id'};
% Make a bunch of storage data structures for plotting. Trim them later.

q_storage = zeros(10000, 7);
dq_storage = zeros(10000, 7);
ddq_storage = zeros(10000, 3); % 3 for 3 Newton's equations
ctrl_storage = zeros(10000, 4);
F_storage = zeros(10000, 4);
bd_storage = zeros(10000, 3);
b_storage = zeros(10000, 3);
contact_storage = zeros(10000, 1);
% --- run mujoco_runner ---
k = 0;
while py.mj_sim.is_running()
    k = k + 1;
    %% Step 1: Access State Data
    q  = double(data.qpos);   % [x_torso, z_torso, torso_pitch, left_hip, left_knee, right_hip, right_knee]
    dq = double(data.qvel);   % corresponding velocities
    ddq = double(data.qacc); % corresponding accelerations
    xpos  = double(data.xpos);
    cvel = double(data.cvel);
    left_foot_pos = xpos(leftfoot_id,  :);    % [x,y,z] world position of left foot
    right_foot_pos = xpos(rightfoot_id,  :);   % [x,y,z] world position of right foot
    left_foot_vel  = cvel(leftfoot_id,  :);    
    right_foot_vel = cvel(rightfoot_id,  :);   

    left_foot_vel  = left_foot_vel(4:6);   % linear velocity [vx,vy,vz] of left foot
    right_foot_vel = right_foot_vel(4:6);  % linear velocity [vx,vy,vz] of right foot

    %% Step 1.5: Check Foot Contact
    contact_state = int32(module.foot_contact_state(model, data));
    %% Step 2: QP Control for Standing
    p1z = q(2) -  left_foot_pos(3); % CoM_y - foot_y
    p1x = q(1) - left_foot_pos(1); % CoM_x - foot_x 
    p2z = q(2) -  right_foot_pos(3); % CoM_y - foot_y
    p2x = q(1) - right_foot_pos(1); % CoM_x - foot_x
  
    A = [ 1 0 1 0;
        0 1 0 1;
        -p1z p1x -p2z p2x];
    alpha = 1;

    % Set up b vector i.e. the REAL current forces and torques
    m = py.mujoco.mj_getTotalmass(model);
    I = .0638; % calculated by hand
    b = [m*ddq(1);
        m*ddq(2);
        I*ddq(3);];

    % Set up PD for desired dynamics on x, z, and theta
    Kp1 = 1; Kp2 = 1; Kp3 = 100; Kd1 = .1; Kd2 = .1; Kd3 = 10;
    xd = 0; zd = .45; thetad = 0;
    bd = [m*(Kp1*(xd-q(1)) + Kd1*(0 - dq(1)));
        m*(Kp2*(zd-q(2)) + Kd2*(0 - dq(2))) + m*9.81;
        I*(Kp3*(thetad-q(3)) + Kd3*(0 - dq(3)))];


    % Set up inequality constraints on the contact forces
    mu = .7;
    Aqp = [0 1 0 0;
        0 0 0 1;
        0 -1 0 0;
        0 0 0 -1;
        1 -mu 0 0;
        0 0 1 -mu;
        -1 -mu 0 0; 
        0 0 -1 -mu];
        
    bqp = [250; 250; -10; -10; 0; 0; 0; 0];

        
    
    % Time to do quadratic programming
    H = 2*(A')*A+2*alpha*eye(4);
    ft = -2*bd'*A; 
    f = ft';
    F = quadprog(H, f, Aqp, bqp, [], [], [], [], []);

    % Force-Torque Mapping
    % TODO: Replace these three lines with foot-contact logic i.e. 
    if k < 100
        ctrl = [0; 0; 0; 0];
    else
  %  ctrl = footForcesToTorques(model,data,-F);
    
    %% Step 3.5: Zargai attempt to do Force- Torque Mapping
    jacp_left = np.zeros(py.tuple({int32(3), int32(7)}));
    jacr_left = np.zeros(py.tuple({int32(3), int32(7)}));
    module.mujoco.mj_jacBody(model, data, jacp_left, jacr_left, int32(4)) % Left Foot Jacobian
    jacp_left = double(jacp_left);
    T_left = -jacp_left(:,4:5)'*[F(1), 0, F(2)]'; % Left hip, left knee torque mapping to F1x, F1y
    jacp_right = np.zeros(py.tuple({int32(3), int32(7)}));
    jacr_right = np.zeros(py.tuple({int32(3), int32(7)}));
    module.mujoco.mj_jacBody(model, data, jacp_right, jacr_right, int32(7)) % Right Foot Jacobian
    jacp_right = double(jacp_right);
    T_right = -jacp_right(:,6:7)'*[F(3), 0, F(4)]'; % Right hip, right knee torque mapping to F2x, F2y
    ctrl = [T_left' T_right']';
    end
    %% Step 3: Set the control and step simulation
    py.mj_sim.set_ctrl(ctrl');
    py.mj_sim.step();

    % Store the data in the storage arrays
    q_storage(k, :) = q;
    dq_storage(k, :) = dq;
    ctrl_storage(k, :) = ctrl;
    F_storage(k,:) = F';
    bd_storage(k,:) = bd;
    b_storage(k,:) = b;
    contact_storage(k) = contact_state;
end


% Trim data structures to match the time spent in simulation

q_storage = q_storage(1:k,:);
dq_storage =  dq_storage(1:k,:);
ctrl_storage =  ctrl_storage(1:k,:);
F_storage =  F_storage(1:k,:);
bd_storage = bd_storage(1:k,:);
b_storage = b_storage(1:k,:);
contact_storage = contact_storage(1:k);

% Plotting

figure;
subplot(2,1,1);
plot(q_storage);
legend('x', 'z', 'pitch', 'left_hip', 'left_knee', 'right_hip', 'right_knee');
title('Joint Positions');
xlabel('Step'); ylabel('Position (rad)');
grid on;

subplot(2,1,2);
plot(ctrl_storage);
legend('Left Hip', 'Left Knee', 'Right Hip', 'Right Knee');
title('Control Inputs');
xlabel('Step'); ylabel('Torque (Nm)');
grid on;

figure
subplot(5,3,1);
plot(q_storage(:,1));
ylabel('Torso x');
subplot(5,3,2);
plot(q_storage(:,2));
ylabel('Torso y');
subplot(5,3,3);
plot(q_storage(:,3));
ylabel('Torso theta');

subplot(5,3,4);
plot(dq_storage(:,1));
ylabel('Torso dx');
subplot(5,3,5);
plot(dq_storage(:,2));
ylabel('Torso dy');
subplot(5,3,6);
plot(dq_storage(:,3));
ylabel('Torso dtheta');

subplot(5,3,7);
plot(bd_storage(:,1));
ylabel('Bd x');
subplot(5,3,8);
plot(bd_storage(:,2));
ylabel('Bd y');
subplot(5,3,9);
plot(bd_storage(:,3));
ylabel('Bd theta');

subplot(5,3,10);
plot(b_storage(:,1));
ylabel('B x');
subplot(5,3,11);
plot(b_storage(:,2));
ylabel('B y');
subplot(5,3,12);
plot(b_storage(:,3));
ylabel('B theta');

subplot(5,3,13);
plot(F_storage(:,1) + F_storage(:,3));
ylabel('Fx');
subplot(5,3,14);
plot(F_storage(:,2) + F_storage(:,4));
ylabel('Fy');
subplot(5,3,15);
plot(contact_storage);
ylabel('Contact State');

py.mj_sim.close();
