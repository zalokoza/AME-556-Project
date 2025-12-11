% mj_sim_runner.m
clear
clc
close all

% --- ensure Python can see mj_sim.py ---
scriptFolder = fileparts(mfilename('fullpath'));
if count(py.sys.path, scriptFolder) == 0
    insert(py.sys.path, int32(0), scriptFolder);
end

% Import module into MATLAB
py.mj_sim.init("G:\My Drive\Courses\AME 556\Final Project MATLAB\hector.xml");
module = py.importlib.import_module('mj_sim');
np = py.importlib.import_module("numpy");
% Load model & data from your XML file
model = module.model;
data  = module.data;

% Define dt from MuJoCo
dt = double(model.opt.timestep); 

% Hardcode useful body IDs. Add 1 because MATLAB indexing. DOUBLE CHECK THESE IDS
% IF NEW BODIES ARE CREATED
% NOTE: NEXT SECTION  CREATES DICTIONARIES FOR GEOM/BODY IDS. USE THAT
% INSTEAD
torso_id     = 2 + 1;
leftfoot_id  = 4 + 1;
rightfoot_id = 7 + 1;

% Gait states 
standing = 0; 
right_support = 1; % right foot is on ground, left foot swings 
left_support = 2; % left foot is on ground, right foot swings 

% Initialize state, time, and step count (standing first then move to walking) 
gait_state = standing; 
phase_time = 0;
step_count = 0; 

% Timing parameters (how long the robot should walk and how long each step
% is) 
standing_duration = 2; % stands for 2 seconds 
step_duration = 2; % each step takes 1.5 seconds 

step_length = 0.12; % x-direction 
step_height = 0.05; % z-direction 


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
AF_storage = zeros(10000, 3);
gait_state_storage = zeros(10000,1); 
% --- run mujoco_runner ---
k = 0;

while py.mj_sim.is_running()
    k = k + 1;

    phase_time = phase_time + dt;

    if gait_state == standing

        if phase_time >= standing_duration
            gait_state = right_support; % Transition to right support phase
            phase_time = 0; % Reset phase time for the next gait phase
            step_count = step_count + 1; 
        end

    elseif gait_state == right_support 

        if phase_time > step_duration 
            gait_state = left_support; % Transition to left support phase 
            phase_time = 0; % Reset phase time for the next gait phase
            step_count = step_count + 1; 
        end

    elseif gait_state == left_support 

        if phase_time > step_duration 
            gait_state = right_support; 
            phase_time = 0; % Reset phase time for the next gait phase 
            step_count = step_count +1; 
        end
    end

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

    if gait_state == standing 

    %% Step 2: QP Control for Standing
    p1z = q(2) -  left_foot_pos(3); % CoM_y - foot_y
    p1x = q(1) - left_foot_pos(1); % CoM_x - foot_x 
    p2z = q(2) -  right_foot_pos(3); % CoM_y - foot_y
    p2x = q(1) - right_foot_pos(1); % CoM_x - foot_x
  
    A = [ 1 0 1 0;
        0 1 0 1;
        -p1z p1x -p2z p2x];
    alpha = 0;

    % Set up b vector i.e. the REAL current forces and torques
    m = py.mujoco.mj_getTotalmass(model);
    I = .0638; % calculated by hand
    b = [m*ddq(1);
        m*ddq(2);
        I*ddq(3);];

    % Set up PD for desired dynamics on x, z, and theta
    Kp1 = 10; Kp2 = 10; Kp3 = 100; 
    Kd1 = 1; Kd2 = 1; Kd3 = 10;
    xd = 0; zd = .45; thetad = 0;
    if k < 500
        zd = .45;
    elseif k < 1000
        zd = .55;
    elseif k < 1500
        zd = .40;
    end
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
    S = diag([1, 1, 1]); % For prioritizing x, y, or theta accelerations
    H = 2*(A')*S*A+2*alpha*eye(4);
    ft = -2*bd'*S*A; 
    f = ft';
    F = quadprog(H, f, Aqp, bqp, [], [], [], [], []);

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

    elseif gait_state == right_support 

        p2x = q(1) - right_foot_pos(1); 
        p2z = q(2) - right_foot_pos(3); 

        A = [1 0;
            0 1;
            -p2z p2x]; 

        xd = right_foot_pos(1) + 0.1; 
        zd = 0.45; 
        thetad = 0; 

        Kp1 = 10; Kp2 = 10; Kp3 = 100; 
        Kd1 = 1; Kd2 = 1; Kd3 = 10; 

        bd = [m*(Kp1*(xd-q(1)) + Kd1*(0 - dq(1)));
        m*(Kp2*(zd-q(2)) + Kd2*(0 - dq(2))) + m*9.81;
        I*(Kp3*(thetad-q(3)) + Kd3*(0 - dq(3)))];

        mu = 0.7; 
        alpha = 0; 

        Aqp = [0 1;
            0 -1;
            1 -mu;
            -1 -mu];

        bqp = [250;-10;0;0]; 

        % Time to do quadratic programming
        S = diag([1, 1, 1]); % For prioritizing x, y, or theta accelerations
        H = 2*(A')*S*A+2*alpha*eye(2);
        ft = -2*bd'*S*A; 
        f = ft';
        F = quadprog(H, f, Aqp, bqp, [], [], [], [], []);


        jacp_right = np.zeros(py.tuple({int32(3), int32(7)}));
        jacr_right = np.zeros(py.tuple({int32(3), int32(7)}));
        module.mujoco.mj_jacBody(model, data, jacp_right, jacr_right, int32(7)) % Left Foot Jacobian
        jacp_right = double(jacp_right);
        T_right = -jacp_right(:,6:7)'*[F(1), 0, F(2)]'; % Left hip, left knee torque mapping to F1x, F1y

        s = phase_time/step_duration; 

        target_foot_x = right_foot_pos(1) + step_length; 
        target_foot_z = step_height*sin(pi*s); 

        rel_x = target_foot_x - q(1); 
        rel_z = target_foot_z - q(2); 

        [q_hip_des, q_knee_des] = leg_inverse_kinematics(rel_x,rel_z); 

        Kp_swing = 100; 
        Kd_swing = 10; 

        q_err = [q_hip_des; q_knee_des] - [q(4);q(5)];
        dq_err = [0; 0] - [dq(4);dq(5)];

        T_left = Kp_swing*q_err+Kd_swing*dq_err; 
        ctrl = [T_left;T_right]; 



        
    else

        p1x = q(1) - left_foot_pos(1); 
        p1z = q(2) - left_foot_pos(3); 

        A = [1 0;
            0 1;
            -p1z p1x]; 

        xd = left_foot_pos(1) + 0.05; 
        zd = 0.45; 
        thetad = 0; 

        Kp1 = 10; Kp2 = 10; Kp3 = 100; 
        Kd1 = 1; Kd2 = 1; Kd3 = 10; 

        bd = [m*(Kp1*(xd-q(1)) + Kd1*(0 - dq(1)));
        m*(Kp2*(zd-q(2)) + Kd2*(0 - dq(2))) + m*9.81;
        I*(Kp3*(thetad-q(3)) + Kd3*(0 - dq(3)))];

        mu = 0.7; 
        alpha = 0; 

        Aqp = [0 1;
            0 -1;
            1 -mu;
            -1 -mu];

        bqp = [250;-10;0;0]; 

        % Time to do quadratic programming
        S = diag([1, 1, 1]); % For prioritizing x, y, or theta accelerations
        H = 2*(A')*S*A+2*alpha*eye(2);
        ft = -2*bd'*S*A; 
        f = ft';
        F = quadprog(H, f, Aqp, bqp, [], [], [], [], []);


        jacp_left = np.zeros(py.tuple({int32(3), int32(7)}));
        jacr_left = np.zeros(py.tuple({int32(3), int32(7)}));
        module.mujoco.mj_jacBody(model, data, jacp_left, jacr_left, int32(4)) % right Foot Jacobian
        jacp_left = double(jacp_left);
        T_left = -jacp_left(:,4:5)'*[F(1), 0, F(2)]'; % right hip, left knee torque mapping to F2x, F2y

        % F = [0;0;0;0]; 
        % bd = [0;0;0]; 
        % ctrl = [0;0;0;0]; 

    end

    %% Step 3: Set the control and step simulation
    py.mj_sim.set_ctrl(ctrl');
    py.mj_sim.step();

    % Store the data in the storage arrays
    % q_storage(k, :) = q;
    % dq_storage(k, :) = dq;
    % ctrl_storage(k, :) = ctrl;
    % F_storage(k,:) = F';
    % bd_storage(k,:) = bd;
    % b_storage(k,:) = b;
    % contact_storage(k) = contact_state;
    % AF_storage(k,:) = (A*F)';
    % gait_state_storage(k) = gait_state; 
end


% Trim data structures to match the time spent in simulation
% 
% q_storage = q_storage(1:k,:);
% dq_storage =  dq_storage(1:k,:);
% ctrl_storage =  ctrl_storage(1:k,:);
% F_storage =  F_storage(1:k,:);
% bd_storage = bd_storage(1:k,:);
% b_storage = b_storage(1:k,:);
% contact_storage = contact_storage(1:k);
% AF_storage = AF_storage(1:k,:);
% gait_state_storage = gait_state_storage(1:k); 
% 
% % Plotting
% 
% figure;
% subplot(2,1,1);
% plot(q_storage);
% legend('x', 'z', 'pitch', 'left_hip', 'left_knee', 'right_hip', 'right_knee');
% title('Joint Positions');
% xlabel('Step'); ylabel('Position (rad)');
% grid on;
% 
% subplot(2,1,2);
% plot(ctrl_storage);
% legend('Left Hip', 'Left Knee', 'Right Hip', 'Right Knee');
% title('Control Inputs');
% xlabel('Step'); ylabel('Torque (Nm)');
% grid on;
% 
% figure
% subplot(6,3,1);
% plot(q_storage(:,1));
% ylabel('Torso x');
% subplot(6,3,2);
% plot(q_storage(:,2));
% ylabel('Torso y');
% subplot(6,3,3);
% plot(q_storage(:,3));
% ylabel('Torso theta');
% 
% subplot(6,3,4);
% plot(dq_storage(:,1));
% ylabel('Torso dx');
% subplot(6,3,5);
% plot(dq_storage(:,2));
% ylabel('Torso dy');
% subplot(6,3,6);
% plot(dq_storage(:,3));
% ylabel('Torso dtheta');
% 
% subplot(6,3,7);
% plot(bd_storage(:,1));
% ylabel('Bd x');
% subplot(6,3,8);
% plot(bd_storage(:,2));
% ylabel('Bd y');
% subplot(6,3,9);
% plot(bd_storage(:,3));
% ylabel('Bd theta');
% 
% subplot(6,3,10);
% plot(b_storage(:,1));
% ylabel('B x');
% subplot(6,3,11);
% plot(b_storage(:,2));
% ylabel('B y');
% subplot(6,3,12);
% plot(b_storage(:,3));
% ylabel('B theta');
% 
% subplot(6,3,13);
% plot(AF_storage(:,1));
% ylabel('AF x');
% subplot(6,3,14);
% plot(AF_storage(:,2));
% ylabel('AF y');
% subplot(6,3,15);
% plot(AF_storage(:,3));
% ylabel('AF theta');
% 
% subplot(6,3,16);
% plot(F_storage(:,1) + F_storage(:,3));
% ylabel('Fx');
% subplot(6,3,17);
% plot(F_storage(:,2) + F_storage(:,4));
% ylabel('Fy');
% subplot(6,3,18);
% plot(contact_storage);
% ylabel('Contact State');
% 
% figure; 
% plot(gait_state_storage, 'LineWidth', 2);
% ylim([-0.5, 2.5]);
% yticks([0, 1, 2]);
% yticklabels({'STANDING', 'LEFT_SUPPORT', 'RIGHT_SUPPORT'});
% ylabel('Gait State');
% xlabel('Simulation Step');
% title('State Machine Transitions');
% grid on;

py.mj_sim.close();


%% Inverse Kinematics for 2-Link Leg
function [q_hip, q_knee] = leg_inverse_kinematics(foot_x, foot_z)
    % Computes hip and knee angles to place foot at (foot_x, foot_z)
    % relative to the hip joint
    
    L1 = 0.22;  % thigh length (meters)
    L2 = 0.22;  % shin length (meters)
    
    % Distance from hip to target foot position
    r = sqrt(foot_x^2 + foot_z^2);
    
    % Clamp to reachable workspace
    r_min = abs(L1 - L2) + 0.01;  % can't fold completely
    r_max = L1 + L2 - 0.01;      
    r = max(min(r, r_max), r_min);
    
    % Knee angle using law of cosines
    cos_knee = (r^2 - L1^2 - L2^2) / (2*L1*L2);
    cos_knee = max(min(cos_knee, 1), -1);  % clamp to valid range
    q_knee = acos(cos_knee);
    
    % Hip angle using geometry
    alpha = atan2(foot_x, -foot_z);  % angle to target
    beta = asin(L2 * sin(q_knee) / r);  
    q_hip = alpha - beta;
end
