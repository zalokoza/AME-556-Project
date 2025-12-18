% mj_sim_runner.m
clear
clc
close all

%% Step 1: Setup

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

% Define dt from MuJoCo
dt = double(model.opt.timestep); 

% Create dictionaries for mapping geom and body names to their IDs
maps = module.get_mujoco_name_maps(model);
body_map = maps{'body_name_to_id'};
geom_map = maps{'geom_name_to_id'};

% Pull feet ID's
left_foot_id = double(body_map{'left_foot'}) + 1;
right_foot_id = double(body_map{'right_foot'}) + 1;

% Define Contact Schedule
right = repmat([0,1], .5/dt/2 , 1); % Right foot on ground for .5s
left = repmat([1,0], .5/dt/2 , 1); % Left foot on ground for .5s
standing = repmat([1,1 ], .5/dt, 1); % Both feet contacted for .5s
%contact_schedule = repmat([right;left], 8, 1); % Right, left, right, left... for 8s
%contact_schedule = ones(4000,2); % Standing schedule
contact_schedule = repmat([standing; right; standing; left;], 8, 1);

% Joint motor saturation limits
lo_vec = [-30; -60; -30; -60];
hi_vec = [30; 60; 30; 60];
k = 1;

%% Plotting set up
x_mpc_storage = zeros(10000,7);
F_storage =  zeros(10000,4);
T_storage = zeros(10000, 4);

%% Step 2: Main Loop

while py.mj_sim.is_running()
    if mod(k-1,10) == 0
        x_mpc = data_extracter(data, body_map, 'mpc_state'); % Generate the MPC-style state vector
        moments = data_extracter(data, body_map, 'moment_arms'); % Get moment arms from CoM to feet
        p1x = moments(1); p1z = moments(2); p2x = moments(3); p2z = moments(4);
        x_qp = horizon(model, x_mpc, contact_schedule, k, p1x, p1z, p2x, p2z); % MPC function. Horizon = 20, discretization = 10*dt
        F_mpc = x_qp(1:4)' ; % Extract U(k) from Xqp
        % Determine contact foot based on schedule
        if contact_schedule(k,1) == 0 && contact_schedule(k,2) == 1
            contact_foot = 'right';
            F_swing = swing_controller(model, data, body_map, contact_schedule, k);
            F = [F_swing F_mpc(3:4)]
            ctrl = force_torque_mapping(module, model, left_foot_id, right_foot_id, np, data, F)
        elseif contact_schedule(k,1) == 1 && contact_schedule(k, 2) == 0
            contact_foot = 'left';
            F_swing = swing_controller(model, data, body_map, contact_schedule, k);
            F = [F_mpc(1:2) F_swing];
            ctrl = force_torque_mapping(module, model, left_foot_id, right_foot_id, np, data, F);
        elseif contact_schedule(k,1) == 1 && contact_schedule(k,2) == 1
            F = F_mpc;
            ctrl = force_torque_mapping(module, model, left_foot_id, right_foot_id, np, data, F);
        end
        ctrl = min(max(ctrl, lo_vec), hi_vec); % clamp maximum torque values
        py.mj_sim.set_ctrl(ctrl');
        py.mj_sim.step();
        x_mpc_storage(k,:) = x_mpc;
        F_storage(k,:) = F;
        T_storage(k,:) = ctrl';
        k = k + 1;
        physical_constraints(data, body_map) % Check if physical constraints are violated or not
    else
        x_mpc = data_extracter(data, body_map, 'mpc_state'); % Generate the MPC-style state vector
        x_mpc_storage(k,:) = x_mpc;
        F_storage(k,:) = F;
        T_storage(k,:) = ctrl';
        ctrl = min(max(ctrl, lo_vec), hi_vec);
        py.mj_sim.set_ctrl(ctrl');
        py.mj_sim.step();
        k = k + 1;
        physical_constraints(data, body_map) % Check if physical constraints are violated or not
    end
end


%% Trim Arrays

x_mpc_storage = x_mpc_storage(1:k-1,:);
F_storage = F_storage(1:k-1,:);
T_storage = T_storage(1:k-1,:);

subplot(3,1,1);
plot(x_mpc_storage);
ylabel('State');
subplot(3,1,2);
plot(F_storage);
ylabel('F');
legend()
subplot(3,1,3);
plot(T_storage);
ylabel('T');
legend()
py.mj_sim.close();

figure
t = (1:k-1)*dt;

subplot(3,4,1);
plot(t, x_mpc_storage(:,(1)));
ylabel('X (m)');
subplot(3,4,2);
plot(t, x_mpc_storage(:,(4)));
ylabel('X Velocity (m/s)');
subplot(3,4,3);
plot(t, x_mpc_storage(:,(2)));
ylabel('Y (m)');
subplot(3,4,4);
plot(t, x_mpc_storage(:,(3)));
ylabel('Theta (rad)');

subplot(3,4,5);
plot(t, F_storage(:,(1)));
ylabel('F1x (N)');
subplot(3,4,6);
plot(t, F_storage(:,(2)));
ylabel('F1y (N)');
subplot(3,4,7);
plot(t, F_storage(:,(3)));
ylabel('F2x (N)');
subplot(3,4,8);
plot(t, F_storage(:,(4)));
ylabel('F2y (N)');

subplot(3,4,9);
plot(t, T_storage(:,(1)));
ylabel('Right Hip Torque (Nm)');
subplot(3,4,10);
plot(t, T_storage(:,(2)));
ylabel('Right Knee Torque (Nm)');
subplot(3,4,11);
plot(t, T_storage(:,(3)));
ylabel('Left Hip Torque (Nm)');
subplot(3,4,12);
plot(t, T_storage(:,(4)));
ylabel('Left Knee Torque (Nm)');