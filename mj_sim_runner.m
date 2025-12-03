% mj_sim_runner.m

% --- ensure Python can see mj_sim.py ---
scriptFolder = fileparts(mfilename('fullpath'));
if count(py.sys.path, scriptFolder) == 0
    insert(py.sys.path, int32(0), scriptFolder);
end

% Import module into MATLAB
py.mj_sim.init("C:\Users\Zargai\Documents\AME_556\.venv\hector.xml");
module = py.importlib.import_module('mj_sim');

% Load model & data from your XML file
model = module.model;
data  = module.data;

% Hardcode useful body IDs. Add 1 because MATLAB indexing. DOUBLE CHECK THESE IDS
% IF NEW BODIES ARE CREATED
torso_id     = 2 + 1;
leftfoot_id  = 4 + 1;
rightfoot_id = 7 + 1;

% Make a bunch of storage data structures for plotting. Trim them later.

q_storage = zeros(10000, 7);
dq_storage = zeros(10000, 7);
ctrl_storage = zeros(10000, 4);

% --- run mujoco_runner ---
k = 0;
while py.mj_sim.is_running()
    k = k + 1;
    %% Step 1: Access State Data
    q  = double(data.qpos);   % [x_torso, z_torso, torso_pitch, left_hip, left_knee, right_hip, right_knee]
    dq = double(data.qvel);   % corresponding velocities

    xpos  = double(data.xpos);
    cvel = double(data.cvel);
    left_foot_pos = xpos(leftfoot_id,  :);    % [x,y,z] world position of left foot
    right_foot_pos = xpos(rightfoot_id,  :);   % [x,y,z] world position of right foot
    left_foot_vel  = cvel(leftfoot_id,  :);    
    right_foot_vel = cvel(rightfoot_id,  :);   

    left_foot_vel  = left_foot_vel(4:6);   % linear velocity [vx,vy,vz] of left foot
    right_foot_vel = right_foot_vel(4:6);  % linear velocity [vx,vy,vz] of right foot


    %%
    % do come ctrl processing 
    ctrl = py.numpy.array([0, 0, 0, 0]);
    py.mj_sim.set_ctrl(ctrl);
    py.mj_sim.step();

    % Store the state data in the storage arrays
    q_storage(k, :) = q;
    dq_storage(k, :) = dq;
    ctrl_storage(k, :) = ctrl;
end


% Trim data structures to match the time spent in simulation

q_storage = q_storage(1:k,:);
dq_storage =  dq_storage(1:k,:);
ctrl_storage =  ctrl_storage(1:k,:);

% Plotting

py.mj_sim.close();
