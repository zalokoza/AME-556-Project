function data = data_extracter(data, body_map, desired_data)

left_foot_id = double(body_map{'left_foot'}) + 1;
right_foot_id = double(body_map{'right_foot'}) + 1;
q = double(data.qpos);
xpos = double(data.xpos);
right_foot_pos = xpos(right_foot_id,  :);
left_foot_pos = xpos(left_foot_id,  :);

if strcmp(desired_data, 'joint_values')
    data = double(data.qpos);
elseif strcmp(desired_data, 'joint_velocities')
    data = double(data.qvel);
elseif strcmp(desired_data, 'joint_accelerations')
    data = double(data.qacc);
elseif strcmp(desired_data, 'body_positions')
    data = double(data.xpos);
elseif strcmp(desired_data, 'body_velocities')
    cvel = double(data.cvel);
    data = cvel;
elseif strcmp(desired_data, 'mpc_state')
    qpos = double(data.qpos);
    qvel = double(data.qvel);
    data = [qpos(1:3)'; qvel(1:3)'; 9.81];
elseif strcmp(desired_data, 'right_foot_location')
    data = right_foot_pos;
elseif strcmp(desired_data, 'left_foot_location')
    data = left_foot_pos;
elseif strcmp(desired_data, 'right_foot_velocity')
    cvel = double(data.cvel);
    data = cvel(right_foot_id,  :);
elseif strcmp(desired_data, 'left_foot_velocity')
    cvel = double(data.cvel);
    data = cvel(left_foot_id,  :);
elseif strcmp(desired_data, 'moment_arms')
    p1z = q(2) -  left_foot_pos(3); % CoM_y - foot_y
    p1x = q(1) - left_foot_pos(1); % CoM_x - foot_x 
    p2z = q(2) -  right_foot_pos(3); % CoM_y - foot_y
    p2x = q(1) - right_foot_pos(1); % CoM_x - foot_x
    data = [p1x, p1z, p2x, p2z];
else
    error('Invalid desired_data input');
end
end
       