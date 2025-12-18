function status = physical_constraints(data, body_map)
joint_pos = data_extracter(data, body_map, 'joint_values')
joint_pos = joint_pos(4:7)*(180/pi); % Select only the joints and convert to degree
joint_vel = data_extracter(data, body_map, 'joint_velocities')
joint_vel = joint_vel(4:7); % Select only the joints. keep in rad/s
if joint_pos(1) < - 120 || joint_pos(1) > 30 || joint_pos(3) < - 120 || joint_pos(3) > 30
    status = 0;
    display(joint_pos);
    error('Hip joint angle limit exceeded')
elseif joint_vel(1) < -30 || joint_vel(1) > 30 || joint_vel(3) < - 30 || joint_vel(3) > 30
    status = 0;
    display(joint_vel)
    error('Hip rotating too fast')
elseif joint_pos(2) < -.1 || joint_pos(2) > 160 || joint_pos(4) < -.1 || joint_pos(4) > 160 % -.1 because -0.0000 would throw an error
    status = 0;
    display(joint_pos);
    error ('Knee joint angle limit exceeded')
elseif joint_vel(2) < -15 || joint_vel(2) > 15 || joint_vel(4) < - 15 || joint_vel(4) > 15
    status = 0;
    display(joint_vel)
    error ('Knee rotating too fast')
else
    status = 1;
end