 function F = swing_controller(model, data, body_map, contact_schedule, k)

alpha = contact_schedule(k,:);
left_foot_pos = data_extracter(data, body_map, 'left_foot_location');
right_foot_pos = data_extracter(data, body_map, 'right_foot_location');
left_foot_velocity = data_extracter(data, body_map, 'left_foot_velocity');
gait_length = .50*.50; % .5 seconds to complete the step. .5 m/s. should travel .125m in that time
Kd1 = 16; Kd2 = 16; % Kd1 is gain for the hip joint, Kd2 is gain for the knee

if alpha(1) == 0
    desired_position = right_foot_pos(1) + .4 %gait_length;
    F = -[Kd1 Kd2]*(desired_position(1)-left_foot_pos(1));
elseif alpha(1) == 1
    desired_position = left_foot_pos(1) + .4 %gait_length;
    F = -[Kd1 Kd2]*(desired_position(1)-left_foot_pos(1));
end
%error('wee')
end