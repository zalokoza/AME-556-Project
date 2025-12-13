function ctrl = force_torque_mapping(module, model, left_foot_id, right_foot_id, np, data, F)
% F is the direct result of horizon.m
    jacp_left = np.zeros(py.tuple({int32(3), int32(7)}));
    jacr_left = np.zeros(py.tuple({int32(3), int32(7)}));
    module.mujoco.mj_jacBody(model, data, jacp_left, jacr_left, py.int(left_foot_id -1)) % Left Foot Jacobian
    jacp_left = double(jacp_left);
    T_left = -jacp_left(:,4:5)'*[F(1), 0, F(2)]'; % Left hip, left knee torque mapping to F1x, F1y
    jacp_right = np.zeros(py.tuple({int32(3), int32(7)}));
    jacr_right = np.zeros(py.tuple({int32(3), int32(7)}));
    module.mujoco.mj_jacBody(model, data, jacp_right, jacr_right, py.int(right_foot_id-1)) % Right Foot Jacobian
    jacp_right = double(jacp_right);
    T_right = -jacp_right(:,6:7)'*[F(3), 0, F(4)]'; % Right hip, right knee torque mapping to F2x, F2y
    ctrl = [T_left' T_right']';
end