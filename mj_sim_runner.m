% mj_sim_runner.m

% --- ensure Python can see mj_sim.py ---
scriptFolder = fileparts(mfilename('fullpath'));
if count(py.sys.path, scriptFolder) == 0
    insert(py.sys.path, int32(0), scriptFolder);
end

% --- run mujoco_runner ---
py.mj_sim.init("put the absolute path to your mjcf file here");

while py.mj_sim.is_running()
    % do come ctrl processing 
    ctrl = py.numpy.array([0.5, 0.5]);
    py.mj_sim.set_ctrl(ctrl);
    py.mj_sim.step();
end

py.mj_sim.close();
