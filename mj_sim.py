# mj_sim.py
import time
import mujoco
import mujoco.viewer

# Global variables to keep model, data, and viewer alive
model = None
data = None
viewer = None
start = None

def init(xml_path="examples/ball_floot.xml"):
    """Initialize the Mujoco model, data, and viewer."""
    global model, data, viewer, start
    model = mujoco.MjModel.from_xml_path(xml_path)
    data = mujoco.MjData(model)
    viewer = mujoco.viewer.launch_passive(model, data)
    start = time.time()
    return True

def step():
    """Perform one simulation step and update the viewer."""
    global model, data, viewer
    if viewer is None:
        return False
    
    if not viewer.is_running():
        close()
        return False

    stepstart = time.time()
    mujoco.mj_step(model, data)

    # Toggle contact points every 2 seconds
    with viewer.lock():
        viewer.opt.flags[mujoco.mjtVisFlag.mjVIS_CONTACTPOINT] = int(data.time % 2)

    viewer.sync()

    # Sleep until next step
    time_until_next_step = model.opt.timestep - (time.time() - stepstart)
    if time_until_next_step > 0:
        time.sleep(time_until_next_step)

    return True

def is_running():
    """Check if viewer is open and within time limit."""
    global viewer, start
    if viewer is None:
        return False
    return viewer.is_running()

def close():
    """Close the viewer cleanly."""
    global viewer
    if viewer is not None:
        viewer.close()
        viewer = None
    return True

def set_ctrl(ctrl):
    """Set control vector"""
    global model, data
    if model is None or data is None:
        return False
    data.ctrl[:] = ctrl
    return True

if __name__ == "__main__":
    if init("examples/03_two_link_robot.xml"):
        while is_running():
            set_ctrl([0.5, 0.0])  # Example control input
            step()
        close()