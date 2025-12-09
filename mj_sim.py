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

    data.qpos[1] = .572
    data.qpos[3] = .3
    data.qpos[4] = 0
    data.qpos[5] = -.3
    data.qpos[6] = 0
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

def get_mujoco_name_maps(model):
    """
    Return two dictionaries mapping names -> IDs:
        body_name_to_id : { 'body_name': int(id), ... }
        geom_name_to_id : { 'geom_name': int(id), ... }

    This is handy for resolving names once and then reusing IDs.
    """
    body_name_to_id = {}
    geom_name_to_id = {}

    # Bodies
    for i in range(model.nbody):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i)
        if name is not None:
            body_name_to_id[name] = int(i)

    # Geoms
    for i in range(model.ngeom):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, i)
        if name is not None:
            geom_name_to_id[name] = int(i)

    return {
        "body_name_to_id": body_name_to_id,
        "geom_name_to_id": geom_name_to_id,
    }


def foot_contact_state(
    model,
    data,
    left_geom_name="left_contact",
    right_geom_name="right_contact",
    floor_geom_name="floor",
):
    """
    Determine which feet are in contact with the ground.

    Returns an int state code:
        0 : neither foot in contact
        1 : left foot in contact
        2 : right foot in contact
        3 : both feet in contact

    Parameters
    ----------
    model : mujoco.MjModel
    data  : mujoco.MjData
    left_geom_name  : str, name of left foot contact geom in the MJCF
    right_geom_name : str, name of right foot contact geom in the MJCF
    floor_geom_name : str, name of floor geom in the MJCF
    """

    maps = get_mujoco_name_maps(model)
    geom_map = maps["geom_name_to_id"]

    try:
        gid_left = geom_map[left_geom_name]
        gid_right = geom_map[right_geom_name]
        gid_floor = geom_map[floor_geom_name]
    except KeyError as e:
        raise ValueError(
            f"Geom name {e.args[0]!r} not found in model. "
            "Check your MJCF names or arguments to foot_contact_state()."
        )

    left_contact = False
    right_contact = False

    # Loop over all active contacts
    for i in range(data.ncon):
        con = data.contact[i]
        g1 = con.geom1
        g2 = con.geom2

        # Left foot vs floor
        if ((g1 == gid_left and g2 == gid_floor) or
            (g2 == gid_left and g1 == gid_floor)):
            left_contact = True

        # Right foot vs floor
        if ((g1 == gid_right and g2 == gid_floor) or
            (g2 == gid_right and g1 == gid_floor)):
            right_contact = True

    state = 0
    if left_contact:
        state += 1
    if right_contact:
        state += 2

    return int(state)

if __name__ == "__main__":
    if init("examples/03_two_link_robot.xml"):
        while is_running():
            set_ctrl([0.5, 0.0])  # Example control input
            step()
        close()