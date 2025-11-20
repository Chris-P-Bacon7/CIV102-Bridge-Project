# -------Defining Cross Section--------

#                            top_flange
#              ***************************************
#              ***************************************
#                     ******    glue    ******
#                     **                    **
#                     **                    **
#                     **                    **
#              web    **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                     ************************
#                             bot_flange
top_flange = 100
top_flange_thickness = 1.27
web = 75 - 1.27
web_thickness = 1.27
glue = 5
glue_thickness = 1.27
bot_flange = 80
bot_flange_thickness = 1.27

# -------General Properties------------
bridge_length = 1200
train_weight = 400
x_pos_train = [52, 228, 392, 568, 732, 908]
p_wheel_train = [train_weight / 6] * 6

# -------Solve for SFD/BMD-------------
import matplotlib.pyplot as plt
import numpy as np


def calc_reaction_forces(x_pos, adjustment, p_wheel_train):
    Ay = 0
    By = 0
    total_load = 0
    x_pos_train = x_pos[:]

    # Adjusts Wheel X Position
    for i in range(len(x_pos_train)):
        x_pos_train[i] += adjustment

    # Find Total Load
    for i in range(len(x_pos_train)):
        if x_pos_train[i] >= 0 and x_pos_train[i] <= 1200:
            total_load += p_wheel_train[i]

    # Take Moment About B to Calculate Ay
    moment = 0
    for i in range(len(x_pos_train)):
        if x_pos_train[i] >= 0 and x_pos_train[i] <= 1200:
            moment += p_wheel_train[i] * (bridge_length - x_pos_train[i])
        Ay = moment / bridge_length

    # Find By
    By = total_load - Ay

    return Ay, By, x_pos_train


def graph_sfd(x_pos, adjustment, p_wheel_train, graph):

    # Calculate Reaction Forces + x_pos_train
    Ay, By, x_pos_train = calc_reaction_forces(x_pos, adjustment, p_wheel_train)

    x_points = [0]
    y_points = [0]

    for i in range(0, 1201):
        y_val = 0

        for wheel in range(len(x_pos_train)):
            if x_pos_train[wheel] <= i and x_pos_train[wheel] >= 0:
                y_val -= p_wheel_train[wheel]

        y_val += Ay

        if i == 1200:
            y_val += By

        x_points.append(i)
        y_points.append(y_val)

    xpoints = np.array(x_points)
    ypoints = np.array(y_points)

    if graph:
        plt.plot(xpoints, ypoints)

    y_min = min(y_points)
    x_min = x_points[y_points.index(y_min)]
    y_max = max(y_points)
    x_max = x_points[y_points.index(y_max)]

    if abs(y_min) > abs(y_max):
        y_abs = abs(y_min)
        x_abs = x_min
    else:
        y_abs = abs(y_max)
        x_abs = x_max

    return x_pos_train, y_abs, y_min, y_max, x_abs, x_min, x_max


def greatest_shear(x_pos_train, p_wheel_train, graph):
    res = [0, 0]
    for i in range(-908, 1149):
        sfd_res = graph_sfd(x_pos_train, i, p_wheel_train, graph)
        if sfd_res[1] > res[1]:
            res[0], res[1] = sfd_res[0], sfd_res[1]

    return res


# -------Create SFD/BMD Envelope-------
def shear_envelope(x_pos_train, p_wheel_train, graph):
    x_points_abs = []
    y_points_abs = []
    x_points_min = []
    y_points_min = []
    x_points_max = []
    y_points_max = []

    for j in range(0, 1201):
        temp_min = 0
        temp_max = 0
        temp_abs = 0
        for i in range(-908, 1149):
            Ay, By, x_pos = calc_reaction_forces(x_pos_train, i, p_wheel_train)

            y_val = 0

            for wheel in range(len(x_pos)):
                if x_pos[wheel] <= j and x_pos[wheel] >= 0:
                    y_val -= p_wheel_train[wheel]

            y_val += Ay

            if i == 1200:
                y_val += By

            temp_max = max(temp_max, y_val)
            temp_min = min(temp_min, y_val)
            temp_abs = max(abs(temp_max), abs(temp_min), temp_abs)
        # print(temp_max, temp_min, temp_abs)
        y_points_abs.append(temp_abs)
        y_points_min.append(temp_min)
        y_points_max.append(temp_max)
        x_points_abs.append(j)
        x_points_min.append(j)
        x_points_max.append(j)

    xpoints_abs = np.array(x_points_abs)
    ypoints_abs = np.array(y_points_abs)
    xpoints_min = np.array(x_points_min)
    ypoints_min = np.array(y_points_min)
    xpoints_max = np.array(x_points_max)
    ypoints_max = np.array(y_points_max)

    plt.plot(xpoints_abs, ypoints_abs, label="Absolute Shear")
    plt.plot(xpoints_max, ypoints_max, linestyle="-.", label="Max Shear")
    plt.plot(xpoints_min, ypoints_min, linestyle="-.", label="Min Shear")
    plt.xlabel("Distance Along Beam (mm)")
    plt.ylabel("Shear Force (N)")
    plt.title("Shear Force Envelope")
    plt.legend(loc="upper right")
    plt.axis([0, 1200, -300, 300])


# print(calc_reaction_forces(x_pos_train, -52, p_wheel_train))
# print(greatest_shear(x_pos_train, p_wheel_train, True))
# print(graph_sfd(x_pos_train, -51.999999999999, p_wheel_train))
shear_envelope(x_pos_train, p_wheel_train, False)
plt.show()

# -------Cross-sectional Properties----
A_top = top_flange * top_flange_thickness
y_top = bot_flange_thickness + web + 0.5 * top_flange_thickness

A_web = web * web_thickness
y_web = bot_flange_thickness + 0.5 * web

A_glue = glue * glue_thickness
y_glue = bot_flange_thickness + web - 0.5 * glue_thickness

A_bot = bot_flange * bot_flange_thickness
y_bot = 0.5 * bot_flange_thickness

y_bar = (A_top * y_top + 2 * A_web * y_web + 2 * A_glue * y_glue + A_bot * y_bot) / (
    A_top + 2 * A_web + 2 * A_glue + A_bot
)

I = (
    (top_flange * top_flange_thickness**3) / 12
    + 2 * (web**3 * web_thickness) / 12
    + (bot_flange * bot_flange_thickness**3) / 12
    + 2 * (glue * glue_thickness**3) / 12
    + (A_top * (y_top - y_bar) ** 2)
    + 2 * (A_web * (y_web - y_bar) ** 2)
    + 2 * (A_glue * (y_glue - y_bar) ** 2)
    + (A_bot * (y_bot - y_bar) ** 2)
)

# -------Calculate Applied Stresses----
# max_tensile_stress = max_moment * y_bar / I
# max_compressive_stress = max_moment * (height-y_bar) / I

# -------Material/Thin Plate Buckling--

# -------FOS---------------------------
# FOS_tension = 30/max_tensile_stress
# FOS_compression = 6/max_compressive_stress
# -------Vfail and Mfail---------------

# -------Output Plots------------------
print("x_pos_train:", x_pos_train, "\np_wheel_train:", p_wheel_train)
print("------------------")
print("y_bar:", y_bar, "\nI:", I)
