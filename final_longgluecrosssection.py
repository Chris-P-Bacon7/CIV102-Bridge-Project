# In this program, the bottom flange is not included in the calculations as it introduces
# an illegal cross-section for the project requirements.

import math
import matplotlib.pyplot as plt
import numpy as np

# -------Defining Cross Section--------

# 0 -> 1/3 and 2/3 -> 3/3
#                            top_flange
#              ***************************************
#              ***************************************
#         glue (20)   *******         ********
#                     **                    **
#                     **                    **
#                     **                    **
#              web    **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                     **                    **

# 1/3 -> 2/3                  top_flange
#              ***************************************
#              ***************************************
#                     ************************
#                     **                    **
#                     **                    **
#                     **                    **
#              web    **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                     **                    **
#                       
top_flange = 100
top_flange_thickness = 2.54 
web = 120
web_thickness = 1.27
glue = 80-2.54
glue_thickness = 1.27
web_separation = 80

# -------General Properties------------
bridge_length = 1200
train_weight = 400
x_pos_train = [52, 228, 392, 568, 732, 908]
p_wheel_train = [train_weight / 6] * 6
E = 4000
mu = 0.2
a = 100
# -------Solve for SFD/BMD-------------

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
    y_max = max(y_points)

    if abs(y_min) > abs(y_max):
        y_abs = abs(y_min)
    else:
        y_abs = abs(y_max)

    return x_pos_train, y_abs, xpoints[1:], ypoints[1:]


def graph_bmd(x_pos, adjustment, p_wheel_train, graph):
    sfd_vals = graph_sfd(x_pos, adjustment, p_wheel_train, False)
    x_pos, xpoints, ypoints = sfd_vals[0], sfd_vals[2], sfd_vals[3]

    dx = np.diff(xpoints)
    avg_height = (ypoints[:-1] + ypoints[1:]) / 2
    trapezoid_areas = avg_height * dx

    bmd_ypoints = np.insert(np.cumsum(trapezoid_areas), 0, 0)
    bmd_xpoints = xpoints

    if graph:
        plt.gca().invert_yaxis()
        plt.plot(bmd_xpoints, bmd_ypoints)

    y_min = min(bmd_ypoints)
    y_max = max(bmd_ypoints)

    if abs(y_min) > abs(y_max):
        y_abs = abs(y_min)
    else:
        y_abs = abs(y_max)

    return x_pos, y_abs, bmd_xpoints, bmd_ypoints.tolist()


def greatest_shear(x_pos_train, p_wheel_train, graph):
    res = [0, 0]
    for i in range(-908, 1149):
        sfd_res = graph_sfd(x_pos_train, i, p_wheel_train, graph)
        if sfd_res[1] > res[1]:
            res[0], res[1] = sfd_res[0], sfd_res[1]

    return res


def greatest_moment(x_pos_train, p_wheel_train, graph):
    res = [0, 0]
    for i in range(-908, 1149):
        bmd_res = graph_bmd(x_pos_train, i, p_wheel_train, graph)
        if bmd_res[1] > res[1]:
            res[0], res[1] = bmd_res[0], bmd_res[1]

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

    plt.figure()
    if graph:
        plt.plot(xpoints_abs, ypoints_abs, label="Absolute Shear")
        plt.plot(xpoints_max, ypoints_max, linestyle="-.", label="Max Shear")
        plt.plot(xpoints_min, ypoints_min, linestyle="-.", label="Min Shear")
    else:
        plt.plot(
            xpoints_abs,
            ypoints_abs,
            linestyle="-.",
            label="Absolute Shear",
            linewidth=2,
        )
    plt.xlabel("Distance Along Beam (mm)")
    plt.ylabel("Shear Force (N)")
    plt.title("Shear Force Envelope")
    plt.legend(loc="upper right")
    return xpoints_abs, ypoints_abs


def moment_envelope(x_pos_train, p_wheel_train):
    x_points_max = []
    for j in range(0, 1201):
        x_points_max.append(j)
    y_points_max = [0] * 1201

    for i in range(-908, 1149):
        bmd_res = graph_bmd(x_pos_train, i, p_wheel_train, False)
        y_temp = bmd_res[3]

        for x in range(len(y_temp)):
            if y_temp[x] > y_points_max[x]:
                y_points_max[x] = y_temp[x]

    xpoints_max = np.array(x_points_max)
    ypoints_max = np.array(y_points_max)

    plt.figure()
    plt.plot(xpoints_max, ypoints_max, linestyle="-.", label="Max Moment", linewidth=2)
    plt.xlabel("Distance Along Beam (mm)")
    plt.ylabel("Moment (Nmm)")
    plt.title("Bending Moment Envelope")
    plt.gca().invert_yaxis()
    return xpoints_max, ypoints_max


# print(graph_bmd(x_pos_train, -51, p_wheel_train, True)[0:2])
print(greatest_moment(x_pos_train, p_wheel_train, True))
# moment_envelope(x_pos_train, p_wheel_train)

# print(calc_reaction_forces(x_pos_train, -52, p_wheel_train))
# print(greatest_shear(x_pos_train, p_wheel_train, True))
print(graph_sfd(x_pos_train, -51.999999999999, p_wheel_train, True))
# shear_envelope(x_pos_train, p_wheel_train)

# -------Cross-sectional Properties----
A_top = top_flange * top_flange_thickness
y_top = web + glue_thickness + 0.5 * top_flange_thickness

A_web = web * web_thickness
y_web = 0.5 * web

A_glue = glue * glue_thickness
y_glue = web + 0.5 * glue_thickness


y_bar = (A_top * y_top + 2 * A_web * y_web + A_glue * y_glue) / (
    A_top + 2 * A_web + A_glue
)

I = (
    (top_flange * top_flange_thickness**3) / 12
    + 2 * (web**3 * web_thickness) / 12
    + (glue * glue_thickness**3) / 12
    + (A_top * (y_top - y_bar) ** 2)
    + 2 * (A_web * (y_web - y_bar) ** 2)
    + (A_glue * (y_glue - y_bar) ** 2)
)

y_tot = web + top_flange_thickness
y_from_top = y_tot - y_bar

b_mat = web_thickness * 2
b_glue = glue

Q_mat = 2 * (web_thickness * (y_bar) 
             * (y_bar) / 2
)
Q_glue = A_top * (y_top - y_bar)

# -------Calculate Applied Stresses----
max_tensile_stress = greatest_moment(x_pos_train, p_wheel_train, False)[1] * y_bar / I
max_compressive_stress = (
    greatest_moment(x_pos_train, p_wheel_train, False)[1] * y_from_top / I
)
max_shear_stress_mat = (
    greatest_shear(x_pos_train, p_wheel_train, False)[1] * Q_mat
) / (I * b_mat)
max_shear_stress_glue = (
    greatest_shear(x_pos_train, p_wheel_train, False)[1] * Q_glue
) / (I * b_glue)

# print("max_tensile_stress", max_tensile_stress)
# print("max_compressive_stress", max_compressive_stress)
# print("max_shear_stress_mat", max_shear_stress_mat)
# print("max_shear_stress_glue", max_shear_stress_glue)

# -------Material/Thin Plate Buckling--
buckling_case1 = (((4 * math.pi**2) * E) / (12 * (1 - mu**2))) * (
    top_flange_thickness / (web_separation-web_thickness)
) ** 2
buckling_case2 = (((0.425 * math.pi**2) * E) / (12 * (1 - mu**2))) * (
    top_flange_thickness / ((top_flange - web_separation) / 2)
) ** 2
buckling_case3 = (((6 * math.pi**2) * E) / (12 * (1 - mu**2))) * (
    (web_thickness) / ((web - y_bar))
) ** 2
buckling_case4 = (((5 * math.pi**2) * E) / (12 * (1 - mu**2))) * (
    ((web_thickness) / (web)) ** 2
    + ((web_thickness) / (a)) ** 2
)

# print("buckling_case1", buckling_case1)
# print("buckling_case2", buckling_case2)
# print("buckling_case3", buckling_case3)
# print("buckling_case4", buckling_case4)
# -------FOS---------------------------
FOS_tension = 30 / max_tensile_stress
FOS_compression = 6 / max_compressive_stress
FOS_shear_mat = 4 / max_shear_stress_mat
FOS_shear_glue = 2 / max_shear_stress_glue
FOS_buck_1 = buckling_case1 / max_compressive_stress
FOS_buck_2 = buckling_case2 / max_compressive_stress
FOS_buck_3 = buckling_case3 / max_compressive_stress
FOS_buck_4 = buckling_case4 / max_shear_stress_mat

print("FOS_tension", FOS_tension)
print("FOS_compression", FOS_compression)
print("FOS_shear_mat", FOS_shear_mat)
print("FOS_shear_glue", FOS_shear_glue)
print("FOS_buck_1", FOS_buck_1)
print("FOS_buck_2", FOS_buck_2)
print("FOS_buck_3", FOS_buck_3)
print("FOS_buck_4", FOS_buck_4)
# -------Vfail and Mfail---------------
moment_env = moment_envelope(x_pos_train, p_wheel_train)
M_fail_tens_y = moment_env[1] * FOS_tension
M_fail_comp_y = moment_env[1] * FOS_compression
M_fail_buck1_y = moment_env[1] * FOS_buck_1
M_fail_buck2_y = moment_env[1] * FOS_buck_2
M_fail_buck3_y = moment_env[1] * FOS_buck_3
plt.plot(moment_env[0], M_fail_tens_y, label="Tension Failure")
plt.plot(moment_env[0], M_fail_comp_y, label="Compression Failure")
plt.plot(moment_env[0], M_fail_buck1_y, label="Buckling 1 Failure")
plt.plot(moment_env[0], M_fail_buck2_y, label="Buckling 2 Failure")
plt.plot(moment_env[0], M_fail_buck3_y, label="Buckling 3 Failure")

plt.legend(loc="upper right")

shear_env = shear_envelope(x_pos_train, p_wheel_train, False)
S_fail_mat_y = shear_env[1] * FOS_shear_mat
S_fail_glue_y = shear_env[1] * FOS_shear_glue
S_fail_buck4_y = shear_env[1] * FOS_buck_4

plt.plot(shear_env[0], S_fail_mat_y, label="Mat Shear Failure")
# plt.plot(shear_env[0], S_fail_glue_y, label="Glue Shear Failure")
plt.plot(shear_env[0], S_fail_buck4_y, label="Buckling 4 Shear Failure")

plt.legend(loc="upper right")

plt.show()
# -------Output Plots------------------
print("x_pos_train:", x_pos_train, "\np_wheel_train:", p_wheel_train)
print("------------------")
print("y_bar:", y_bar, "\nI:", I)