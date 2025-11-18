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

# -------Solve for SFD/BMD-------------

# -------Create SFD/BMD Envelope-------

# -------Cross-sectional Properties----
y_bar = (
    top_flange
    * top_flange_thickness
    * (bot_flange_thickness + web + 0.5 * top_flange_thickness)
    + 2 * (web * web_thickness * (bot_flange_thickness + 0.5 * web))
    + 2 * (glue * glue_thickness * (bot_flange_thickness + web - 0.5 * glue_thickness))
    + bot_flange * bot_flange_thickness * 0.5 * bot_flange_thickness
) / (
    top_flange * top_flange_thickness
    + 2 * (web * web_thickness)
    + 2 * (glue * glue_thickness)
    + bot_flange * bot_flange_thickness
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
