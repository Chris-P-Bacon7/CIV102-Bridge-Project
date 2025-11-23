import math
import matplotlib.pyplot as plt
import numpy as np

# ==========================================
# PART 1: DEFINE GEOMETRY VARIABLES
# ==========================================

# --- Dimensions Constants ---
# These stay the same for the whole bridge
top_flange_width = 104
top_flange_thickness = 2.54 # Two layers glued? Or just thick board?
web_height = 120
web_thickness = 1.27
web_separation = 80 # Distance between inside faces of webs

# --- Variable Dimensions (Optimization) ---
# ENDS (0-400, 800-1200): We remove the extra strip, leaving only small glue tabs
glue_width_ends = 20  # Small tabs just to hold the web
glue_thick_ends = 1.27

# MIDDLE (400-800): We keep the full strip (Doubler) for strength
glue_width_mid = 80 - 2.54 # Full width between webs (approx)
glue_thick_mid = 1.27

# ==========================================
# PART 2: CALCULATE PROPERTIES FUNCTION
# ==========================================

def get_section_properties(glue_w, glue_t):
    """
    Calculates I, y_bar, and Q for a given glue/doubler width.
    """
    # 1. Top Flange (Sits on top of web)
    # Centroid is web_height + glue_t + half_flange
    # (Assuming the glue/doubler sits BETWEEN web top and top flange)
    A_tf = top_flange_width * top_flange_thickness
    y_tf = web_height + 0.5 * top_flange_thickness 
    # Note: If glue is INSIDE the web span, adjust y_tf to just 'web_height + 0.5*tf'
    # Let's assume standard Pi-Beam: Web Top -> Glue/Doubler -> Top Flange
    
    # 2. Webs (2 of them)
    A_web = 2 * (web_height * web_thickness)
    y_web = 0.5 * web_height
    
    # 3. Glue/Doubler Layer (Sits at top of web, under flange)
    # We model this as 2 strips if it's tabs, or 1 wide strip if it's a doubler.
    # Since 'glue_w' represents total width of glue material:
    A_glue = glue_w * glue_t
    y_glue = web_height - 0.5 * glue_t # Sits inside the web depth (gusset/tab style)
    # OR does it sit on top?
    # Your previous code implied it sits inside/under. Let's stick to "Inside top of web".
    
    # --- Centroid (y_bar) ---
    total_area = A_tf + A_web + A_glue
    y_bar = (A_tf * y_tf + A_web * y_web + A_glue * y_glue) / total_area
    
    # --- Moment of Inertia (I) ---
    # Parallel Axis Theorem: I = I_own + A * d^2
    
    # Top Flange
    I_tf = (top_flange_width * top_flange_thickness**3)/12 + A_tf * (y_tf - y_bar)**2
    
    # Webs (x2)
    I_web = 2 * ((web_thickness * web_height**3)/12 + (web_height * web_thickness) * (y_web - y_bar)**2)
    
    # Glue/Doubler
    I_glue = (glue_w * glue_t**3)/12 + A_glue * (y_glue - y_bar)**2
    
    I_total = I_tf + I_web + I_glue
    
    # --- Q values for Shear ---
    # Q_mat (at Neutral Axis) - Area of webs below NA
    # Distance from NA to center of the web-part-below-NA
    # Height of web below NA = y_bar
    Q_mat = 2 * (web_thickness * y_bar * (y_bar / 2))
    
    # Q_glue (at the joint) - Area of top flange * distance to NA
    Q_glue = A_tf * (y_tf - y_bar)
    
    return {
        "I": I_total,
        "y_bar": y_bar,
        "y_top_dist": (web_height + top_flange_thickness) - y_bar,
        "Q_mat": Q_mat,
        "Q_glue": Q_glue
    }

# Calculate properties for BOTH sections
props_ends = get_section_properties(glue_width_ends, glue_thick_ends)
props_mid = get_section_properties(glue_width_mid, glue_thick_mid)

print(f"Ends Stiffness (I):   {props_ends['I']:.0f} mm^4")
print(f"Middle Stiffness (I): {props_mid['I']:.0f} mm^4 ({(props_mid['I']/props_ends['I'] - 1)*100:.1f}% stiffer)")

# ==========================================
# PART 3: LOAD ANALYSIS (SFD/BMD)
# ==========================================
# (Standard Load Logic - Unchanged)
bridge_length = 1200
train_weight = 400
p_wheel = [train_weight / 6] * 6
# Relative wheel positions
wheel_spacing = [0, 176, 164, 176, 164, 176] 
# Convert relative to cumulative
cum_pos = 0
x_pos_train = []
for s in wheel_spacing:
    cum_pos += s
    x_pos_train.append(cum_pos)

def get_max_forces():
    # Scans the train across the bridge
    max_V = 0
    max_M = 0
    
    # 5mm steps
    for offset in range(-1000, 1300, 5):
        # reaction calc
        current_pos = [x + offset for x in x_pos_train]
        load_on_bridge = sum([p for p, x in zip(p_wheel, current_pos) if 0 <= x <= 1200])
        moment_sum = sum([p * (1200 - x) for p, x in zip(p_wheel, current_pos) if 0 <= x <= 1200])
        
        if load_on_bridge == 0: continue
            
        Ay = moment_sum / 1200
        By = load_on_bridge - Ay
        
        # Max Shear (Approx at supports)
        max_V = max(max_V, abs(Ay), abs(By))
        
        # Max Moment (Under wheels)
        for i, w_pos in enumerate(current_pos):
            if 0 <= w_pos <= 1200:
                # Moment at w_pos
                # Cut left
                m = Ay * w_pos
                for j, w_sub in enumerate(current_pos):
                    if 0 <= w_sub < w_pos:
                        m -= p_wheel[j] * (w_pos - w_sub)
                max_M = max(max_M, m)
                
    return max_V, max_M

MAX_V, MAX_M = get_max_forces()
print(f"Max Shear: {MAX_V:.2f} N")
print(f"Max Moment: {MAX_M:.2f} Nmm")

# ==========================================
# PART 4: FACTOR OF SAFETY CHECKS
# ==========================================

# --- 1. Check SHEAR at the ENDS ---
# We use 'props_ends' because max shear happens at the ends where the bridge is lighter
tau_mat = (MAX_V * props_ends['Q_mat']) / (props_ends['I'] * (2 * web_thickness))
# Note: Glue shear width is the tab width
tau_glue = (MAX_V * props_ends['Q_glue']) / (props_ends['I'] * (2 * glue_width_ends)) 

# --- 2. Check MOMENT at the MIDDLE ---
# We use 'props_mid' because max moment happens in the middle
sigma_tens = (MAX_M * props_mid['y_bar']) / props_mid['I']
sigma_comp = (MAX_M * props_mid['y_top_dist']) / props_mid['I']

# --- 3. Check BUCKLING ---
E = 4000
mu = 0.2
a = 100 # Diaphragm spacing

# Case 1: Deck Buckling (Use MID properties - critical comp zone)
# If we have the wide doubler, the unsupported gap is smaller!
# Gap = (Web_Sep - Doubler_Width) / 2
gap_case1 = (web_separation - glue_width_mid)/2
if gap_case1 < 1: gap_case1 = 1 # Prevent div/0 if doubler fills gap

buck_1 = (((4 * math.pi**2) * E) / (12 * (1 - mu**2))) * (top_flange_thickness / gap_case1)**2

# Case 2: Wing Buckling (Use MID properties)
# Overhang = (Top_Flange_Width - Web_Sep)/2
overhang = (top_flange_width - web_separation)/2
buck_2 = (((0.425 * math.pi**2) * E) / (12 * (1 - mu**2))) * (top_flange_thickness / overhang)**2

# Case 3: Web Flexural (Use MID properties - y_bar changes!)
h_comp = web_height - props_mid['y_bar']
buck_3 = (((6 * math.pi**2) * E) / (12 * (1 - mu**2))) * (web_thickness / h_comp)**2

# Case 4: Web Shear (Use END properties - critical shear zone)
# Use full web height
buck_4 = (((5 * math.pi**2) * E) / (12 * (1 - mu**2))) * ( (web_thickness/web_height)**2 + (web_thickness/a)**2 )


# ==========================================
# PART 5: OUTPUT RESULTS
# ==========================================

# Limits
LIMIT_TENS = 30
LIMIT_COMP = 6
LIMIT_SHEAR_MAT = 4
LIMIT_SHEAR_GLUE = 2

print("\n--- FACTOR OF SAFETY REPORT ---")
print(f"FOS Tension (Mid):      {LIMIT_TENS / sigma_tens:.2f}")
print(f"FOS Compression (Mid):  {LIMIT_COMP / sigma_comp:.2f}")
print(f"FOS Shear Mat (End):    {LIMIT_SHEAR_MAT / tau_mat:.2f}")
print(f"FOS Shear Glue (End):   {LIMIT_SHEAR_GLUE / tau_glue:.2f}")
print("-" * 30)
print(f"FOS Buckling 1 (Deck):  {buck_1 / sigma_comp:.2f}")
print(f"FOS Buckling 2 (Wing):  {buck_2 / sigma_comp:.2f}")
print(f"FOS Buckling 3 (Web):   {buck_3 / sigma_comp:.2f}")
print(f"FOS Buckling 4 (Shear): {buck_4 / tau_mat:.2f}")

print("\nNOTE: FOS Buckling 1 is likely very high because the doubler strip reduces the gap.")