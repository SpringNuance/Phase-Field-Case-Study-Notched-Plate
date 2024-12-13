import os 
import shutil 
# Going to current directory
os.chdir(os.getcwd())

working_CAE_path = 'process_input_file_to_UEL_millimeter'
output_simulation_path = 'simulation_results_millimeter'
inp_file_name = 'phase_field_3D_millimeter'

########################################
# STEP 0: Deleting previous sim files  #
########################################

# Now we would delete everything in output_simulation_path 
files_to_keep = []

# List of files to delete
files_to_delete = [f for f in os.listdir(output_simulation_path) if f not in files_to_keep]

# Loop through each file in the list of files to delete
for file in files_to_delete:
    # Construct the full path to the file
    file_path = os.path.join(output_simulation_path, file)
    # Check if the file exists
    if os.path.exists(file_path):
        # Remove the file
        os.remove(file_path)
print(f"Successfully deleting all files in previous simulation_results")

################################
# STEP 1: Copying the inp file #
################################

combined_CAE_inp_path = f'{working_CAE_path}/{inp_file_name}.inp'
combined_original_inp_path = f'{output_simulation_path}/{inp_file_name}.inp'
combined_UEL_inp_path = f'{output_simulation_path}/{inp_file_name}_UEL.inp'

shutil.copyfile(combined_CAE_inp_path, combined_original_inp_path)

#####################################
# STEP 2: Combining the subroutines #
#####################################

# List of Fortran file names in correct order
file_names = ['srt_part0_description.f90',
              'srt_part1_precision.f90', 
              'srt_part2_common_block.f90',
              'srt_part3_uexternaldb.f90',
              'srt_part4_shape_grad.f90',
              'srt_part5_uel.f90',
              'srt_part6_umat_elastic.f90',
              #'srt_part7_uhard.f90',
              #'srt_part8_umatht_diffusion.f90',
              'srt_part9_umat_visual.f90',
              ]  

# Define the output file
combined_subroutine_output_path = f'{output_simulation_path}/phase_field_combined.f90'

# Open the output file for writing
with open(combined_subroutine_output_path, 'w', encoding="UTF-8") as outfile:
    # Loop through each file in the provided list
    for filename in file_names:
        # Open and read the contents of each .f90 file
        with open(filename, 'r', encoding="UTF-8") as infile:
            # Write the contents to the output file
            outfile.write(infile.read() + '\n\n')

print("All files have been successfully concatenated into:", combined_subroutine_output_path)

###########################################
# STEP 3: Modifying normal inp to UEL inp #
###########################################

import numpy as np

def return_UEL_property(description_properties_dict): 
    mechanical_properties_list = list(description_properties_dict["mechanical_properties"].values())
    mechanical_description_list = list(description_properties_dict["mechanical_properties"].keys())
    flow_curve_true_strain = description_properties_dict["flow_curve_properties"]["equivalent_plastic_strain"]
    flow_curve_true_stress = description_properties_dict["flow_curve_properties"]["equivalent_plastic_stress"]
    
    flow_curve_zipped = []
    for stress, strain in zip(flow_curve_true_stress, flow_curve_true_strain):
        flow_curve_zipped.append(stress)
        flow_curve_zipped.append(strain)
        
    
    phase_field_properties_list = list(description_properties_dict["phase_field_damage_properties"].values())
    phase_field_description_list = list(description_properties_dict["phase_field_damage_properties"].keys())
    hydrogen_diffusion_properties_list = list(description_properties_dict["hydrogen_diffusion_properties"].values())
    hydrogen_diffusion_description_list = list(description_properties_dict["hydrogen_diffusion_properties"].keys())

    # Abaqus needs to define 8 properties each line
    mech_prop_num_lines = int(np.ceil(len(mechanical_properties_list)/8))
    mech_prop_num_properties = int(mech_prop_num_lines*8)
    phase_field_prop_num_lines = int(np.ceil(len(phase_field_properties_list)/8))
    phase_field_prop_num_properties = int(phase_field_prop_num_lines*8)
    hydrogen_diffusion_prop_num_lines = int(np.ceil(len(hydrogen_diffusion_properties_list)/8))
    hydrogen_diffusion_prop_num_properties = int(hydrogen_diffusion_prop_num_lines*8)
    flow_curve_num_lines = int(np.ceil(len(flow_curve_zipped)/8))
    flow_curve_num_properties = int(flow_curve_num_lines*8)

    total_num_properties = mech_prop_num_properties + phase_field_prop_num_properties +\
                           hydrogen_diffusion_prop_num_properties + flow_curve_num_properties

    UEL_property = [
        "*******************************************************",
        "*UEL PROPERTY, ELSET=SOLID                             ",
    ]
    
    # The last line would be padded with 0.0 and their corresponding description would be "none"
    # If the number of properties is not a multiple of 8

    # For mechanical properties
    UEL_property.append("**")
    UEL_property.append("** =====================")
    UEL_property.append("**")
    UEL_property.append("** MECHANICAL PROPERTIES")
    UEL_property.append("**")
    for line_index in range(mech_prop_num_lines):
        if line_index != mech_prop_num_lines - 1:
            subset_properties = mechanical_properties_list[line_index*8:(line_index+1)*8]
            subset_description = mechanical_description_list[line_index*8:(line_index+1)*8]
            UEL_property.append(", ".join(subset_properties))
            UEL_property.append("** " + ", ".join(subset_description[0:4]))
            UEL_property.append("** " + ", ".join(subset_description[4:8]))
        else:
            subset_properties = mechanical_properties_list[line_index*8:] + ["0.0"]*(8-len(mechanical_properties_list[line_index*8:]))
            subset_description = mechanical_description_list[line_index*8:] + ["none"]*(8-len(mechanical_description_list[line_index*8:]))
            UEL_property.append(", ".join(subset_properties))
            UEL_property.append("** " + ", ".join(subset_description[0:4]))
            UEL_property.append("** " + ", ".join(subset_description[4:8]))
    
    UEL_property.append("**")
    UEL_property.append("** ======================")
    # For phase field properties
    UEL_property.append("**")
    UEL_property.append("** PHASE FIELD PROPERTIES")
    UEL_property.append("**")
    for line_index in range(phase_field_prop_num_lines):
        if line_index != phase_field_prop_num_lines - 1:
            subset_properties = phase_field_properties_list[line_index*8:(line_index+1)*8]
            subset_description = phase_field_description_list[line_index*8:(line_index+1)*8]
            UEL_property.append(", ".join(subset_properties))
            UEL_property.append("** " + ", ".join(subset_description[0:4]))
            UEL_property.append("** " + ", ".join(subset_description[4:8]))
        else:
            subset_properties = phase_field_properties_list[line_index*8:] + ["0.0"]*(8-len(phase_field_properties_list[line_index*8:]))
            subset_description = phase_field_description_list[line_index*8:] + ["none"]*(8-len(phase_field_description_list[line_index*8:]))
            UEL_property.append(", ".join(subset_properties))
            UEL_property.append("** " + ", ".join(subset_description[0:4]))
            UEL_property.append("** " + ", ".join(subset_description[4:8]))

    # For hydrogen diffusion properties
    UEL_property.append("**")
    UEL_property.append("** =============================")
    UEL_property.append("**")
    UEL_property.append("** HYDROGEN DIFFUSION PROPERTIES")
    UEL_property.append("**")

    for line_index in range(hydrogen_diffusion_prop_num_lines):
        if line_index != hydrogen_diffusion_prop_num_lines - 1:
            subset_properties = hydrogen_diffusion_properties_list[line_index*8:(line_index+1)*8]
            subset_description = hydrogen_diffusion_description_list[line_index*8:(line_index+1)*8]
            UEL_property.append(", ".join(subset_properties))
            UEL_property.append("** " + ", ".join(subset_description[0:4]))
            UEL_property.append("** " + ", ".join(subset_description[4:8]))
        else:
            subset_properties = hydrogen_diffusion_properties_list[line_index*8:] + ["0.0"]*(8-len(hydrogen_diffusion_properties_list[line_index*8:]))
            subset_description = hydrogen_diffusion_description_list[line_index*8:] + ["none"]*(8-len(hydrogen_diffusion_description_list[line_index*8:]))
            UEL_property.append(", ".join(subset_properties))
            UEL_property.append("** " + ", ".join(subset_description[0:4]))
            UEL_property.append("** " + ", ".join(subset_description[4:8]))

    # For flow curve properties
    # Important: DO NOT PAD THIS TIME FOR FLOW CURVE
    UEL_property.append("**")
    UEL_property.append("** =====================")
    UEL_property.append("**")
    UEL_property.append("** FLOW CURVE PROPERTIES")
    UEL_property.append("**")
    
    UEL_property.append("** True stress (Pa) - PEEQ (dimless) value pairs")
    for line_index in range(flow_curve_num_lines):
        if line_index != flow_curve_num_lines - 1:
            str_values = [str(value) for value in flow_curve_zipped[line_index*8:(line_index+1)*8]]
            UEL_property.append(", ".join(str_values))
        else:
            str_values = [str(value) for value in flow_curve_zipped[line_index*8:]]
            UEL_property.append(", ".join(str_values))
    
    UEL_property.append("**")
    UEL_property.append("*******************************************************")

    return UEL_property, total_num_properties

# 210000., 0.3, 0.05, 2.7, 0.0127  
description_properties_dict = {
    "mechanical_properties": {
        "Young's modulus E (MPa)": "210000.0d0", # correct
        "Poisson ratio nu (dimless)": "0.3", # correct
    },
    # 0.05, 2.7, 0.0127, Gibbs free 3.0d7 N.mm/mol   
    # Gc(0) = 2.7 MPa / mm = 2700 Pa / m = 2700 J/m^2
    # Note: J = N*m
    "phase_field_damage_properties": {
        "length scale lc (mm)": "0.05d0", # correct
        "critical energy release rate Gc (N*mm/mm^2)": "2.7e0", # correct
        "well-conditioning parameter xkap (dimless)": "1.0e-7", # correct
        "fitting slope to DFT data X (chi) (dimless)": "0.89",  # correct
        # For iron, chi = 0.89
        # For nickel, chi = 0.41
        # For aluminum, chi = 0.67
        "Gibbs free energy difference between the decohering interface and the surrounding material delta_g_b0 (N*mm/mol)": "3.0e7", # correct
    },
    # R = 8314.5 N.mm/(mol.K), T = 300 K, VH = 2000 mm^3/mol, DL = 0.0127 mm^2/s
    "hydrogen_diffusion_properties": {
        "universal gas constant R (N*mm)/(mol*K)": "8314.5e0", # correct
        "temperature T (K)": "300e0", # correct
        "partial molar volume VH (mm^3/mol)": "2000.0e0", # correct
        "diffusion coefficient DL (mm^2/s)": "0.0127e0", # correct
        "Avogadro constant NA (1/mol)": "6.023e23",
        "number of solvent atoms N_L (1/m^3)": "9.24e28",
        "number of interstitial sites per solvent atom - beta (dimless)": "6.0",
        "number of H per trap site of dislocation - alpha_dis (dimless)": "1.0",
        "number of H per trap site of grain boundaries - alpha_gb (dimless)": "1.0",
        "number of H per trap site of carbides - alpha_carb (dimless)": "1.0",
        "lattice parameter a_lattice (m)": "2.86e-10",
        "binding energy of H to dislocations trap WB_dis (N*m/mol)": "-20.2e3",
        "binding energy of H to grain boundaries trap WB_gb (N*m/mol)": "-58.6e3",
        "binding energy of H to carbides trap WB_carb (N*m/mol)": "-11.5e3",
        "number of trap type of dislocations NT_dis (1/m^2)": "0",
        "number of trap type of grain boundaries NT_gb (1/m^2)": "8.464e22",
        "number of trap type of carbides NT_carb (1/m^2)": "8.464e26",
        "gamma fitting parameter (1/m^2)": "2e16",
        "dislocation density for annealed material rho_d0 (1/m^2)": "1e10"
    },
    "flow_curve_properties": {
        "equivalent_plastic_strain": 
            [0.00e+00, 2.00e-04, 4.00e-04, 6.00e-04, 8.00e-04, 1.00e-03, 2.00e-03, 3.00e-03,
                                4.00e-03, 5.00e-03, 6.00e-03, 7.00e-03, 8.00e-03, 9.00e-03, 1.00e-02, 2.00e-02,
                                3.00e-02, 4.00e-02, 5.00e-02, 6.00e-02, 7.00e-02, 8.00e-02, 9.00e-02, 1.00e-01,
                                1.25e-01, 1.50e-01, 1.75e-01, 2.00e-01, 2.25e-01, 2.50e-01, 2.75e-01, 3.00e-01,
                                3.25e-01, 3.50e-01, 3.75e-01, 4.00e-01, 4.25e-01, 4.50e-01, 4.75e-01, 5.00e-01,
                                5.25e-01, 5.50e-01, 5.75e-01, 6.00e-01, 6.25e-01, 6.50e-01, 6.75e-01, 7.00e-01,
                                7.25e-01, 7.50e-01, 7.75e-01, 8.00e-01, 8.25e-01, 8.50e-01, 8.75e-01, 9.00e-01,
                                9.25e-01, 9.50e-01, 9.75e-01, 1.00e+00, 1.05e+00, 1.10e+00, 1.15e+00, 1.20e+00,
                                1.25e+00, 1.30e+00, 1.35e+00, 1.40e+00, 1.45e+00, 1.50e+00, 1.55e+00, 1.60e+00,
                                1.65e+00, 1.70e+00, 1.75e+00, 1.80e+00, 1.85e+00, 1.90e+00, 1.95e+00, 2.00e+00,
                                2.05e+00, 2.10e+00, 2.15e+00, 2.20e+00, 2.25e+00, 2.30e+00, 2.35e+00, 2.40e+00,
                                2.45e+00, 2.50e+00, 2.55e+00, 2.60e+00, 2.65e+00, 2.70e+00, 2.75e+00, 2.80e+00,
                                2.85e+00, 2.90e+00, 2.95e+00, 3.00e+00],
        "equivalent_plastic_stress": 
            [418501696.8, 419158989.7, 419816282.5, 420473122.5, 421129509.7, 421785444.7, 425058345.2, 428319983.9,
                                431570399.5, 434809630.8, 438037716.1, 441254693.8, 444460602.2, 447655479.2, 450839363, 482081826.2,
                                512265755.3, 541427014.8, 569600254, 596818948.3, 623115438.8, 648520971, 673065731.8, 696778885,
                                752606048.5, 803823955.2, 850813159.2, 893922794.7, 933473170.4, 969758149.4, 1003047333, 1033588062,
                                1061607258, 1087313105, 1110896601, 1132532974, 1152382983, 1170594117, 1187301685, 1202629827, 1216692431,
                                1229593985, 1241430349, 1252289467, 1262252024, 1271392043, 1279777434, 1287470503, 1294528409, 1301003593,
                                1306944167, 1312394269, 1317394395, 1321981695, 1326190254, 1330051343, 1333593648, 1336843490, 1339825015,
                                1342560377, 1347372227, 1351422317, 1354831240, 1357700499, 1360115527, 1362148234, 1363859144, 1365299202,
                                1366511285, 1367531484, 1368390176, 1369112929, 1369721263, 1370233293, 1370664263, 1371027007, 1371332325,
                                1371589308, 1371805609, 1371987667, 1372140904, 1372269882, 1372378441, 1372469815, 1372546723, 1372611456,
                                1372665941, 1372711800, 1372750400, 1372782889, 1372810235, 1372833251, 1372852624, 1372868930, 1372882654,
                                1372894206, 1372903929, 1372912113, 1372919001, 1372924799],

    }

}

UEL_PROPERTY, total_num_properties = return_UEL_property(description_properties_dict)


ndim = 3
nnodes = 8
nstatev = 15
nsvars = nnodes * nstatev

def return_user_element(total_num_properties, ndim, nnodes, nsvars):
    USER_ELEMENT = [
        "*************************************************",
       f"*User element, nodes={nnodes}, type=U1, properties={total_num_properties}, coordinates={ndim}, variables={nsvars}",
        "1, 2, 3",
        "1, 4",
        "1, 11",
        "*************************************************"
    ]
    return USER_ELEMENT

# The user element variables can be defined so as to order the degrees of freedom on the element 
# in any arbitrary fashion. You specify a list of degrees of freedom for the first node on the element. 
# All nodes with a nodal connectivity number that is less than the next connectivity number for which a 
# list of degrees of freedom is specified will have the first list of degrees of freedom. The second list 
# of degrees of freedom will be used for all nodes until a new list is defined, etc. If a new list of degrees 
# of freedom is encountered with a nodal connectivity number that is less than or equal to that given with 
# the previous list, the previous list's degrees of freedom will be assigned through the last node of 
# the element. This generation of degrees of freedom can be stopped before the last node on the element 
# by specifying a nodal connectivity number with an empty (blank) list of degrees of freedom.

USER_ELEMENT = return_user_element(total_num_properties, ndim, nnodes, nsvars)

def return_depvar(nstatev):
    # DEPVAR = [
    #     "*Depvar       ",
    #     f"  {nstatev},      ",       
    #     "1, AR1_sig11, AR1_sig11   ",
    #     "2, AR2_sig22, AR2_sig22   ",
    #     "3, AR3_sig33, AR3_sig33   ",
    #     "4, AR4_sig12, AR4_sig12   ",
    #     "5, AR5_sig13, AR5_sig13   ",
    #     "6, AR6_sig23, AR6_sig23   ",
    #     "7, AR7_eps11, AR7_eps11   ",
    #     "8, AR8_eps22, AR8_eps22   ",
    #     "9, AR9_eps33, AR9_eps33   ",
    #     "10, AR10_eps12, AR10_eps12   ",
    #     "11, AR11_eps13, AR11_eps13   ",
    #     "12, AR12_eps23, AR12_eps23   ",
    #     "13, AR13_eqplas, AR13_eqplas   ",
    #     "14, AR14_sig_H, AR14_sig_H   ",
    #     "15, AR15_phi, AR15_phi   ", 
    #     "16, AR16_history, AR16_history   ",
    #     "17, AR17_rho_d, AR17_rho_d   ",
    #     "18, AR18_C, AR18_C   ",	
    #     "19, AR19_CL, AR19_CL   ",   
    #     "20, AR20_CT, AR20_CT   "   
    # ]


    DEPVAR = [
        "*Depvar       ",
        f"  {nstatev},      ",       
        "1, AR1_sig11, AR1_sig11   ",
        "2, AR2_sig22, AR2_sig22   ",
        "3, AR3_sig33, AR3_sig33   ",
        "4, AR4_sig12, AR4_sig12   ",
        "5, AR5_sig13, AR5_sig13   ",
        "6, AR6_sig23, AR6_sig23   ",
        "7, AR7_eps11, AR7_eps11   ",
        "8, AR8_eps22, AR8_eps22   ",
        "9, AR9_eps33, AR9_eps33   ",
        "10, AR10_eps12, AR10_eps12   ",
        "11, AR11_eps13, AR11_eps13   ",
        "12, AR12_eps23, AR12_eps23   ",
        "13, AR13_phi, AR13_phi      ",
        "14, AR14_sig_H, AR14_sig_H   ",
        "15, AR15_CL, AR15_CL   ", 
    ]
    return DEPVAR

DEPVAR = return_depvar(nstatev)

def constructing_visualization_mesh(input_file_path):
    """
    Visualization mesh is simply a copy of the original mesh
    Except that it has different element ID and a standard Abaqus element type
    in order to store values of the user element at integration points
    We cannot view the user element directly in Abaqus/Viewer since UEL subroutine
    only has knowledge to the stiffness matrix amatrx and residual right hand side rhs
    """

    # Open the file
    with open(input_file_path, 'r') as fid:
        flines = fid.readlines()

    # Process the lines
    flines = [line.strip() for line in flines]
    flines_upper = [line.upper() for line in flines]
    start_elements = [i for i, line in enumerate(flines_upper) if '*ELEMENT' in line and '*ELEMENT OUTPUT' not in line]
    start_element_index = start_elements[0]
    element_indices = [] # list of element index
    element_connvectivity = [] # list of of list of nodes that make up the element

    #print("The starting element index is: ", start_element_index)
    #print(start_element_index)

    mesh_index = start_element_index + 1

    while flines_upper[mesh_index][0] != "*" and flines_upper[mesh_index][0] != " ":
       
        # remove all empty spaces and split by comma
        # each line look like this: 1,    35,     2,    36,  2503,  5502,  5503,  5504,  5505
        split_line = flines_upper[mesh_index].replace(" ", "").split(",")
        
        element_indices.append(split_line[0])
        element_connvectivity.append(split_line[1:])
        mesh_index += 1
    
    end_element_index = mesh_index

    # print("The element indices are: ", element_indices)
    # print("The element connectivity are: ", element_connvectivity)

    # Now we would count the number of elements 
    num_elements = len(element_indices)
    #print("The number of elements is: ", num_elements)  

    # We would start to reconstruct an identical mesh, except that we replace node indices with new indices
    # New indices should start from a power of 10 that is at least greater than the number of elements
    # For example, if num_elements = 59, then the new indices should start from 101 to 159
    # If num_elements = 128, then the new indices should start from 1001 to 1128
    # If num_elements = 3456, then the new indices should start from 10001 to 13456

    # We would find the order of the number of elements
    order = int(np.floor(np.log10(num_elements)))
    offset = 10 * 10**order + 1
    #print("The offset is: ", offset)

    new_element_indices = [i + offset for i in range(num_elements)]

    #print("The new element indices are: ", new_element_indices)

    visualization_mesh = [
        "*ELEMENT, TYPE=C3D8, ELSET=VISUALIZATION",
    ]

    for i in range(num_elements):
        reconstructed_line = ", ".join([str(value) for value in [new_element_indices[i]] + element_connvectivity[i]])
        visualization_mesh.append(reconstructed_line)
    return visualization_mesh, start_element_index, end_element_index

VISUALIZATION_MESH, start_element_index, end_element_index =\
    constructing_visualization_mesh(combined_original_inp_path)

# Open the file
with open(combined_original_inp_path, 'r') as fid:
    flines = fid.readlines()

# Process the lines
flines = [line.strip() for line in flines]

# Now, we would reconstruct the input file as follows

# 1. Replace the original *ELEMENT section line with this line
#    *Element, type=U1, elset=SOLID

# 2. Add the USER ELEMENT section just above the *ELEMENT section in step 1

# 3. Add the UEL property just below the element connectivity section of the original mesh *ELEMENT above

# 4. Add the visualization mesh just below the USER ELEMENT section

# 5. Finally, in **SECTION, we change to this 
# ** Section: Section-1
# *Solid Section, elset=VISUALIZATION, material=(whatever material you define here)

# 6. We would also modify the *Depvar section to include the key descriptions

# do step 1

import copy
flines_new = copy.deepcopy(flines)
flines_new[start_element_index] = "*ELEMENT, TYPE=U1, ELSET=SOLID"

# do step 2 to 4

flines_new = flines_new[:start_element_index] + USER_ELEMENT\
            + flines_new[start_element_index:end_element_index] + UEL_PROPERTY\
            + VISUALIZATION_MESH + flines_new[end_element_index:]

# do step 5
solid_section_index = [i for i, line in enumerate(flines_new) if '*SOLID SECTION' in line.upper()][0]

# we should change this line
# *Solid Section, elset=<whatever name>, material=<whatever name>
# to *Solid Section, elset=VISUALIZATION, material=<whatever name>

# find the index where the word material is found
starting_index_string = flines_new[solid_section_index].find("material")
# print("The starting index of the word material is: ", starting_index_string)
flines_new[solid_section_index] = "*Solid Section, elset=VISUALIZATION, " + flines_new[solid_section_index][starting_index_string:]

# do step 6

# find the index of the *Depvar section
depvar_index = [i for i, line in enumerate(flines_new) if '*DEPVAR' in line.upper()][0]
flines_new = flines_new[:depvar_index] + DEPVAR + flines_new[depvar_index+2:]
# write this to see how it looks like
with open(combined_UEL_inp_path, 'w') as fid:
    for line in flines_new:
        fid.write(line + "\n")

        