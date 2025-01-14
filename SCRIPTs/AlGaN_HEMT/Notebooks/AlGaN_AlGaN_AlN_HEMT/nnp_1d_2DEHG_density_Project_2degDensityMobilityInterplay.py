#!/usr/bin/env python
# coding: utf-8

# # Project: AlGaN/AlGaN/AlN HEMT device simulation using nexnano++ and nextnanopy

# ## 1. General settings

# In[41]:


submit_cluster = 0 # Submit the job to the cluster or not.


# ### 1.1 Import modules

# In[42]:


if not submit_cluster:
    get_ipython().run_line_magic('reload_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')


# ### 1.1.1 Adding local module path to python module search path

# In[43]:


from pathlib import Path
import sys, os
current_script_path = Path().absolute()
module_path = os.path.abspath(os.path.join(current_script_path,'../..'))
sys.path.append(module_path)


# #### 1.1.2 Import global modules

# In[44]:


import shutil
import json
import nextnanopy as nn
from nextnanopy.utils.misc import mkdir_if_not_exist
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib.tri as tri 
import pandas as pd
from matplotlib.widgets import Slider


# #### 1.1.2 Import local defined modules

# In[45]:


from src.PlotFunctions import general_plot_functions, Plot1DFuns, PlotQuasi3DFuns
lpltgen = general_plot_functions()
lplt1d = Plot1DFuns()
lpltq3d = PlotQuasi3DFuns()


# ### 1.2 Matplotlib settings

# In[46]:


params = {'figure.figsize': (8, 6), 'legend.fontsize': 18, 'axes.labelsize': 24, 'axes.titlesize': 24,
          'xtick.labelsize':24, 'xtick.major.width':2, 'xtick.major.size':5, 'ytick.labelsize': 24,
          'ytick.major.width':2, 'ytick.major.size':5, 'xtick.minor.width':2, 'xtick.minor.size':3,
          'ytick.minor.width':2, 'ytick.minor.size':3, 'errorbar.capsize':2}
plt.rcParams.update(params)
plt.rc('font', size=24)


# ### 1.3 nextnanopy settings

# In[47]:


#%% ===========================================================================
#++++++++++++++++++++++++++++++++++++++++++++++
# Software
#++++++++++++++++++++++++++++++++++++++++++++++
software_ = 'nextnano++'

#++++++++++++++++++++++++++++++++++++++++++++++
# NextNANO configuration file
#++++++++++++++++++++++++++++++++++++++++++++++
if submit_cluster:
    software_path = "/sfihome/local/rocks7-temp/codes/nextnano"
    software_path = "/sfihome/badal.mondal/local/Softwares/nextnano"
    software_version = 'v1.21.24'
    software_executable = os.path.join(software_path, f"RHEL/{software_version}/nextnano++_gcc_RHEL76_nextgen.exe")
    software_license = os.path.join(software_path, "License/license.txt") 
    software_database = os.path.join(software_path, "Database/database_nnp.in")
    
    nn.config.set(software_,'exe',software_executable)
    nn.config.set(software_,'runmode','--old')
    nn.config.set(software_,'license',software_license)
    nn.config.set(software_,'database',software_database)
    nn.config.set(software_,'outputdirectory','')
    
    nn.config.save() #save permanently

print(f'-nextnano config: {nn.config}')


# ### 1.4 Set tasks to perform

# In[48]:


run_sim = 0 # run single simulations
run_sweep = 0 # run sweep simulations
run_sim_specific_sample = False # run single simulation for specific sample device
run_sweep_specific_sample = False # run sweep simulation for specific sample device
create_data_sweep = 1 # create required data sheet from sweep simulation results
append_new_sweep_data_2_data_base = 0 # append the new sweep data to existing data files. if False, rewrite the complete excel data file.
do_plot = 1  # plot results of single simulation
do_plot_sweep = 1 # plot results of sweep simulation
savefigure = True # save the figures generated


# ### 1.5 Input and output directories/files

# In[49]:


#++++++++++++++++++++++++++++++++++++++++++++++
# Parent simulation folder
#++++++++++++++++++++++++++++++++++++++++++++++
sim_dir = str(current_script_path).split('SCRIPTs')[0]

#++++++++++++++++++++++++++++++++++++++++++++++
# Set input file 
#++++++++++++++++++++++++++++++++++++++++++++++
# extension of the input file
FileExtension = '.in' 
# project ID to track the simulations project-wise. 
# Note: During 1st set of simulations we did not have the project id. 
# The simulations were already done with out project id.
my_project_id_old = ''
my_project_id = '2DEG_DensityMobilityInterplay'
upgrade_figs_folder = not bool(len(my_project_id_old.strip()))
replace_figs_path = ('nnp', f'nnp/{my_project_id}')
# inner location of input file within INPUTs folder
folder_input_lower_level = f'AlGaN_HEMT/AlGaN_AlGaN_AlN_HEMT/nnp/{my_project_id}/2DEHG_density/OneDim'
folder_input_lower_level_old = f'AlGaN_HEMT/AlGaN_AlGaN_AlN_HEMT/nnp/{my_project_id_old}/2DEHG_density/OneDim'
# complete path of input file folder
folder_input_ = os.path.join(sim_dir, 'INPUTs', folder_input_lower_level_old)
# list of name of input files
input_filename_list = [r'sim.in']
# list of all input files complete path
input_files_dest = [os.path.join(folder_input_, fname) for fname in input_filename_list]

#++++++++++++++++++++++++++++++++++++++++++++++
# Create input file for specific sample
#++++++++++++++++++++++++++++++++++++++++++++++
# settings for specific sample(s)
specific_samples = {'Sample1': {'AlContentBarrier': 0.75, 
                                'AlContentChannel': 0.55, 
                                'ThicknessAlGaNBarrier': 42}
                    }
if run_sim_specific_sample or run_sweep_specific_sample:
    # base input file to use for creating specific sample
    base_input_file = input_filename_list[0]
    base_input_file_name = base_input_file.replace(FileExtension, '')
    input_path_ = os.path.join(folder_input_, base_input_file)
    # read base input file
    input_file = nn.InputFile(input_path_)
    if (not run_sim) and (not run_sweep): 
        input_files_dest = []
    # set-update parameters for specific sample(s)
    # temporarily save generated input files
    for specific_sample_params in specific_samples.values():
        add_ftext = ''
        for key, val in specific_sample_params.items():
            input_file.set_variable(key, value=val)
            add_ftext += f'_{key}_{val:.2f}'
        input_file_specific = f'specific_sample_{base_input_file_name}_{add_ftext}.in'
        input_path_specific = os.path.join(folder_input_, 'tmp', input_file_specific)
        input_file.save(input_path_specific, overwrite=True)
        input_files_dest.append(input_path_specific)
            
#++++++++++++++++++++++++++++++++++++++++++++++
# Simulation output path
#++++++++++++++++++++++++++++++++++++++++++++++
folder_output_ = os.path.join(sim_dir, 'OUTPUTs', folder_input_lower_level_old)
mkdir_if_not_exist(folder_output_)

#++++++++++++++++++++++++++++++++++++++++++++++
# Postprocess data path
#++++++++++++++++++++++++++++++++++++++++++++++
folder_post_data_ = os.path.join(sim_dir, 'DATAs', folder_input_lower_level)
mkdir_if_not_exist(folder_post_data_)

#++++++++++++++++++++++++++++++++++++++++++++++
# Figure folder path
#++++++++++++++++++++++++++++++++++++++++++++++
folder_figs_ = os.path.join(sim_dir, 'FIGs', folder_input_lower_level) 
mkdir_if_not_exist(folder_figs_)

#++++++++++++++++++++++++++++++++++++++++++++++
# Specify output image format: 
# ['.pdf','.svg','.jpg','.png']
#++++++++++++++++++++++++++++++++++++++++++++++
FigFormat = '.png'
FigFormat4Paper = '.eps'
FigDpi = 300
color_map = 'viridis'


# ### 1.5 Sweep parameters

# In[50]:


# Specify sweep variables:
# SweepVariablesSet = {name_of_set: {sweep_variable: [start value, end value, number of points]}}
# name_of_set == name of the set. choosen by user. arbitrary. 
# sweep_variable == name of sweep variable that you want to sweep. name must be in the input file.
SweepVariablesSet = { 
    #----------------------- Project: 2DEG_DensityMobilityInterplay -----------------------
    'NeumannBCEndDevice' : 
                    {'ThicknessAlNSub'          : [300, 600, 1200, 1500, 1800, 2000]},
    'SchottkyBarrierEndDevice' : 
                    {'ThicknessAlNSub'          : [300, 600, 1200, 1500, 1800, 2000]},
    'SchottkyContactScan'       : 
                    {'SchottkyBarrierHeight'    : np.linspace( 0.0,   4.0, num= 9)},
    'BandgapBowingScan'         :
                    {'AlGaNbandgapBowing'       : [0.5, 0.7, 0.9, 1.1, 1.3]},
    'PyroelectricBowingScan'    : 
                    {'AlGaNpyroelectricBowing'  : [-0.1, -0.021, 0.0]},
    'BarrierThicknessScan' : 
                    {'ThicknessAlGaNBarrier'    : [10, 25, 50, 75, 100, 150, 200, 250, 300]},
    'ChannelThicknessScan'      : 
                    {'ThicknessAlGaNChannel'    : np.linspace(50.0, 600.0, num=23)},
    'Al85Const2DEGReverseEng' : 
                    {'ThicknessAlGaNBarrier'    : [50, 75, 100, 150, 200, 250, 300, 350]},  
    'TemperatureScan'           : 
                    {'Temperature'              : np.linspace(50.0, 600.0, num=23)},
    'GaNChannelStudy' : 
                    {'Temperature'              : np.linspace(50.0, 600.0, num=23)},
    'CompositionBThicknessScan' : 
                    {'AlContentBarrier'         : np.linspace( 0.5,   1.0, num=11), 
                     'AlContentChannel'         : np.linspace( 0.5,   1.0, num=11),
                     'ThicknessAlGaNBarrier'    : np.linspace( 5.0,  50.0, num=10)
                     }
}

# For specific samples only sweep channel thickness and composition gradient
SpecificSamplesCase = ['BarrierThicknessScan', 'ChannelThicknessScan','TemperatureScan']

# Note: negative IntentionalDopingBCLength means doping in Barrier, positive means in Channel
# SwitchKeyVal: for IntentionalDopingConcScan is length of doping region in nm
TemporaryInputFiles4 = {'SchottkyBarrierEndDevice': {'SwitchKey': ['end_device_air'], 
                                                     'SwitchKeyVal': [1],
                                                     'input_file_name': 'end_bc_Schottky'},
                       'GaNChannelStudy'              : {'SwitchKey': ['AlContentBarrier','AlContentChannel',
                                                                       'AlContentSub', 'ThicknessAlNSub', 
                                                                       'ThicknessAlGaNBarrier','DHG_QR_width'], 
                                                         'SwitchKeyVal': [0.25, 0.0, 0.0, 5, 25, 4],
                                                         'input_file_name': 'GaN_C'},
                       'Al85Const2DEGReverseEng'      : {'SwitchKey': ['AlContentBarrier','AlContentChannel'], 
                                                         'SwitchKeyVal': [1.0, 0.85],
                                                         'input_file_name': 'Al85_C'}
                       }

# Here we create mapping between the long data sheet name to short one.
# Windows system can not handle sheet name > 31 charachters
MappingShortDataSheetName = {
                             'sim_sweep__AlGaNpyroelectricBowing':
                             {'abbr':'Psp_bowing_scan', 
                              'description': 'AlGaN pyroelectric bowing parameter variation simulations for AlN(25nm)/Al75Ga25N(300nm)/AlN(300nm)'
                             },
                             'sim_sweep__AlGaNbandgapBowing': 
                            {'abbr':'Eg_bowing_scan', 
                              'description': 'AlGaN bandgap bowing parameter variation simulations for AlN(25nm)/Al75Ga25N(300nm)/AlN(300nm)'
                             },
                            'sim_sweep__SchottkyBarrierHeight':
                             {'abbr':'Hschottky_scan', 
                              'description': 'Schottky barrier height variation simulations for AlN(25nm)/Al75Ga25N(300nm)/AlN(300nm)'
                             },
                            'GaN_C_sweep__Temperature':
                             {'abbr':'GaN_T_scan', 
                              'description': 'Temperature variation simulations for Al75Ga25N(25nm)/GaN(300nm)/GaN'
                             },
                             'sim_sweep__Temperature':
                             {'abbr':'T_scan', 
                              'description': 'Temperature variation simulations for AlN(25nm)/Al75Ga25N(300nm)/AlN(300nm)'
                             },
                            'Al85_C_sweep__ThicknessAlGaNBarrier': 
                             {'abbr':'Al85GaN_Lb_scan', 
                              'description': 'AlGaN barrier thickness variation simulations for AlN(Lb)/Al85Ga15N(300nm)/AlN(300nm)'
                             },
                            'sim_sweep__ThicknessAlNSub':
                             {'abbr':'Lsub_scan', 
                              'description': 'Substrate thickness variation simulations with Neumann end contact for AlN(25nm)/Al75Ga25N(300nm)/AlN(Lsub)'
                             },
                             'end_bc_Schottky_sweep__ThicknessAlNSub':
                             {'abbr':'end_scottky_Lsub_scan', 
                              'description': 'Substrate thickness variation simulations with Schottky end contact for AlN(25nm)/Al75Ga25N(300nm)/AlN(Lsub)'
                             },
                             'sim_sweep__ThicknessAlGaNBarrier':
                             {'abbr':'Lb_scan', 
                              'description': 'AlGaN barrier thickness variation simulations for AlN(Lb)/Al75Ga25N(300nm)/AlN(300nm)'
                             },
                             'sim_sweep__ThicknessAlGaNChannel':
                             {'abbr':'Lc_scan', 
                              'description': 'AlGaN channel thickness variation simulations for AlN(25nm)/Al75Ga25N(Lc)/AlN(300nm)'
                             },
                             'sim_sweep__AlContentBarrier__AlContentChannel__ThicknessAlGaNBarrier':
                             {'abbr':'x_y_Lb_scan', 
                              'description': 'AlGaN composition and barrier thickness variation simulations for AlyGa1-yN(Lb)/AlxGa1-xN(300nm)/AlN(300nm)'
                             }
                            }
# Creating helper.txt of this mapping in the DATAs folder
with open(os.path.join(folder_post_data_,'helper.txt'),'w') as helper_file:
    header_txt = f'''Simulation Software: Nextnano++(v1.21.24 RHEL compilation)
Simulation setup: 1D Schrodinger-Poisson, Cation-face growth direction
Sample: AlGaN/AlGaN/AlN HEMT
\n
'''
    
    sheet_column_name = f''' 
### {'='*72}
# Sheet column name: description
### {'='*72}
0th column    : Row indices
1st column(s) : The parameter(s) varied in simulations
2DEG_device   : Intergrated 2D electron gas density over whole device (in carriers/cm^2)
2DEG_BC       : Intergrated 2D electron gas density around the barrier-channel interface (within quantum region) (in carriers/cm^2)
2DEG_SC       : Intergrated 2D electron gas density around the Substrate/Buffer-channel interface (within quantum region) (in carriers/cm^2)
2DHG_device   : Intergrated 2D hole gas density over whole device (in carriers/cm^2)
2DHG_BC       : Intergrated 2D hole gas density around the barrier-channel interface (within quantum region) (in carriers/cm^2)
2DHG_SC       : Intergrated 2D hole gas density around the Substrate/Buffer-channel interface (within quantum region) (in carriers/cm^2)
\n### {'='*72}
'''
    helper_file.write(header_txt)
    helper_file.write(f"### {'='*72}\n")
    helper_file.write('# Sheet name: description\n')
    helper_file.write(f"### {'='*72}")
    
    for sim_name_full in MappingShortDataSheetName:
        save_text = f"""
{MappingShortDataSheetName[sim_name_full]['abbr']:<31} : {MappingShortDataSheetName[sim_name_full]['description']}
"""
        helper_file.write(save_text)
    helper_file.write(f"\n")
    helper_file.write(sheet_column_name)
    


# In[51]:


##=============================================================================
## Functions to screen conditional sweep variables set
# Note: Nextnanopy still did not include this functionality. I locally modified my installation of
# nextnanopy to accept this function. 
# Refer to: https://github.com/nextnanopy/nextnanopy/issues/12
def condition_screen_fn(combination):
    '''
    Returns variables combination where AlContentBarrier >  AlContentChannel.
    '''
    SweepVariablesKeys = list(SweepVariables.keys())
    i_index = SweepVariablesKeys.index("AlContentBarrier")
    j_index = SweepVariablesKeys.index("AlContentChannel")
    if combination[i_index] > combination[j_index]:
        return True
    return False

def create_tmp_input_file_4_sweep(ScanVariableName, base_input_path, mapps_, FileExtension='.in'):
    # read old input file
    input_file = nn.InputFile(base_input_path)
    old_input_path = base_input_path
    # turn on composition gradient keyword
    key_listt = mapps_[ScanVariableName]['SwitchKey']
    val_listt = mapps_[ScanVariableName]['SwitchKeyVal']
    assert len(key_listt) ==  len(val_listt), 'unmatched length of SwitchKey and SwitchKeyVal'
    for ii in range(len(key_listt)):
        input_file.set_variable(key_listt[ii], value=val_listt[ii])
    input_path = os.path.join(folder_input_, 'tmp', f'{mapps_[ScanVariableName]["input_file_name"]}{FileExtension}')
    input_file.save(input_path, overwrite=True)
    return old_input_path, input_path


# ## 2. Perform simulations

# In[52]:


for input_path in input_files_dest:
    print(f"{'*'*72}")
    #==============================================================================
    # run nextnano simulation for single set of variables
    if run_sim or run_sim_specific_sample:
        print('* Running single parameter set simulation...')
        print(f'- Input file: {input_path}', flush=True)
        input_file = nn.InputFile(input_path)
        # print(f"List of variables: {input_file.variables}")
        # for var in input_file.variables.values():
        #     print(f'{var.text}') 
        input_file.execute(outputdirectory=folder_output_, convergenceCheckbool=True)
    
    #==============================================================================    
    # run nextnano sweep simulation
    if run_sweep or run_sweep_specific_sample:
        print('* Running sweep simulation...')
        print(f'- Input file for Sweeping: {input_path}', flush=True)
        
        #***************************************************************************
        # loop over sweep sets defined above
        for ScanVariableName, SweepVariables in SweepVariablesSet.items():
            if 'specific_sample_' in input_path and ScanVariableName not in SpecificSamplesCase:
                continue
            if ScanVariableName in TemporaryInputFiles4:
                old_input_path, input_path = create_tmp_input_file_4_sweep(ScanVariableName, input_path, 
                                                                           TemporaryInputFiles4, 
                                                                           FileExtension=FileExtension)
            print(f"{'='*52}")
            print(f'- Sweeping set: {ScanVariableName}...', flush=True)
            print(f'- Sweeping variables: {list(SweepVariables.keys())}')
            sweep = nn.Sweep(SweepVariables, input_path)
            if ScanVariableName == 'CompositionBThicknessScan':
                sweep.save_sweep(delete_old_files = True, variables_comb_screen_fns=condition_screen_fn)
            else:
                sweep.save_sweep(delete_old_files = True)
            sweep.execute_sweep(delete_input_files=True, outputdirectory=folder_output_,
                                overwrite=True, convergenceCheckbool=True, show_log=False)   

            if ScanVariableName in TemporaryInputFiles4:
                input_path = old_input_path
    #==============================================================================
    # remove temporary folder
    shutil.rmtree(os.path.join(folder_input_, 'tmp'), ignore_errors=True) 
    #==============================================================================


# ## 3. Create post-processed data sheet from sweep simulations

# In[53]:


what_to_plots = ['2DEG', '2DHG']
plot_data_files = ['integrated_density_electron.dat', 'integrated_density_hole.dat']


# In[ ]:


if create_data_sweep:
    for input_path in input_files_dest:
        print(f"{'*'*72}", flush=True)
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '') 
        
        #==============================================================================
        # Set folder for post-processed data to save 
        post_processed_data_save = folder_output_.replace('OUTPUTs','DATAs')
        mkdir_if_not_exist(post_processed_data_save) 
        print(f'- Post-processed Data directory: {post_processed_data_save}')
        if append_new_sweep_data_2_data_base:
            fopen_mode = 'a'  
            if_sheet_exists='replace'
        else:
            fopen_mode = 'w'
            if_sheet_exists = None

        with pd.ExcelWriter(f'{folder_post_data_}/{input_file_name_variable}_post_process_data.xlsx', 
                            if_sheet_exists=if_sheet_exists, mode=fopen_mode,
                            engine="openpyxl") as writer:
            #==============================================================================
            # loop over sweep sets defined above
            for ScanVariableName, SweepVariables in SweepVariablesSet.items():
                if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                    continue
                
                #***************************************************************************
                # Set the output folder path
                SweepVariablesKeys = list(SweepVariables.keys())
                SweepVariablesFolder ='__'.join(SweepVariablesKeys)
                
                if ScanVariableName in TemporaryInputFiles4:
                    tmp_txt = TemporaryInputFiles4[ScanVariableName]['input_file_name']
                    fname_variable_text = f'{tmp_txt}_sweep__{SweepVariablesFolder}'
                else:
                    fname_variable_text =  f'{input_file_name_variable}_sweep__{SweepVariablesFolder}'
    
                sweep_folder_path = os.path.join(folder_output_, fname_variable_text)
                print(f'- Output data folder: {sweep_folder_path}')
                print(f'- Sweep variables: {SweepVariablesKeys}')
                
                #***************************************************************************
                # Read the sweep information from the output folder
                data_folder_ = nn.DataFolder(sweep_folder_path)
                sweep_info_file = data_folder_.file('sweep_infodict.json')
                
                with open(sweep_info_file) as json_file:
                    sweep_info_data = json.load(json_file)    
                #print(f'Sweep information: {sweep_info_data}')
                
                #=============================================================================================
                df = pd.DataFrame()
                for sweep_var_file in sweep_info_data:
                    df2_dict = sweep_info_data[sweep_var_file]
                    sweep_var_value = [sweep_info_data[sweep_var_file][sweep_variable_] for sweep_variable_ in SweepVariables]
                    #s_data_folder = nn.DataFolder(sweep_var_file)
                    s_data_folder = nn.DataFolder(sweep_var_file)
                    for ii, plot_file in enumerate(plot_data_files):
                        get_file_loc = nn.DataFile(s_data_folder.file(plot_file), product=software_)
                        # integrated 2DEHG for the whole device
                        df2_dict[f'{what_to_plots[ii]}_device'] = get_file_loc.variables['Whole_Device'].value[0]
                        try:
                            # integrated 2DEHG only in the barrier-channel interface
                            df2_dict[f'{what_to_plots[ii]}_BC'] = get_file_loc.variables['Qregion_2EG'].value[0]
                            # integrated 2DEHG only in the substrate-channel interface
                            df2_dict[f'{what_to_plots[ii]}_SC'] = get_file_loc.variables['Qregion_2HG'].value[0]
                        except:
                            pass
                        
                    df = pd.concat([df, pd.DataFrame([df2_dict])], ignore_index=True)

                print(f"- Saving data sheet: {fname_variable_text} => {MappingShortDataSheetName[fname_variable_text]['abbr']}")
                df.to_excel(writer, sheet_name=MappingShortDataSheetName[fname_variable_text]['abbr'])
                


# ## 4. Plottings

# ### 4.1 Plot band diagram from single simulation results (** Require original simulation results)

# In[54]:


# Define the region of band digram you want to zoom in
# [[xmin, xmax], [which_band(s)], [shift_yr, y_left_major_locator_distance]]
# zoom_band_diagram_regions = [[['QRegion_Left_2DEG', 'QRegion_Right_2DEG'], ['Gamma_', 'electron_Fermi_level_'], [-1, 1]],
#                              [['QRegion_Left_2DHG', 'QRegion_Right_2DHG'], ['HH_', 'LH_', 'SO_', 'electron_Fermi_level_'], [-1, 0.2]]]
zoom_band_diagram_regions = [[['EndAlGaNBarrier', 10], ['Gamma_', 'electron_Fermi_level_'], [-1, 0.4]],
                             [['EndAlGaNChannel', 10], ['HH_', 'LH_', 'SO_', 'electron_Fermi_level_'], [-1, 0.4]]]


# In[182]:


if do_plot:
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')
        data_folder_ = nn.DataFolder(f'{folder_output_}/{input_file_name_variable}')
        output_figs = data_folder_.fullpath.replace('OUTPUTs','FIGs')
        if upgrade_figs_folder: output_figs = output_figs.replace(replace_figs_path[0], replace_figs_path[1])
        print('- Output data folder:', data_folder_.fullpath)
        print('- Figure folder:', output_figs)
        #==========================================================================
        # Plot band edges + e-h desity
        lplt1d.PlotBandDiagrams(data_folder_, figs_path=output_figs, software_=software_, 
                                filename='band_edges', 
                                savefigure=savefigure, FigDpi=FigDpi, FigFormat=FigFormat4Paper,
                                density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                                plot_density=True, show_twin_yaxis_labels=False)
        #==========================================================================
        # Plot band edges + device structure
        lplt1d.PlotBandDiagrams(data_folder_, show_doping=False, show_Qregion=True, 
                                figs_path=output_figs, FigDpi=FigDpi, filename='band_edges_device_qr',
                                savefigure=savefigure, software_=software_, FigFormat=FigFormat, 
                                density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                                plot_density=True, plot_device_sketch=True, show_twin_yaxis_labels=True)
        #==========================================================================
        # Plot band edges zoomed in the 2DEG and 2DHG regions
        i = 0
        
        for zoom_region in zoom_band_diagram_regions:
            density_list=[('Electron_density', 'r'), ('Hole_density', 'b')]
            if i == 0: density_list=density_list[-1::-1]
            fig, ax, ax0, ax2 =\
            lplt1d.PlotBandDiagrams(data_folder_, figs_path=output_figs, software_=software_, 
                                    filename=f'band_edges_zoom_r{i}', savefigure=savefigure,
                                    FigDpi=FigDpi, FigFormat=FigFormat,xaxis_n_locator=None,
                                    x_zoom_region=zoom_region[0], plot_bands=zoom_region[1],
                                    x_zoom_2nd_no_shift=True, right_yaxis_shift=zoom_region[2][0],
                                    y_left_major_locator_distance=zoom_region[2][1],
                                    density_list=density_list, plot_density=True, 
                                    show_twin_yaxis_labels=True, align_left_right_yaxis=True)
            i+=1


# In[184]:


if do_plot:
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')
        data_folder_ = nn.DataFolder(f'{folder_output_}/{input_file_name_variable}')
        output_figs = data_folder_.fullpath.replace('OUTPUTs','FIGs')
        if upgrade_figs_folder: output_figs = output_figs.replace(replace_figs_path[0], replace_figs_path[1])
        print('- Output data folder:', data_folder_.fullpath)
        print('- Figure folder:', output_figs)
        # Plot band edges zoomed in the 2DEG and 2DHG regions
        kindex = 0
        plot_band_indices = {1:'tab:red',2:'tab:blue',3:'cyan'} #band_index, color 
        band_plot = {'Gamma_': 'quantum_2DEG_Gamma', 
                     'HH_': 'quantum_2DHG_kp6'} #'quantum_2DHG_HH'}
        for zoom_region in zoom_band_diagram_regions:
            i = 0
            for band_index, color_ in plot_band_indices.items():
                density_list=[('PsiSqare', color_)]
                if i == 0: 
                    fig=None; ax=None; ax0=None; ax2=None
                fig, ax, ax0, ax2 = \
                    lplt1d.PlotBandDiagrams(data_folder_, fig=fig, ax=ax, ax0=ax0, ax2=ax2, 
                                            software_=software_, savefigure=False, xaxis_n_locator=None,
                                            x_zoom_region=zoom_region[0], plot_bands=zoom_region[1],
                                            x_zoom_2nd_no_shift=True, right_yaxis_shift=None,
                                            y_left_major_locator_distance=zoom_region[2][1],
                                            band_file=band_plot[zoom_region[1][0]], 
                                            band_index=band_index, kindex=kindex,
                                            subband_energy_level=True, plot_density_on_left_axis=True,
                                            density_list=density_list, plot_density=True, ylabel_twin=None,
                                            show_twin_yaxis_labels=True)
                i += 1
                    
            lpltgen.save_figs(fig, filename=f'Psi_sqr_{band_plot[zoom_region[1][0]]}', 
                              figs_path=output_figs, savefigure=savefigure,
                              FigFormat=FigFormat, FigDpi=FigDpi)


# In[ ]:


if do_plot:
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')
        data_folder_ = nn.DataFolder(f'{folder_output_}/{input_file_name_variable}')
        output_figs = data_folder_.fullpath.replace('OUTPUTs','FIGs')
        if upgrade_figs_folder: output_figs = output_figs.replace(replace_figs_path[0], replace_figs_path[1])
        print('- Output data folder:', data_folder_.fullpath)
        print('- Figure folder:', output_figs)
        #==========================================================================
        # Plot band edges + e-h desity
        lplt1d.PlotBandDiagrams(data_folder_, figs_path=output_figs, software_=software_, 
                                filename='band_edges_eh', savefigure=savefigure,
                                FigDpi=FigDpi, FigFormat=FigFormat,
                                density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                                plot_density=True, show_twin_yaxis_labels=True)
        #==========================================================================
        # Plot band edges + polarization charge
        lplt1d.PlotBandDiagrams(data_folder_, figs_path=output_figs, software_=software_,
                                filename='band_edges_charges', savefigure=savefigure, 
                                density_list=[('Polarization_density', 'r')],
                                right_yaxis_shift=None,
                                plot_density=True, FigDpi=FigDpi, FigFormat=FigFormat)
        #==========================================================================
        # Plot band edges + device structure
        lplt1d.PlotBandDiagrams(data_folder_, show_doping=False, show_Qregion=False,
                                figs_path=output_figs, FigDpi=FigDpi, filename='band_edges_device',
                                savefigure=savefigure, software_=software_, FigFormat=FigFormat4Paper,
                                density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                                plot_density=True, plot_device_sketch=True, show_twin_yaxis_labels=True)


# ### 4.2 Plot band diagram from sweep results (** Require original simulation results)

# In[ ]:


if do_plot_sweep:
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')
        
        for ScanVariableName, SweepVariables in SweepVariablesSet.items():
            if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                continue
            
            #=============================================================================================
            SweepVariablesKeys = list(SweepVariables.keys())  
            SweepVariablesFolder ='__'.join(SweepVariablesKeys)

            if ScanVariableName in TemporaryInputFiles4:
                tmp_txt = TemporaryInputFiles4[ScanVariableName]['input_file_name']
                sweep_folder_path = os.path.join(folder_output_, f'{tmp_txt}_sweep__{SweepVariablesFolder}')
            else:
                sweep_folder_path = os.path.join(folder_output_, f'{input_file_name_variable}_sweep__{SweepVariablesFolder}')
            
            print(f'- Output data folder: {sweep_folder_path}')
            print(f'- Sweep variables: {SweepVariablesKeys}')
            
            #=============================================================================================
            data_folder_ = nn.DataFolder(sweep_folder_path)
            sweep_info_file = data_folder_.file('sweep_infodict.json')
            
            with open(sweep_info_file) as json_file:
                sweep_info_data = json.load(json_file)   
                for sweep_folder_path_ in sweep_info_data:
                    data_folder_sweep = nn.DataFolder(sweep_folder_path_)
                    output_figs_sweep = data_folder_sweep.fullpath.replace('OUTPUTs','FIGs')
                    if upgrade_figs_folder: output_figs_sweep = output_figs_sweep.replace(replace_figs_path[0], replace_figs_path[1])
                    print(f'-- Plotting: {data_folder_sweep.fullpath}')
                    #==========================================================================
                    # Plot band edges + carrier density + device structure
                    lplt1d.PlotBandDiagrams(data_folder_sweep, show_doping=False, show_Qregion=False,
                                            figs_path=output_figs_sweep, FigDpi=FigDpi, filename='band_edges_device',
                                            savefigure=savefigure, software_=software_, FigFormat=FigFormat, 
                                            density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                                            plot_density=True, plot_device_sketch=True, 
                                            show_twin_yaxis_labels=True)
                    
                    if SweepVariablesKeys[0] in ['ThicknessAlNSub']:
                        # Plot band edges + potential + device structure
                        fig, ax, ax0, ax2 = \
                        lplt1d.PlotBandDiagrams(data_folder_sweep, show_doping=False, show_Qregion=False,
                                                software_=software_, FigFormat=FigFormat, line_alpha=0.4,
                                                density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                                                plot_density=True, plot_device_sketch=True, show_twin_yaxis_labels=True)
                        # Create a Rectangle patch
                        input_variables_file_ = data_folder_sweep.file('variables_input.txt')
                        df_input_variables = nn.DataFile(input_variables_file_, product=software_)
                        end_device_pos = df_input_variables['EndDevice'].value
                        rect_ = patches.Rectangle((0.97,-0.05), 0.05, 1.1, 
                                  linewidth=1, edgecolor='r', linestyle='--',fill=False, transform=ax0.transAxes, clip_on=False)
                        # Add the patch to the Axes
                        ax0.add_patch(rect_)
                        lplt1d.PlotBandDiagrams(data_folder_sweep, fig=fig, ax=ax, 
                                                show_doping=False, show_Qregion=False,
                                                figs_path=output_figs_sweep, FigDpi=FigDpi, filename='potential',
                                                savefigure=savefigure, software_=software_, FigFormat=FigFormat4Paper, 
                                                density_list=[('Potential', 'c')], ylabel_twin=None,
                                                plot_density=True, plot_density_on_left_axis=True,
                                                plot_device_sketch=False, right_yaxis_shift=None,
                                                show_twin_yaxis_labels=True)


# #### 4.2.1 Set folder path for Sweep simulations run (** Require only post-processed data file)

# ##### 4.2.1.1 Set mapping of x-axis labels and x-ticks locator for different sweep plots

# In[55]:


rescale_2deg_fact = 1e13  # Rescalings 2DEG in 10^13 unit


# In[56]:


## Map of some variables Sweep variables for plotting
# {SweepVariablesKeys: {x_label_text: ..., ticks_multiplicator: 
# [xticks major multiplicator, xticks minor multiplicator, 
# yticks major multiplicator,yticks minor multiplicator]}}
# ticks_multiplicator_plot1: only for 2DEG plot
# ticks_multiplicator_plot2: Plot both 2DEG and 2DHG in same plot
mappp_ = {'SchottkyBarrierHeight':{'x_label_text': 'Schottky barrier height (eV)', 
                                   'ticks_multiplicator_plot1': [1.0, 0.5, 0.1, 0.05],
                                   'ticks_multiplicator_plot2': [1.0, 0.5, 0.1, 0.05]},
          'ThicknessAlNSub':{'x_label_text': 'Buffer thickness (nm)', 
                                   'ticks_multiplicator_plot1': [300, 100, None, None],
                                   'ticks_multiplicator_plot2': [300, 100, None, None]},
          'ThicknessAlGaNChannel':{'x_label_text': r'Channel thickness, L$_\mathrm{C}$ (nm)', 
                                   'ticks_multiplicator_plot1': [100, 50, 0.1, 0.05],
                                   'ticks_multiplicator_plot2': [100, 50, 0.2, 0.1]},
          'ThicknessAlGaNBarrier':{'x_label_text': r'Barrier thickness, L$_\mathrm{B}$ (nm)', 
                                   'ticks_multiplicator_plot1': [50, 25, None, None],
                                   'ticks_multiplicator_plot2': [50, 25, None, None]},
          'Temperature':{'x_label_text': 'Temperature (K)', 
                                   'ticks_multiplicator_plot1': [100, 50, None, None],
                                   'ticks_multiplicator_plot2': [100, 50, None, None]},
          'AlGaNbandgapBowing':{'x_label_text': 'Bandgap bowing (eV)', 
                                   'ticks_multiplicator_plot1': [0.2, 0.1, None, None],
                                   'ticks_multiplicator_plot2': [0.2, 0.1, 0.1, 0.05]},
         'AlGaNpyroelectricBowing':{'x_label_text': 'Pyroelectric bowing', 
                                   'ticks_multiplicator_plot1': [0.1, 0.1, 0.2, 0.1],
                                   'ticks_multiplicator_plot2': [0.1, 0.1, 0.5, 0.25]}
         }


# ##### 4.2.1.2 Plot 1D sweep variables vs property (e.g. 2DEG density)

# In[ ]:


if do_plot_sweep:
    y_label_text1 = '2DEG density, n$_\\mathrm{2D}$ ($10^{13}$ $\\mathrm{cm}^{-2}$)' #2DEG density
    y_label_text1_H = '2DHG density ($10^{13}$ $\\mathrm{cm}^{-2}$)'
    y_label_text2 = '2DCG density ($10^{13}$ $\\mathrm{cm}^{-2}$)'
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')  
        excel_file = f'{folder_post_data_}/{input_file_name_variable}_post_process_data.xlsx'        
        #==============================================================================
        # loop over sweep sets defined above
        for ScanVariableName, SweepVariables in SweepVariablesSet.items():
            if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                continue
            if ScanVariableName in ['CompositionBThicknessScan', 'GaNChannelStudy']: continue
            #***************************************************************************
            # Set the output folder path
            SweepVariablesKeys = list(SweepVariables.keys())
            SweepVariablesFolder ='__'.join(SweepVariablesKeys) 

            if ScanVariableName in TemporaryInputFiles4:
                tmp_txt = TemporaryInputFiles4[ScanVariableName]['input_file_name']
                data_sheet_name_map = f'{tmp_txt}_sweep__{SweepVariablesFolder}'
            else:
                data_sheet_name_map =  f'{input_file_name_variable}_sweep__{SweepVariablesFolder}'
            data_sheet_name = MappingShortDataSheetName[data_sheet_name_map]['abbr']  
            print(f'- Reading data sheet: {data_sheet_name}')
            #=============================================================================================
            # Set figures directory
            output_figs_sweep = os.path.join(folder_output_,data_sheet_name_map).replace('OUTPUTs','FIGs')
            if upgrade_figs_folder: output_figs_sweep = output_figs_sweep.replace(replace_figs_path[0], replace_figs_path[1])
            mkdir_if_not_exist(output_figs_sweep) 
            print(f'- Figs directory: {output_figs_sweep}')
            # #=============================================================================================
            df = pd.read_excel(excel_file, sheet_name=data_sheet_name, index_col=0)
            #=============================================================================================
            x_label_text = mappp_.get(SweepVariablesKeys[0])['x_label_text']
            ticks_multiplicator_plot1 = mappp_.get(SweepVariablesKeys[0])['ticks_multiplicator_plot1']
            ticks_multiplicator_plot2 = mappp_.get(SweepVariablesKeys[0])['ticks_multiplicator_plot2']
            #=============================================================================================
            # device: whole region
            which_regions_to_plot = ['device'] 
            x_log_scale = False
            fname_save_fig_extra = ''
            #=============================================================================================
            XX = df[SweepVariablesKeys[0]]
            for JJJ in which_regions_to_plot:
                YY_2DEG =  df[f'2DEG_{JJJ}'].copy()/rescale_2deg_fact
                YY_2DHG =  df[f'2DHG_{JJJ}'].copy()/rescale_2deg_fact
                #=============================================================================================
                # Plot only 2DEG
                fig, ax = lplt1d.PlotSweepsData(XX, YY_2DEG, x_label=x_label_text, y_label=y_label_text1,
                                                tick_multiplicator=ticks_multiplicator_plot1,
                                                FigDpi=FigDpi, FigFormat=FigFormat4Paper,
                                                figs_path=output_figs_sweep, filename=f'2DEG_{JJJ}', 
                                                savefigure=0, x_log_scale=x_log_scale)
                if SweepVariablesKeys[0] in ['ThicknessAlGaNBarrier', 'ThicknessAlGaNChannel']: 
                    ax.axhline(y=1, c='k', ls='--')
                lplt1d.save_figs(fig, filename=f'2DEG_{JJJ}', figs_path=output_figs_sweep, 
                                 savefigure=savefigure, FigDpi=FigDpi, FigFormat=FigFormat4Paper)
                #=============================================================================================
                if SweepVariablesKeys[0] in ['ThicknessAlNSub']:
                    # Plot only 2DHG
                    lplt1d.PlotSweepsData(XX, YY_2DHG, x_label=x_label_text, y_label=y_label_text1_H,
                                          tick_multiplicator=ticks_multiplicator_plot1,
                                          FigDpi=FigDpi, FigFormat=FigFormat4Paper,color='b',
                                          figs_path=output_figs_sweep, filename=f'2DHG_{JJJ}', 
                                          savefigure=savefigure, x_log_scale=x_log_scale)
                #=============================================================================================
                # Plot 2DEG and 2DHG
                fig, ax = lplt1d.PlotSweepsData(XX, YY_2DEG)
                fig, ax = lplt1d.PlotSweepsData(XX, YY_2DHG,
                                                fig=fig, ax=ax, x_label=x_label_text, y_label=y_label_text2,
                                                tick_multiplicator=ticks_multiplicator_plot2,
                                                FigDpi=FigDpi, FigFormat=FigFormat4Paper, color='b',
                                                figs_path=output_figs_sweep, filename=f'2DEHG_{JJJ}',
                                                savefigure=True, x_log_scale=x_log_scale)
                #=============================================================================================


# In[ ]:


if do_plot_sweep:
    deg_label_text = r'2DEG density, n$_\mathrm{2D}$' #'2DEG density' 
    y_label_text1 = f'{deg_label_text} ($10^{{13}}$ $\\mathrm{{cm}}^{{-2}}$)'
    y_label_text1_H = f'2DHG density ($10^{{13}}$ $\\mathrm{{cm}}^{{-2}}$)'
    y_label_text2 = f'{deg_label_text} ($10^{{13}}$ $\\mathrm{{cm}}^{{-2}}$)'
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')  
        excel_file = f'{folder_post_data_}/{input_file_name_variable}_post_process_data.xlsx'        
        #==============================================================================
        # loop over sweep sets defined above
        ii = 0; tmp_savefig=False
        for ScanVariableName in ['SchottkyBarrierEndDevice', 'NeumannBCEndDevice']:
            if ScanVariableName not in SweepVariablesSet: continue 
            if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                continue
            SweepVariables = SweepVariablesSet[ScanVariableName]
            #***************************************************************************
            # Set the output folder path
            SweepVariablesKeys = list(SweepVariables.keys())
            SweepVariablesFolder ='__'.join(SweepVariablesKeys) 

            if ScanVariableName in TemporaryInputFiles4:
                tmp_txt = TemporaryInputFiles4[ScanVariableName]['input_file_name']
                data_sheet_name_map = f'{tmp_txt}_sweep__{SweepVariablesFolder}'
            else:
                data_sheet_name_map =  f'{input_file_name_variable}_sweep__{SweepVariablesFolder}'
            data_sheet_name = MappingShortDataSheetName[data_sheet_name_map]['abbr']   
            print(f'- Reading data sheet: {data_sheet_name}')
            #=============================================================================================
            # Set figures directory
            output_figs_sweep = os.path.join(folder_output_,data_sheet_name_map).replace('OUTPUTs','FIGs')
            if upgrade_figs_folder: output_figs_sweep = output_figs_sweep.replace(replace_figs_path[0], replace_figs_path[1])
            mkdir_if_not_exist(output_figs_sweep) 
            print(f'- Figs directory: {output_figs_sweep}')
            # #=============================================================================================
            df = pd.read_excel(excel_file, sheet_name=data_sheet_name, index_col=0)
            #=============================================================================================
            x_label_text = mappp_.get(SweepVariablesKeys[0])['x_label_text']
            ticks_multiplicator_plot1 = mappp_.get(SweepVariablesKeys[0])['ticks_multiplicator_plot1']
            ticks_multiplicator_plot2 = mappp_.get(SweepVariablesKeys[0])['ticks_multiplicator_plot2']
            #=============================================================================================
            which_regions_to_plot = ['device'] 
            x_log_scale = False
            fname_save_fig_extra = ''
            #=============================================================================================
            XX = df[SweepVariablesKeys[0]]
            for JJJ in which_regions_to_plot:
                YY_2DHG =  df[f'2DHG_{JJJ}'].copy()/rescale_2deg_fact
                #=============================================================================================
                # Plot 2DEG and 2DHG
                if ii in [0]: 
                    fig, ax = lplt1d.PlotSweepsData(XX, YY_2DHG)
                    line_style = '-'
                else:
                    tmp_savefig = True
                    line_style = '--'
                fig, ax = lplt1d.PlotSweepsData(XX, YY_2DHG, fig=fig, ax=ax,
                                                x_label=x_label_text, y_label=y_label_text1_H,
                                                tick_multiplicator=ticks_multiplicator_plot1,
                                                FigDpi=FigDpi, FigFormat=FigFormat4Paper,color='b',
                                                line_style=line_style,
                                                figs_path=output_figs_sweep, filename=f'2DHG_compare_{JJJ}',
                                                savefigure=tmp_savefig, x_log_scale=x_log_scale)
                #=============================================================================================
            ii += 1


# ##### 4.2.1.3 Plot 3D sweep variables vs property (e.g. 2DEG density)

# In[ ]:


if do_plot_sweep:
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')  
        excel_file = f'{folder_post_data_}/{input_file_name_variable}_post_process_data.xlsx'
        #==============================================================================
        # loop over sweep sets defined above
        ScanVariableNamesList = ['CompositionBThicknessScan']
        for ScanVariableName in ScanVariableNamesList:
            if ScanVariableName not in SweepVariablesSet: continue 
            SweepVariables = SweepVariablesSet[ScanVariableName]
            if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                continue
            #***************************************************************************
            # Set the output folder path
            SweepVariablesKeys = list(SweepVariables.keys())
            SweepVariablesFolder ='__'.join(SweepVariablesKeys)
            data_sheet_name_map =  f'{input_file_name_variable}_sweep__{SweepVariablesFolder}'
            data_sheet_name = MappingShortDataSheetName[data_sheet_name_map]['abbr']
            print(f'- Reading data sheet: {data_sheet_name}')
            #=============================================================================================
            # Set figures directory
            output_figs_sweep = os.path.join(folder_output_,data_sheet_name_map).replace('OUTPUTs','FIGs')
            if upgrade_figs_folder: output_figs_sweep = output_figs_sweep.replace(replace_figs_path[0], replace_figs_path[1])
            mkdir_if_not_exist(output_figs_sweep) 
            print(f'- Figs directory: {output_figs_sweep}')
        
            #=============================================================================================
            df = pd.read_excel(excel_file, sheet_name=data_sheet_name, index_col=0)
            df['2DEG_device_rescale'] =  df['2DEG_device'].copy()/rescale_2deg_fact
            df['2DHG_device_rescale'] =  df['2DHG_device'].copy()/rescale_2deg_fact
            
            df_cap_thickness_group = df.groupby(['ThicknessAlGaNBarrier']) # Group by thickness
            x_label_text = 'Al composition contrast, $\\Delta_{yx}$'
            z_label_text = r'Barrier thickness, L$_\mathrm{C}$ (nm)'
            deg_label_text = r'2DEG density, n$_\mathrm{2D}$' #'2DEG density'
            y_label_text = [f'{deg_label_text} ($10^{{13}}$ $\\mathrm{{cm}}^{{-2}}$)',
                            '2DHG density ($10^{13}$ $\\mathrm{cm}^{-2}$)']
            
            output_dataframe = {} 
            fig1, ax1 = plt.subplots(figsize=(8,6))
            fig2, ax2 = plt.subplots()
            fig = [fig1, fig2]
            ax = [ax1, ax2]
            for name, group in df_cap_thickness_group:
                output_dataframe[f'{name[0]:.2f}'] = {}
                for iii, what_to_plot in enumerate(what_to_plots):
                    ##### Generate data 
                    if '2DEG' in what_to_plot:
                        group['AlContentContrast'] = (group['AlContentBarrier'] - group['AlContentChannel']) #/ group['AlContentChannel']
                        x_label_text = 'Al composition contrast, y-x'
                    elif '2DHG' in what_to_plot:
                        group['AlContentContrast'] = (1 - group['AlContentChannel']) #/ group['AlContentSubstrate']
                        x_label_text = 'Al composition contrast, 1-x'
                    else:
                        raise AttributeError('Noting to plot')

                    #=====================================================================================
                    ax[iii].scatter(group['AlContentContrast'], group[f'{what_to_plot}_device_rescale'], label=f'{name[0]:.2f} nm')
                    
                    #=====================================================================================
                    tmp_df = group[group[f'{what_to_plot}_device_rescale']>0.01]
                    A = np.vstack([tmp_df['AlContentContrast'], np.ones(len(tmp_df))]).T
                    m, c = np.linalg.lstsq(A, tmp_df[f'{what_to_plot}_device_rescale'])[0]
                    y_eq_zero_point = -c/m
                    output_dataframe[f'{name[0]:.2f}'][f'{what_to_plot}_slope']       = f'{m:.3f}' 
                    output_dataframe[f'{name[0]:.2f}'][f'{what_to_plot}_y_intersect'] = f'{c:.3f}'
                    output_dataframe[f'{name[0]:.2f}'][f'{what_to_plot}_x_intersect'] = f'{y_eq_zero_point:.3f}'
                    new_x = np.insert(np.array(tmp_df['AlContentContrast']),0,y_eq_zero_point)
                    ax[iii].plot(new_x, m*new_x + c, 'gray')
                    
                    ax[iii].legend(ncols=1, labelspacing=0.001, columnspacing=0.1, handletextpad=0.01)
                    ax[iii].set_xlabel(x_label_text)
                    ax[iii].set_ylabel(y_label_text[iii], size=20)
                    ax[iii].axhline(color='k',ls='--')
            #=============================================================================================
            for iii, what_to_plot in enumerate(what_to_plots):
                lpltgen.save_figs(fig[iii], filename=f'Comp_contrast_thickness_{what_to_plot}', 
                                  figs_path=f'{output_figs_sweep}/Others', savefigure=savefigure, 
                                  FigDpi=FigDpi, FigFormat=FigFormat)          

            #=============================================================================================
            output_data = pd.DataFrame.from_dict(output_dataframe, orient='index', dtype=float)
            ## Plot cut-off composition contrast for each barrier thickness
            tick_multiplicator = [10,5,0.1, 0.05]
            XX = np.array(output_data.index, dtype=float)
            YY = np.array(output_data['2DEG_x_intersect'], dtype=float)
            _ = lplt1d.PlotSweepsData(XX, YY, x_label='Barrier thickness, L$_\\mathrm{B}$ (nm)', y_label='Al composition contrast, $\\Delta_{yx}$',
                                      tick_multiplicator=tick_multiplicator,
                                      FigDpi=FigDpi, FigFormat=FigFormat, color='k', marker='o',
                                      figs_path=f'{output_figs_sweep}/Others', filename=f'Critical_comp_contrast_2DEG',
                                      savefigure=savefigure)


# In[ ]:


if do_plot_sweep:
    for input_path in input_files_dest:
        print(f"{'*'*72}")
        input_filename = input_path.split('/')[-1]
        input_file_name_variable = input_filename.replace(FileExtension, '')  
        excel_file = f'{folder_post_data_}/{input_file_name_variable}_post_process_data.xlsx'
        #==============================================================================
        # loop over sweep sets defined above
        ScanVariableNamesList = ['CompositionBThicknessScan']
        for ScanVariableName in ScanVariableNamesList:
            if ScanVariableName not in SweepVariablesSet: continue
            SweepVariables = SweepVariablesSet[ScanVariableName]
            if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                continue
            #***************************************************************************
            # Set the output folder path
            SweepVariablesKeys = list(SweepVariables.keys())
            SweepVariablesFolder ='__'.join(SweepVariablesKeys)
            data_sheet_name_map =  f'{input_file_name_variable}_sweep__{SweepVariablesFolder}'
            data_sheet_name = MappingShortDataSheetName[data_sheet_name_map]['abbr']
            print(f'- Reading data sheet: {data_sheet_name}')
            #=============================================================================================
            # Set figures directory
            output_figs_sweep = os.path.join(folder_output_,data_sheet_name_map).replace('OUTPUTs','FIGs')
            if upgrade_figs_folder: output_figs_sweep = output_figs_sweep.replace(replace_figs_path[0], replace_figs_path[1])
            mkdir_if_not_exist(output_figs_sweep) 
            print(f'- Figs directory: {output_figs_sweep}')
        
            #=============================================================================================
            df = pd.read_excel(excel_file, sheet_name=data_sheet_name, index_col=0)
            df['2DEG_device_rescale'] =  df['2DEG_device'].copy()/rescale_2deg_fact
            df['2DHG_device_rescale'] =  df['2DHG_device'].copy()/rescale_2deg_fact
            
            df_cap_thickness_group = df.groupby(['ThicknessAlGaNBarrier']) # Group by thickness
            x_label_text = 'Al composition channel, x'
            y_label_text = 'Al composition barrier, y'
            z_label_unit = r'($10^{13}$ $\mathrm{cm}^{-2}$)'
            for iii, what_to_plot in enumerate(what_to_plots):
                vmin, vmax = df[f'{what_to_plot}_device_rescale'].min(), df[f'{what_to_plot}_device_rescale'].max()
                norm, cbar_mapable = lpltq3d.CreateColorbarMapableObject(vmin=vmin, vmax=vmax, color_map=color_map)
                for name, group in df_cap_thickness_group:
                    print(f'- Plotting barrier thickness - {what_to_plot} = {name[0]:.2f} nm')
                    ##### Generate data with interpolation
                    xi, yi, zi = lpltq3d.InterPolation(group['AlContentChannel'], 
                                                       group['AlContentBarrier'],
                                                       group[f'{what_to_plot}_device_rescale'],
                                                       method='linear', interpolation_points=20)
    
                    ##### Plot composition map vs 2DEG
                    output_figs_sweep_tmp = os.path.join(output_figs_sweep, what_to_plot)
                    color_bar_text = f'{what_to_plot} density {z_label_unit}'
                    if what_to_plot == '2DEG':
                        color_bar_text = f'{what_to_plot} density, n$_\\mathrm{{2D}}$ {z_label_unit}'
                    fig, ax = lpltq3d.PlotContour(xi, yi, zi, vmin=vmin, vmax=vmax, 
                                                  x_label=x_label_text, y_label=y_label_text,
                                                  tick_multiplicator=[0.1, 0.05, 0.1, 0.05],
                                                  title_label=None,#f"Barrier thickness = {name[0]:.2f} nm", 
                                                  cbar_mappable=cbar_mapable, norm=norm,
                                                  show_contour_lines=True, show_colorbar=True,
                                                  cbar_text=color_bar_text,
                                                  filename= f'Barrier_{name[0]:.2f}',
                                                  FigDpi=FigDpi, FigFormat=FigFormat, 
                                                  figs_path=output_figs_sweep_tmp,
                                                  savefigure=savefigure)


# ##### 4.2.1.4 Plot 2DEG distributions for selected sweep samples (* require original simulation)¶

# In[176]:


fname_lists = [['sim__AlContentBarrier_1.0_AlContentChannel_0.5_ThicknessAlGaNBarrier_50.0_',
              'sim__AlContentBarrier_1.0_AlContentChannel_0.75_ThicknessAlGaNBarrier_50.0_',
              'sim__AlContentBarrier_1.0_AlContentChannel_0.9_ThicknessAlGaNBarrier_50.0_'],
              ['sim__AlContentBarrier_1.0_AlContentChannel_0.5_ThicknessAlGaNBarrier_50.0_',
              'sim__AlContentBarrier_1.0_AlContentChannel_0.5_ThicknessAlGaNBarrier_25.0_',
              'sim__AlContentBarrier_1.0_AlContentChannel_0.5_ThicknessAlGaNBarrier_10.0_']]
scale_x_axis_ = [[0, 0, 0], [40, 15, 0]]
clsp = ['r', 'c', 'm', 'g', 'y']
zoom_band_diagram_regions = [[['EndAlGaNBarrier', 10], ['Gamma_', 'electron_Fermi_level_'], [-160, 0.4]],
                             [['EndAlGaNChannel', 10], ['HH_', 'LH_', 'SO_', 'electron_Fermi_level_'], [-160, 0.4]]]


# In[181]:


if do_plot_sweep:
    for input_path in input_files_dest:
            print(f"{'*'*72}")
            input_filename = input_path.split('/')[-1]
            input_file_name_variable = input_filename.replace(FileExtension, '')
            ScanVariableNamesList = ['CompositionBThicknessScan']
            for ScanVariableName in ScanVariableNamesList:
                if ScanVariableName not in SweepVariablesSet: continue
                SweepVariables = SweepVariablesSet[ScanVariableName]
                if 'specific_sample_' in input_filename and ScanVariableName not in SpecificSamplesCase:
                    continue
                #***************************************************************************
                # Set the output folder path
                SweepVariablesKeys = list(SweepVariables.keys())
                SweepVariablesFolder ='__'.join(SweepVariablesKeys)
                sweep_folder_path = os.path.join(folder_output_, f'{input_file_name_variable}_sweep__{SweepVariablesFolder}')
                
                print(f'- Output data folder: {sweep_folder_path}')
                print(f'- Sweep variables: {SweepVariablesKeys}')
                output_figs_sweep = f'{sweep_folder_path}/Others'.replace('OUTPUTs','FIGs')
                if upgrade_figs_folder: output_figs_sweep = output_figs_sweep.replace(replace_figs_path[0], replace_figs_path[1])
                mkdir_if_not_exist(output_figs_sweep) 
                print(f'- Figs directory: {output_figs_sweep}')

                for kk, fname_list in enumerate(fname_lists):
                    for zoom_region in zoom_band_diagram_regions:
                        i = 0
                        legend_txtx = []
                        for lll, sweep_folder_path_ in enumerate(fname_list):
                            bandedge_characterstics = {'Gamma_':('Gamma', clsp[i]),
                                        'HH_':('heavy hole', 'y'),
                                        'LH_':('light hole', 'tab:blue'),
                                        'SO_':('crystal-field hole', 'g'),
                                        'electron_Fermi_level_':('Fermi level', 'gray')
                                        }
                            data_folder_sweep = nn.DataFolder(os.path.join(sweep_folder_path, sweep_folder_path_))
                            print(f'-- Plotting: {data_folder_sweep.fullpath}')
                            density_list=[('Electron_density', clsp[i])]
                            if i == 0: 
                                #density_list=density_list[-1::-1]
                                fig=None; ax=None; ax0=None; ax2=None
                            fig, ax, ax0, ax2 =\
                            lplt1d.PlotBandDiagrams(data_folder_sweep, fig=fig, ax=ax, ax0=ax0, ax2=ax2,
                                                    figs_path=output_figs_sweep, software_=software_, 
                                                    savefigure=False,bands_characters=bandedge_characterstics,
                                                    band_edge_ls='--', scale_x_axis = scale_x_axis_[kk][lll],
                                                    FigDpi=FigDpi, FigFormat=FigFormat4Paper,xaxis_n_locator=None,
                                                    x_zoom_region=zoom_region[0], plot_bands=zoom_region[1],
                                                    x_zoom_2nd_no_shift=True, 
                                                    right_yaxis_shift= -160,#zoom_region[2][0],
                                                    y_left_major_locator_distance=zoom_region[2][1],
                                                    density_list=density_list, plot_density=True, 
                                                    show_twin_yaxis_labels=1, align_left_right_yaxis=False)
                            i+=1
                            ttt = sweep_folder_path_.split('_')
                            legend_txtx.append(f'AlN({ttt[-2]})/Al$_{{{float(ttt[5]):.2}}}$Ga$_{{{1.0-float(ttt[5]):.2}}}$N') 
                        ax.set_ylim(ymax=1, ymin=-0.8)
                        ax2.set_yticks([0,50,100,150,200])
                        ax2.set_yticklabels([0,50,100,150,200])
                        ax.axhline(y=0,c='gray')
                        #ax.legend(legend_txtx)
                        lpltgen.save_figs(fig, filename=f'{zoom_region[1][0]}{kk}', figs_path=output_figs_sweep, 
                                          savefigure=savefigure, FigFormat=FigFormat4Paper, FigDpi=FigDpi)


# In[ ]:




