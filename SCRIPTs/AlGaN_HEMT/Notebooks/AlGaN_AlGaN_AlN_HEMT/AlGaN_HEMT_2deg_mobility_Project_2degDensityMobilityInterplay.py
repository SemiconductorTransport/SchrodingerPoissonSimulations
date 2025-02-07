#!/usr/bin/env python
# coding: utf-8

# # Project: AlGaN/AlGaN/AlN HEMT 2DEG mobility calculation using mobilitypy
# 
# The mobility models are implemented using the following refs:
# 
# Ref-1: J. Bassaler, J. Mehta, I. Abid, L. Konczewicz, S. Juillaguet, S. Contreras, S. Rennesson, S. Tamariz, M. Nemoz, F. Semond, J. Pernot, F. Medjdoub, Y. Cordier, P. Ferrandis, Al-Rich AlGaN Channel High Electron Mobility Transistors on Silicon: A Relevant Approach for High Temperature Stability of Electron Mobility. Adv. Electron. Mater. 2024, 2400069. https://doi.org/10.1002/aelm.202400069
# 
# Ref-2: Zhang, J., Hao, Y., Zhang, J. et al. The mobility of two-dimensional electron gas in AlGaN/GaN heterostructures with varied Al content. Sci. China Ser. F-Inf. Sci. 51, 780–789 (2008). https://doi.org/10.1007/s11432-008-0056-7
# 
# Ref-3: This publication
# 

# # 1. Settings

# ## 1.1 Import modules

# In[ ]:


import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colors, cm
import matplotlib.ticker as ticker
from matplotlib.patches import Ellipse


# In[ ]:


#%load_ext autoreload
#%autoreload 2


# In[ ]:


import sys, os
from pathlib import Path
from mobilitypy import AlloyParams, Mobility2DCarrier
from mobilitypy import Plottings, PlotQuasi3DFuns


# ## 1.2 Set some parameters

# In[ ]:


redo_mobility_cal = True # re-do mobility calculation
resave_mobilities = True # re-save the calculated mobility values. Make sure you use redo_mobility_cal=True.


# In[ ]:


savefigure = 0
FigFormat = 'png'
FigFormatPaper = 'eps'
fig_save_dpi = 300
color_map_plot = 'viridis'
print_align_space = 20


# ## 1.3 Set folders path

# In[ ]:


current_script_path = Path().absolute()
# Simulation directory
sim_dir = str(current_script_path).split('SCRIPTs')[0]
# Project id
my_project_id = '2DEG_DensityMobilityInterplay'
# inner location of input file within INPUTs folder
folder_input_lower_level = f'AlGaN_HEMT/AlGaN_AlGaN_AlN_HEMT/nnp/{my_project_id}'
# complete path of data file folder
folder_2deg_data_parent = os.path.join(sim_dir, 'DATAs', folder_input_lower_level)
folder_2deg_density_file = os.path.join(folder_2deg_data_parent, '2DEHG_density/OneDim', 'sim_post_process_data.xlsx')
folder_2deg_mobility_folder = os.path.join(folder_2deg_data_parent, '2DEHG_mobility/OneDim')
folder_2deg_mobility_file = os.path.join(folder_2deg_mobility_folder,'sim_post_process_data.xlsx') 
save_figure_dir = os.path.join(sim_dir, 'FIGs', folder_input_lower_level, '2DEHG_mobility/OneDim')

Path(folder_2deg_mobility_folder).mkdir(parents=True, exist_ok=True)
Path(save_figure_dir).mkdir(parents=True, exist_ok=True)


# In[ ]:


print(f'* 2DEG density data file: {folder_2deg_density_file}')
print(f'* 2DEG mobility data file: {folder_2deg_mobility_file}')
print(f'* 2DEG mobility figure file: {folder_2deg_mobility_file}')


# ## 1.4 Set physical constants

# In[ ]:


'''
n_2d => in nm^-2
rms_roughness => nm^-1
corr_len => nm^-1
n_dis => nm^-2
f_dis => unit less
'''
alloy_disordered_effect=1
interface_roughness_effect=1
dislocation_effect=1
deformation_potential_effect=1 
piezoelectric_effect=1
acoustic_phonon_effect=1
polar_optical_phonon_effect=1
total_mobility=1
mobility_model='Bassaler'
density_2deg = 0.1 # nm^-2
irf_rms_roughness = 0.3 # nm
irf_corr_length = 3.0 # nm
dislocation_density = 1e-4 # nm^-2
occup_dislocation = 0.3
T_sweep = [10, 100, 200, 300, 400, 500, 600, 700, 800] # in K , sweep temperature
T_corect_bandgap_in_LFOM = False


# In[ ]:


plot_mobilities = ['IFR', 'DIS', 'AD', 'AP', 'POP', 'TOT', 'LFOM', 'LFOMnorm', '2DEG']
plots_list_individual = ['TOT', 'LFOM', 'LFOMnorm', '2DEG', '2DHG']
len_plots = len(plot_mobilities)
rescale_2deg_fact=1e13

x_label_text = 'Al composition barrier, y'
x_label_text_2d = 'Barrier thickness, L$_\\mathrm{B}$ (nm)'
x_p_label_text = 'Al composition channel, x'
y_label_text = 'Barrier thickness, L$_\\mathrm{B}$(nm)'
z_label = {'General': r'2DEG mobility, $\mu$ ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)', 
           'TOT': r'2DEG mobility, $\mu$ ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)',
           'LFOM': r'LFOM (MW/cm$^2$)',
           'LFOMnorm': r'LFOM$_\mathrm{(Al,Ga)N}$/LFOM$_\mathrm{highest}$', 
           'LFOMnorm_2': r'LFOM$_\mathrm{(Al,Ga)N}$/LFOM$_\mathrm{GaN}$',
           '2DEG': r'2DEG density, n$_{\mathrm{2D}}$ ($\mathrm{10}^{13}\mathrm{cm}^{-2}$)', 
           '2DHG': r'2DHG density ($\mathrm{10}^{13}\mathrm{cm}^{-2}$)'} #Sheet resistance


# In[ ]:


params = {'axes.labelsize': 20,
          'axes.titlesize': 20,
          'xtick.labelsize':20,
          'xtick.major.width':2,
          'xtick.major.size':5,
          'xtick.minor.width':2,
          'xtick.minor.size':3,
          'ytick.labelsize': 20,
          'ytick.major.width':2,
          'ytick.major.size':5,
          'ytick.minor.width':2,
          'ytick.minor.size':3,
          'errorbar.capsize':2}
plt.rcParams.update(params)
plt.rc('font', size=20)


# In[ ]:


rescale_2deg_fact=1e13


# # 2. Calculate mobilities

# ## 2.1 Al0.25Ga0.75N(25nm)/GaN(305nm) HEMT LFOM for normalization

# ### 2.1.1 Read 2DEG data from nextnano simulation

# In[ ]:


if redo_mobility_cal:
    data_sheet_name =  'GaN_T_scan' #'GaN_C_sweep__Temperature'
    n_2d_df_ref = pd.read_excel(folder_2deg_density_file, sheet_name=data_sheet_name, index_col=0)
    density_2deg_comp_ref = np.array(n_2d_df_ref['2DEG_device'], dtype=float)/1e+14 # cm^-2 to nm^-2 conversion


# In[ ]:


if redo_mobility_cal:
    fig, ax = plt.subplots(figsize=(9,6), constrained_layout=True)
    XX = np.array(n_2d_df_ref['Temperature'], dtype=float)
    ax.plot(XX, n_2d_df_ref['2DEG_device']/rescale_2deg_fact, 'o-', c='k')
    ax.set_xlabel('Al composition channel, x')
    ax.set_ylabel(r'2DEG density, n$_\mathrm{2D}$ ($\mathrm{10}^{13}\mathrm{cm}^{-2}$)')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))


# Conclusion: 2DEG density negligible depend on temperature

# ### 2.1.2 Calculate Al0.25Ga0.75N(25nm)/GaN(305nm) mobilities

# In[ ]:


if redo_mobility_cal:
    channel_comp_ref = 0 # Channel is pure GaN
    density_2deg_comp_ref = 0.0625 # in nm^-2
    mu2deg = Mobility2DCarrier(compositions=channel_comp_ref, binaries=['AlN', 'GaN'], alloy='AlGaN', system='ternary', 
                               print_log=None, eps_n_2d=1e-10)


# In[ ]:


if redo_mobility_cal:
    mobility_ref_dff = pd.DataFrame()
    
    for T in T_sweep:
        # Calculate mobilities
        mobility_ref = mu2deg.calculate_sheet_mobility(n_2d=density_2deg_comp_ref, rms_roughness=irf_rms_roughness,
                                                       corr_len=irf_corr_length, n_dis=dislocation_density,
                                                       f_dis=occup_dislocation, T=T,
                                                       alloy_disordered_effect=False,
                                                       interface_roughness_effect=interface_roughness_effect,
                                                       dislocation_effect=dislocation_effect,
                                                       deformation_potential_effect=deformation_potential_effect,
                                                       piezoelectric_effect=piezoelectric_effect,
                                                       acoustic_phonon_effect=acoustic_phonon_effect,
                                                       polar_optical_phonon_effect=polar_optical_phonon_effect,
                                                       total_mobility=total_mobility,
                                                       #calculate_total_mobility_only=True,
                                                       mobility_model=mobility_model)
        
        ref_mu_total = np.array(mobility_ref['TOT'], dtype=float) # cm^2.V^-1.s^-1
        # Add temperature column
        mobility_ref['T'] = T # K
        mobility_ref['2DEG_device'] = density_2deg_comp_ref*1e14 #cm^-2
        # Calculate LFOM
        mobility_ref['LFOM'] = mu2deg.calculate_figure_of_merit(density_2deg_comp_ref, ref_mu_total,
                                                                temp=T, T_corect_bandgap=T_corect_bandgap_in_LFOM) # LFOM: MW/cm^2
        # Concatenate dataframes
        mobility_ref_dff = pd.concat([mobility_ref_dff, mobility_ref], ignore_index=True)


# ### 2.1.3 Save Al0.25Ga0.75N(25nm)/GaN(305nm) mobilities

# In[ ]:


if resave_mobilities:
    # Create the file if it does not exist
    if not os.path.exists(folder_2deg_mobility_file):
        with pd.ExcelWriter(folder_2deg_mobility_file, engine='openpyxl', mode="w") as writer:  
            mobility_ref_dff.to_excel(writer, sheet_name='ref_Al25Ga75N_GaN')
    else:
        with pd.ExcelWriter(folder_2deg_mobility_file, engine='openpyxl', mode="a", if_sheet_exists="replace") as writer:  
            mobility_ref_dff.to_excel(writer, sheet_name='ref_Al25Ga75N_GaN')


# ## 2.3 AlGaN/AlGaN mobilities and figure-of-merit

# ### 2.3.1 Read 2DEG data from nextnano simulation

# In[ ]:


if redo_mobility_cal:
    data_sheet_name = 'x_y_Lb_scan' #'sim_sweep__AlContentBarrier__AlContentChannel__ThicknessAlGaNBarrier'
    n_2d_df = pd.read_excel(folder_2deg_density_file, sheet_name=data_sheet_name, index_col=0)
    density_2deg_comp = np.array(n_2d_df['2DEG_device'], dtype=float)/1e+14 # cm^-2 to nm^-2 conversion
    channel_comp = np.array(n_2d_df['AlContentChannel'], dtype=float)


# ### 2.3.2 Calculate mobilities

# In[ ]:


if redo_mobility_cal:
    mu2deg = Mobility2DCarrier(compositions=channel_comp, binaries=['AlN', 'GaN'], alloy='AlGaN', system='ternary')


# In[ ]:


if redo_mobility_cal:
    mobility_dff = pd.DataFrame()
    
    for T in T_sweep:
        # Calculate mobilities
        mobility_df = mu2deg.calculate_sheet_mobility(n_2d=density_2deg_comp, rms_roughness=irf_rms_roughness, 
                                                      corr_len=irf_corr_length, n_dis=dislocation_density, 
                                                      f_dis=occup_dislocation, T=T,
                                                      alloy_disordered_effect=alloy_disordered_effect,
                                                      interface_roughness_effect=interface_roughness_effect,
                                                      dislocation_effect=dislocation_effect,
                                                      deformation_potential_effect=deformation_potential_effect,
                                                      piezoelectric_effect=piezoelectric_effect,
                                                      acoustic_phonon_effect=acoustic_phonon_effect,
                                                      polar_optical_phonon_effect=polar_optical_phonon_effect,
                                                      total_mobility=total_mobility,
                                                      mobility_model=mobility_model)
        mu_total = np.array(mobility_df['TOT'], dtype=float) # cm^2.V^-1.s^-1
        # Add temperature column
        mobility_df['T'] = T # K
        # Calculate LFOM
        mobility_df['LFOM'] = mu2deg.calculate_figure_of_merit(density_2deg_comp, mu_total, 
                                                               temp=T, T_corect_bandgap=T_corect_bandgap_in_LFOM) # LFOM: MW/cm^2
        mobility_df['LFOM_2DEG1e13'] = mu2deg.calculate_figure_of_merit(0.1, mu_total, temp=T, 
                                                                        T_corect_bandgap=T_corect_bandgap_in_LFOM) # LFOM: MW/cm^2
        # Concatenate dataframes
        Combined_dff = pd.concat([n_2d_df, mobility_df], axis=1)
        mobility_dff = pd.concat([mobility_dff, Combined_dff], ignore_index=True)


# ### 2.3.3 Save mobilities

# In[ ]:


if resave_mobilities:
    with pd.ExcelWriter(folder_2deg_mobility_file, engine='openpyxl', mode="a", if_sheet_exists="replace") as writer:  
        mobility_dff.to_excel(writer, sheet_name='AlyGaN_AlxGaN')


# ### 2.3.4 Create helper for the data folder

# In[ ]:


if resave_mobilities:
    # Creating helper.txt of this mapping in the DATAs folder
    with open(os.path.join(folder_2deg_mobility_folder,'helper.txt'),'w') as helper_file:
        header_txt = f'''Simulation Software: mobilitypy
2DEG density simulation setup: 1D Schrodinger-Poisson from Nextnano++, Cation-face growth direction
Sample: AlGaN/AlGaN/AlN HEMT
\n
'''
        sheet_column_name = f''' 
### {'='*72}
# Sheet column name: description
### {'='*72}
0th column           : Row indices
comp                 : Al content in AlGaN channel (in mole fraction)
AlContentBarrier     : Al content in AlGaN barrier (in mole fraction)
AlContentChannel     : Al content in AlGaN channel (in mole fraction)
ThicknessAlGaNBarrier: AlGaN barrier thickness (in nm)
2DEG_device          : Intergrated 2D electron gas density over whole device (in carriers/cm^2)
2DEG_BC              : Intergrated 2D electron gas density around the barrier-channel interface (within quantum region) (in carriers/cm^2)
2DEG_SC              : Intergrated 2D electron gas density around the Substrate/Buffer-channel interface (within quantum region) (in carriers/cm^2)
2DHG_device          : Intergrated 2D hole gas density over whole device (in carriers/cm^2)
2DHG_BC              : Intergrated 2D hole gas density around the barrier-channel interface (within quantum region) (in carriers/cm^2)
2DHG_SC              : Intergrated 2D hole gas density around the Substrate/Buffer-channel interface (within quantum region) (in carriers/cm^2)
IFR                  : Interface roughness limited mobility contribution (in cm^2.V^-1.s^-1)
AD                   : Alloy disordered limited mobility contribution (in cm^2.V^-1.s^-1)
DIS                  : Dislocation limited mobility contribution (in cm^2.V^-1.s^-1)
DP                   : Deformation potential limited mobility contribution (in cm^2.V^-1.s^-1)
PE                   : Piezoelectric limited mobility contribution (in cm^2.V^-1.s^-1)
AP                   : Acoustic phonon limited mobility contribution (in cm^2.V^-1.s^-1)
POP                  : Polar optical phonon limited mobility contribution (in cm^2.V^-1.s^-1)
TOT                  : Total mobility (in cm^2.V^-1.s^-1)
T                    : Temperature (in K)
LFOM                 : Lateral figure-of-merit (in MW/cm^2)
LFOM_2DEG1e13        : LFOM when considering constant 2DEG density of 1x10^13 carriers/cm^2
\n### {'='*72}
'''
        helper_file.write(header_txt)
        helper_file.write(f"### {'='*72}\n")
        helper_file.write('# Sheet name: description\n')
        helper_file.write(f"### {'='*72}")
        save_text = f"""
ref_Al25Ga75N_GaN : Mobilities for reference Al25Ga75N(25nm)/GaN HEMT
AlyGaN_AlxGaN : Mobilities for AlyGa1-xN(Lb)/AlxGa1-xN(300nm)/AlN(300nm) HEMT
"""
        helper_file.write(save_text)
        helper_file.write(f"\n")
        helper_file.write(sheet_column_name)


# # 3. Read mobilities, FOMs from saved database

# ## 3.1 Al0.25Ga0.75N(25nm)/GaN(305nm) HEMT LFOM for normalization

# In[ ]:


mobility_ref_dff = pd.read_excel(folder_2deg_mobility_file, index_col=0, sheet_name='ref_Al25Ga75N_GaN') 


# In[ ]:


mobility_ref_df = mobility_ref_dff[mobility_ref_dff['T'] == 300].copy() # Reference data at 300 K
ref_LFOM = mobility_ref_df['LFOM'].iloc[0]


# Conclusion: 2DEG density negligible depend on temperature

# In[ ]:


print(f'{"- Reference sample":<{print_align_space}}: Al0.25Ga0.75N(25nm)/GaN(305nm)')
print(f'{"-- 2DEG density":<{print_align_space}}: {mobility_ref_df['2DEG_device'].iloc[0]:.5e} cm^-2')
print(f'{"-- 2DEG mobility":<{print_align_space}}: {mobility_ref_df['TOT'].iloc[0]:.5e} cm^2.V^-1.s^-1')
print(f'{"-- LFOM":<{print_align_space}}: {ref_LFOM:.5e} MW/cm^2')


# ## 3.2 AlGaN/AlGaN mobilities and figure-of-merit¶

# In[ ]:


mobility_dff = pd.read_excel(folder_2deg_mobility_file, index_col=0, sheet_name='AlyGaN_AlxGaN') #, dtype=float) 
## Not sure why 'IFR' column is read as object
mobility_dff['IFR'] = pd.to_numeric(mobility_dff['IFR'], errors='coerce')
# Collect data for 300 K temperature only
mobility_dff_300k = mobility_dff[mobility_dff['T'] == 300].copy()


# In[ ]:


mobility_dff_300k['2DEG'] =  mobility_dff_300k['2DEG_device']/rescale_2deg_fact
mobility_dff_300k['2DHG'] =  mobility_dff_300k['2DHG_device']/rescale_2deg_fact
mobility_dff_300k['LFOMnorm'] = mobility_dff_300k['LFOM']/ref_LFOM


# In[ ]:


mobility_df_300k = mobility_dff_300k.copy()

clean_up_mobilities = ['IFR', 'DIS', 'AD', 'AP', 'POP', 'TOT']
XXX_ = mobility_df_300k['2DEG'] <1e-2 # Put np.nan for 2DEG densities < 1e-2
for which_column in ['2DEG', 'LFOM', 'LFOMnorm', 'LFOM_2DEG1e13']+clean_up_mobilities:
    mobility_df_300k.loc[XXX_, which_column] = np.nan 


# ## 3.3 Plot AlGaN/AlGaN HEMT mobilities

# In[ ]:


tick_multiplicator = [0.1, 0.05, 0.1, 0.05]


# In[ ]:


df_cap_thickness_group = mobility_df_300k.groupby(['ThicknessAlGaNBarrier']) # Group by thickness


# ### 3.3.1 All contributions together

# In[ ]:


lpltq3d = PlotQuasi3DFuns(save_figure_dir=f'{save_figure_dir}/Together')


# In[ ]:


fig_ncols, fig_nrows = 3, 3
for name, group in df_cap_thickness_group:
    barrier_thickness = f'{name[0]:.2f}'
    print(f'- Plotting barrier thickness = {barrier_thickness} nm')
    fig, axs = plt.subplots(fig_ncols, fig_nrows, figsize=(16.,12.5), 
                            sharey=True,constrained_layout=True)    
    for ii in range(fig_ncols*fig_nrows):
        if ii < len(plot_mobilities):
            which_mobility = plot_mobilities[ii]
        else:
            #axs[ii//3][ii%3].axis("off") 
            continue
        print(f'\t-- Plotting {which_mobility}...')
        xx, yy, zz = group['AlContentChannel'], group['AlContentBarrier'], group[which_mobility]
        vmax = min(mobility_df_300k[which_mobility].max(numeric_only=True), 1e10)
        vmin = mobility_df_300k[which_mobility].min(numeric_only=True) if which_mobility in ['2DEG', 'LFOM', 'LFOMnorm'] else 20
        #vmin, vmax = zz.min()+1e-20, zz.max()
        norm = colors.Normalize(vmin=vmin, vmax=vmax) if which_mobility in ['LFOMnorm','2DEG'] else colors.LogNorm(vmin=vmin, vmax=vmax)
        cbar_mapable = cm.ScalarMappable(norm=norm, cmap='viridis')
        x_label_text_ = x_p_label_text #if ii-6>=0 else ''
        y_label_text_ = x_label_text #if ii%3==0 else ''
        z_label_ = ''
        #if ii%3-2==0: 
        z_label_ = f'$\\mu_{{\\mathrm{{{which_mobility}}}}}$ ($\\mathrm{{cm}}^2\\mathrm{{V}}^{{-1}}\\mathrm{{s}}^{{-1}}$)'
        if which_mobility in ['2DEG']: z_label_= r'n$_{\mathrm{2D}}$ ($\mathrm{10}^{13}\mathrm{cm}^{-2}$)'
        if which_mobility in ['LFOM']: z_label_= z_label[which_mobility]
        if which_mobility == 'LFOMnorm': z_label_= z_label[f'{which_mobility}_2'] 
        

        ##### Plot composition map 
        fig, axt = lpltq3d.Plotq3D(xx,yy,zz, fig=fig, ax=axs[ii//3][ii%3],
                                 xmin=0.475, xmax=0.975, ymin=0.525, ymax=1.025,
                                 x_label='', y_label=y_label_text_,
                                 interpolation_method='linear',
                                 interpolation_points = 20,
                                 tick_multiplicator=tick_multiplicator,
                                 title_label=None, #which_mobility,
                                 cbar_mappable=cbar_mapable, norm=norm,
                                 show_contour_lines=False, marker='s', marker_size=24**2,
                                 cbar_text=z_label_,show_colorbar=False,
                                 plot_controur=0, plot_scatter=1,
                                 savefigure=False, show_plot=False)
        #plt.suptitle(f"Barrier thickness = {barrier_thickness} nm", size=24)
        axt.set_xlabel(x_label_text_, size=18)
        axt.set_ylabel(y_label_text_, size=18)
        cbar = fig.colorbar(cbar_mapable, ax=axt)         
        cbar.ax.set_ylabel(z_label_, size=18)
        #if ii%3 != 0: axs[ii//3][ii%3].get_yaxis().set_visible(False)
    lpltq3d.save_figure(f'Barrier_{barrier_thickness}.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)
    print('- Done', flush=True)
    #break


# ### 3.3.2 Individual plots

# In[ ]:


lpltq3d = PlotQuasi3DFuns(save_figure_dir=f'{save_figure_dir}/Individual')


# In[ ]:


for name, group in df_cap_thickness_group:
    #if name[0] < 49: continue
    barrier_thickness = f'{name[0]:.2f}'
    print(f'- Plotting barrier thickness = {barrier_thickness} nm')
    for which_mobility in plots_list_individual:
        print(f'\t-- Plotting {which_mobility}...')
        xx, yy, zz = group['AlContentChannel'], group['AlContentBarrier'], group[which_mobility]
        vmax = min(mobility_df_300k[which_mobility].max(numeric_only=True), 1e10)
        vmin = mobility_df_300k[which_mobility].min(numeric_only=True) #if which_mobility in ['2DEG', 'LFOMnorm', 'R'] else 20 
        norm = colors.Normalize(vmin=vmin, vmax=vmax) if which_mobility in ['2DEG', '2DHG', 'LFOMnorm', 'TOT'] else colors.LogNorm(vmin=vmin, vmax=vmax)
        cbar_mapable = cm.ScalarMappable(norm=norm, cmap=color_map_plot)
        z_label_= z_label[which_mobility]
        if which_mobility == 'LFOMnorm': z_label_= z_label[f'{which_mobility}_2'] 
        FigFormat_tmp = FigFormatPaper if which_mobility in ['2DEG','2DHG'] else FigFormat
        ##### Plot composition map 
        _ = lpltq3d.Plotq3D(xx, yy, zz, x_label=x_p_label_text, y_label=x_label_text,
                            xmin=0.475, xmax=0.975, ymin=0.525, ymax=1.025,
                            interpolation_method='linear', interpolation_points = 20,
                            tick_multiplicator=tick_multiplicator, norm=norm,color_map=color_map_plot,
                            cbar_mappable=cbar_mapable, show_contour_lines=1,
                            marker='s', marker_size=38**2, cbar_text=z_label_,
                            show_colorbar=True, plot_controur=0, plot_scatter=1,
                            save_file_name=f'Barrier_{barrier_thickness}_{which_mobility}.{FigFormat_tmp}',
                            savefigure=savefigure, show_plot=False, dpi=fig_save_dpi)
    print('- Done', flush=True)
    #break


# In[ ]:


#### This is specific for plottings for paper
for name, group in df_cap_thickness_group:
    barrier_thickness = f'{name[0]:.2f}'
    if name[0] < 49: continue
    for which_mobility in ['2DEG', 'TOT', 'LFOMnorm', 'LFOM']: #
        xx, yy, zz = group['AlContentChannel'], group['AlContentBarrier'], group[which_mobility]
        vmax = min(zz.max(numeric_only=True), 1e10)
        vmin = zz.min(numeric_only=True) #if which_mobility in ['2DEG', 'LFOMnorm', 'R'] else 20 
        norm = colors.Normalize(vmin=vmin, vmax=vmax) if which_mobility in ['2DEG', 'LFOMnorm', 'LFOM', 'TOT'] else colors.LogNorm(vmin=vmin, vmax=vmax)
        cbar_mapable = cm.ScalarMappable(norm=norm, cmap=color_map_plot)
        z_label_= z_label[which_mobility]
        if which_mobility == 'LFOMnorm': z_label_= z_label[f'{which_mobility}_2'] 
        show_colorbar_ = True
        savefigure_ = savefigure
        if which_mobility in ['LFOM','R']:
            show_colorbar_=False
            savefigure_ = False
        figsize=(7.9,6)
        if which_mobility in ['TOT']: figsize=(8.1,6)
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ##### Plot composition map 
        fig, ax = lpltq3d.Plotq3D(xx, yy, zz, fig=fig, ax=ax, x_label=x_p_label_text, y_label=x_label_text,
                            xmin=0.475, xmax=0.975, ymin=0.525, ymax=1.025,
                            interpolation_method='linear', interpolation_points = 100,
                            tick_multiplicator=tick_multiplicator, norm=norm,color_map=color_map_plot,
                            cbar_mappable=cbar_mapable, show_contour_lines=1,
                            marker='s', marker_size=38**2, cbar_text=z_label_,
                            show_colorbar=show_colorbar_, plot_controur=0, plot_scatter=1,
                            save_file_name=f'Barrier_{barrier_thickness}_{which_mobility}_Paper.{FigFormatPaper}',
                            savefigure=savefigure_, show_plot=False, dpi=fig_save_dpi)
        if which_mobility in ['LFOM']:
            cbar = fig.colorbar(cbar_mapable, ax=ax)         
            cbar.ax.set_ylabel(z_label_)
            cbar.formatter.set_powerlimits((0, 0))
            # to get 10^3 instead of 1e3
            cbar.formatter.set_useMathText(True)
            cbar.ax.yaxis.set_offset_position('left')
            cbar.update_ticks()

        lpltq3d.save_figure(f'Barrier_{barrier_thickness}_{which_mobility}_Paper.{FigFormatPaper}', 
                            savefig=not savefigure_, show_plot=False, fig=fig, CountFig=None, dpi=fig_save_dpi) 
        
    print('- Done', flush=True)
    #break


# In[ ]:


plt.close('all')


# In[ ]:


plt2deg = Plottings(save_figure_dir=f'{save_figure_dir}/Others')


# In[ ]:


output_dataframe = {}
fig, ax = plt.subplots(figsize=(8,7))
tick_multiplicator = [0.1, 0.05, 1, 0.5]
markers_ = ['o', '*', 'd', '^', 's', 'p']*3
ii=0
for name, group in mobility_dff_300k.groupby(['ThicknessAlGaNBarrier']):
    output_dataframe[f'{name[0]:.2f}'] = {}
    group['AlContentContrast'] = (group['AlContentBarrier'] - group['AlContentChannel']) #/ group['AlContentChannel']
    
    #=====================================================================================
    if name[0] in [5,10,20,50]:
        ax.scatter(group['AlContentContrast'], group['2DEG'], marker=markers_[ii],label=f'{name[0]:.2f} nm',s=100)
        ii+=1
    
    #=====================================================================================
    tmp_df = group[group[f'2DEG']>0.01]
    A = np.vstack([tmp_df['AlContentContrast'], np.ones(len(tmp_df))]).T
    m, c = np.linalg.lstsq(A, tmp_df[f'2DEG'])[0]
    y_eq_zero_point = -c/m
    output_dataframe[f'{name[0]:.2f}'][f'2DEG_slope']       = f'{m:.3f}' 
    output_dataframe[f'{name[0]:.2f}'][f'2DEG_y_intersect'] = f'{c:.3f}'
    output_dataframe[f'{name[0]:.2f}'][f'2DEG_x_intersect'] = f'{y_eq_zero_point:.3f}'
    new_x = np.insert(np.array(tmp_df['AlContentContrast']),0,y_eq_zero_point)
    #ax.plot(new_x, m*new_x + c, 'gray')
    if name[0] in [5,10,20,50]: 
        ax.plot(new_x, m*new_x + c, 'gray')
    
    ax.legend(ncols=1, labelspacing=0.001, columnspacing=0.1, handletextpad=0.01)
    ax.set_xlabel('Al composition contrast, $\\Delta_{yx}$')
    ax.set_ylabel(z_label['2DEG'])
    ax.axhline(color='k',ls='--')
ax.axhline(y=1, c='c', ls='--')
ax.set_xlim(0.05,0.5)
plt2deg.set_tickers(ax, tick_multiplicator)
circle = Ellipse((0.35,1.6), width=0.025, height=0.35, 
                        edgecolor='r', fc='None', lw=2)
ax.add_patch(circle)
plt2deg.save_figure(f'Comp_contrast_thickness_2DEG.{FigFormatPaper}', savefig=savefigure, show_plot=True,
                    fig=fig, CountFig=None, dpi=fig_save_dpi)      


# In[ ]:


#=============================================================================================
output_data = pd.DataFrame.from_dict(output_dataframe, orient='index', dtype=float)
## Plot cut-off composition contrast for each barrier thickness
tick_multiplicator = [0.1, 0.05,10,5]
XX = np.array(output_data.index, dtype=float)
YY = np.array(output_data['2DEG_x_intersect'], dtype=float)
fig, axs = plt.subplots(constrained_layout=True)
axs.plot(YY, XX, 'ko-', ms=12)
#axs.set_yscale('log')
#axs.set_ylabel('Critical barrier thickness (nm)')
#axs.set_xlabel('Critical composition contrast')
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left', title='Al-content') 
axs.set_xlim(0.05, 0.54)
plt2deg.set_tickers(axs, tick_multiplicator)
 

fit_intercept = np.array(output_data['2DEG_y_intersect'], dtype=float)
fit_slope = np.array(output_data['2DEG_slope'], dtype=float)
YY = (1 - fit_intercept)/fit_slope
axs.plot(YY, XX, 'cs-', ms=12)
#axs.set_yscale('log')
axs.set_ylabel('Critical barrier thickness (nm)')
#axs.set_ylabel('Barrier thickness, L$_\\mathrm{B}$ (nm)')
axs.set_xlabel('Al composition contrast, $\\Delta_{yx}$')

plt2deg.save_figure(f'Critical_comp_contrast_2DEG.{FigFormatPaper}', savefig=savefigure, show_plot=True,
                    fig=fig, CountFig=None, dpi=fig_save_dpi) 

# _ = lplt1d.PlotSweepsData(XX, YY, x_label=, y_label=,
#                           tick_multiplicator=tick_multiplicator,
#                           FigDpi=FigDpi, FigFormat=FigFormat, color='k', marker='o',
#                           figs_path=f'{output_figs_sweep}/Others', filename=f'Critical_comp_contrast_2DEG',
#                           savefigure=savefigure)


# ## 3.4 Plots for AlN/AlGaN

# In[ ]:


plt2deg = Plottings(save_figure_dir=f'{save_figure_dir}/Others')


# ### 3.4.1 Mobility plot

# In[ ]:


mobility_df_300k_ = mobility_dff_300k[mobility_dff_300k['AlContentBarrier']>0.99]
df_channel_comp_group = mobility_df_300k_.groupby(['AlContentChannel']) # Group by channel composition


# In[ ]:


x_label = 'Barrier thickness, L$_\\mathrm{B}$ (nm)'
y_label = z_label['TOT'] #r'Electron mobility ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)'
save_file_name_ = f'mu_Tb_AlN_AlGaN.{FigFormatPaper}'
tick_multiplicator=[10, 5, 50, 25]
plot_selected_compositions = ['0.50', '0.65', '0.75', '0.85', '0.90', '0.95']
markers_ = ['o', '*', 's', '^', 'd', 'p']*3


# In[ ]:


fig, axs = plt.subplots(constrained_layout=True)
ii = 0
for name, group in df_channel_comp_group:
    xx, zz = group['ThicknessAlGaNBarrier'],  group['TOT']
    if f'{name[0]:0.2f}' in plot_selected_compositions: 
        axs.plot(xx, zz, '-', marker=markers_[ii],label=f'x={name[0]:.2f}', ms=12)
        ii+=1
for name, group in df_channel_comp_group:
    xx, zz = group['ThicknessAlGaNBarrier'],  group['TOT']
    if f'{name[0]:0.2f}' in plot_selected_compositions: 
        if any(zz<1e-5):
            undefined_mobilities = np.where(zz<1)[0]
            axs.plot(xx.iloc[undefined_mobilities], zz.iloc[undefined_mobilities], '-', color='gray',
                     marker=markers_[ii], mec='k', mfc='white', ms=12)
        ii+=1
#axs.set_yscale('log')
axs.set_xlabel(x_label)
axs.set_ylabel(y_label)
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left', title='Al-content') 
axs.set_xlim(5,50)
plt2deg.set_tickers(axs, tick_multiplicator)
plt2deg.save_figure(save_file_name_, fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# In[ ]:


save_file_name_ = f'channel_comp_legends.{FigFormatPaper}'
# Now create a new image with legends only
# adjust the figure size as necessary
fig_leg = plt.figure(figsize=(1,1))
ax_leg = fig_leg.add_subplot(111)
# add the legend from the previous axes
ax_leg.legend(*axs.get_legend_handles_labels(), ncol=3, loc='center')
# hide the axes frame and the x/y labels
ax_leg.axis('off')
#ax_leg.set_frame_on(False)
plt2deg.save_figure(save_file_name_, fig=fig_leg, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# In[ ]:


fig, axs = plt.subplots(constrained_layout=True) #figsize=(9.5,6)
ii = 0
for name, group in df_channel_comp_group:
    if name[0]< 0.74 or name[0]> 0.91: 
        axs.plot([None], [None], label=None)
        continue
    xx, zz = group['ThicknessAlGaNBarrier'],  group['TOT']
    axs.plot(xx, zz, '-', marker=markers_[ii],label=f'{name[0]:.2f}', ms=12)
    ii+=1
#axs.set_yscale('log')
axs.set_xlabel(x_label)
axs.set_ylabel(y_label)
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left') 
plt2deg.set_tickers(axs, tick_multiplicator)
axs.set_xlim(5,50)
plt2deg.save_figure(f'zoom_{save_file_name_}', fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# ### 3.4.2 Individual mobility contributions for AlN/Al0.9Ga0.1N

# In[ ]:


Al90Ga10N_mu_df_300k_ = mobility_df_300k_[(mobility_df_300k_['AlContentChannel']>0.89) & (mobility_df_300k_['AlContentChannel']<0.91)]


# In[ ]:


fig, axs = plt.subplots(figsize=(9.5,6),constrained_layout=True)
ii=0
for name in plot_mobilities[:-3]:
    axs.plot(Al90Ga10N_mu_df_300k_['ThicknessAlGaNBarrier'], Al90Ga10N_mu_df_300k_[name], '-', marker=markers_[ii],label=f'{name}', ms=12)
    ii+=1
axs.set_yscale('log')
axs.set_xlabel(x_label)
axs.set_ylabel(y_label)
axs.legend(bbox_to_anchor=(1, 0.5),loc='center left') 
axs.set_ylim(100, 100000)
axs.set_xlim(5,50)
plt2deg.set_tickers(axs, [10,5,None, None])
plt2deg.save_figure(f'mu_contibs_Tb_AlN_Al0.9Ga0.1N.{FigFormat}', fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# ### 3.4.3 FOM plot

# In[ ]:


x_label = x_label_text_2d
y_label = z_label['LFOMnorm_2']
save_file_name_ = f'lfom_Tb_AlN_AlGaN.{FigFormatPaper}'
tick_multiplicator=[10, 5, 0.5, 0.25]


# In[ ]:


fig, axs = plt.subplots(constrained_layout=True)
ii=0
for name, group in df_channel_comp_group:
    xx, zz = group['ThicknessAlGaNBarrier'],  group['LFOMnorm']
    if f'{name[0]:0.2f}' in plot_selected_compositions:
        axs.plot(xx, zz, '-', marker=markers_[ii],label=f'{name[0]:.2f}', ms=12)
        ii+=1
for name, group in df_channel_comp_group:
    xx, zz = group['ThicknessAlGaNBarrier'],  group['LFOMnorm']
    if f'{name[0]:0.2f}' in plot_selected_compositions: 
        if any(zz<1e-5):
            undefined_mobilities = np.where(zz<1e-5)[0]
            axs.plot(xx.iloc[undefined_mobilities], zz.iloc[undefined_mobilities], '-', color='gray',
                     marker=markers_[ii], mec='k', mfc='white', ms=12)
        ii+=1
#axs.set_yscale('log')
axs.set_xlabel(x_label)
axs.set_ylabel(y_label)
axs.axhline(y=1.0, c='k', ls='--')
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left') 
plt2deg.set_tickers(axs, tick_multiplicator)
axs.set_xlim(5,50)
plt2deg.save_figure(save_file_name_, fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# In[ ]:


fig, axs = plt.subplots(constrained_layout=True)
ii=0
for name, group in df_channel_comp_group:
    if name[0]< 0.74 or name[0]> 0.91: 
        axs.plot([None], [None])
        continue
    xx, zz = group['ThicknessAlGaNBarrier'],  group['LFOMnorm']
    axs.plot(xx, zz, '-', marker=markers_[ii],label=f'{name[0]:.2f}', ms=12)
    ii+=1
#axs.set_yscale('log')
axs.axhline(y=1.0, c='k', ls='--')
axs.set_xlabel(x_label)
axs.set_ylabel(y_label)
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left') 
plt2deg.set_tickers(axs, tick_multiplicator)
axs.set_xlim(5,50)
plt2deg.save_figure(f'zoom_{save_file_name_}', fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# ### 3.4.4 2DEG density plot

# In[ ]:


x_label = x_label_text_2d
y_label = z_label['2DEG']
save_file_name_ = f'deg_Tb_AlN_AlGaN.{FigFormatPaper}'
tick_multiplicator=[10, 5, None, None]


# In[ ]:


fig, axs = plt.subplots(figsize=(8,6),constrained_layout=True)
ii=0
for name, group in df_channel_comp_group:
    xx, zz = group['ThicknessAlGaNBarrier'],  group['2DEG']
    if f'{name[0]:0.2f}' in plot_selected_compositions:
        axs.plot(xx, zz, '-', marker=markers_[ii],label=f'{name[0]:.2f}', ms=12)
        ii+=1
#axs.set_yscale('log')
axs.set_xlabel(x_label)
axs.set_ylabel(y_label)
axs.axhline(y=1.0, c='k', ls='--')
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left') 
plt2deg.set_tickers(axs, tick_multiplicator)
axs.set_xlim(5,50)
plt2deg.save_figure(save_file_name_, fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# In[ ]:


fig, axs = plt.subplots(figsize=(8,6),constrained_layout=True)
ii=0
for name, group in df_channel_comp_group:
    xx, zz = group['ThicknessAlGaNBarrier'],  group['2DEG']
    if f'{name[0]:0.2f}' in ['0.50']: 
        axs.plot(xx, zz, 'r-', marker=markers_[ii],label=f'{name[0]:.2f}', ms=12)
        ii+=1
axs.set_xlabel(x_label_text_2d)
axs.set_ylabel(y_label)
axs.axhline(y=1.0, c='k', ls='--')
#axs.legend(bbox_to_anchor=(1, 0.5),loc='center left') 
plt2deg.set_tickers(axs, [10, 5, None, None])
#axs.set_xlim(5,50)
plt2deg.save_figure(f'two_deg_Tb_AlN_Al0.50GaN.{FigFormatPaper}', fig=fig, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# ### 3.4.5 Mobility contributions for AlN(50nm)/AlGaN

# In[ ]:


aln_algan_300k = mobility_dff_300k[(mobility_dff_300k['AlContentBarrier']>0.99) & (mobility_dff_300k['ThicknessAlGaNBarrier']>49)]
aln_algan_dff = aln_algan_300k[['AlContentChannel','IFR','AD','DIS','DP','PE','AP','POP','TOT']].copy()
aln_algan_dff.rename(columns={'AlContentChannel':'comp'}, inplace=True)

x_label = x_p_label_text
y_label = z_label['TOT']
save_file_name_ = f'mu_AlN50_AlGaN.{FigFormatPaper}'
tick_multiplicator=[0.1, 0.05, None, None]


# In[ ]:


# Mobility contributions at 300K
save_file_name_ = f'Mobility_contribs_AlN50AlGaN_300K.{FigFormatPaper}'
x_label = x_p_label_text
y_label = z_label['TOT']
fig, ax,_ = plt2deg.plot_2d_carrier_mobilities(aln_algan_dff, save_file_name=save_file_name_,
                                               ymin=5e1, ymax=2e5, xmax=0.9, xmin=0.5, y_scale_log=True, 
                                               annotate_pos=(2,2), annotatetextoffset=(0,-20),show_right_ticks=True,
                                               mode='2d_carrier_mobility', yaxis_label=y_label, xaxis_label=x_label,
                                               color=None, color_map='viridis', savefig=0, show_plot=False)
plt2deg.set_tickers(ax, tick_multiplicator)
plt2deg.save_figure(save_file_name_, fig=fig, savefig= savefigure, dpi=fig_save_dpi, show_plot=False)


# In[ ]:


fig, ax = plt.subplots(figsize=(9,6), constrained_layout=True)
XX = np.array(aln_algan_300k['AlContentChannel'], dtype=float)
ax.plot(XX, aln_algan_300k['2DEG_device']/rescale_2deg_fact, 'o-', c='k')
ax.set_xlim(0.5,0.9)
ax.set_xlabel(x_p_label_text)
ax.set_ylabel(z_label['2DEG'])
plt2deg.set_tickers(ax, tick_multiplicator)

# Plotting mobility
ax2 = ax.twinx() 
ax2.plot(XX, aln_algan_300k['TOT'], 'o-', c='r')

# Plotting bandgap
bandgap_ = (6.20*XX + 3.43*(1-XX) + 0.7*XX*(1-XX))**5 
ax2.plot(XX, bandgap_, 'o-', c='m')

ax2.set_ylabel(z_label['TOT'], c='r')
ax2.set_yscale('log')
ax2.set_ylim(ymin=70, ymax=1e4)
ax2.annotate(r'E$_\mathrm{g}^{5}$', xy=(XX[3], bandgap_[3]+1e3), color='m')

# Plotting LFOM
LFOM_300k = np.array(aln_algan_300k['LFOM'], dtype=float)/ref_LFOM 
ax.plot(XX, LFOM_300k, 'o-', c='b')
ax.annotate(r'LFOM$_{\mathrm{norm}}$', xy=(XX[3], LFOM_300k[3]-0.2), color='b')

plt2deg.save_figure(f'AlN50AlGaN_300K_mu_FOM.{FigFormat}', savefig=savefigure, show_plot=True,
                    fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
XX = np.array(aln_algan_300k['AlContentChannel'], dtype=float)
yy_2deg_density = np.array(aln_algan_300k['2DEG_device']/rescale_2deg_fact, dtype=float)
ax.plot(XX, yy_2deg_density, 's-.', c='g', ms=12)
ax.axhline(y=1, c='k', ls='--')
ax.annotate(z_label['2DEG'][14:], xy=(XX[3], yy_2deg_density[3]+0.08), color='g')
circle0 = Ellipse((XX[2], yy_2deg_density[2]), width=0.015, height=0.2, 
                        edgecolor='g', fc='None', lw=2)
ax.add_patch(circle0)
ax.annotate("", xy=(XX[2]-0.05, yy_2deg_density[2]-0.1), xytext=(XX[2], yy_2deg_density[2]-0.1),
             arrowprops=dict(arrowstyle="->",color='g', linewidth=2))
ax.set_xlim(0.5,0.9)
ax.set_xlabel(x_p_label_text)
ax.set_ylabel('')
plt2deg.set_tickers(ax, tick_multiplicator)

#=========================================================================
ax2 = ax.twinx() 
ax2.set_ylabel('', c='r')
ax2.set_yscale('log')
ax2.set_ylim(ymin=50, ymax=1e4)

# Plotting mobility
yy_mobility = np.array(aln_algan_300k['TOT'], dtype=float)
ax2.plot(XX, yy_mobility, '*-.', c='c', ms=12)
ax2.annotate(z_label['TOT'][15:], xy=(XX[3], yy_mobility[3]-25), color='c')
circle1 = Ellipse((XX[6], yy_mobility[6]), width=0.014, height=30, 
                        edgecolor='c', fc='None', lw=2)
ax2.add_patch(circle1)
ax2.annotate("", xy=(XX[6]+0.05, yy_mobility[6]-15), xytext=(XX[6], yy_mobility[6]-15),
             arrowprops=dict(arrowstyle="->",color='c', linewidth=2))

# Plotting bandgap
bandgap_ = (6.25*XX + 3.51*(1-XX) + 0.7*XX*(1-XX))**5 
ax2.plot(XX, bandgap_, '^-.', c='m', ms=12)
ax2.annotate(r'E$_\mathrm{g}^{5}$ (eV$^{5}$)', xy=(XX[3], bandgap_[3]+1.2e3), color='m')
circle2 = Ellipse((XX[6], bandgap_[6]), width=0.015, height=1.6e3, 
                        edgecolor='m', fc='None', lw=2)
ax2.add_patch(circle2)
ax2.annotate("", xy=(XX[6]+0.05, bandgap_[3]+1.e3), xytext=(XX[6], bandgap_[3]+1.e3),
             arrowprops=dict(arrowstyle="->",color='m', linewidth=2))

# Plotting LFOM
LFOM_300k = np.array(aln_algan_300k['LFOM'], dtype=float)/ref_LFOM 
ax.plot(XX, LFOM_300k, 'o-', c='r', ms=12)
ax.annotate(r'LFOM$^{\mathrm{B}}_{\mathrm{norm}}$', xy=(XX[1], LFOM_300k[1]+0.3), color='r')
circle3 = Ellipse((XX[2], LFOM_300k[2]), width=0.015, height=0.2, 
                        edgecolor='r', fc='None', lw=2)
ax.add_patch(circle3)
ax.annotate("", xy=(XX[2]-0.05, LFOM_300k[2]+0.1), xytext=(XX[2], LFOM_300k[2]+0.1),
             arrowprops=dict(arrowstyle="->",color='r', linewidth=2))

LFOM_300k_cnst = np.array(aln_algan_300k['LFOM_2DEG1e13'], dtype=float)/ref_LFOM 
ax.plot(XX, LFOM_300k_cnst, 'd-', c='b', ms=12)
ax.annotate(r'LFOM$^{\mathrm{A}}_{\mathrm{norm}}$', xy=(XX[1], LFOM_300k_cnst[1]+0.3), color='b')
circle4 = Ellipse((XX[2], LFOM_300k_cnst[2]), width=0.015, height=0.2, 
                        edgecolor='b', fc='None', lw=2)
ax.add_patch(circle4)
ax.annotate("", xy=(XX[2]-0.05, LFOM_300k_cnst[2]+0.1), xytext=(XX[2], LFOM_300k_cnst[2]+0.1),
             arrowprops=dict(arrowstyle="->",color='b', linewidth=2))

plt2deg.set_tickers(ax, [0.1,0.05,0.5,0.25])

plt2deg.save_figure(f'AlN50AlGaN_300K_mu_FOM_paper.{FigFormatPaper}', savefig=savefigure, show_plot=True,
                    fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


ii=0
for ii in range(5):
    fig, ax = plt.subplots(figsize=(8,6.), constrained_layout=True)
    XX = np.array(aln_algan_300k['AlContentChannel'], dtype=float)

    if ii >2:
        yy_2deg_density = np.array(aln_algan_300k['2DEG_device']/rescale_2deg_fact, dtype=float)
        ax.plot(XX, yy_2deg_density, 's-', c='k', ms=12)
        ax.annotate(z_label['2DEG'][14:], xy=(XX[3], yy_2deg_density[3]+0.08), color='k')
        circle0 = Ellipse((XX[2], yy_2deg_density[2]), width=0.015, height=0.2, 
                                edgecolor='k', fc='None', lw=2)
        ax.add_patch(circle0)
        ax.annotate("", xy=(XX[2]-0.05, yy_2deg_density[2]-0.1), xytext=(XX[2], yy_2deg_density[2]-0.1),
                     arrowprops=dict(arrowstyle="->",color='k', linewidth=2))
        
    ax.axhline(y=1, c='k', ls='--')
    ax.set_xlim(0.5,0.9)
    ax.set_ylim(-0.12,3.55)
    ax.set_xlabel(x_p_label_text)
    ax.set_ylabel('')
    plt2deg.set_tickers(ax, [0.1,0.05,0.5,0.25])
    
    #=========================================================================
    ax2 = ax.twinx() 
    ax2.set_ylabel('', c='r')
    ax2.set_yscale('log')
    ax2.set_ylim(ymin=50, ymax=1e4)
    
    # Plotting mobility
    if ii >1:
        yy_mobility = np.array(aln_algan_300k['TOT'], dtype=float)
        ax2.plot(XX, yy_mobility, '*-', c='c', ms=12)
        ax2.annotate(z_label['TOT'][15:], xy=(XX[3], yy_mobility[3]-25), color='c')
        circle1 = Ellipse((XX[6], yy_mobility[6]), width=0.014, height=30, 
                                edgecolor='c', fc='None', lw=2)
        ax2.add_patch(circle1)
        ax2.annotate("", xy=(XX[6]+0.05, yy_mobility[6]-15), xytext=(XX[6], yy_mobility[6]-15),
                     arrowprops=dict(arrowstyle="->",color='c', linewidth=2))
    if ii >0:
        # Plotting bandgap
        bandgap_ = (6.25*XX + 3.51*(1-XX) + 0.7*XX*(1-XX))**5 
        ax2.plot(XX, bandgap_, '^-', c='m', ms=12)
        ax2.annotate(r'E$_\mathrm{g}^{5}$ (eV$^{5}$)', xy=(XX[3], bandgap_[3]+1.2e3), color='m')
        circle2 = Ellipse((XX[6], bandgap_[6]), width=0.015, height=1.6e3, 
                                edgecolor='m', fc='None', lw=2)
        ax2.add_patch(circle2)
        ax2.annotate("", xy=(XX[6]+0.05, bandgap_[3]+1.e3), xytext=(XX[6], bandgap_[3]+1.e3),
                     arrowprops=dict(arrowstyle="->",color='m', linewidth=2))
    
    # Plotting LFOM
    LFOM_300k = np.array(aln_algan_300k['LFOM'], dtype=float)/ref_LFOM 
    ax.plot(XX, LFOM_300k, 'o-', c='r', ms=12)
    ax.annotate(r'LFOM$_{\mathrm{norm}}$', xy=(XX[1], LFOM_300k[1]+0.3), color='r')
    circle3 = Ellipse((XX[2], LFOM_300k[2]), width=0.015, height=0.2, 
                            edgecolor='r', fc='None', lw=2)
    ax.add_patch(circle3)
    ax.annotate("", xy=(XX[2]-0.05, LFOM_300k[2]+0.1), xytext=(XX[2], LFOM_300k[2]+0.1),
                 arrowprops=dict(arrowstyle="->",color='r', linewidth=2))
    if ii >3:
        LFOM_300k_cnst = np.array(aln_algan_300k['LFOM_2DEG1e13'], dtype=float)/ref_LFOM 
        ax.plot(XX, LFOM_300k_cnst, 'd-', c='b', ms=12)
        ax.annotate(r'LFOM$^*_{\mathrm{norm}}$', xy=(XX[1], LFOM_300k_cnst[1]+0.3), color='b')
        circle4 = Ellipse((XX[2], LFOM_300k_cnst[2]), width=0.015, height=0.2, 
                                edgecolor='b', fc='None', lw=2)
        ax.add_patch(circle4)
        ax.annotate("", xy=(XX[2]-0.05, LFOM_300k_cnst[2]+0.1), xytext=(XX[2], LFOM_300k_cnst[2]+0.1),
                     arrowprops=dict(arrowstyle="->",color='b', linewidth=2))
    
    plt2deg.save_figure(f'{ii}_AlN50AlGaN_300K_mu_FOM_paper.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)
    ii+=1
    #break


# In[ ]:


fig, ax = plt.subplots(constrained_layout=True)
XX = np.array(aln_algan_300k['AlContentChannel'], dtype=float)
#ax.plot(XX, aln_algan_300k['2DEG_device']/rescale_2deg_fact, 'o-', c='k')
ax.set_xlim(0.5,0.9)
ax.set_xlabel(x_p_label_text)
ax.set_ylabel(z_label['2DEG'])
plt2deg.set_tickers(ax, tick_multiplicator)

# # Plotting mobility
# ax2 = ax.twinx() 
# ax2.plot(XX, aln_algan_300k['TOT'], 'o-', c='r')

# # Plotting bandgap
# bandgap_ = (6.20*XX + 3.43*(1-XX) + 0.7*XX*(1-XX))**5 
# ax2.plot(XX, bandgap_, 'o-', c='m')

# ax2.set_ylabel(z_label['TOT'], c='r')
# ax2.set_yscale('log')
# ax2.set_ylim(ymin=70, ymax=1e4)
# ax2.annotate(r'E$_\mathrm{g}^{5}$', xy=(XX[3], bandgap_[3]+1e3), color='m')

# # Plotting sheet resistance
# R_300k = np.array(aln_algan_300k['R'], dtype=float) #/20000 #1/20000 is to rescale
# ax2.plot(XX, R_300k, 'o-', c='g')
# ax2.annotate(r'R', xy=(XX[3], R_300k[3]-1e3), color='g')

# Plotting LFOM
LFOM_300k = np.array(aln_algan_300k['LFOM'], dtype=float)/ref_LFOM 
ax.plot(XX, LFOM_300k, 'o-', c='r')
LFOM_300k = np.array(aln_algan_300k['LFOM_2DEG1e13'], dtype=float)/ref_LFOM 
ax.plot(XX, LFOM_300k, 'o-', c='b')
ax.set_ylabel(z_label['LFOMnorm_2'])
ax.axhline(y=1, c='k', ls='--')
#ax.annotate(r'LFOM$_{\mathrm{norm}}$', xy=(XX[3], LFOM_300k[3]-0.2), color='b')

plt2deg.save_figure(f'AlN50AlGaN_300K_mu_FOM_paper1.{FigFormat}', savefig=savefigure, show_plot=True,
                    fig=fig, CountFig=None, dpi=fig_save_dpi)


# # 4. 2DEG mobility plottings - temperature variation

# ## 4.1 Read mobilities, FOMs from saved database

# 2DEG density negligibly depend on Temperature for both GaN- and AlGaN-channel HEMT

# ### 4.1.1 Al0.25Ga0.75N(25nm)/GaN(305nm) HEMT LFOM for normalization

# In[ ]:


mobility_ref_dff = pd.read_excel(folder_2deg_mobility_file, index_col=0, sheet_name='ref_Al25Ga75N_GaN') 


# ### 4.1.2 AlGaN/AlGaN mobilities and figure-of-merit

# In[ ]:


mobility_dff = pd.read_excel(folder_2deg_mobility_file, index_col=0, sheet_name='AlyGaN_AlxGaN') 


# ### 4.1.3 Collect AlN/AlGaN mobility and lateral figure-of-merit
# 
# Note: Highest LFOM is found for AlN in Barrier and barrier thickness of 50 nm always
# 
# Note: AlN/AlGaN-channel/AlN normalized LFOM are normalized w.r.t Al0.25Ga0.75N(25nm)/GaN(305nm) GaN-channel LFOM at the same temperature

# In[ ]:


mobility_dfff = {}
for tb in [50, 25]:
    if tb == 50:
        mobility_df = mobility_dff[(mobility_dff['AlContentBarrier']>0.99) & (mobility_dff['ThicknessAlGaNBarrier']>49)].copy()
    elif tb == 25:
        mobility_df = mobility_dff[(mobility_dff['AlContentBarrier']>0.99) & (mobility_dff['ThicknessAlGaNBarrier']>24) & (mobility_dff['ThicknessAlGaNBarrier']<26)].copy()
    else:
        break
    mobility_df.reset_index(drop=True, inplace=True)
    mobility_df['2DEG'] =  mobility_df['2DEG_device']/rescale_2deg_fact
    mobility_df['2DHG'] =  mobility_df['2DHG_device']/rescale_2deg_fact
    mobility_df['LFOMnorm'] = np.nan

    mobility_dfff[tb] = mobility_df


# In[ ]:


plots_data = {}
for tb, mobility_df in mobility_dfff.items():
    for index, row in mobility_df.iterrows():
        ref_LFOM = mobility_ref_dff.loc[mobility_ref_dff['T']==row['T'], 'LFOM'].iloc[0]
        mobility_df.loc[index,'LFOMnorm'] = row['LFOM']/ref_LFOM

    plots_data[tb] = {}
    plots_data[tb]['comp_'] = np.array(mobility_df['AlContentChannel'], dtype=float)
    plots_data[tb]['temperature_'] = np.array(mobility_df['T'], dtype=float)
    plots_data[tb]['lfom_'] = np.array(mobility_df['LFOM'], dtype=float)
    plots_data[tb]['lfom_norm_'] = np.array(mobility_df['LFOMnorm'], dtype=float)
    tot_mu = np.array(mobility_df['TOT'], dtype=float)
    tot_mu_cp = tot_mu.copy()
    tot_mu_cp[tot_mu_cp <1] = np.nan # put to 1 below 1 cm^2V^-1s^-1 mobilities. For better plottings
    plots_data[tb]['tot_mu'] = tot_mu_cp


# In[ ]:


plots_data_maxms = {}
for tb, mobility_df in mobility_dfff.items():
    XX, YY = [], [] # where norm_lfom > 1
    XXX, YYY = [], [] # where norm_lform is highest at each T
    df_group_T = mobility_df.groupby(['T'])
    for name, group in df_group_T:
        pp_ = group['LFOM'].argmax()
        XXX.append(group.iloc[pp_]['AlContentChannel'])
        YYY.append(name[0])
        tmp_ = group[group['LFOMnorm']>1]
        if len(tmp_) > 0:
            XX.append(np.array(tmp_['AlContentChannel']))
            YY.append(np.array(tmp_['T']))
    plots_data_maxms[tb] = (np.concatenate(XX), np.concatenate(YY), XXX, YYY)


# ## 4.2 Plot AlN/AlGaN HEMT mobilities and FOM

# In[ ]:


tick_multiplicator = [0.1, 0.05, 200, 100]


# In[ ]:


lpltq3d = PlotQuasi3DFuns(save_figure_dir=f'{save_figure_dir}/Temperature')


# In[ ]:


for tb in mobility_dfff:
    comp_,temperature_,tot_mu = plots_data[tb]['comp_'], plots_data[tb]['temperature_'], plots_data[tb]['tot_mu']
    fig, ax = plt.subplots(figsize=(8,5.5), constrained_layout=True)
    fig, _ = lpltq3d.Plotq3D(comp_,temperature_,tot_mu, fig=fig, ax=ax,
                             xmin=0.475, xmax=0.975, ymin=-40, ymax=850,
                             x_label=x_p_label_text, y_label= 'T (K)',
                             tick_multiplicator=tick_multiplicator,
                             show_contour_lines=False, marker='s', marker_size=38**2,
                             cbar_text=z_label['TOT'],show_colorbar=True,
                             plot_controur=0, plot_scatter=1, show_plot=False)
    lpltq3d.save_figure(f'Highest_FOM_mobility_Tb{tb}.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


# Plot mobilities in log scale
for tb in mobility_dfff:
    comp_,temperature_,tot_mu = plots_data[tb]['comp_'], plots_data[tb]['temperature_'], plots_data[tb]['tot_mu']
    fig, ax = plt.subplots(figsize=(8,5.5), constrained_layout=True)
    vmin = np.nanmin(tot_mu)
    vmax = np.nanmax(tot_mu)
    norm = colors.LogNorm(vmin=vmin, vmax=vmax) 
    cbar_mapable = cm.ScalarMappable(norm=norm, cmap='viridis')
    fig, _ = lpltq3d.Plotq3D(comp_,temperature_,tot_mu, fig=fig, ax=ax,
                             xmin=0.475, xmax=0.975, ymin=-40, ymax=850,
                             x_label=x_p_label_text, y_label= 'T (K)',
                             tick_multiplicator=tick_multiplicator,cbar_mappable=cbar_mapable,
                             norm=norm, vmin=vmin, vmax=vmax,
                             show_contour_lines=False, marker='s', marker_size=38**2,
                             cbar_text=z_label['TOT'],show_colorbar=True,
                             plot_controur=0, plot_scatter=1, show_plot=False)
    lpltq3d.save_figure(f'Highest_FOM_log_mobility_Tb{tb}.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


for tb in mobility_dfff:
    comp_,temperature_,lfom_norm_ = plots_data[tb]['comp_'], plots_data[tb]['temperature_'], plots_data[tb]['lfom_norm_']
    XX, YY, XXX, YYY = plots_data_maxms[tb]
    fig, ax = plt.subplots(figsize=(7.6,5.5), constrained_layout=True)
    fig, _ = lpltq3d.Plotq3D(comp_,temperature_,lfom_norm_, fig=fig, ax=ax,
                             xmin=0.475, xmax=0.975, ymin=-40, ymax=850,
                             x_label=x_p_label_text, y_label= 'T (K)',
                             tick_multiplicator=tick_multiplicator,
                             show_contour_lines=False, marker='s', marker_size=38**2,
                             cbar_text=z_label['LFOMnorm_2'],show_colorbar=True,
                             plot_controur=0, plot_scatter=1, show_plot=False)
    ax.scatter(XX, YY, marker='x', s=20**2, c='k')
    ax.scatter(XXX, YYY, marker='o', s=10**2, c='k')
    lpltq3d.save_figure(f'Highest_FOM_norm_Tb{tb}.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


markers_ = ['^', '*', 'd', 'o', 's', 'p']*3
markers_ = ['.','>','D','o','*','h','p','s','^', '*', 'd', 'o', 's', 'p']*3
colors_ = ['black','darkred', 'indianred','red', 'chocolate', 'peru', 'orange', 'gold', 'yellow']
for tb, mobility_df in mobility_dfff.items():
    ii = 0
    df_group_T = mobility_df.groupby(['T'])
    fig, ax = plt.subplots(constrained_layout=True)
    for name, group in df_group_T:
        #if name[0] in [100,200,300,400,600,800]:
        markeredgecolor='darkgoldenrod' if name[0]>400 else 'none'
        ax.plot(group['AlContentChannel'], group['LFOMnorm'], '-',marker=markers_[ii],
                color=colors_[ii],label=f'T={name[0]} K', ms=12, mec=markeredgecolor,mew=0.5)
        ii+=1
    #ax.set_yscale('log')
    ax.set_xlim(0.5,0.9)
    ax.set_xlabel(x_p_label_text)
    ax.set_ylabel(z_label['LFOMnorm_2'])
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5)) 
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    #ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.axhline(y=1,c='k', ls='--')
    lpltq3d.save_figure(f'Highest_FOM_norm_line_Tb{tb}.{FigFormatPaper}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


save_file_name_ = f'temperature_legends.{FigFormatPaper}'
# Now create a new image with legends only
# adjust the figure size as necessary
fig_leg = plt.figure(figsize=(1,1))
ax_leg = fig_leg.add_subplot(111)
# add the legend from the previous axes
ax_leg.legend(*ax.get_legend_handles_labels(), ncol=3, loc='center')
# hide the axes frame and the x/y labels
ax_leg.axis('off')
#ax_leg.set_frame_on(False)
lpltq3d.save_figure(save_file_name_, fig=fig_leg, savefig=savefigure, dpi=fig_save_dpi, show_plot=False)


# In[ ]:


for tb, mobility_df in mobility_dfff.items():
    ii=0
    df_group_T = mobility_df.groupby(['T'])
    fig, ax = plt.subplots(figsize=(9.5,6), constrained_layout=True)
    ax2 = ax.twinx()
    for name, group in df_group_T:
        #if name[0] in [100,200,300,400,600,800]:
        ax2.plot(group['AlContentChannel'], group['LFOM'], '-',marker=markers_[ii], color=colors_[ii], label=name[0], ms=12)
        ref_LFOM = mobility_ref_dff.loc[mobility_ref_dff['T']==name[0], 'LFOM'].iloc[0]
        ax.scatter(0.47,ref_LFOM, marker=markers_[ii], color=colors_[ii], s=100)
        ii+=1
    ax.set_yscale('log')
    ax.set_xlabel(x_p_label_text)
    ax.set_ylabel(z_label['LFOM'])
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    
    ax.axvline(x=0.5, c='k', ls='--')
    #ax2.scatter([0.48]*len(lfom_r),lfom_r, marker='s')
    ax2.set_yscale('log')
    if tb == 50: 
        yminn = 7e3
    elif tb==25:
        yminn = 2e3
    ax2.set_ylim(ymin=yminn, ymax=5e4)
    #ax2.yaxis.set_major_locator(ticker.LogLocator(numticks=2))
    #ax2.yaxis.set_minor_formatter(ticker.NullFormatter())
    #ax2.minorticks_off()
    
    tmp = ax.get_xticks()
    tmp_txt = [f'{xx:.1f}' for xx in tmp]
    ax.set_xticks(list(tmp)+[0.47], tmp_txt+['0'])
    ax.set_xlim(0.45,0.9)
    lpltq3d.save_figure(f'Highest_FOM_line_Tb{tb}.{FigFormatPaper}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


for tb, mobility_df in mobility_dfff.items():
    ii=0
    df_group_T = mobility_df.groupby(['T'])
    fig, ax = plt.subplots(figsize=(9,6), constrained_layout=True)
    ax2 = ax.twinx()
    for name, group in df_group_T:
        ax2.plot(group['AlContentChannel'], group['LFOM'], '-',marker=markers_[ii], color=colors_[ii], label=name[0], ms=12)
        ref_LFOM = mobility_ref_dff.loc[mobility_ref_dff['T']==name[0], 'LFOM'].iloc[0]
        ax.scatter(0.47,ref_LFOM, color=colors_[ii], marker=markers_[ii],s=100)
        ii+=1
    ax.set_yscale('log')
    ax.set_xlabel(x_p_label_text)
    ax.set_ylabel(z_label['LFOM'])
    
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=3)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    
    ax.axvline(x=0.5, c='k', ls='--')
    ax2.set_yscale('log')
    ax2.set_ylim(ax.get_ylim())
    tmp = ax.get_xticks()
    tmp_txt = [f'{xx:.1f}' for xx in tmp]
    ax.set_xticks(list(tmp)+[0.47], tmp_txt+['0'])
    ax.set_xlim(0.45,0.9)
    lpltq3d.save_figure(f'Highest_FOM_line_zoom_Tb{tb}.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


for tb, mobility_df in mobility_dfff.items():
    XX, YY, XXX, YYY = plots_data_maxms[tb]
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_xlabel('T (K)')
    ax.set_ylabel('Al composition channel, x') #Highest LFOM (MW/cm$^2$)
    ax.plot(YYY, XXX, 'ko-', ms=12)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.025))
    ax.set_xlim(0,800)
    lpltq3d.save_figure(f'Highest_FOM_T_comp_Tb{tb}.{FigFormatPaper}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:


for tb in mobility_dfff:
    comp_,temperature_,lfom_ = plots_data[tb]['comp_'], plots_data[tb]['temperature_'], plots_data[tb]['lfom_']
    XX, YY, XXX, YYY = plots_data_maxms[tb]
    vmin = np.nanmin(lfom_) #lfom_[lfom_>1].min()
    vmax = np.nanmax(lfom_)
    norm = colors.LogNorm(vmin=vmin, vmax=vmax) 
    cbar_mapable = cm.ScalarMappable(norm=norm, cmap='viridis')
    
    fig, ax = plt.subplots(figsize=(8,4.5), constrained_layout=True)
    
    fig, _ = lpltq3d.Plotq3D(comp_,temperature_,lfom_, fig=fig, ax=ax,
                             xmin=0.475, xmax=0.975, ymin=-40, ymax=850,
                             vmin=vmin, vmax=vmax,cbar_mappable=cbar_mapable,
                             x_label=x_p_label_text, y_label= 'T (K)',
                             tick_multiplicator=tick_multiplicator, norm='log',
                             show_contour_lines=False, marker='s', marker_size=38**2,
                             cbar_text=z_label['LFOM'],show_colorbar=True,
                             plot_controur=0, plot_scatter=1, show_plot=False)
    ax.scatter(XX, YY, marker='x', s=20**2, c='k')
    ax.scatter(XXX, YYY, marker='o', s=10**2, c='k')
    lpltq3d.save_figure(f'Highest_FOM_Tb{tb}.{FigFormat}', savefig=savefigure, show_plot=True,
                        fig=fig, CountFig=None, dpi=fig_save_dpi)


# In[ ]:




