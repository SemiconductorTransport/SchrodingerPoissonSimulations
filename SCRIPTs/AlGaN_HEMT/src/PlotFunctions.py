#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 08:56:43 2024

@author: badal.mondal
"""
import os
import numpy as np
import pandas as pd
import nextnanopy as nn
from nextnanopy.utils.misc import mkdir_if_not_exist
import matplotlib as mpl
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib.tri as tri 

## ============================================================================
params = {'figure.figsize': (8, 6), 'legend.fontsize': 18, 'axes.labelsize': 24, 'axes.titlesize': 24,
          'xtick.labelsize':24, 'xtick.major.width':2, 'xtick.major.size':5, 'ytick.labelsize': 24,
          'ytick.major.width':2, 'ytick.major.size':5, 'xtick.minor.width':2, 'xtick.minor.size':3,
          'ytick.minor.width':2, 'ytick.minor.size':3, 'errorbar.capsize':2}
plt.rcParams.update(params)
plt.rc('font', size=24)

## ============================================================================
class _general_plot_functions:
    def __init__(self):
        pass
    
    @classmethod
    def _save_figs(cls, fig, filename:str='test', figs_path='.', savefigure:bool=False, 
                   FigFormat='.png', FigDpi:int=75):
        if savefigure:
            filename_ = f'{filename}{FigFormat}'
            fig.savefig(os.path.join(figs_path, filename_), bbox_inches='tight', dpi=FigDpi) 
        return
    

class Plot1DFuns(_general_plot_functions):
    def __init__(self):
        pass

    def PlotBandEdges(self, ax, df_band_edge, density_list, show_right_yaxis:bool=False, xlabel:str="Distance", 
                      ylabel:str='Energy', yright_label:str='Electron concentration',
                      show_legend:bool=False, set_ymin:bool=True, xaxis_n_locator:int=6):
        XX = df_band_edge.coords['x'].value
        xunit = f'{df_band_edge.coords['x'].unit}'
        xlimit = [XX.min(), XX.max()]
        device_length = xlimit[1] - xlimit[0]
        ttm = device_length//xaxis_n_locator
        x_major_locator_division = ttm - (ttm%10)
        ax.plot(XX, df_band_edge.variables['Gamma_'].value, label='Gamma', color='k')
        ax.plot(XX, df_band_edge.variables['HH_'].value, label='heavy hole', color='y')
        ax.plot(XX, df_band_edge.variables['LH_'].value, label='light hole', color='tab:blue')
        ax.plot(XX, df_band_edge.variables['SO_'].value, label='crystal-field hole', color='g')
        ax.plot(XX, df_band_edge.variables['electron_Fermi_level_'].value, color='gray')
        
        yunit = f'{df_band_edge.variables['Gamma_'].unit}' # eV
        ax.set_ylabel(f"{ylabel} ({yunit})")
        ax.set_xlabel(f'{xlabel} ({xunit})')
        # xticks_location =[KK.value for KK in df_input_variables if KK.name.startswith('End')]
        # ax.set_xticks(xticks_location)
        # ax.xaxis.set_tick_params(rotation=45)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=x_major_locator_division))
        ax.xaxis.set_major_formatter('{x:.0f}')
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=x_major_locator_division//2))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=2))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=1))
        ax.set_xlim(xlimit)
        
        ax2 = ax.twinx()  
        ax2_unit = '10$^{{18}}$ cm$^{{-3}}$' #f'{df_e_density.variables['Electron_density'].unit}' # 10$^{{18}}$ cm$^{{-3}}$
        
            
        for key, color_, val in density_list:
            ax2.plot(val.coords['x'].value, val.variables[key].value , color=color_)
            
        if show_right_yaxis:
            ax2.set_ylabel(f'{yright_label} ({ax2_unit})', size=18) 
            ax2.tick_params(axis='y', labelcolor='r')
        else:
            ax2.set_yticks([])
        
        if set_ymin: ax2.set_ylim(ymin=-1)
            
        if show_legend: ax.legend()
        return ax, ax2
    
    def PlotDeviceSketch(self, ax0, df_input_variables, df_composition, 
                         show_doping:bool=False, show_Qregion:bool=False):
        #--------------------
        # Create the 2D grid for pcolormesh
        XX = df_composition.coords['x'].value
        X_grid, Y_grid = np.meshgrid(XX, [0, 1])
        # Create the 2D Z grid for pcolormesh
        ZZ = df_composition.variables['Alloy_x_'].value
        Z_grid = np.tile(ZZ, (2, 1))
        # Create the continuous plot using pcolormesh
        c = ax0.pcolormesh(X_grid, Y_grid, Z_grid, cmap='ocean',vmin=0, vmax=1.02)
        ax0.set_axis_off()
        
        # -------------------
        if show_doping and df_input_variables['Doping'].value:
            doping_lines = [KK.value for KK in df_input_variables if KK.name.endswith('DopingRegion')]
            ax0.fill_between(doping_lines, -0.49, 1.5, color='gray')
        #--------------------
        if show_Qregion:
            QR_lines = [KK.value for KK in df_input_variables if KK.name.startswith('QRegion_')]
            for pp in QR_lines:
                ax0.axvline(x=pp, linewidth=1, color='r', linestyle='--')
        
        rect1 = patches.Rectangle((df_input_variables['StartContact'].value, -0.5), 
                                  df_input_variables['ThicknessContact'].value, 2, 
                                  linewidth=1, edgecolor=None, linestyle='--',facecolor='yellow')
        rect2 = patches.Rectangle((df_input_variables['StartContact'].value, -0.5), 
                                  df_input_variables['EndDevice'].value-1, 2, 
                                  linewidth=1, edgecolor='k', linestyle='-',fill=False)
        ax0.add_patch(rect1)
        ax0.add_patch(rect2)
        ax0.set_ylim(-0.55, 1.5)
        return ax0
    
    def PlotBandDiagrams(self, data_folder_, output_figs, FigDpi:int=300, 
                         show_doping:bool=False, show_Qregion:bool=True, 
                         savefigure:bool=True, software_='nextnano++',
                         FigFormat:str='.png'):
        #==============================================================================
        print('- Output data folder:', data_folder_.fullpath)
        bias_folder = data_folder_.folders['bias_00000']
        strain_folder = data_folder_.folders['Strain']
        #quantum_data_folder = data_folder_.go_to('bias_00000', 'Quantum')
        structure_data_folder = data_folder_.folders['Structure']
        
        #==============================================================================
        input_variables_file = data_folder_.file('variables_input.txt')
        composition_file = structure_data_folder.file('alloy_composition.dat')
        band_edges_file = bias_folder.file('bandedges.dat')
        e_density_file = bias_folder.file('density_electron.dat')
        h_density_file = bias_folder.file('density_hole.dat')
        polarization_file = strain_folder.file('density_polarization_charge.dat')
        #piezoelectric_file = strain_folder.file('density_piezoelectric_charge.dat')
        #pyroelectric_file = strain_folder.file('density_pyroelectric_charge.dat')
        
        #==============================================================================
        df_input_variables = nn.DataFile(input_variables_file, product=software_)
        df_composition = nn.DataFile(composition_file, product=software_)
        df_band_edge = nn.DataFile(band_edges_file, product=software_)
        df_e_density = nn.DataFile(e_density_file, product=software_)
        df_h_density = nn.DataFile(h_density_file, product=software_)
        df_pol_density = nn.DataFile(polarization_file, product=software_)
        #df_piezo_density = nn.DataFile(piezoelectric_file, product=software_)
        #df_pyro_density = nn.DataFile(pyroelectric_file, product=software_)
        
        #==============================================================================
        mkdir_if_not_exist(output_figs)
    
        #%% ***************************************************************************
        ##### 3.1.2.1 Plot band edges
        ### ***************************************************************************
        fig, ax = plt.subplots(1)
        ax, ax2 = self.PlotBandEdges(ax, df_band_edge, density_list=[('Electron_density', 'r', df_e_density),
                                                                     ('Hole_density', 'b', df_h_density)])
        #fig.tight_layout()
        self._save_figs(fig, filename='band_edges', figs_path=output_figs, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
    
        #%% ***************************************************************************
        ##### 3.1.2.2 Plot band edges with polarization charge
        ### ***************************************************************************
        fig, ax = plt.subplots(1)
        # ax, ax2 = PlotBandEdges(ax, df_band_edge, density_list=[('Density', 'r', df_piezo_density),
        #                                                         ('Density', 'b', df_pyro_density)],
        #                         show_right_yaxis=True, yright_label='Charge density', set_ymin=False)
        ax, ax2 = self.PlotBandEdges(ax, df_band_edge, density_list=[('Density', 'r', df_pol_density)],
                                     show_right_yaxis=True, yright_label='Pol. charge density', set_ymin=False)
        ax.set_ylim(ymin=None)
        
        self._save_figs(fig, filename='band_edges_charges', figs_path=output_figs, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
    
        #%% ***************************************************************************
        ##### 3.1.2.3 Plot band edges with device structure
        ### ***************************************************************************
        fig, (ax0, ax) = plt.subplots(2,1, height_ratios=[1, 6], sharex=True, constrained_layout=True)
        ax0 = self.PlotDeviceSketch(ax0, df_input_variables, df_composition, 
                                    show_doping=show_doping, show_Qregion=show_Qregion)
        #--------------------
        ax, ax2 = self.PlotBandEdges(ax, df_band_edge, density_list=[('Electron_density', 'r', df_e_density), 
                                                                     ('Hole_density', 'b', df_h_density)], 
                                     show_right_yaxis=True)
        
        self._save_figs(fig, filename='band_edges_full', figs_path=output_figs, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return   
        
    def PlotSweepsData(self, XX, YY, fig=None, ax=None, x_label:str='', y_label:str='', x_log_scale:bool=False, 
                   tick_multiplicator:list=[None, None, None, None],
                   line_style='-', marker='.', color='r', FigDpi:int=75, FigFormat='.png',
                   figs_path='.', filename:str='test', savefigure:bool=False):
        if fig is None and ax is None:
            fig, ax = plt.subplots(1, constrained_layout=True)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.plot(XX, YY, ls=line_style, marker=marker, color=color)

        if all(tick_multiplicator[:2]):
            ax.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))
        
        #ax.set_xlim(xmin=45)
        #ax.set_ylim(ymin=0.835)

        if x_log_scale: ax.set_xscale('log')

        self._save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, ax