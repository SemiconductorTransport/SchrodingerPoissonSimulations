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
            mkdir_if_not_exist(figs_path)
            filename_ = f'{filename}{FigFormat}'
            fig.savefig(os.path.join(figs_path, filename_), bbox_inches='tight', dpi=FigDpi) 
            #plt.close()
        return
    

class Plot1DFuns(_general_plot_functions):
    def __init__(self):
        pass

    def _plot_band_edges(self, ax, df_band_edge, xlabel:str="Distance", 
                         ylabel:str='Energy', xaxis_n_locator:int=6):
        XX = df_band_edge.coords['x'].value
        xunit = f"{df_band_edge.coords['x'].unit}"
        xlimit = [XX.min(), XX.max()]
        device_length = xlimit[1] - xlimit[0]
        ttm = device_length//xaxis_n_locator
        x_major_locator_division = ttm - (ttm%10)
        #======================================================================
        ax.plot(XX, df_band_edge.variables['Gamma_'].value, label='Gamma', color='k')
        ax.plot(XX, df_band_edge.variables['HH_'].value, label='heavy hole', color='y')
        ax.plot(XX, df_band_edge.variables['LH_'].value, label='light hole', color='tab:blue')
        ax.plot(XX, df_band_edge.variables['SO_'].value, label='crystal-field hole', color='g')
        ax.plot(XX, df_band_edge.variables['electron_Fermi_level_'].value, color='gray')
        #======================================================================
        yunit = f"{df_band_edge.variables['Gamma_'].unit}" # eV
        ax.set_ylabel(f"{ylabel} ({yunit})")
        ax.set_xlabel(f'{xlabel} ({xunit})')
        # xticks_location =[KK.value for KK in df_input_variables if KK.name.startswith('End')]
        # ax.set_xticks(xticks_location)
        # ax.xaxis.set_tick_params(rotation=45)
        #======================================================================
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=x_major_locator_division))
        ax.xaxis.set_major_formatter('{x:.0f}')
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=x_major_locator_division//2))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=2))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=1))
        #======================================================================
        ax.set_xlim(xlimit)
        return ax
        
    def _plot_right_twin_band_edges(self, ax, density_list, y_label:str='Electron concentration', 
                                    show_twin_yaxis_labels:bool=False, set_ymin:bool=True):
        ax2 = ax.twinx()  
        ax2_unit = '10$^{{18}}$ cm$^{{-3}}$' #f'{df_e_density.variables['Electron_density'].unit}' # 10$^{{18}}$ cm$^{{-3}}$
        #======================================================================
        for key, color_, val in density_list:
            ax2.plot(val.coords['x'].value, val.variables[key].value , color=color_)
        #======================================================================
        if show_twin_yaxis_labels:
            ax2.set_ylabel(f'{y_label} ({ax2_unit})', size=18) 
            ax2.tick_params(axis='y', labelcolor='r')
        else:
            ax2.set_yticks([])
        #======================================================================
        if set_ymin: ax2.set_ylim(ymin=-1)
        return ax, ax2
    
    def _plot_device_sketch(self, ax0, df_input_variables, df_composition, 
                            show_doping:bool=False, show_Qregion:bool=False,
                            device_cmap='ocean'):
        #======================================================================
        # Create the 2D grid for pcolormesh
        XX = df_composition.coords['x'].value
        X_grid, Y_grid = np.meshgrid(XX, [0, 1])
        # Create the 2D Z grid for pcolormesh
        ZZ = df_composition.variables['Alloy_x_'].value
        Z_grid = np.tile(ZZ, (2, 1))
        # Create the continuous plot using pcolormesh
        ax0.pcolormesh(X_grid, Y_grid, Z_grid, cmap=device_cmap,vmin=0, vmax=1.02)
        ax0.set_axis_off()
        
        #======================================================================
        if show_doping and df_input_variables['Doping'].value:
            doping_lines = [KK.value for KK in df_input_variables if KK.name.endswith('DopingRegion')]
            ax0.fill_between(doping_lines, -0.49, 1.5, color='gray')
        #======================================================================
        if show_Qregion:
            QR_lines = [KK.value for KK in df_input_variables if KK.name.startswith('QRegion_')]
            for pp in QR_lines:
                ax0.axvline(x=pp, linewidth=1, color='r', linestyle='--')
        #======================================================================
        rect1 = patches.Rectangle((df_input_variables['StartContact'].value, -0.5), 
                                  df_input_variables['ThicknessContact'].value, 2, 
                                  linewidth=1, edgecolor=None, linestyle='--',facecolor='yellow')
        rect2 = patches.Rectangle((df_input_variables['StartContact'].value, -0.5), 
                                  df_input_variables['EndDevice'].value-1, 2, 
                                  linewidth=1, edgecolor='k', linestyle='-',fill=False)
        #======================================================================
        ax0.add_patch(rect1)
        ax0.add_patch(rect2)
        ax0.set_ylim(-0.55, 1.5)
        return ax0
    
    def PlotBandDiagrams(self, data_folder_, show_doping:bool=False, show_Qregion:bool=False, 
                         figs_path='.', FigDpi:int=300, filename:str='test', device_cmap='ocean',
                         savefigure:bool=False, software_='nextnano++',
                         FigFormat:str='.png', plot_device_sketch:bool=False,
                         plot_pol_charge:bool=False, plot_eh_density:bool=False,
                         show_twin_yaxis_labels:bool=False, show_legend:bool=False):
        # default is band edges only. extra plots will be drawn on the twin axis
        # plot_eh_density: plot electron-hole density
        # plot_pol_charge: plot polarization change density
        #======================================================================
        bias_folder = data_folder_.folders['bias_00000']
        strain_folder = data_folder_.folders['Strain']
        #quantum_data_folder = data_folder_.go_to('bias_00000', 'Quantum')
        structure_data_folder = data_folder_.folders['Structure']
        #======================================================================
        ax, ax0, ax2 = None, None, None
        if plot_device_sketch:
            fig, (ax0, ax) = plt.subplots(2,1, height_ratios=[1, 6], sharex=True, constrained_layout=True)
        else:
            fig, ax = plt.subplots(1)
        
        # *********************************************************************
        ##### Plot band edges
        # *********************************************************************
        band_edges_file = bias_folder.file('bandedges.dat')
        df_band_edge = nn.DataFile(band_edges_file, product=software_)
        ax = self._plot_band_edges(ax, df_band_edge)
        
        # *********************************************************************
        ##### Plot electron and hole density
        # *********************************************************************
        if plot_eh_density:
            #==================================================================
            e_density_file = bias_folder.file('density_electron.dat')
            h_density_file = bias_folder.file('density_hole.dat')
            #==================================================================
            df_e_density = nn.DataFile(e_density_file, product=software_)
            df_h_density = nn.DataFile(h_density_file, product=software_)
            #==================================================================
            ax, ax2 = self._plot_right_twin_band_edges(ax, density_list=[('Electron_density', 'r', df_e_density),
                                                                         ('Hole_density', 'b', df_h_density)], 
                                                       y_label='Electron concentration', set_ymin=True,
                                                       show_twin_yaxis_labels=show_twin_yaxis_labels)
            #==================================================================
    
        # *********************************************************************
        ##### Plot band edges with polarization charge
        # *********************************************************************
        if plot_pol_charge:
            #==================================================================
            polarization_file = strain_folder.file('density_polarization_charge.dat')
            #piezoelectric_file = strain_folder.file('density_piezoelectric_charge.dat')
            #pyroelectric_file = strain_folder.file('density_pyroelectric_charge.dat')
            #==================================================================
            df_pol_density = nn.DataFile(polarization_file, product=software_)
            #df_piezo_density = nn.DataFile(piezoelectric_file, product=software_)
            #df_pyro_density = nn.DataFile(pyroelectric_file, product=software_)
            #==================================================================
            # ax, ax2 = PlotRightTwinBandEdges(ax, density_list=[('Density', 'r', df_piezo_density),
            #                                                         ('Density', 'b', df_pyro_density)],
            #                         show_twin_yaxis_labels=True, y_label='Charge density', set_ymin=False)
            ax, ax2 = self._plot_right_twin_band_edges(ax, density_list=[('Density', 'r', df_pol_density)], 
                                                       y_label='Pol. charge density', set_ymin=False,
                                                       show_twin_yaxis_labels=show_twin_yaxis_labels)
            ax.set_ylim(ymin=None)
            #==================================================================
    
        # *********************************************************************
        ##### Plot band edges with device structure
        # *********************************************************************
        if plot_device_sketch:
            #==================================================================
            input_variables_file = data_folder_.file('variables_input.txt')
            composition_file = structure_data_folder.file('alloy_composition.dat')
            #==================================================================
            df_input_variables = nn.DataFile(input_variables_file, product=software_)
            df_composition = nn.DataFile(composition_file, product=software_)
            #==================================================================
            ax0 = self._plot_device_sketch(ax0, df_input_variables, df_composition, device_cmap=device_cmap,
                                           show_doping=show_doping, show_Qregion=show_Qregion)
            #==================================================================
        
        if show_legend: ax.legend()
        
        self._save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, ax, ax0, ax2
        
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
    
    
class PlotQuasi3DFuns(_general_plot_functions):
    def __init__(self):
        pass    
    
    @classmethod
    def _triangulation(cls, x_values, y_values):
        return tri.Triangulation(x_values, y_values)
    
    @classmethod
    def _meshgrid(cls, x_values, y_values, npoints:int=20): 
        xi, yi = np.meshgrid(np.linspace(x_values.min(), x_values.max(), npoints), 
                             np.linspace(y_values.min(), y_values.max(), npoints))
        return xi, yi
    
    @classmethod                    
    def _linear_interpolation(cls, triangles_, x_values, y_values, z_values):
        interp_lin = tri.LinearTriInterpolator(triangles_, z_values)
        return interp_lin(x_values, y_values)
    
    def InterPolation(self, x_values, y_values, z_values, method:str='linear',
                      interpolation_points:int=20):
        assert method in ['linear'], 'Requested interpolation method not implemented yet.'
        triang = self._triangulation(x_values, y_values)
        xi, yi = self._meshgrid(x_values, y_values, npoints=interpolation_points)
        zi = None
        if method == 'linear':
            zi = self._linear_interpolation(triang, xi, yi, z_values)
        return xi, yi, zi
    
    def PlotContour(self, xi, yi, zi, vmin=0, vmax=3.2, x_label:str='', y_label:str='',
                    title_label:str=None, z_label:str='', 
                    tick_multiplicator:list=[None, None, None, None],
                    FigDpi:int=75, FigFormat='.png', figs_path='.', 
                    filename:str='test', savefigure:bool=False):
        # Set up the figure
        fig, axs = plt.subplots(constrained_layout=True)

        # Plot linear interpolation to quad grid.
        pcm0 = axs.contourf(xi, yi, zi, vmin = 0, vmax = 3.2)
        #axs.contour(cs, colors='k', linewidths=1)
        # Plot the triangulation.
        #cs = axs.tricontourf(triang, ZZ, vmin = vmin, vmax = vmax)
        #axs.tricontour(cs, colors='k', linewidths=1)
        axs.set_ylabel(y_label)
        axs.set_xlabel(x_label)
        if title_label is not None: axs.set_title(title_label)
        
        print(pcm0)
        cbar = fig.colorbar(pcm0, ax=axs, location='right', label=z_label,
                           norm=mpl.colors.Normalize(vmin=0, vmax=3.2), boundaries=np.linspace(0, 3.2, 6))
        #cbar.set_clim(0, 3.2)
        #cbar =fig.colorbar(m, boundaries=np.linspace(0, 3., 6))
        if all(tick_multiplicator[:2]):
            axs.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            axs.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            axs.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            axs.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))
            
        self._save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, axs