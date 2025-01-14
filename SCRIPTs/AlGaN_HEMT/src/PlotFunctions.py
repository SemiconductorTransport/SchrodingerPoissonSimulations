#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 08:56:43 2024

@author: badal.mondal
"""
import os
import numpy as np
import nextnanopy as nn
from nextnanopy.utils.misc import mkdir_if_not_exist
from matplotlib import colors, cm
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
class general_plot_functions:
    def __init__(self):
        pass
    
    @classmethod
    def save_figs(cls, fig, filename:str='test', figs_path='.', savefigure:bool=False, 
                   FigFormat='.png', FigDpi:int=75):
        if savefigure and fig is not None:
            mkdir_if_not_exist(figs_path)
            filename_ = f'{filename}{FigFormat}'
            fig.savefig(os.path.join(figs_path, filename_), bbox_inches='tight', dpi=FigDpi) 
            plt.close()
        return
    
    def _align_yaxis_np(self, axes):
        """
        Align zeros of the two axes, zooming them out by same ratio
        Source: https://stackoverflow.com/a/46901839
        """
        axes = np.array(axes)
        extrema = np.array([ax.get_ylim() for ax in axes])
    
        # reset for divide by zero issues
        for i in range(len(extrema)):
            if np.isclose(extrema[i, 0], 0.0):
                extrema[i, 0] = -1
            if np.isclose(extrema[i, 1], 0.0):
                extrema[i, 1] = 1
    
        # upper and lower limits
        lowers = extrema[:, 0]
        uppers = extrema[:, 1]
        #print(lowers, uppers)
        # if all pos or all neg, don't scale
        all_positive = False
        all_negative = False
        if lowers.min() > 0.0:
            all_positive = True
    
        if uppers.max() < 0.0:
            all_negative = True
    
        if all_negative or all_positive:
            # don't scale
            return
    
        # pick "most centered" axis
        res = abs(uppers+lowers)
        min_index = np.argmin(res)
    
        # scale positive or negative part
        multiplier1 = abs(uppers[min_index]/lowers[min_index])
        multiplier2 = abs(lowers[min_index]/uppers[min_index])
    
        for i in range(len(extrema)):
            # scale positive or negative part based on which induces valid
            if i != min_index:
                lower_change = extrema[i, 1] * -1*multiplier2
                upper_change = extrema[i, 0] * -1*multiplier1
                if upper_change < extrema[i, 1]:
                    extrema[i, 0] = lower_change
                else:
                    extrema[i, 1] = upper_change
    
            # bump by 10% for a margin
            extrema[i, 0] *= 1.1
            extrema[i, 1] *= 1.1
    
        # set axes limits
        [axes[i].set_ylim(*extrema[i]) for i in range(len(extrema))]

    @classmethod
    def _set_colorbar(cls, axs, fig, CS=None, cbar_mappable=None, cbar_text:str=None):
        if cbar_mappable is None:
                cbar = fig.colorbar(CS, ax=axs)
        else:
            cbar = fig.colorbar(cbar_mappable, ax=axs)         
        cbar.ax.set_ylabel(cbar_text)
        return cbar

    @classmethod
    def _set_labels(cls, axs, x_label:str='', y_label:str='', title_label:str=None):
        axs.set_ylabel(y_label)
        axs.set_xlabel(x_label)
        if title_label is not None: axs.set_title(title_label)

    @classmethod
    def _set_tickers(cls, axs, tick_multiplicator=[None]*4):
        if all(tick_multiplicator[:2]):
            axs.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            axs.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            axs.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            axs.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))

class Plot1DFuns(general_plot_functions):
    def __init__(self):
        pass
    
    @staticmethod
    def _bands_characterstics_plot(bands_characterstics=None):
        # band_name: (band_label_name_to_pu_in_plot_label, color_to_use)
        defau_characterstics = {'Gamma_':('Gamma', 'k'),
                                'HH_':('heavy hole', 'y'),
                                'LH_':('light hole', 'tab:blue'),
                                'SO_':('crystal-field hole', 'g'),
                                'electron_Fermi_level_':('Fermi level', 'gray')
                                }
        if bands_characterstics is None:
            return defau_characterstics
        else:
            tmp = defau_characterstics.copy()
            for chrt in defau_characterstics:
                if bands_characterstics.get(chrt):
                    tmp[chrt] = bands_characterstics.get(chrt)
            return tmp

    def _plot_band_edges(self, ax, df_band_edge, xlabel:str="Distance", 
                         ylabel:str='Energy', xaxis_n_locator:int=6,
                         scale_x_axis:float = None,
                         y_major_locator_distance:float=2,
                         x_zoom_region:list=[None, None],
                         plot_bands:list=['Gamma_','HH_','LH_','SO_',
                                          'electron_Fermi_level_'],
                         bands_characters:dict = None, band_edge_ls='-'):
        if plot_bands is None or len(plot_bands) < 1:
            print('Skipping band edge plot. Nothing to plot.')
            return ax
        XX = df_band_edge.coords['x'].value
        if scale_x_axis:
            XX -= scale_x_axis
        xunit = f"{df_band_edge.coords['x'].unit}"
        yunit = f"{df_band_edge.variables['Gamma_'].unit}" # eV

        xlimit = [XX.min(), XX.max()]
        xlimit = [x_zoom_ if x_zoom_ is not None else xlimit[ii] for ii, x_zoom_ in enumerate(x_zoom_region)]
        pos_ = [np.argmax(XX > xlimit[0]), np.argwhere(XX< xlimit[1])[-1][0]]
               
        bands_characters = self._bands_characterstics_plot(bands_characterstics=bands_characters)
        #======================================================================
        for band_plot in plot_bands:
            YY = df_band_edge.variables[band_plot].value
            YY_ = YY[pos_[0]:pos_[1]]
            ax.plot(XX[pos_[0]:pos_[1]], YY_, linestyle=band_edge_ls,
                    label=bands_characters.get(band_plot)[0],
                    color=bands_characters.get(band_plot)[1])
        #======================================================================
        ax.set_ylabel(f"{ylabel} ({yunit})")
        ax.set_xlabel(f'{xlabel} ({xunit})')
        # xticks_location =[KK.value for KK in df_input_variables if KK.name.startswith('End')]
        # ax.set_xticks(xticks_location)
        # ax.xaxis.set_tick_params(rotation=45)
        #======================================================================
        if xaxis_n_locator is not None:
            ttm = (xlimit[1] - xlimit[0])//xaxis_n_locator
            x_major_locator_division = ttm - (ttm%10)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(base=x_major_locator_division))
            ax.xaxis.set_major_formatter('{x:.0f}')
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=x_major_locator_division/2))
            
        if y_major_locator_distance:
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base=y_major_locator_distance))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=y_major_locator_distance/2))
        #======================================================================
        ax.set_xlim(xlimit)
        return ax
        
    def _plot_density_plots(self, ax2, density_list, ax2_unit:str=None,
                            label_size=18, scale_x_axis:float=None,
                            y_label:str='Carrier concentration', line_alpha=1.0,
                            show_twin_yaxis_labels:bool=False, set_ymin:float=None):
        if ax2_unit is None:
            ax2_unit = '10$^{{18}}$ cm$^{{-3}}$' #f'{df_e_density.variables['Electron_density'].unit}' # 10$^{{18}}$ cm$^{{-3}}$
        #======================================================================
        for key, color_, val in density_list:
            x_vals_plt = val.coords['x'].value
            if scale_x_axis:
                x_vals_plt -= scale_x_axis
            ax2.plot(x_vals_plt, val.variables[key].value , color=color_, alpha=line_alpha)
        #======================================================================
        if show_twin_yaxis_labels:
            if y_label and label_size: 
                ax2.set_ylabel(f'{y_label} ({ax2_unit})', size=label_size) 
            ax2.tick_params(axis='y', labelcolor='k')
        else:
            ax2.set_yticks([])
        #======================================================================
        if set_ymin: ax2.set_ylim(ymin=set_ymin)
        return ax2
    
    def _plot_device_sketch(self, ax0, df_input_variables, df_composition, 
                            show_doping:bool=False, show_Qregion:bool=False,
                            device_cmap='ocean',scale_x_axis:float=None):
        #======================================================================
        # Create the 2D grid for pcolormesh
        XX = df_composition.coords['x'].value
        if scale_x_axis: XX -= scale_x_axis
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
                                  df_input_variables['EndDevice'].value, 2, 
                                  linewidth=1, edgecolor='k', linestyle='-',fill=False)
        #======================================================================
        ax0.add_patch(rect1)
        ax0.add_patch(rect2)
        ax0.set_ylim(-0.55, 1.5)
        ax0.set_xlim(xmax=df_input_variables['EndDevice'].value+1)
        return ax0
    
    def PlotBandDiagrams(self, data_folder_, fig=None, ax=None, ax0=None, ax2=None,
                         show_doping:bool=False, show_Qregion:bool=False, line_alpha=1.0,
                         xlabel:str="Distance", ylabel:str='Energy', ylabel_twin:str='Carrier concentration',
                         figs_path='.', FigDpi:int=300, filename:str='test', device_cmap='ocean',
                         xaxis_n_locator:int=6, y_left_major_locator_distance:float=2,
                         scale_x_axis:float = None, 
                         savefigure:bool=False, software_='nextnano++',
                         FigFormat:str='.png', plot_device_sketch:bool=False,
                         density_list=[('Electron_density', 'r'), ('Hole_density', 'b')],
                         band_file:str=None, band_index:int=0, kindex:int=0, 
                         subband_energy_level:bool=True, plot_density:bool=False, 
                         show_twin_yaxis_labels:bool=False, twin_yaxis_unit:str=None,
                         right_yaxis_shift:float=-1, show_legend:bool=False,
                         x_zoom_region:list=[None, None], x_zoom_2nd_no_shift:bool=False,
                         bands_characters:dict = None, band_edge_ls='-', 
                         plot_density_on_left_axis:bool=False,
                         plot_bands:list=['Gamma_','HH_','LH_','SO_', 'electron_Fermi_level_'], 
                         align_left_right_yaxis:bool=False):
        # default is band edges only. extra plots will be drawn on the twin axis
        # plot_eh_density: plot electron-hole density
        # plot_pol_charge: plot polarization change density
        #======================================================================
        bias_folder = data_folder_.folders['bias_00000']
        strain_folder = data_folder_.folders['Strain']
        quantum_data_folder = data_folder_.go_to('bias_00000', 'Quantum')
        structure_data_folder = data_folder_.folders['Structure']
        #======================================================================
        if any([ax, ax0, ax2]):
            pass
        else:
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
        
        if any(x_zoom_region): 
            input_variables_file = data_folder_.file('variables_input.txt')
            df_input_variables = nn.DataFile(input_variables_file, product=software_)
            x_zoom = [df_input_variables.variables[jj].value if isinstance(jj, str) else jj for jj in x_zoom_region]
            if x_zoom_2nd_no_shift:
                x_zoom = [x_zoom[0]-x_zoom[1], x_zoom[0]+x_zoom[1]]
        else:
            x_zoom = x_zoom_region[:]
        if scale_x_axis: 
                x_zoom = [x_zoom_-scale_x_axis for x_zoom_ in x_zoom]
        
        ax = self._plot_band_edges(ax, df_band_edge, xlabel=xlabel, ylabel=ylabel, 
                                   xaxis_n_locator=xaxis_n_locator, scale_x_axis=scale_x_axis, 
                                   y_major_locator_distance=y_left_major_locator_distance,
                                   x_zoom_region=x_zoom, plot_bands=plot_bands,
                                   bands_characters=bands_characters,
                                   band_edge_ls=band_edge_ls)
        
        # *********************************************************************
        ##### Plot electron and hole density
        # *********************************************************************
        if plot_density:
            new_density_list = []
            for ddensity_details in density_list:
                if ddensity_details[0] == 'Electron_density':
                    density_file = bias_folder.file('density_electron.dat')
                    col_name = ['Electron_density']
                elif ddensity_details[0] == 'Hole_density':
                    density_file = bias_folder.file('density_hole.dat')
                    col_name = ['Hole_density']
                elif ddensity_details[0] == 'Potential':
                    density_file = bias_folder.file('potential.dat')
                    col_name = ['Potential']
                elif ddensity_details[0] == 'PsiSqare':
                    assert band_file is not None, 'Provide the which band file to plot.'
                    psi_sqr_data_folder = quantum_data_folder.go_to(band_file, f'kIndex_{kindex:05d}')
                    if 'kp' in band_file:
                        tf_name = f'probability_shift_{band_file}_{kindex:05d}_{band_index:04d}.dat'
                    else:
                        tf_name = f'probability_shift_{band_file}_{band_index:04d}.dat'
                    density_file = psi_sqr_data_folder.file(tf_name)
                    col_name = [f'Psi^2_{band_index}_']
                    if subband_energy_level: col_name.append(f'E_{band_index}_')
                elif ddensity_details[0] == 'Polarization_density':
                    #density_file = strain_folder.file('density_polarization_charge.dat')
                    density_file = strain_folder.file('polarization_charge_density_total.dat')
                    col_name = ['Density']
                elif  ddensity_details[0] == 'Pizoelectric_density':
                    density_file = strain_folder.file('density_piezoelectric_charge.dat')
                    col_name = ['Density']
                elif ddensity_details[0] == 'Pyroelectric_density':
                    density_file = strain_folder.file('density_pyroelectric_charge.dat')
                    col_name = ['Density']
                else:
                    raise ValueError('Requested density plot is not implemented yet.')
                
                for col_nam in col_name:
                    new_density_list.append((col_nam, ddensity_details[1], 
                                             nn.DataFile(density_file, product=software_)))
            #==================================================================
            if plot_density_on_left_axis:
                ax = self._plot_density_plots(ax, density_list=new_density_list,
                                              scale_x_axis=scale_x_axis, 
                                              y_label=ylabel_twin, set_ymin=right_yaxis_shift,
                                              show_twin_yaxis_labels=show_twin_yaxis_labels,
                                              ax2_unit=twin_yaxis_unit, line_alpha=line_alpha)
            else:
                if ax2 is None: ax2 = ax.twinx() 
                ax2 = self._plot_density_plots(ax2, density_list=new_density_list,line_alpha=line_alpha,
                                               scale_x_axis=scale_x_axis, 
                                               y_label=ylabel_twin, set_ymin=right_yaxis_shift,
                                               show_twin_yaxis_labels=show_twin_yaxis_labels)
                
                if align_left_right_yaxis: 
                    self._align_yaxis_np([ax, ax2])
    
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
            ax0 = self._plot_device_sketch(ax0, df_input_variables, df_composition, 
                                           scale_x_axis=scale_x_axis, device_cmap=device_cmap,
                                           show_doping=show_doping, show_Qregion=show_Qregion)
            #==================================================================
        
        if show_legend: ax.legend()
        
        
        self.save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, ax, ax0, ax2
        
    def PlotSweepsData(self, XX, YY, fig=None, ax=None, x_label:str='', y_label:str='', x_log_scale:bool=False, 
                       tick_multiplicator:list=[None, None, None, None], title_label:str=None,
                       line_style='-', marker='o', color='r', FigDpi:int=75, FigFormat='.png',
                       figs_path='.', filename:str='test', savefigure:bool=False):
        if fig is None and ax is None:
            fig, ax = plt.subplots(1, constrained_layout=True)
        ax.plot(XX, YY, ls=line_style, marker=marker, color=color, ms=12)
        self._set_labels(ax, x_label=x_label, y_label=y_label, title_label=title_label)
        self._set_tickers(ax, tick_multiplicator=tick_multiplicator)

        if x_log_scale: ax.set_xscale('log')

        self.save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, ax
    
    
class PlotQuasi3DFuns(general_plot_functions):
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
    
    def CreateColorbarMapableObject(self, vmin=None, vmax=None, color_map='viridis'):
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        mappable = cm.ScalarMappable(norm=norm, cmap=color_map)
        return norm, mappable
    
    def PlotContour(self, xi, yi, zi, fig=None, axs=None,
                    x_label:str='', y_label:str='',
                    title_label:str=None, z_label:str='', 
                    tick_multiplicator:list=[None, None, None, None],
                    FigDpi:int=75, FigFormat='.png', figs_path='.', 
                    vmin=None, vmax=None, cbar_mappable=None, norm=None,
                    color_map='viridis', show_contour_lines:bool=False, 
                    cbar_text:str=None,show_colorbar:bool=False,
                    filename:str='test', savefigure:bool=False):

        self.fig, axs = self._set_figure(fig=fig, axs=axs)

        CS = axs.contourf(xi, yi, zi, cmap=color_map, norm=norm)
        
        if show_contour_lines: 
            CS2 = axs.contour(CS, levels=CS.levels, colors='k')
        if show_colorbar:
            cbar=self._set_colorbar(axs, self.fig, CS=CS, cbar_mappable=cbar_mappable, 
                                    cbar_text=cbar_text)        
        self._set_labels(axs, x_label=x_label, y_label=y_label, title_label=title_label)
        self._set_tickers(axs, tick_multiplicator=tick_multiplicator)
        self.save_figs(self.fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return self.fig, axs

    def _set_figure(self, fig=None, axs=None):
        if axs is None: 
            self.fig, axs = plt.subplots(constrained_layout=True)
        else:
            self.fig = fig 
        return self.fig, axs
