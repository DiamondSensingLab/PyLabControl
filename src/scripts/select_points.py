
# This file is part of PyLabControl, software for laboratory equipment control for scientific experiments.
# Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell
#
#
# PyLabControl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyLabControl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyLabControl.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
import pandas as pd
import trackpy as tp
import scipy.spatial
import time
from matplotlib import patches

from PyLabControl.src.core import Script, Parameter

class SelectPoints(Script):
    """
Script to select points on an image. The selected points are saved and can be used in a superscript to iterate over.
    """
    _DEFAULT_SETTINGS = [
        Parameter('patch_size', 0.01),
        Parameter('nv_selection_method', 'manual',['manual', 'automatic', 'loaded'],'Identify NVs by mouse clicks or by automatic trackpy routine'),
        Parameter('nv_size_est', 5, int, 'estimate of NV size [pixels]'),
        Parameter('nv_mass_min', 0, float, 'maximum NV integrated brightness expected [pixels]'),
        Parameter('coords_file', 'C:\Users\sensing\Documents\RAW_DATA\___.csv', str, 'file to load NV coordinates from for "loaded" selection'),
        Parameter('coord_offset_x', 0, float, 'x offset btw current glavo image and coordinate list'),
        Parameter('coord_offset_y', 0, float, 'y offset btw current glavo image and coordinate list')
    ]

    _INSTRUMENTS = {}
    _SCRIPTS = {}

    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Select points by clicking on an image
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

        self.patches = []
        self.plot_settings = {}

    def _function(self):
        """
        Waits until stopped to keep script live. Gui must handle calling of Toggle_NV function on mouse click.
        """

        self.data = {'nv_locations': [], 'nv_masses': [], 'nv_sizes':[], 'image_data': None, 'extent': None}

        self.progress = 50
        self.updateProgress.emit(self.progress)
        # keep script alive while NVs are selected
        while not self._abort:
            time.sleep(1)


    def plot(self, figure_list):
        '''
        Plots a dot on top of each selected NV, with a corresponding number denoting the order in which the NVs are
        listed.
        Precondition: must have an existing image in figure_list[0] to plot over
        Args:
            figure_list:
        '''
        # if there is not image data get it from the current plot
        print('here we are')
        if not self.data == {} and self.data['image_data'] is  None:
            axes = figure_list[0].axes[0]
            if len(axes.images)>0:
                self.data['image_data'] = np.array(axes.images[0].get_array())
                self.data['extent'] = np.array(axes.images[0].get_extent())
                self.plot_settings['cmap'] = axes.images[0].get_cmap().name
                self.plot_settings['xlabel'] = axes.get_xlabel()
                self.plot_settings['ylabel'] = axes.get_ylabel()
                self.plot_settings['title'] = axes.get_title()
                self.plot_settings['interpol'] = axes.images[0].get_interpolation()

        if self.data['nv_locations']==[]:
            if self.settings['nv_selection_method']=='automatic':
                nvsauto = tp.locate(self.data['image_data'],self.settings['nv_size_est'], minmass = self.settings['nv_mass_min'], invert=False)
                nvcoords = [self.pixel_to_voltage(p, self.data['extent'], np.shape(self.data['image_data'])) for p in nvsauto[['x', 'y']].as_matrix()]
                self.data['nv_locations'] = nvcoords
                # self.data['nv_masses'] = nvsauto['mass']
                # self.data['nv_sizes'] = nvcoords['size']
            elif self.settings['nv_selection_method']=='loaded':
                coordstmp = np.loadtxt(self.settings['coords_file'], delimiter=',')
                coordstmp[:,0] = coordstmp[:,0] + self.settings['coord_offset_x']
                coordstmp[:,1] = coordstmp[:,1] + self.settings['coord_offset_y']
                coordstmp = coordstmp.tolist()
                coordstmp.pop(0)
                self.data['nv_locations'] = coordstmp


        Script.plot(self, figure_list)

    #must be passed figure with galvo plot on first axis
    def _plot(self, axes_list):
        '''
        Plots a dot on top of each selected NV, with a corresponding number denoting the order in which the NVs are
        listed.
        Precondition: must have an existing image in figure_list[0] to plot over
        Args:
            figure_list:
        '''

        axes = axes_list[0]

        if self.plot_settings:
            axes.imshow(self.data['image_data'], cmap=self.plot_settings['cmap'], interpolation=self.plot_settings['interpol'], extent=self.data['extent'])
            axes.set_xlabel(self.plot_settings['xlabel'])
            axes.set_ylabel(self.plot_settings['ylabel'])
            axes.set_title(self.plot_settings['title'])

        self._update(axes_list)

    def _update(self, axes_list):

        axes = axes_list[0]

        patch_size = self.settings['patch_size']

        #first clear all old patches (circles and numbers), then redraw all
        if not self.patches == []:
            try: #catch case where plot has been cleared, so old patches no longer exist. Then skip clearing step.
                for patch in self.patches:
                    patch.remove()
            except ValueError:
                pass

        self.patches = []

        for index, pt in enumerate(self.data['nv_locations']):
            # axes.plot(pt, fc='b')

            circ = patches.Circle((pt[0], pt[1]), patch_size, fc='b')
            axes.add_patch(circ)
            self.patches.append(circ)

            text = axes.text(pt[0], pt[1], '{:d}'.format(index),
                    horizontalalignment='center',
                    verticalalignment='center',
                    color='white'
                    )
            self.patches.append(text)

        # histaxes = axes_list[1]
        # # histaxes.subplot(121)
        # histaxes.hold(False)
        # histaxes.hist(self.data['nv_masses'],bins=20)
        # histaxes.subplot(122)
        # histaxes.hold(False)
        # histaxes.hist(self.data['nv_sizes'],bins=20)

    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        the default creates a single axes object on each figure
        This can/should be overwritten in a child script if more axes objects are needed
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """

        # only pick the first figure from the figure list, this avoids that get_axes_layout clears all the figures
        return super(SelectPoints, self).get_axes_layout([figure_list[0]])

    def toggle_NV(self, pt):
        '''
        If there is not currently a selected NV within self.settings[patch_size] of pt, adds it to the selected list. If
        there is, removes that point from the selected list.
        Args:
            pt: the point to add or remove from the selected list

        Poststate: updates selected list

        '''
        if self.data['nv_locations']==[]: #if self.data is empty so this is the first point
            self.data['nv_locations'].append(pt)
            self.data['image_data'] = None # clear image data

        else:
            # use KDTree to find NV closest to mouse click
            tree = scipy.spatial.KDTree(self.data['nv_locations'])
            #does a search with k=1, that is a search for the nearest neighbor, within distance_upper_bound
            d, i = tree.query(pt,k = 1, distance_upper_bound = self.settings['patch_size'])

            # removes NV if previously selected
            if d is not np.inf:
                self.data['nv_locations'].pop(i)
                # if not self.data['nv_masses']==[]:
                #     self.data['nv_masses'].pop(i)
                # if not self.data['nv_sizes']==[]:
                #     self.data['nv_sizes'].pop(i)
            # adds NV if not previously selected
            else:
                self.data['nv_locations'].append(pt)

            print('chosen [x,y]:', pt)

    @staticmethod
    def pixel_to_voltage(pt, extent, image_dimensions):
        """"
        pt: point in pixels
        extent: [xVmin, Vmax, Vmax, yVmin] in volts
        image_dimensions: dimensions of image in pixels

        Returns: point in volts
        """

        image_x_len, image_y_len = image_dimensions
        image_x_min, image_x_max, image_y_max, image_y_min = extent

        assert image_x_max > image_x_min
        assert image_y_max > image_y_min

        volt_per_px_x = (image_x_max - image_x_min) / (image_x_len - 1)
        volt_per_px_y = (image_y_max - image_y_min) / (image_y_len - 1)

        V_x = volt_per_px_x * pt[0] + image_x_min
        V_y = volt_per_px_y * pt[1] + image_y_min

        return [V_x, V_y]

if __name__ == '__main__':


    script, failed, instr = Script.load_and_append({'SelectPoints':'SelectPoints'})

    print(script)
    print(failed)
    print(instr)