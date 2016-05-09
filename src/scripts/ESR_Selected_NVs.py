from src.core import Script, Parameter
from src.scripts import StanfordResearch_ESR
from src.instruments.NIDAQ import DAQ
import time
import scipy.spatial
import numpy as np
from matplotlib import patches
from PySide.QtCore import Signal, QThread
from src.core.plotting import plot_fluorescence




class ESR_Selected_NVs(Script, QThread):
    updateProgress = Signal(int)

    _DEFAULT_SETTINGS = Parameter([
        Parameter('nv_coor_path', 'Z:/Lab/Cantilever/Measurements/__test_data_for_coding/', str, 'path for NV coordinates'),
        Parameter('nv_coor_tag', 'dummy_tag', str, 'tag for NV coordinates'),
        Parameter('image_path', 'Z:/Lab/Cantilever/Measurements/__test_data_for_coding/', str, 'path for data'),
        Parameter('image_tag', 'some_name', str, 'some_name'),
        Parameter('path', 'Z:/Lab/Cantilever/Measurements/', str, 'path for data'),
        Parameter('tag', 'dummy_tag', str, 'tag for data'),
        Parameter('save', True, bool, 'save data on/off')
    ])

    _INSTRUMENTS = {'daq':  DAQ}
    _SCRIPTS = {'StanfordResearch_ESR': StanfordResearch_ESR}

    #updateProgress = Signal(int)

    #This is the signal that will be emitted during the processing.
    #By including int as an argument, it lets the signal know to expect
    #an integer argument when emitting.

    def __init__(self, instruments = None, scripts = None, name = None, settings = None,  log_output = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        self._abort = False

        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_output = log_output)

        QThread.__init__(self)

        self._plot_type = 2

        self.index = 0

        self.scripts['StanfordResearch_ESR'].updateProgress.connect(self._receive_signal)

        self.image_data = None
        self.nv_locs = None

    def _receive_signal(self, progress_sub_script):
        # calculate progress of this script based on progress in subscript

        progress = int((float(self.index)/self.num_esrs)*100 + float(progress_sub_script)/self.num_esrs)
        self.updateProgress.emit(progress)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        self._abort = False

        self.load_coors_and_image()

        esr_settings = self.scripts['StanfordResearch_ESR'].settings
        esr_instruments = self.scripts['StanfordResearch_ESR'].instruments

        self.data = {'image_data': self.image_data, 'nv_coors': self.nv_locs, 'ESR_freqs': [], 'ESR_data': np.zeros((len(self.nv_locs), esr_settings['freq_points']))}

        for index, pt in enumerate(self.nv_locs):
            if not self._abort:
                self.index = index
                self.cur_pt = pt
                self.log('Taking ESR of point ' + str(index + 1) + ' of ' + str(len(self.nv_locs)))
                #set laser to new point
                self.instruments['daq']['instance'].AO_init(["ao0","ao1"], pt)
                self.instruments['daq']['instance'].AO_run()
                self.instruments['daq']['instance'].AO_waitToFinish()
                self.instruments['daq']['instance'].AO_stop()

                #run the ESR
                # self.scripts['StanfordResearch_ESR']['instance'].tag = self.scripts['StanfordResearch_ESR']['instance'].tag + '_NV_no_' + index
                self.scripts['StanfordResearch_ESR'].run()
                self.scripts['StanfordResearch_ESR'].wait() #wait for previous ESR thread to complete

                self.data['ESR_freqs'] = self.scripts['StanfordResearch_ESR'].data[-1]['frequency']
                self.data['ESR_data'][index] = self.scripts['StanfordResearch_ESR'].data[-1]['data']

                #can't call run twice on the same object, so start a new, identical ESR script for the next run
                self.scripts['StanfordResearch_ESR'] = StanfordResearch_ESR(esr_instruments, settings = esr_settings)
                self.scripts['StanfordResearch_ESR'].updateProgress.connect(self._receive_signal)



        self.updateProgress.emit(100)
        if self.settings['save']:
            self.save()
            self.save_data()

    def stop(self):
        self._abort = True
        self.scripts['StanfordResearch_ESR'].stop()

    def plot(self, axes_Image, axes_ESR = None):
        if self.image_data is None or self.nv_locs is None:
            self.load_coors_and_image()
            print(self.nv_locs)
        image = self.image_data
        extend = [self.x_min, self.x_max, self.y_max, self.y_min]
        plot_fluorescence(image, extend, axes_Image)

        if self.isRunning():
            patch = patches.Circle((self.cur_pt[0], self.cur_pt[1]), .005, fc='b')
            axes_Image.add_patch(patch)
        else:
            for pt in self.nv_locs:
                patch = patches.Circle((pt[0], pt[1]), .005, fc='b')
                axes_Image.add_patch(patch)

        if axes_ESR is not None:
            self.scripts['StanfordResearch_ESR'].plot(axes_ESR)

    def load_coors_and_image(self):
        self.image_data = Script.load_data(self.settings['image_path'], data_name_in='image_data')
        bounds = Script.load_data(self.settings['image_path'], data_name_in='bounds')
        self.x_min = bounds[0][0]
        self.x_max = bounds[1][0]
        self.y_min = bounds[2][0]
        self.y_max = bounds[3][0]

        print('settings', self.settings)
        self.nv_locs = Script.load_data(self.settings['nv_coor_path'], data_name_in='nv_locations')
        self.num_esrs = len(self.nv_locs)


if __name__ == '__main__':
    from src.core import Instrument

    instruments, instruments_failed = Instrument.load_and_append({'daq':  'DAQ'})

    script, failed, instruments = Script.load_and_append(script_dict={'ESR_Selected_NVs':'ESR_Selected_NVs'}, instruments = instruments)

    # print(script)
    # print(failed)
    print(instruments)