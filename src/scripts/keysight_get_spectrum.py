from src.core import Script, Parameter
from PySide.QtCore import Signal, QThread
from src.instruments import SpectrumAnalyzer, MicrowaveGenerator, CryoStation
from collections import deque
import time
import numpy as np

class KeysightGetSpectrum(Script):

    # NOTE THAT THE ORDER OF Script and QThread IS IMPORTANT!!
    _DEFAULT_SETTINGS = Parameter([
        Parameter('path', 'Z:/Lab/Cantilever/Measurements/----data_tmp_default----', str, 'path for data'),
        Parameter('tag', 'dummy_tag', str, 'tag for data'),
        Parameter('save', True, bool, 'save data on/off'),
        Parameter('start_frequency', 2.7e9, float, 'start frequency of spectrum'),
        Parameter('stop_frequency', 3e9, float, 'end frequency of spectrum'),
        Parameter('output_power',0.0, float, 'output power (dBm)'),
        Parameter('output_on',True, bool, 'enable output'),
    ])

    _INSTRUMENTS = {
        'spectrum_analyzer' : SpectrumAnalyzer
    }

    _SCRIPTS = {}

    def __init__(self, instruments = None, name = None, settings = None,  log_output = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, log_output = log_output)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        def setup_instrument():
            print('self.settings',self.settings)
            inst = self.instruments['spectrum_analyzer']
            if inst.settings['start_frequency'] != self.settings['start_frequency']:
                inst.start_frequency = self.settings['start_frequency']

            if inst.settings['stop_frequency'] != self.settings['stop_frequency']:
                inst.stop_frequency = self.settings['stop_frequency']

            if self.settings['output_on']:
                if inst.settings['mode'] != 'TrackingGenerator':
                    inst.mode = 'TrackingGenerator'
                if inst.settings['output_power'] != self.settings['output_power']:
                    inst.output_power = self.settings['output_power']
                if inst.settings['output_on'] != self.settings['output_on']:
                    inst.output_on = self.settings['output_on']

        setup_instrument()

        trace = self.instruments['spectrum_analyzer'].trace

        self.data = {
            'spectrum' : [item[1] for item in trace],
            'frequency' : [item[0] for item in trace]
        }

        self.save()



    def plot(self, axes):

        spectrum = self.data['spectrum']
        freq = self.data['frequency']

        axes.plot(freq, spectrum)