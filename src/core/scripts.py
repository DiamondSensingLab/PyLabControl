import datetime
import time
from abc import ABCMeta, abstractmethod, abstractproperty
from copy import deepcopy
from src.core import Parameter, Instrument
from PyQt4 import QtCore
from collections import deque
import os
import pandas as pd
import glob
import json as json
from PySide.QtCore import Signal, QThread
from src.core.read_write_functions import save_b26_file
from src.core.loading import instantiate_instruments
class Script(object):
    # __metaclass__ = ABCMeta

    # ========================================================================================
    # ======= Following old_functions have to be customized for each instrument subclass =========
    # ========================================================================================

    # '''
    # returns the default settings of the script
    # settings contain Parameters, Instruments and Scripts
    # :return:
    # '''
    # _DEFAULT_SETTINGS = Parameter([
    #     Parameter('parameter', 1),
    #     Parameter('file_path', './some/path')
    # ])
    #

    # ========================================================================================
    # ======= Following old_functions are generic ================================================
    # ========================================================================================

    def __init__(self, name = None, settings = None, instruments = None, scripts = None, log_output = None):
        """
        executes scripts and stores script parameters and settings
        Args:
            name (optinal):  name of script, if not provided take name of function
            settings (optinal): a Parameter object that contains all the information needed in the script
            instruments (optinal): instruments used in the script
            scripts (optinal):  sub_scripts used in the script
            log_output(optinal): function reference that takes a string
        """
        if name is None:
            name = self.__class__.__name__
        self.name = name

        self._instruments = {}
        if instruments is None:
            instruments = {}
        else:
            assert isinstance(instruments, dict)
            assert set(instruments.keys()) == set(self._INSTRUMENTS.keys())
        self.instruments = instruments

        self._scripts = {}
        if scripts is None:
            scripts = {}
        self.scripts = scripts

        # set end time to be before start time if script hasn't been excecuted this tells us
        self.start_time = datetime.datetime.now()
        self.end_time = self.start_time - datetime.timedelta(seconds=1)

        self._settings = deepcopy(self._DEFAULT_SETTINGS)

        if settings is not None:
            self.update(settings)
        self._abort = False
        # data hold the data generated by the script,
        # this should either be a dictionary or a deque of dictionaries
        self.data = {}

        # a log for status outputs
        self.log_data = deque()
        # this can be overwritten
        self.log_output = log_output


    # @abstractmethod
    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        # some generic function
        raise NotImplementedError

    def log(self, string):
        """
        appends input string to log file and sends it to log function (self.log_output)
        Returns:

        """

        self.log_data.append(string)
        if self.log_output is None:
            print(string)
        else:
            self.log_output(string)


    @property
    def _DEFAULT_SETTINGS(self):
        """
        returns the default parameter_list of the script this function should be over written in any subclass
        """
        raise NotImplementedError("Subclass did not implement _DEFAULT_SETTINGS")
    @property
    def _INSTRUMENTS(self):
        """

        Returns: a dictionary of the instruments, where key is the instrument name and value is the instrument class
        if there is not instrument it should return an empty dict

        """
        raise NotImplementedError("Subclass did not implement _INSTRUMENTS")
    @property
    def _SCRIPTS(self):
        """

        Returns: a dictionary of the instruments, were key is the instrument name and value is the instrument class
        if there is not instrument it should return an empty dict

        """
        raise NotImplementedError("Subclass did not implement _SCRIPTS")

    def __str__(self):

        output_string = '{:s} (class type: {:s})\n'.format(self.name, self.__class__.__name__)

        output_string += 'settings:\n'
        for key, value in self.settings.iteritems():
            output_string += "{:s} : {:s}\n".format(key, str(value))
        return output_string
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        if isinstance(value, unicode):
            value = str(value)
        assert isinstance(value, str)
        self._name = value

    @property
    def instrumets(self):
        return self._instruments
    @instrumets.setter
    def instruments(self, instrument_dict):
        assert isinstance(instrument_dict, dict)
        assert set(instrument_dict.keys()) == set(self._INSTRUMENTS.keys()), "keys in{:s}\nkeys expected{:s}".format(str(instrument_dict.keys()), str( self._INSTRUMENTS.keys()))
        for key, value in self._INSTRUMENTS.iteritems():
            self._instruments.update({key: instrument_dict[key]})

    @property
    def scripts(self):
        return self._scripts
    @scripts.setter
    def scripts(self, script_dict):
        assert isinstance(script_dict, dict)
        assert script_dict.keys() == self._SCRIPTS.keys(), "keys in{:s}\nkeys expected{:s}".format(str(script_dict.keys()), str( self._SCRIPTS.keys()))

        for key, value in self._SCRIPTS.iteritems():
            assert isinstance(script_dict[key], self._SCRIPTS[key])
            self._scripts.update({key: script_dict[key]})

    @property
    def settings(self):
        '''
        :return: returns the settings of the script
        settings contain Parameters, Instruments and Scripts
        '''
        return self._settings

    def update(self, settings):
        '''
        updates the internal dictionary and sends changed values to instrument
        Args:
            settings: parameters to be set
        # mabe in the future:
        # Returns: boolean that is true if update successful

        '''
        self._settings.update(settings)

    @property
    def end_time(self):
        '''
        time when script execution started
        :return:
        '''
        return self._time_stop
    @end_time.setter
    def end_time(self, value):
        assert isinstance(value, datetime.datetime)
        self._time_stop = value

    @property
    def start_time(self):
        '''
        time when script execution started
        :return:
        '''
        return self._time_start
    @start_time.setter
    def start_time(self, value):
        assert isinstance(value, datetime.datetime)
        self._time_start = value

    @property
    def excecution_time(self):
        '''
        :return: script excecition time as time_delta object to get time in seconds use .total_seconds()
        '''
        return self.end_time - self.start_time

    def run(self):
        '''
        executes the script
        :return: boolean if execution of script finished succesfully
        '''
        self.is_running = True
        self.start_time  = datetime.datetime.now()
        self.log('starting script {:s} at {:s} on {:s}'.format(self.name, self.start_time.strftime('%H:%M:%S'),self.start_time.strftime('%d/%m/%y')))
        self._function()

        self.end_time  = datetime.datetime.now()
        self.log('script {:s} finished at {:s} on {:s}'.format(self.name, self.end_time.strftime('%H:%M:%S'),self.end_time.strftime('%d/%m/%y')))
        success = self._abort == False
        self.is_running = False
        return success

    def stop(self):
        self._abort == True

    def filename(self, appendix):
        """
        creates a filename based
        Args:
            appendix: appendix for file

        Returns: filename

        """
        path = self.settings['path']
        tag = self.settings['tag']

        if os.path.exists(path) == False:
            os.makedirs(path)

        filename = "{:s}\\{:s}_{:s}{:s}".format(
            path,
            self.start_time.strftime('%y%m%d-%H_%M_%S'),
            tag,
            appendix
        )
        return filename

    def save_data(self, filename = None):
        """
        saves the script data to a file: filename is filename is not provided, it is created from internal function

        Returns:

        """

        if filename is None:
            filename = self.filename('-data.csv')

        # if deque object, take the last dataset, which is the most recent
        if isinstance(self.data, deque):
            data = self.data[-1]
        elif isinstance(self.data, dict):
            data = self.data
        else:
            raise TypeError("unknown datatype!")

        if len(set([len(v) for k, v in data.iteritems()])) == 1:
            # if all entries of the dictionary are the same length we can write the data into a single file

            df = pd.DataFrame(data)
            df.to_csv(filename, index=False)

        else:
            # otherwise, we write each entry into a separate file into a subfolder data
            for key, value in self.data.iteritems():
                df = pd.DataFrame(value)
                df.to_csv(filename.replace('-data.csv', '-{:s}.csv'.format(key)), index=False)

    def save_log(self, filename = None):
        """
        save log to file
        Returns:

        """

        if filename is None:
            filename = self.filename('-info.txt')

        with open(filename, 'w') as outfile:
            for item in self.log_data:
                outfile.write("%s\n" % item)

    def save(self, filename=None):
        """
        saves the script settings to a file: filename is filename is not provided, it is created from internal function
        """
        if filename is None:
            filename = self.filename('.b26')

        save_b26_file(filename, scripts=self.to_dict())

    def plot(self, axes):
        """
        plots the data contained in self.data, which should be a dictionary or a deque of dictionaries
        for the latter use the last entry

        """

        if isinstance(self.data, deque):
            data = self.data[-1]
        elif isinstance(self.data, dict):
            data = self.data
        else:
            data = {}

        assert isinstance(data, dict)

        if data == {}:
            self.log("warning, not data found that can be plotted")
        else:
            for key, value in data.iteritems():
                axes.plot(value)

    def to_dict(self):
        """

        Returns: itself as a dictionary

        """

        dictator = {self.name: {'class' : self.__class__.__name__}}

        if self.scripts != {}:
            dictator[self.name].update({'scripts': {} })
            for subscript_name, subscript in self.scripts.iteritems():
                dictator[self.name]['scripts'].update(subscript.to_dict() )

        if self.instruments != {}:
            # dictator[self.name].update({'instruments': self.instruments})
            # dictator[self.name].update({'instruments': {} })
            # for instrument_name, instrument in self.instruments.iteritems():
            #     dictator[self.name]['instruments'].update(instrument.to_dict())

            dictator[self.name].update({'instruments': {
                instrument_name: {'class': instrument['instance'].__class__.__name__, 'settings':instrument['settings']}
                for instrument_name, instrument in self.instruments.iteritems()
            }})

        dictator[self.name]['settings'] = self.settings

        print('GGG', dictator)
        return dictator

    @staticmethod
    def load_data(path, time_tag_in = None):
        """
        loads the data that has been save with Script.save
        Args:
            path: path to data
            time_tag: (optional) time tag of data if None, returns all data sets that have been found
        Returns:
            a dictionary with the data of form
            data = {time_tag:data_set}
            with
            data_set = {'varible1':data, 'variable2' : data, etc.}

            if there is only a single time tag, return data_set
        """
        data = {}
        if 'data' in os.listdir(path):
            data_files = os.listdir(path + '\data')
            data_names = set([f.split('.')[1] for f  in  data_files])
            time_tags = set(['_'.join(f.split('_')[0:3]) for f  in  data_files])


            for time_tag in time_tags:
                sub_data = {}
                for data_name in data_names:
                    file_path = "{:s}\\data\\{:s}*.{:s}".format(DIR, time_tag, data_name)
                    df = pd.read_csv(glob.glob(file_path)[0])
                    sub_data.update({data_name:list(df[list(df)[0]])})
                data.update({time_tag: sub_data})

        else:
            data_files = glob.glob(path + '*.dat')

            for data_file in data_files:
                time_tag = '_'.join(data_file.split('\\')[-1].split('_')[0:3])
                df = pd.read_csv(data_file)
                data_names = list(df)
                data.update({time_tag:
                            {data_name:list(df[data_name]) for data_name in data_names}
                            })


        if time_tag_in != None:
            # if user specifies a time_tag we only return that specific data set
            data = data[time_tag_in]
        elif len(data) == 1:
            # if there is only one data_set, we strip of the time_tag level
            data = data[data.keys()[0]]

        return data


    @staticmethod
    def load_and_append(script_dict, scripts = {}, instruments = {}, log_function = None):
        """
        load script from script_dict and append to scripts, if additional instruments are required create them and add them to instruments

        Args:
            script_dict: dictionary of form

                script_dict = {
                name_of_script_1 :
                    {"settings" : settings_dictionary, "class" : name_of_class}
                name_of_instrument_2 :
                    {"settings" : settings_dictionary, "class" : name_of_class}
                ...
                }

            or

                script_dict = {
                name_of_script_1 : name_of_class,
                name_of_script_2 : name_of_class
                ...
                }

            where name_of_class is either a class or the name of a class

            scripts: dictionary of form

                scripts = {
                name_of_script_1 : instance_of_instrument_1,
                name_of_script_2 : instance_of_instrument_2,
                ...
                }

            instruments: dictionary of form

                instruments = {
                name_of_instrument_1 : instance_of_instrument_1,
                name_of_instrument_2 : instance_of_instrument_2,
                ...
                }

        Returns:
                list of form

                loaded_failed = [name_of_script_1, name_of_script_2, ....]

        """
        loaded_failed = []
        updated_scripts = {}
        updated_scripts.update(scripts)
        updated_instruments = {}
        updated_instruments.update(instruments)

        def get_script_information(script_information):
            """

            Args:
                script_information: information of the script. This can be
                    - a dictionary
                    - a Script instance
                    - name of Script class

            Returns: module_path, script_class_name, script_settings, script_instruments, script_sub_scripts

            """
            script_settings = None
            script_instruments = None
            script_sub_scripts = None
            script_class_name = None

            print('script_information', script_information,  type(script_information))

            if isinstance(script_information, dict):
                script_settings = script_information['settings']
                script_class_name = str(script_information['class'])
                if 'instruments' in script_information:
                    script_instruments = script_information['instruments']
                if 'scripts' in script_information:
                    script_sub_scripts = script_information['scripts']
            elif isinstance(script_information, str):
                script_class_name = script_information
            elif issubclass(script_information, Script):
                # watch out when testing this code from __main__, then classes might not be identified correctly because the path is different
                # to avoid this problem call from src.core import Script (otherwise the path to Script is __main__.Script)
                script_class_name = script_information.__name__
            if len(script_class_name.split('.')) == 1:
                module_path = 'src.scripts'
            else:
                module_path = 'src.scripts.' + '.'.join(script_class_name.split('.')[0:-1])
                script_class_name = script_information.split('.')[-1]

            return module_path, script_class_name, script_settings, script_instruments, script_sub_scripts

        def get_instruments(class_of_script, script_instruments, instruments):
            """

            creates the dictionary with the instruments needed for the script and update the instrument dictionary if new instruments are required

            Args:
                class_of_script: the class of the script
                instruments: the instruments that have been loaded already

            Returns: dictionary with the instruments that the script needs and the updated instruments dictionary

            """
            default_instruments = getattr(class_of_script, '_INSTRUMENTS')
            instrument_dict = {}
            instruments_updated = {}
            instruments_updated.update(instruments)
            # check if instruments needed by script already exist, if not create an instance
            for instrument_name, instrument_class in default_instruments.iteritems():
                # check if instruments needed by script already exist
                instrument = [instance for name, instance in instruments_updated.iteritems() if
                              isinstance(instance, instrument_class) and name == instrument_name]
                if len(instrument) == 0:
                    # create new instance of instrument
                    instruments_updated, __ = Instrument.load_and_append({instrument_name: instrument_class.__name__}, instruments_updated)

                if script_instruments is not None and instrument_name in script_instruments:
                    instrument_settings = script_instruments[instrument_name]['settings']
                else:
                    instrument_settings = instruments_updated[instrument_name].settings

                instrument_instance = instruments_updated[instrument_name]

                instrument_dict.update({instrument_name: {"instance":instrument_instance, "settings":instrument_settings}})


            return instrument_dict, instruments_updated

        def get_sub_scripts(class_of_script, instruments):
            """

            creates the dictionary with the sub scripts needed by the script and updates the instrument dictionary if new instruments are required

            Args:
                class_of_script: the class of the script
                instruments: the instruments that have been loaded already

            Returns: dictionary with the sub scripts that the script needs

            """
            default_scripts = getattr(class_of_script, '_SCRIPTS')
            # create instruments that script needs
            sub_scripts = {}
            sub_scripts, scripts_failed, instruments_updated = Script.load_and_append(default_scripts, sub_scripts, instruments)

            # todo: update the scripts with the settings if there are any
            # .....
            return sub_scripts, instruments_updated

        for script_name, script_class_name in script_dict.iteritems():

            # check if script already exists
            if script_name in scripts.keys():
                print('WARNING: script {:s} already exists. Did not load!'.format(script_name))
                loaded_failed.append(script_name)
            else:
                module_path, script_class_name, script_settings, script_instruments, script_sub_scripts = get_script_information(script_class_name)
                script_instance = None

                # try to import the script
                module = __import__(module_path, fromlist=[script_class_name])
                # this returns the name of the module that was imported.
                class_of_script = getattr(module, script_class_name)

                #  ========= get the instruments that are needed by the script =========
                script_instruments, updated_instruments = get_instruments(class_of_script, script_instruments, updated_instruments)
                #  ========= create the scripts that are needed by the script =========
                sub_scripts, updated_instruments = get_sub_scripts(class_of_script, updated_instruments)

                class_creation_string = ''
                if script_instruments != {}:
                    class_creation_string += ', instruments = script_instruments'
                if sub_scripts != {}:
                    class_creation_string += ', scripts = sub_scripts'
                if script_settings != {}:
                    class_creation_string += ', settings = script_settings'
                if log_function is not None:
                    class_creation_string += ', log_output = log_function'
                class_creation_string = 'class_of_script(name=script_name{:s})'.format(class_creation_string)

                script_instance = eval(class_creation_string)
                updated_scripts.update({script_name :script_instance})

        return updated_scripts, loaded_failed, updated_instruments

class QThreadWrapper(QThread):


    updateProgress = Signal(int)

    def __init__(self, script):
        """
        This is a wrapper for scripts that are not QThread, to execute them on a different thread than the gui
        Args:
            script: script to be executed

        """
        self.script = script
        QThread.__init__(self)

    def run(self):

        self._running = True
        self.updateProgress.emit(1)
        self.script.run()
        self.updateProgress.emit(100)


#
# if __name__ == '__main__':
#     from src.core.read_write_functions import load_b26_file
#     from src.instruments import DummyInstrument
#     filename = "Z:\Lab\Cantilever\Measurements\\__tmp\\XX.b26"
#     data = load_b26_file(filename)
#     scripts = {}
#     instruments = {"dummy_instrument": DummyInstrument("dummy_instrument")}
#     scripts_failed = Script.load_and_append(data['scripts'], scripts, instruments)
#     print('loaded', scripts)
#     print('failed', scripts_failed)
#
#
#     print(scripts['dummy script with inst'].settings)
#     print(scripts['dummy script with inst'].instruments['dummy_instrument'])
#     print(scripts['dummy script with inst'].scripts)



if __name__ == '__main__':
    from src.core.read_write_functions import load_b26_file
    from src.instruments import DummyInstrument
    from src.scripts import ScriptDummyWithInstrument
    path = "Z:\Lab\Cantilever\Measurements\\__tmp\\"
    from src.core import Script
    # create script


    scripts, loaded_failed, instruments = Script.load_and_append({"script with inst": 'ScriptDummyWithInstrument', "script with suscript": 'ScriptDummyWithSubScript'})
    print(instruments)
    print(scripts)

    for name, script in scripts.iteritems():
        script.save('{:s}{:s}.b26'.format(path, name.replace(' ', '_')))

    #
    #
    # data = load_b26_file(filename)
    # scripts = {}
    # instruments = {}
    # scripts, scripts_failed, instruments_2 = Script.load_and_append(data['scripts'], scripts, instruments)
    # print('loaded', scripts)
    # print('failed', scripts_failed)
    #
    # print(instruments)
    # print(instruments_2)