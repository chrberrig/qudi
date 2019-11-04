# -*- coding: utf-8 -*-
"""
Use Swabian Instruments PulseStreamer8/2 as a pulse generator.

Protobuf (pb2) and grpc files generated from pulse_streamer.proto
file available at https://www.swabianinstruments.com/static/documentation/PulseStreamer/sections/interface.html#grpc-interface.

Regenerate files for an update proto file using the following:
python3 -m grpc_tools.protoc -I=./ --python_out=. --grpc_python_out=. ./pulse_streamer.proto

Qudi is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Qudi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Qudi. If not, see <http://www.gnu.org/licenses/>.

Copyright (c) the Qudi Developers. See the COPYRIGHT.txt file at the
top-level directory of this distribution and at <https://github.com/Ulm-IQO/qudi/>
"""

from core.module import Base
from core.configoption import ConfigOption
from core.util.modules import get_home_dir
from core.util.helpers import natural_sort
from interface.pulser_interface import PulserInterface, PulserConstraints
from collections import OrderedDict

# import numpy as np
# import time
import pickle
import grpc
import os
import hardware.swabian_instruments.pulse_streamer_pb2 as pulse_streamer_pb2
#import dill


class PulseStreamer(Base, PulserInterface):
    """ Methods to control PulseStreamer.

    Example config for copy-paste:

    pulse_streamer:
        module.Class: 'swabian_instruments.pulse_streamer.PulseStreamer'
        pulsestreamer_ip: '192.168.1.100'
        laser_channel: 0
        uw_x_channel: 2

    """

    _pulsestreamer_ip = ConfigOption('pulsestreamer_ip', '10.54.10.64', missing='warn')
    _laser_channel = ConfigOption('laser_channel', 0, missing='warn')
    _uw_x_channel = ConfigOption('uw_x_channel', 2, missing='warn')



    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        if 'pulsed_file_dir' in config.keys():
            self.pulsed_file_dir = config['pulsed_file_dir']

            if not os.path.exists(self.pulsed_file_dir):
                homedir = get_home_dir()
                self.pulsed_file_dir = os.path.join(homedir, 'pulsed_files')
                self.log.warning('The directory defined in parameter '
                            '"pulsed_file_dir" in the config for '
                            'PulseStreamer does not exist!\n'
                            'The default home directory\n{0}\n will be taken '
                            'instead.'.format(self.pulsed_file_dir))
        else:
            homedir = get_home_dir()
            self.pulsed_file_dir = os.path.join(homedir, 'pulsed_files')
            self.log.warning('No parameter "pulsed_file_dir" was specified in the config for '
                             'PulseStreamer as directory for the pulsed files!\nThe default home '
                             'directory\n{0}\nwill be taken instead.'.format(self.pulsed_file_dir))

        # self.host_waveform_directory = self._get_dir_for_name('sampled_hardware_files')

        self.current_status = -1
        self.sample_rate = 1e9
        self.current_loaded_asset = {}, None

        self._channel = grpc.insecure_channel(self._pulsestreamer_ip + ':50051')

    def on_activate(self):
        """ Establish connection to pulse streamer and tell it to cancel all operations """
        self.pulse_streamer = pulse_streamer_pb2.PulseStreamerStub(self._channel)
        self.pulser_off()
        self.current_status = 0

    def on_deactivate(self):
        del self.pulse_streamer

    def get_constraints(self):
        """ Retrieve the hardware constrains from the Pulsing device.

        @return dict: dict with constraints for the sequence generation and GUI

        Provides all the constraints (e.g. sample_rate, amplitude,
        total_length_bins, channel_config, ...) related to the pulse generator
        hardware to the caller.
        The keys of the returned dictionary are the str name for the constraints
        (which are set in this method). No other keys should be invented. If you
        are not sure about the meaning, look in other hardware files to get an
        impression. If still additional constraints are needed, then they have
        to be add to all files containing this interface.
        The items of the keys are again dictionaries which have the generic
        dictionary form:
            {'min': <value>,
             'max': <value>,
             'step': <value>,
             'unit': '<value>'}

        Only the keys 'activation_config' and differs, since it contain the
        channel configuration/activation information.

        If the constraints cannot be set in the pulsing hardware (because it
        might e.g. has no sequence mode) then write just zero to each generic
        dict. Note that there is a difference between float input (0.0) and
        integer input (0).
        ALL THE PRESENT KEYS OF THE CONSTRAINTS DICT MUST BE ASSIGNED!
        """
        constraints = PulserConstraints()


        constraints.waveform_format = ['pstream']
        constraints.sequence_format = []

        constraints.sample_rate.min = 1e9
        constraints.sample_rate.max = 1e9
        constraints.sample_rate.step = 0
        constraints.sample_rate.default = 1e9

        constraints.d_ch_low.min = 0.0
        constraints.d_ch_low.max = 0.0
        constraints.d_ch_low.step = 0.0
        constraints.d_ch_low.default = 0.0

        constraints.d_ch_high.min = 3.3
        constraints.d_ch_high.max = 3.3
        constraints.d_ch_high.step = 0.0
        constraints.d_ch_high.default = 3.3

        # sample file length max is not well-defined for PulseStreamer, which collates sequential identical pulses into
        # one. Total number of not-sequentially-identical pulses which can be stored: 1 M.
        constraints.waveform_length.min = 1
        constraints.waveform_length.max = 134217728
        constraints.waveform_length.step = 1
        constraints.waveform_length.default = 1

        # the name a_ch<num> and d_ch<num> are generic names, which describe UNAMBIGUOUSLY the
        # channels. Here all possible channel configurations are stated, where only the generic
        # names should be used. The names for the different configurations can be customary chosen.
        activation_config = OrderedDict()
        activation_config['digital'] = frozenset({'d_ch1', 'd_ch2', 'd_ch3', 'd_ch4', 'd_ch5', 'd_ch6', 'd_ch7',
                                    'd_ch8'})
        activation_config['analogue'] = frozenset({'a_ch1', 'a_ch2'})
        activation_config['all'] = frozenset({'d_ch1', 'd_ch2', 'd_ch3', 'd_ch4', 'd_ch5', 'd_ch6', 'd_ch7',
                                    'd_ch8', 'a_ch1', 'a_ch2'})
        constraints.activation_config = activation_config

        return constraints

        # # The file formats are hardware specific.
        # constraints['waveform_format'] = 'pstream'
        # constraints['sequence_format'] = None
        #
        # constraints['sample_rate'] = {
        #     'min': 1e9,
        #     'max': 1e9,
        #     'step': 0,
        #     'default': 1e9}
        # constraints['d_ch_low'] = {
        #     'min': 0.0,
        #     'max': 0.0,
        #     'step': 0.0,
        #     'default': 0.0}
        # constraints['d_ch_high'] = {
        #     'min': 3.3,
        #     'max': 3.3,
        #     'step': 0.0,
        #     'default': 3.3}
        # # sample file length max is not well-defined for PulseStreamer, which collates sequential identical pulses into
        # # one. Total number of not-sequentially-identical pulses which can be stored: 1 M.
        # constraints['waveform_length'] = {
        #     'min': 1,
        #     'max': 134217728,
        #     'step': 1,
        #     'default': 1}
        #
        # # the name a_ch<num> and d_ch<num> are generic names, which describe UNAMBIGUOUSLY the
        # # channels. Here all possible channel configurations are stated, where only the generic
        # # names should be used. The names for the different configurations can be customary chosen.
        # activation_config = OrderedDict()
        # activation_config['all'] = ['d_ch1', 'd_ch2', 'd_ch3', 'd_ch4', 'd_ch5', 'd_ch6', 'd_ch7',
        #                             'd_ch8']
        # constraints['activation_config'] = activation_config
        # return constraints


    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:OK, -1:error)
        """
        # start the pulse sequence
        self.pulse_streamer.stream(self._sequence)
        self.log.info('Asset uploaded to PulseStreamer')
        self.pulse_streamer.startNow(pulse_streamer_pb2.VoidMessage())
        self.current_status = 1
        return 0

    def pulser_off(self):
        """ Switches the pulsing device off.

        @return int: error code (0:OK, -1:error)
        """
        # stop the pulse sequence
        channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        self.pulse_streamer.constant(pulse_streamer_pb2.PulseMessage(ticks=0, digi=channels, ao0=0, ao1=0))
        self.current_status = 0
        return 0

    def get_loaded_assets(self):
        """
        Retrieve the currently loaded asset names for each active channel of the device.
        The returned dictionary will have the channel numbers as keys.
        In case of loaded waveforms the dictionary values will be the waveform names.
        In case of a loaded sequence the values will be the sequence name appended by a suffix
        representing the track loaded to the respective channel (i.e. '<sequence_name>_1').

        @return (dict, str):    Dictionary with keys being the channel number and values being the
                                respective asset loaded into the channel,
                                string describing the asset type ('waveform' or 'sequence')
         """
        return self.current_loaded_asset

    def delete_waveform(self, waveform_name):
        """ Delete the waveform with name "waveform_name" from the device memory.

        @param str waveform_name: The name of the waveform to be deleted
                                  Optionally a list of waveform names can be passed.

        @return list: a list of deleted waveform names.
        Unused for digital pulse generators without sequence storage capability
        (PulseBlaster, FPGA).
        """
        pass

    def delete_sequence(self, sequence_name):
        """ Delete the sequence with name "sequence_name" from the device memory.

        @param str sequence_name: The name of the sequence to be deleted
                                  Optionally a list of sequence names can be passed.

        @return list: a list of deleted sequence names.
        Unused for digital pulse generators without sequence storage capability
        (PulseBlaster, FPGA).
        """
        pass

    def write_waveform(self, name, analog_samples, digital_samples, is_first_chunk, is_last_chunk,
                       total_number_of_samples):
        """
        Write a new waveform or append samples to an already existing waveform on the device memory.
        The flags is_first_chunk and is_last_chunk can be used as indicator if a new waveform should
        be created or if the write process to a waveform should be terminated.

        NOTE: All sample arrays in analog_samples and digital_samples must be of equal length!

        @param str name: the name of the waveform to be created/append to
        @param dict analog_samples: keys are the generic analog channel names (i.e. 'a_ch1') and
                                    values are 1D numpy arrays of type float32 containing the
                                    voltage samples normalized to half Vpp (between -1 and 1).
        @param dict digital_samples: keys are the generic digital channel names (i.e. 'd_ch1') and
                                     values are 1D numpy arrays of type bool containing the marker
                                     states.
        @param bool is_first_chunk: Flag indicating if it is the first chunk to write.
                                    If True this method will create a new empty waveform.
                                    If False the samples are appended to the existing waveform.
        @param bool is_last_chunk:  Flag indicating if it is the last chunk to write.
                                    Some devices may need to know when to close the appending wfm.
        @param int total_number_of_samples: The number of sample points for the entire waveform
                                            (not only the currently written chunk)

        @return (int, list): Number of samples written (-1 indicates failed process) and list of
                             created waveform names
        """
        #TODO make a print statement in pulsesequencegeneration logic module to see what kind of object is expected as input.
        #TODO Also, what are these markerstates that is talked about in the doc-string?

    def write_sequence(self, name, sequence_parameters):
        """
        Write a new sequence on the device memory.

        @param str name: the name of the waveform to be created/append to
        @param list sequence_parameters: List containing tuples of length 2. Each tuple represents
                                         a sequence step. The first entry of the tuple is a list of
                                         waveform names (str); one for each channel. The second
                                         tuple element is a SequenceStep instance containing the
                                         sequencing parameters for this step.

        @return: int, number of sequence steps written (-1 indicates failed process)
        """
        pass

    def get_waveform_names(self):
        """ Retrieve the names of all uploaded waveforms on the device.

        @return list: List of all uploaded waveform name strings in the device workspace.
        """
        return []

    def get_sequence_names(self):
        """ Retrieve the names of all uploaded sequence on the device.

        @return list: List of all uploaded sequence name strings in the device workspace.
        """
        return []

    def load_waveform(self, load_dict):
        """ Loads a waveform to the specified channel of the pulsing device.
        For devices that have a workspace (i.e. AWG) this will load the waveform from the device
        workspace into the channel.
        For a device without mass memory this will make the waveform/pattern that has been
        previously written with self.write_waveform ready to play.

        @param load_dict:  dict|list, a dictionary with keys being one of the available channel
                                      index and values being the name of the already written
                                      waveform to load into the channel.
                                      Examples:   {1: rabi_ch1, 2: rabi_ch2} or
                                                  {1: rabi_ch2, 2: rabi_ch1}
                                      If just a list of waveform names if given, the channel
                                      association will be invoked from the channel
                                      suffix '_ch1', '_ch2' etc.

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel, string describing the asset
                             type ('waveform' or 'sequence')
        """
        # ignore if no asset_name is given
        if load_dict is None:
            self.log.warning('"load_asset" called with asset_name = None.')
            return {}

        # check if asset exists
        saved_assets = self.get_waveform_names()
        if load_dict not in saved_assets:
            self.log.error('No waveform with name "{0}" found for PulseStreamer.\n'
                                '"load_waveform" call ignored.'.format(load_dict))
            return -1

        # get samples from file
        # filepath = os.path.join(self.pulsed_file_dir, load_dict + '.ensemble')
        # pulse_waveform_raw = pickle.load(open(filepath, 'rb'))

        pulse_waveform = []
        for pulse in pulse_waveform_raw:
            pulse_waveform.append(pulse_streamer_pb2.PulseMessage(ticks=pulse[0], digi=pulse[1], ao0=0, ao1=1))

        blank_pulse = pulse_streamer_pb2.PulseMessage(ticks=0, digi=0, ao0=0, ao1=0)
        laser_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=self._convert_to_bitmask([self._laser_channel]), ao0=0, ao1=0)
        laser_and_uw_channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        laser_and_uw_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=laser_and_uw_channels, ao0=0, ao1=0)
        self._waveform = pulse_streamer_pb2.SequenceMessage(pulse=pulse_waveform, n_runs=0, initial=laser_on,
            final=laser_and_uw_on, underflow=blank_pulse, start=1)

        self.current_loaded_asset = load_dict, 'waveform'

        return load_dict

    def load_sequence(self, sequence_name):
        """ Loads a sequence to the channels of the device in order to be ready for playback.
        For devices that have a workspace (i.e. AWG) this will load the sequence from the device
        workspace into the channels.

        @param sequence_name:  str, name of the sequence to load

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel, string describing the asset
                             type ('waveform' or 'sequence')
        """

        # ignore if no asset_name is given
        if sequence_name is None:
            self.log.warning('"load_asset" called with asset_name = None.')
            return {}

        # check if asset exists
        saved_assets = self.get_sequence_names()
        if sequence_name not in saved_assets:
            self.log.error('No sequence with name "{0}" found for PulseStreamer.\n'
                           '"load_sequence" call ignored.'.format(sequence_name))
            return -1

        # get samples from file
        filepath = os.path.join(self.pulsed_file_dir, sequence_name + '.pstream')
        pulse_sequence_raw = pickle.load(open(filepath, 'rb'))
        # pulse_sequence_raw = samples_write_methods._write_pstream(data)

        pulse_sequence = []
        for pulse in pulse_sequence_raw:
            pulse_sequence.append(pulse_streamer_pb2.PulseMessage(ticks=pulse[0], digi=pulse[1], ao0=0, ao1=1))

        blank_pulse = pulse_streamer_pb2.PulseMessage(ticks=0, digi=0, ao0=0, ao1=0)
        laser_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=self._convert_to_bitmask([self._laser_channel]), ao0=0,
                                                   ao1=0)
        laser_and_uw_channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        laser_and_uw_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=laser_and_uw_channels, ao0=0, ao1=0)
        self._sequence = pulse_streamer_pb2.SequenceMessage(pulse=pulse_sequence, n_runs=0, initial=laser_on,
                                                            final=laser_and_uw_on, underflow=blank_pulse, start=1)
        self.current_loaded_asset = sequence_name, 'sequence'

        return sequence_name

    def clear_all(self):
        """ Clears all loaded waveforms from the pulse generators RAM.

        @return int: error code (0:OK, -1:error)

        Unused for digital pulse generators without storage capability
        (PulseBlaster, FPGA).
        """
        return 0

    def get_status(self):
        """ Retrieves the status of the pulsing hardware

        @return (int, dict): tuple with an interger value of the current status
                             and a corresponding dictionary containing status
                             description for all the possible status variables
                             of the pulse generator hardware.
        """
        status_dic = dict()
        status_dic[-1] = 'Failed Request or Failed Communication with device.'
        status_dic[0] = 'Device has stopped, but can receive commands.'
        status_dic[1] = 'Device is active and running.'

        return self.current_status, status_dic

    def get_sample_rate(self):
        """ Get the sample rate of the pulse generator hardware

        @return float: The current sample rate of the device (in Hz)

        Do not return a saved sample rate in a class variable, but instead
        retrieve the current sample rate directly from the device.
        """
        return self.sample_rate

    def set_sample_rate(self, sample_rate):
        """ Set the sample rate of the pulse generator hardware.

        @param float sample_rate: The sampling rate to be set (in Hz)

        @return float: the sample rate returned from the device.

        Note: After setting the sampling rate of the device, retrieve it again
              for obtaining the actual set value and use that information for
              further processing.
        """
        self.log.debug('PulseStreamer sample rate cannot be configured')
        return self.sample_rate

    def get_analog_level(self, amplitude=None, offset=None):
        """ Retrieve the analog amplitude and offset of the provided channels.

        @param list amplitude: optional, if a specific amplitude value (in Volt
                               peak to peak, i.e. the full amplitude) of a
                               channel is desired.
        @param list offset: optional, if a specific high value (in Volt) of a
                            channel is desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel
                               number and items being the values for those
                               channels. Amplitude is always denoted in
                               Volt-peak-to-peak and Offset in (absolute)
                               Voltage.
        """
        return {}, {}

    def set_analog_level(self, amplitude=None, offset=None):
        """ Set amplitude and/or offset value of the provided analog channel.

        @param dict amplitude: dictionary, with key being the channel and items
                               being the amplitude values (in Volt peak to peak,
                               i.e. the full amplitude) for the desired channel.
        @param dict offset: dictionary, with key being the channel and items
                            being the offset values (in absolute volt) for the
                            desired channel.

        @return (dict, dict): tuple of two dicts with the actual set values for
                              amplitude and offset.

        If nothing is passed then the command will return two empty dicts.
        """
        return {}, {}

    def get_digital_level(self, low=None, high=None):
        """ Retrieve the digital low and high level of the provided channels.

        @param list low: optional, if a specific low value (in Volt) of a
                         channel is desired.
        @param list high: optional, if a specific high value (in Volt) of a
                          channel is desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel
                               number and items being the values for those
                               channels. Both low and high value of a channel is
                               denoted in (absolute) Voltage.

        Note: Do not return a saved low and/or high value but instead retrieve
              the current low and/or high value directly from the device.

        If no entries provided then the levels of all channels where simply
        returned. If no digital channels provided, return just an empty dict.

        Example of a possible input:
            low = [1,4]
        to obtain the low voltage values of digital channel 1 an 4. A possible
        answer might be
            {1: -0.5, 4: 2.0} {}
        since no high request was performed.

        The major difference to analog signals is that digital signals are
        either ON or OFF, whereas analog channels have a varying amplitude
        range. In contrast to analog output levels, digital output levels are
        defined by a voltage, which corresponds to the ON status and a voltage
        which corresponds to the OFF status (both denoted in (absolute) voltage)

        In general there is no bijective correspondence between
        (amplitude, offset) and (value high, value low)!
        """
        if low is None:
            low = []
        if high is None:
            high = []
        low_dict = {}
        high_dict = {}
        if low is [] and high is []:
            for channel in range(8):
                low_dict[channel] = 0.0
                high_dict[channel] = 3.3
        else:
            for channel in low:
                low_dict[channel] = 0.0
            for channel in high:
                high_dict[channel] = 3.3
        return low_dict, high_dict

    def set_digital_level(self, low=None, high=None):
        """ Set low and/or high value of the provided digital channel.

        @param dict low: dictionary, with key being the channel and items being
                         the low values (in volt) for the desired channel.
        @param dict high: dictionary, with key being the channel and items being
                         the high values (in volt) for the desired channel.

        @return (dict, dict): tuple of two dicts where first dict denotes the
                              current low value and the second dict the high
                              value.

        If nothing is passed then the command will return two empty dicts.

        Note: After setting the high and/or low values of the device, retrieve
              them again for obtaining the actual set value(s) and use that
              information for further processing.

        The major difference to analog signals is that digital signals are
        either ON or OFF, whereas analog channels have a varying amplitude
        range. In contrast to analog output levels, digital output levels are
        defined by a voltage, which corresponds to the ON status and a voltage
        which corresponds to the OFF status (both denoted in (absolute) voltage)

        In general there is no bijective correspondence between
        (amplitude, offset) and (value high, value low)!
        """
        if low is None:
            low = {}
        if high is None:
            high = {}
        self.log.warning('PulseStreamer logic level cannot be adjusted!')
        return 0

    def get_active_channels(self, ch=None):
        if ch is None:
            ch = {}
        d_ch_dict = {}
        if len(ch) < 1:
            for chnl in range(1, 9):
                d_ch_dict['d_ch{0}'.format(chnl)] = True
        else:
            for channel in ch:
                d_ch_dict[channel] = True
        return d_ch_dict

    def set_active_channels(self, ch=None):
        if ch is None:
            ch = {}
        d_ch_dict = {
            'd_ch1': True,
            'd_ch2': True,
            'd_ch3': True,
            'd_ch4': True,
            'd_ch5': True,
            'd_ch6': True,
            'd_ch7': True,
            'd_ch8': True}
        return d_ch_dict

    def get_interleave(self):
        """ Check whether Interleave is ON or OFF in AWG.

        @return bool: True: ON, False: OFF

        Unused for pulse generator hardware other than an AWG.
        """
        return False

    def set_interleave(self, state=False):
        """ Turns the interleave of an AWG on or off.

        @param bool state: The state the interleave should be set to
                           (True: ON, False: OFF)

        @return bool: actual interleave status (True: ON, False: OFF)

        Note: After setting the interleave of the device, retrieve the
              interleave again and use that information for further processing.

        Unused for pulse generator hardware other than an AWG.
        """
        return False

    def reset(self):
        """ Reset the device.

        @return int: error code (0:OK, -1:error)
        """
        channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        self.pulse_streamer.constant(pulse_streamer_pb2.PulseMessage(ticks=0, digi=channels, ao0=0, ao1=0))
        self.pulse_streamer.constant(laser_on)
        return 0

#  ====== Internal functions ======

    def _get_dir_for_name(self, name):
        """ Get the path to the pulsed sub-directory 'name'.

        @param name: string, name of the folder
        @return: string, absolute path to the directory with folder 'name'.
        """
        path = os.path.join(self.pulsed_file_dir, name)
        if not os.path.exists(path):
            os.makedirs(os.path.abspath(path))
        return os.path.abspath(path)


    def _convert_to_bitmask(self, active_channels):
        """ Convert a list of channels into a bitmask.
        @param numpy.array active_channels: the list of active channels like
                            e.g. [0,4,7]. Note that the channels start from 0.
        @return int: The channel-list is converted into a bitmask (an sequence
                     of 1 and 0). The returned integer corresponds to such a
                     bitmask.
        Note that you can get a binary representation of an integer in python
        if you use the command bin(<integer-value>). All higher unneeded digits
        will be dropped, i.e. 0b00100 is turned into 0b100. Examples are
            bin(0) =    0b0
            bin(1) =    0b1
            bin(8) = 0b1000
        Each bit value (read from right to left) corresponds to the fact that a
        channel is on or off. I.e. if you have
            0b001011
        then it would mean that only channel 0, 1 and 3 are switched to on, the
        others are off.
        Helper method for write_pulse_form.
        """
        bits = 0     # that corresponds to: 0b0
        for channel in active_channels:
            # go through each list element and create the digital word out of
            # 0 and 1 that represents the channel configuration. In order to do
            # that a bitwise shift to the left (<< operator) is performed and
            # the current channel configuration is compared with a bitwise OR
            # to check whether the bit was already set. E.g.:
            #   0b1001 | 0b0110: compare elementwise:
            #           1 | 0 => 1
            #           0 | 1 => 1
            #           0 | 1 => 1
            #           1 | 1 => 1
            #                   => 0b1111
            bits = bits | (1<< channel)
        return bits


    def _load_block_from_file(self, block_name):
        """
        De-serializes a PulseBlock instance from file.

        @param str block_name: The name of the PulseBlock instance to de-serialize
        @return PulseBlock: The de-serialized PulseBlock instance
        """
        block = None
        filepath = os.path.join(self.pulsed_file_dir, '{0}.block'.format(block_name))
        if os.path.exists(filepath):
            try:
                with open(filepath, 'rb') as file:
                    block = pickle.load(file)
            except pickle.UnpicklingError:
                self.log.error('Failed to de-serialize PulseBlock "{0}" from file.'
                               ''.format(block_name))
                os.remove(filepath)
            except ModuleNotFoundError:
                self.log.error('Failed to de-serialize PulseBlock "{0}" from file because of missing dependencies.\n'
                               'For better debugging I dumped the traceback to debug.'.format(block_name))
                # self.log.debug('{0!s}'.format(traceback.format_exc()))
        return block

    def _load_ensemble_from_file(self, ensemble_name):
        """
        De-serializes a PulseBlockEnsemble instance from file.

        @param str ensemble_name: The name of the PulseBlockEnsemble instance to de-serialize
        @return PulseBlockEnsemble: The de-serialized PulseBlockEnsemble instance
        """
        ensemble = None
        filepath = os.path.join(self.pulsed_file_dir, '{0}.ensemble'.format(ensemble_name))
        if os.path.exists(filepath):
            try:
                with open(filepath, 'rb') as file:
                    ensemble = pickle.load(file)
            except pickle.UnpicklingError:
                self.log.error('Failed to de-serialize PulseBlockEnsemble "{0}" from file. '
                               'Deleting broken file.'.format(ensemble_name))
                os.remove(filepath)
        return ensemble

    def _load_sequence_from_file(self, sequence_name):
        """
        De-serializes a PulseSequence instance from file.

        @param str sequence_name: The name of the PulseSequence instance to de-serialize
        @return PulseSequence: The de-serialized PulseSequence instance
        """
        filepath = os.path.join(self.pulsed_file_dir, '{0}.sequence'.format(sequence_name))
        if os.path.exists(filepath):
            try:
                with open(filepath, 'rb') as file:
                    sequence = pickle.load(file)
                # FIXME: Due to the pickling the dict namespace merging gets lost on the way.
                # Restored it here but a better way needs to be found.
                for step in range(len(sequence)):
                    sequence[step].__dict__ = sequence[step]
            except pickle.UnpicklingError:
                self.log.error('Failed to de-serialize PulseSequence "{0}" from file.'
                               ''.format(sequence_name))
                os.remove(filepath)
                return None

        # Conversion for backwards compatibility
        if len(sequence) > 0 and not isinstance(sequence[0].flag_high, list):
            self.log.warning('Loading deprecated PulseSequence instances from disk. '
                             'Attempting conversion to new format.\nIf you keep getting this '
                             'message after reloading SequenceGeneratorLogic or restarting qudi, '
                             'please regenerate the affected PulseSequence "{0}".'
                             ''.format(sequence_name))
            for step_no, step_params in enumerate(sequence):
                # Try to convert "flag_high" step parameter
                if isinstance(step_params.flag_high, str):
                    if step_params.flag_high.upper() == 'OFF':
                        sequence[step_no].flag_high = list()
                    else:
                        sequence[step_no].flag_high = [step_params.flag_high]
                elif isinstance(step_params.flag_high, dict):
                    sequence[step_no].flag_high = [flag for flag, state in
                                                   step_params.flag_high.items() if state]
                else:
                    self.log.error('Failed to de-serialize PulseSequence "{0}" from file.'
                                   '"flag_high" step parameter is of unknown type'
                                   ''.format(sequence_name))
                    os.remove(filepath)
                    return None

                # Try to convert "flag_trigger" step parameter
                if isinstance(step_params.flag_trigger, str):
                    if step_params.flag_trigger.upper() == 'OFF':
                        sequence[step_no].flag_trigger = list()
                    else:
                        sequence[step_no].flag_trigger = [step_params.flag_trigger]
                elif isinstance(step_params.flag_trigger, dict):
                    sequence[step_no].flag_trigger = [flag for flag, state in
                                                      step_params.flag_trigger.items() if state]
                else:
                    self.log.error('Failed to de-serialize PulseSequence "{0}" from file.'
                                   '"flag_trigger" step parameter is of unknown type'
                                   ''.format(sequence_name))
                    os.remove(filepath)
                    return None
            self._save_sequence_to_file(sequence)
        return sequence

    def _channel_to_index(self,ch):
        """
        Inputs channel string and outputs corresponding channel index integer

        @param str ch: string corresponding to generic channel number i.e. either 'd_chX' or 'a_chX',
                        where X corresponds to genericv channel number
        @return int: integer corresponding to physical channel index
        """
        # checks for string being digital and analogue channel name,
        # and gives the index of the physical channel, by subtracting 1
        if ch.startswith('d_ch'):
            return int(ch[-1])-1
        elif ch.startswith('a_ch'):
            return int(ch[-1])-1
        else:
            self.log.error('Channel not given or ill-defined')

    def _seq_to_block_list(self, seq):
        """
        Make sequence (from file) into list of blocknames

        @param str seq: name of sequence.
        @return list str: list of strings corresponding to names of blocks in the ensembles that the sequence is made of.
        """
        block_list = []
        seq = self._load_sequence_from_file(seq)
        for ens in seq.ensemble_list:
            for i in range(ens['repetitions'] + 1):
                for block in self._ens_to_block_list(ens['ensemble']):
                    block_list.append(block)
        return block_list

    def _ens_to_block_list(self, ens):
        """
        Make ensemble (from file) into list of blocknames

        @param str ens: name of ensemble.
        @return list str: list of strings corresponding to names of blocks in the ensemble.
        """
        block_list = []
        ens = self._load_ensemble_from_file(ens)
        for block, rep in ens.block_list:
            for i in range(rep + 1):
                block_list.append(block)
        return block_list

    def _block_list_to_msg_list(self, block_list):
        """
        Converts list of blocknames to list of elements that the blocks are constituted of.

        @param list str block_list: list of strings corresponding to names of blocks to load.
        @return list PulseBlockElement: list of elements of type "PulseBlockElements" (internal to qudi)
        """
        msg_list = []
        for item in block_list:
            block = self._load_block_from_file(item)
            for elem in block.element_list:
                msg_list.append(elem)
        return msg_list

    def _msg_list_to_ps_msg_list(self, msg_list):
        """
        Converts list of PulseBlockElements to list of PulseMessage methods.

        @param list PulseBlockElement msg_list: list of PulseBlockElements to be converted
        @return PulseMessage: list of PulseMessages
        """
        return [self._msg_to_ps_msg(elem) for elem in msg_list]
#        ps_list = []
#        for elem in msg_list:
#            # print(self._msg_to_ps_msg(elem))
#            ps_list.append(self._msg_to_ps_msg(elem))
#        return ps_list

    def _msg_to_ps_msg(self, elem):
        """
        Converts list of PulseBlockElements to list of PulseMessage methods.

        @param PulseBlockElement elem: PulseBlockElement to be converted
        @return PulseMessage: PulseMessage to be read into PulseStreamer
        """
        num_ticks = int(elem.init_length_s * 1e9)
        digi_chnls = []
        for key, value in elem.digital_high.items():
            if value:
                digi_chnls.append(self._channel_to_index(key))
        digi_chnls = self._convert_to_bitmask(digi_chnls)
        return pulse_streamer_pb2.PulseMessage(ticks=num_ticks, digi=digi_chnls, ao0=0, ao1=0)

    def _seq_to_ps_msg_list(self, seq):
        return self._msg_list_to_ps_msg_list(self._block_list_to_msg_list(self._seq_to_block_list(seq)))

    def _ens_to_ps_msg_list(self, ens):
        return self._msg_list_to_ps_msg_list(self._block_list_to_msg_list(self._ens_to_block_list(ens)))

    def _msg_list_to_ps_sequences(self, msg_list):
        """

        @param list msg_list: list of messages to be translated to the sequences that are PulseStreamer compliant.
        @return dict: Dictionary containing:
                            Keys str:       strings containing generic channel names
                            Values list:    list of tuples containing the states which are to be sent to pulsestreamer.
                                                Tuples consists of: (duration_time(ns) , digi_state(0/1))
        """
        #TODO need to do rounding of duration lengths correct... still problems regarding sequences...
        chnl_dict = {}
        for chnl_name in list(self.get_active_channels()):
            #initializasion:
            chnl_seq = []
            duration = msg_list[0].init_length_s
            #iterate over all but 0'th index, and compare boolean w. previous value:
            for index in range(1, len(msg_list)):
                if msg_list[index].digital_high[chnl_name] != msg_list[index - 1].digital_high[chnl_name]:
                    chnl_seq.append((int(round(duration, 9) * 1e9), msg_list[index - 1].digital_high[chnl_name]))
                    duration = msg_list[index].init_length_s
                else:
                    duration = duration + msg_list[index].init_length_s
                if index == (len(msg_list)-1):
                    chnl_seq.append((int(round(duration, 9) * 1e9), msg_list[index].digital_high[chnl_name])) # time in units ns
#                    chnl_seq.append((round(duration, 9), msg_list[index].digital_high[chnl_name])) # time in units s
            # cleaning up...
            if (len(chnl_seq) > 1) or ((len(chnl_seq) == 1) and chnl_seq[0][1]):
                chnl_dict[chnl_name] = chnl_seq
            # do rounding correction here instead of inside body of function?
        return chnl_dict

    def _ens_to_ps_sequences(self, ens):
        return self._msg_list_to_ps_sequences(self._block_list_to_msg_list(self._ens_to_block_list(ens)))

    def _seq_to_ps_sequences(self, seq):
        return self._msg_list_to_ps_sequences(self._block_list_to_msg_list(self._seq_to_block_list(seq)))







#        # cleaning up...
#        if (len(chnl_dict[chnl_name]) == 1) and (chnl_dict[chnl_name][0][1] == False):
#            del chnl_dict[chnl_name]





#    def _msg_list_to_ps_msg_list(self, msg_list):
#        ps_list = []
#        for elem in msg_list:
#            num_ticks = elem.init_length_s *1e9
#            digi_chnls = []
#            for key, value in elem.digital_high
#                if value = True:
#                    digi_chnls.append(_channel_to_index(key))
#            digi_chnls.sort()
#            ps_list.append(pulse_streamer_pb2.PulseMessage(ticks=num_ticks, digi=digi_chnls, ao0=0, ao1=0))
#        return ps_list











#    def _cat_ens_to_psseq(self, ens_name):
#        psseq = None
#        ens = self._load_ensemble_from_file(ens_name)
#        for block ,rep in ens.block_list:
#            for i in range(rep):
#                self._append_block(psseq,block)
#        return psseq
#
#    def _cat_blocks_to

#  ====== Old functions ======

    # def upload_asset(self, asset_name=None):
    #     """ Upload an already hardware conform file to the device.
    #         Does NOT load it into channels.
    #
    #     @param name: string, name of the ensemble/seqeunce to be uploaded
    #
    #     @return int: error code (0:OK, -1:error)
    #     """
    #     self.log.debug('PulseStreamer has no own storage capability.\n"upload_asset" call ignored.')
    #     return 0

    # def load_asset(self, asset_name, load_dict=None):
    #     """ Loads a sequence or waveform to the specified channel of the pulsing
    #         device.
    #
    #     @param str asset_name: The name of the asset to be loaded
    #
    #     @param dict load_dict:  a dictionary with keys being one of the
    #                             available channel numbers and items being the
    #                             name of the already sampled
    #                             waveform/sequence files.
    #                             Examples:   {1: rabi_Ch1, 2: rabi_Ch2}
    #                                         {1: rabi_Ch2, 2: rabi_Ch1}
    #                             This parameter is optional. If none is given
    #                             then the channel association is invoked from
    #                             the sequence generation,
    #                             i.e. the filename appendix (_Ch1, _Ch2 etc.)
    #
    #     @return int: error code (0:OK, -1:error)
    #     """
    #     # ignore if no asset_name is given
    #     if asset_name is None:
    #         self.log.warning('"load_asset" called with asset_name = None.')
    #         return 0
    #
    #     # check if asset exists
    #     saved_assets = self.get_saved_asset_names()
    #     if asset_name not in saved_assets:
    #         self.log.error('No asset with name "{0}" found for PulseStreamer.\n'
    #                        '"load_asset" call ignored.'.format(asset_name))
    #         return -1
    #
    #     # get samples from file
    #     filepath = os.path.join(self.host_waveform_directory, asset_name + '.pstream')
    #     pulse_sequence_raw = dill.load(open(filepath, 'rb'))
    #     # pulse_sequence_raw = samples_write_methods._write_pstream(data)
    #
    #     pulse_sequence = []
    #     for pulse in pulse_sequence_raw:
    #         pulse_sequence.append(pulse_streamer_pb2.PulseMessage(ticks=pulse[0], digi=pulse[1], ao0=0, ao1=1))
    #
    #     blank_pulse = pulse_streamer_pb2.PulseMessage(ticks=0, digi=0, ao0=0, ao1=0)
    #     laser_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=self._convert_to_bitmask([self._laser_channel]), ao0=0, ao1=0)
    #     laser_and_uw_channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
    #     laser_and_uw_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=laser_and_uw_channels, ao0=0, ao1=0)
    #     self._sequence = pulse_streamer_pb2.SequenceMessage(pulse=pulse_sequence, n_runs=0, initial=laser_on,
    #         final=laser_and_uw_on, underflow=blank_pulse, start=1)
    #
    #     self.current_loaded_asset = asset_name
    #     return 0

#     def has_sequence_mode(self):
#         """ Asks the pulse generator whether sequence mode exists.
#
#         @return: bool, True for yes, False for no.
#         """
#         return False

#     def tell(self, command):
#         """ Sends a command string to the device.
#
#         @param string command: string containing the command
#
#         @return int: error code (0:OK, -1:error)
#         """
#         return 0
#
#     def ask(self, question):
#         """ Asks the device a 'question' and receive and return an answer from it.
#
#         @param string question: string containing the command
#
#         @return string: the answer of the device to the 'question' in a string
#         """
#         return ''


    # def set_asset_dir_on_device(self, dir_path):
    #     """ Change the directory where the assets are stored on the device.
    #
    #     @param str dir_path: The target directory
    #
    #     @return int: error code (0:OK, -1:error)
    #
    #     Unused for digital pulse generators without changeable file structure
    #     (PulseBlaster, FPGA).
    #     """
    #     return 0

    # def get_asset_dir_on_device(self):
    #     """ Ask for the directory where the hardware conform files are stored on
    #         the device.
    #
    #     @return str: The current file directory
    #
    #     Unused for digital pulse generators without changeable file structure
    #     (PulseBlaster, FPGA).
    #     """
    #     return ''

    # def get_uploaded_asset_names(self):
    #     """ Retrieve the names of all uploaded assets on the device.
    #
    #     @return list: List of all uploaded asset name strings in the current
    #                   device directory. This is no list of the file names.
    #
    #     Unused for digital pulse generators without sequence storage capability
    #     (PulseBlaster, FPGA).
    #     """
    #     names = []
    #     return names

    # def get_saved_asset_names(self):
    #     """ Retrieve the names of all sampled and saved assets on the host PC.
    #     This is no list of the file names.
    #
    #     @return list: List of all saved asset name strings in the current
    #                   directory of the host PC.
    #     """
    #     file_list = self._get_filenames_on_host()
    #
    #     saved_assets = []
    #     for filename in file_list:
    #         if filename.endswith('.pstream'):
    #             asset_name = filename.rsplit('.', 1)[0]
    #             if asset_name not in saved_assets:
    #                 saved_assets.append(asset_name)
    #     return saved_assets

    # def delete_asset(self, asset_name):
    #     """ Delete all files associated with an asset with the passed asset_name from the device memory.
    #
    #     @param str asset_name: The name of the asset to be deleted
    #                            Optionally a list of asset names can be passed.
    #
    #     @return int: error code (0:OK, -1:error)
    #
    #     Unused for digital pulse generators without sequence storage capability
    #     (PulseBlaster, FPGA).
    #     """
    #     return 0

    #    def _get_filenames_on_host(self):
    #        """ Get the full filenames of all assets saved on the host PC.
    #
    #        @return: list, The full filenames of all assets saved on the host PC.
    #        """
    #        filename_list = [f for f in os.listdir(self.host_waveform_directory) if f.endswith('.pstream')]
    #        return filename_list