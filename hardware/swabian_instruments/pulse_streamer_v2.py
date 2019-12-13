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

# import JSON-RPC Pulse Streamer wrapper class, to use Google-RPC import from pulsestreamer.grpc
#from pulsestreamer import PulseStreamer
import pulsestreamer
# import enum types
from pulsestreamer import TriggerStart, TriggerRearm
# import class Sequence and OutputState for advanced sequence building
from pulsestreamer import Sequence, OutputState

import numpy as np
# import time
import pickle
# import grpc
import os
#import hardware.swabian_instruments.pulse_streamer_pb2 as pulse_streamer_pb2
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

        self.current_status = -1
        self.sample_rate = 1e9
        self.to_be_streamed = None # sequence to be streamed like: self.pulsestreamer.stream(self.to_be_streamed)

        # TODO: This can probably be done MUCH cleaner!!!
        self.waveform_names = []
        self.waveform_channel_names = []
        self.ps_waveforms_dict = {} # Dict containing human-readable elements of the waveform eg. {(dummy_ens_ch1:[(1,10), (0,10), ...]), ...}

        self.sequence_names = []
        self.sequence_channel_names = []
        self.ps_sequence_dict = {} # Dict containing human-readable elements of the sequence eg. {(dummy_ens_ch1:[(1,10), (0,10), ...]), ...}
        self.seq_master_dict = {} # dict with keys being sequencenames, and values being
                                # dict of corresponding sequence channelnames as keys with values being the humanreadable sequence.
        # ============================

        self.to_be_streamed = None # current uploaded waveform/sequence in ps object format.
        self.current_loaded_asset = {}, None

        self.pulse_streamer = pulsestreamer.PulseStreamer(self._pulsestreamer_ip)

    def on_activate(self):
        """ Establish connection to pulse streamer and tell it to cancel all operations """
#        self.pulse_streamer = pulse_streamer_pb2.PulseStreamerStub(self._channel)
        self.pulse_streamer = pulsestreamer.PulseStreamer(self._pulsestreamer_ip)
        start = TriggerStart.SOFTWARE
        rearm = TriggerRearm.AUTO
        self.pulse_streamer.setTrigger(start=start, rearm=rearm)

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
        activation_config['digital'] = frozenset({'d_ch1', 'd_ch2', 'd_ch3', 'd_ch4',
                                                  'd_ch5', 'd_ch6', 'd_ch7', 'd_ch8'})
        activation_config['analogue'] = frozenset({'a_ch1', 'a_ch2'})
        activation_config['all'] = frozenset({'d_ch1', 'd_ch2', 'd_ch3', 'd_ch4',
                                              'd_ch5', 'd_ch6', 'd_ch7', 'd_ch8',
                                              'a_ch1', 'a_ch2'})
        constraints.activation_config = activation_config

        return constraints


    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:OK, -1:error)
        """
        # start the pulse sequence
        # self.pulse_streamer.stream(self.to_be_streamed)
        self.log.info('Asset uploaded to PulseStreamer')
        self.pulse_streamer.startNow()
        self.current_status = 1
        return 0

    def pulser_off(self):
        """ Switches the pulsing device off.

        @return int: error code (0:OK, -1:error)
        """
        # stop the pulse sequence
        channels = [self._laser_channel, self._uw_x_channel]
        self.pulse_streamer.constant(OutputState(channels, 0, 0))
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

    def get_waveform_names(self):
        """ Retrieve the names of all uploaded waveforms on the device.

        @return list: List of all uploaded waveform name strings in the device workspace.
        """
        return self.waveform_channel_names

    def get_sequence_names(self):
        """ Retrieve the names of all uploaded sequence on the device.

        @return list: List of all uploaded sequence name strings in the device workspace.
        """
        return self.sequence_names

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

        print(name, digital_samples, total_number_of_samples)

        waveforms = list()

        # Sanity checks
        if len(analog_samples) > 0:
            number_of_samples = len(analog_samples[list(analog_samples)[0]])
        elif len(digital_samples) > 0:
            number_of_samples = len(digital_samples[list(digital_samples)[0]])
        else:
            self.log.error('No analog or digital samples passed to write_waveform method in pulsestreamer.')
            return -1, waveforms
        for chnl, samples in analog_samples.items():
            if len(samples) != number_of_samples:
                self.log.error('Unequal length of analogue sample arrays for different channels in pulsestreamer.')
                return -1, waveforms
        for chnl, samples in digital_samples.items():
            if len(samples) != number_of_samples:
                self.log.error('Unequal length of digital sample arrays for different channels in pulsestreamer.')
                return -1, waveforms

        for chnl, samples in digital_samples.items():
            waveform_name = name + chnl[1:]
            waveforms.append(waveform_name)
            if waveform_name not in self.ps_waveforms_dict.keys():
                self.ps_waveforms_dict[waveform_name] = self._array_to_ps_sequence(samples)
            if waveform_name not in self.waveform_channel_names:
                self.waveform_channel_names.append(waveform_name)

        if name not in self.waveform_names:
            self.waveform_names.append(name)

        # setDigital method here!! Nope. Only add it the the "script memory",
        # since the upload to the channels happens in "load_dict"
        temp_dict = {}
        for wf_name, ps_seq in self.ps_waveforms_dict.items():
            temp_dict[self._channel_to_index(wf_name)+1] = wf_name

        # self.sequence_names = []
        self.current_loaded_asset = temp_dict, 'waveform'

        return number_of_samples, waveforms


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
        # print(name, sequence_parameters)
        # constructing the ps sequence elements (human-readable to be loaded into the channels)
        init = True
        for channels, dict in sequence_parameters:
            # initializes ps_sequence_dict
            if init:
                init = False
                for chnl in channels:
                    sequence_channel_name = name + '_' + str(self._channel_to_index(chnl) + 1)
                    self.ps_sequence_dict[sequence_channel_name] = []
            ensemble = dict['ensemble']
            #Checks that ensemble has been read in as waveform.
            if ensemble not in self.waveform_names:
                self.log.error(
                    'Waveform {0} not written to device. Write waveforms present in sequence to device, before proceeding'.format(
                        ensemble))
                return -1
            reps = dict['repetitions']
            for chnl in channels:
                sequence_channel_name = name + '_' + str(self._channel_to_index(chnl) + 1)
                for i in range(reps + 1):
                    self.ps_sequence_dict[sequence_channel_name] = self.ps_sequence_dict[sequence_channel_name] + self.ps_waveforms_dict[chnl]
        self.seq_master_dict[name] = self.ps_sequence_dict

        # updating variable wrt. what is "loaded" into PS.
        temp_dict = {}
        for seq_name, ps_seq in self.ps_sequence_dict.items():
            temp_dict[self._channel_to_index(seq_name) + 1] = seq_name

        # self.waveform_channel_names = []
        self.sequence_names.append(name)
        self.sequence_channel_names = [s for s in self.ps_sequence_dict.keys()]
        self.current_loaded_asset = temp_dict, 'sequence'

        return len(sequence_parameters)


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

        if not isinstance(load_dict, (list,dict)):
            self.log.error('Input object of type either Dict or List as input.')
            return -1

        # check if type of load_dict is list, and correct to dict format.
        if type(load_dict) == list:
            temp_dict = {}
            for i in load_dict:
                temp_dict[self._channel_to_index(i)+1] = i
            load_dict = temp_dict

        # check if asset exists, and write to Pulsestreamer if it does.
        self.to_be_streamed = self.pulse_streamer.createSequence()
        for channel, waveform in load_dict.items():
            if waveform not in self.get_waveform_names():
                self.log.error('No waveform with name "{0}" found for PulseStreamer.\n'
                                '"load_waveform" call ignored.'.format(load_dict))
                return -1
            self.to_be_streamed.setDigital(channel-1, self.ps_waveforms_dict[waveform])

        # load process goes here:
        # defining pulse structures and respective channels and actually loads sequence into channels.
        blank_pulse = OutputState.ZERO()
        laser_on = OutputState([self._laser_channel], 0, 0)
        laser_and_uw_on = OutputState([self._laser_channel, self._uw_x_channel], 0, 0)
        self.pulse_streamer.stream(self.to_be_streamed, n_runs=0, final=laser_and_uw_on)

        if self.pulse_streamer.hasSequence():
            # self.to_be_streamed.plot()
            self.current_loaded_asset = load_dict, 'waveform'
            return self.current_loaded_asset
        else:
            self.log.error('No waveforms uploaded to PulseStreamer.\n'
                           '"load_waveform" call ignored.')


    def load_sequence(self, sequence_name):
        """ Loads a sequence to the channels of the device in order to be ready for playback.
        For devices that have a workspace (i.e. AWG) this will load the sequence from the device
        workspace into the channels.

        @param sequence_name:  str, name of the sequence to load

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel, string describing the asset
                             type ('waveform' or 'sequence')
        """
        # Sanity checks:
        # ignore if no asset_name is given
        if sequence_name is None:
            self.log.warning('"load_sequence" called with sequence_name = None.')
            return {}
        # check if asset exists
        saved_assets = self.get_sequence_names()
        if sequence_name not in saved_assets:
            self.log.error('No sequence with name "{0}" found for PulseStreamer.\n'
                           '"load_sequence" call ignored.'.format(sequence_name))
            return -1

        # Check is sequence already exists on PulseStreamer Hardware, and clears it.
        if self.pulse_streamer.hasSequence():
            self.reset()
        self.to_be_streamed = self.pulse_streamer.createSequence()
        return_dict = {}
        for seq_chnl_name, ps_seq in self.seq_master_dict[sequence_name].items():
            self.to_be_streamed.setDigital(self._channel_to_index(seq_chnl_name), ps_seq)
            return_dict[self._channel_to_index(seq_chnl_name)+1] = seq_chnl_name

        # defining pulse structures and respective channels and actually loads sequence into channels.
        blank_pulse = OutputState.ZERO()
        laser_on = OutputState([self._laser_channel], 0, 0)
        laser_and_uw_on = OutputState([self._laser_channel, self._uw_x_channel], 0, 0)
        self.pulse_streamer.stream(self.to_be_streamed, n_runs=0, final=laser_and_uw_on)

        # update, confirmation and return
        if self.pulse_streamer.hasSequence():
            # self.to_be_streamed.plot()
            self.current_loaded_asset = return_dict, 'sequence'
            return self.current_loaded_asset
        else:
            self.log.error('No sequence with name "{0}" uploaded to PulseStreamer.\n'
                           '"load_sequence" call ignored.'.format(sequence_name))


    def clear_all(self):
        """ Clears all loaded waveforms from the pulse generators RAM.

        @return int: error code (0:OK, -1:error)

        Unused for digital pulse generators without storage capability
        (PulseBlaster, FPGA).
        """
        return 0
        # if not self.pulse_streamer.hasSequence():
        #     return -1
        # else:
        #     return 0

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
#        channels = [self._laser_channel, self._uw_x_channel]
        self.pulse_streamer.reset()
        start = TriggerStart.SOFTWARE
        rearm = TriggerRearm.AUTO
        self.pulse_streamer.setTrigger(start=start, rearm=rearm)
        # self.pulse_streamer.constant(OutputState([self._laser_channel], 0, 0))
        return 0

#  ====== Internal functions ======

    def _channel_to_index(self, ch):
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
        elif ch[-4:-1] == '_ch':
            return int(ch[-1]) - 1
        elif ch in self.ps_sequence_dict.keys():
            return int(ch[-1]) - 1
        else:
            self.log.error('Channel not given or ill-defined')


    def _array_to_ps_sequence(self, array):
        """
        Transform np array into pulse streamer sequence format

        @param array: nparray containing the state of the digital channel of pullse streamer pr ns.
        @return list: list consisting of tuples:
                        first element in tuple is duration in nsec
                        second element in tuple being the (digital) state of the channel
        """
        seq = []
        state = [array[0]]
        for ind, i in enumerate(array[1:]):
            if i == state[-1]:
                state.append(i)
            else:
                seq.append((len(state), state[0]))
                state = [i]
            if ind == len(array[1:]) - 1:
                seq.append((len(state), state[0]))
        _seq = []
        for dur, bool in seq:
            _seq.append((dur, int(bool)))
        return _seq

    # def _ps_seq_to_array(self, ps_seq):
    #     """
    #     Transform np array into pulse streamer sequence format
    #
    #     @param list ps_seq: list consisting of tuples:
    #                         first element in tuple is duration in nsec
    #                         second element in tuple being the (digital) state of the channel
    #
    #     @return array:  nparray containing the state of the digital channel of pullse streamer pr ns.
    #     """
    #     ret_list = []
    #     for duration, digi_state in ps_seq:
    #         ret_list = ret_list + [digi_state]*duration
    #
    #     return np.asarray(ret_list)


    def _set_up_osc(self, chnl_num, freq):
        """

        :param int chnl_num:
        :param float freq:

        :return int: error code (0:OK, -1:error)
        """

        # convert freq. to a "ON/OFF" single cycle in the PS list-tuple format, and upload it to corresponding channel.
        cycle_duration = int(1e9/freq) # in picosec
        pattern = [(cycle_duration/2, 1), (cycle_duration - cycle_duration/2, 0)]
        wave = self.pulse_streamer.createSequence()
        wave.setDigital(chnl_num, pattern)
        self.pulse_streamer.stream(wave, n_runs=self.pulse_streamer.REPEAT_INFINITELY) #, final=laser_and_uw_on)

        return 0

    def set_constatnt_state(self, channels, a1=0, a2=0):
        """ Switches the pulsing device off.

        @return int: error code (0:OK, -1:error)
        """
        # stop the pulse sequence
        self.pulse_streamer.constant(OutputState(channels, a1, a2))
        self.current_status = 0
        return 0