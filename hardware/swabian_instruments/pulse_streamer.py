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

import numpy as np
import time
import pickle
import grpc
import os
import hardware.swabian_instruments.pulse_streamer_pb2 as pulse_streamer_pb2
#from pulsestreamer import PulseStreamer
#import pulsestreamer
#from pulsestreamer import TriggerStart, TriggerRearm
from pulsestreamer import Sequence, OutputState
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

        # if 'pulsed_file_dir' in config.keys():
        #     self.pulsed_file_dir = config['pulsed_file_dir']
        #
        #     if not os.path.exists(self.pulsed_file_dir):
        #         homedir = get_home_dir()
        #         self.pulsed_file_dir = os.path.join(homedir, 'pulsed_files')
        #         self.log.warning('The directory defined in parameter '
        #                     '"pulsed_file_dir" in the config for '
        #                     'PulseStreamer does not exist!\n'
        #                     'The default home directory\n{0}\n will be taken '
        #                     'instead.'.format(self.pulsed_file_dir))
        # else:
        #     homedir = get_home_dir()
        #     self.pulsed_file_dir = os.path.join(homedir, 'pulsed_files')
        #     self.log.warning('No parameter "pulsed_file_dir" was specified in the config for '
        #                      'PulseStreamer as directory for the pulsed files!\nThe default home '
        #                      'directory\n{0}\nwill be taken instead.'.format(self.pulsed_file_dir))

        # self.host_waveform_directory = self._get_dir_for_name('sampled_hardware_files')

        self.current_status = -1
        self.sample_rate = 1e9
        self.waveform_names = []
        self.ps_waveforms_dict = {}
        self.to_be_streamed = None
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
        self.pulse_streamer.stream(self.to_be_streamed)
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

        waveforms = list()

        # Sanity checks
        if len(analog_samples) > 0:
            number_of_samples = len(analog_samples[list(analog_samples)[0]])
        elif len(digital_samples) > 0:
            number_of_samples = len(digital_samples[list(digital_samples)[0]])
        else:
            self.log.error('No analog or digital samples passed to write_waveform method in dummy '
                           'pulser.')
            return -1, waveforms

        for chnl, samples in analog_samples.items():
            if len(samples) != number_of_samples:
                self.log.error('Unequal length of sample arrays for different channels in dummy '
                               'pulser.')
                return -1, waveforms
        for chnl, samples in digital_samples.items():
            if len(samples) != number_of_samples:
                self.log.error('Unequal length of sample arrays for different channels in dummy '
                               'pulser.')
                return -1, waveforms

#         self._convert_digi_samples_to_quick_blocks(samples) # loads humanreadable ps compliant sequence into ps_waveform_dict
#         self._array_to_ps_sequence(samples)

        for chnl, samples in digital_samples.items():
            waveform_name = name + chnl[1:]
            waveforms.append(waveform_name)
            if waveform_name not in self.ps_waveforms_dict.keys():
                self.ps_waveforms_dict[waveform_name] = samples #samples are still in raw format
            # if waveform_name not in self.waveform_channel_names:
            #     self.waveform_channel_names.append(waveform_name) # waveform_channel_names not used?

        if name not in self.waveform_names:
            self.waveform_names.append(name)

        temp_dict = {}
        for wf_name, ps_seq in self.ps_waveforms_dict.items():
            temp_dict[self._channel_to_index(wf_name) + 1] = wf_name

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
        # print("name: " + str(name))
        # print("sequence_parameters: " + str(sequence_parameters))

        sequence = []
        for channels, dict in sequence_parameters:
            ensemble = dict['ensemble']
            # print(ensemble)
            if ensemble not in self.ps_waveforms_dict.keys():
                self.log.error('Waveform {0} not written to device. Write waveforms present in sequence to device, before proceeding'.format(ensemble))
                return -1
            reps = dict['repetitions']
            for i in range(reps + 1):
                for j in self.ps_waveform[ensemble]:
                    sequence.append(j)
                self.ps_sequence = sequence

        # self.loaded_waveforms.append(name)
        return len(sequence_parameters)

        # name: dummy_seq
        # sequence_parameters: [(('dummy_ens_ch1', 'dummy_ens_ch2', 'dummy_ens_ch3', 'dummy_ens_ch4', 'dummy_ens_ch5',
        #                         'dummy_ens_ch6', 'dummy_ens_ch7', 'dummy_ens_ch8'),
        #                        {'ensemble': 'dummy_ens', 'repetitions': 1, 'go_to': -1, 'event_jump_to': -1,
        #                         'event_trigger': 'OFF', 'wait_for': 'OFF', 'flag_trigger': [], 'flag_high': []}), ((
        #                                                                                                            'dummy_ens_ch1',
        #                                                                                                            'dummy_ens_ch2',
        #                                                                                                            'dummy_ens_ch3',
        #                                                                                                            'dummy_ens_ch4',
        #                                                                                                            'dummy_ens_ch5',
        #                                                                                                            'dummy_ens_ch6',
        #                                                                                                            'dummy_ens_ch7',
        #                                                                                                            'dummy_ens_ch8'),
        #                                                                                                            {
        #                                                                                                                'ensemble': 'dummy_ens',
        #                                                                                                                'repetitions': 0,
        #                                                                                                                'go_to': -1,
        #                                                                                                                'event_jump_to': -1,
        #                                                                                                                'event_trigger': 'OFF',
        #                                                                                                                'wait_for': 'OFF',
        #                                                                                                                'flag_trigger': [],
        #                                                                                                                'flag_high': []})]

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
        print(load_dict)
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

        # converting dictionary to fir format of internal functions:
        formatted_load_dict = {}
        for chnl_num, chnl_wf_name in load_dict:
            formatted_load_dict['d_ch' + chnl_num] = self.ps_waveforms_dict[chnl_wf_name]

        # writing list of pulse_messages
        msgs = []
        for duration, bitmask in self._convert_digi_samples_to_quick_blocks(formatted_load_dict):
            msgs.append(pulse_streamer_pb2.PulseMessage(ticks=duration, digi=bitmask, ao0=0, ao1=0))

        laser_channel = self._convert_to_bitmask([self._laser_channel])
        laser_and_uw_channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])

        blank_pulse = pulse_streamer_pb2.PulseMessage(ticks=0, digi=0, ao0=0, ao1=0)
        laser_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=laser_channel, ao0=0, ao1=0)
        laser_and_uw_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=laser_and_uw_channels, ao0=0, ao1=0)
        # actually loads waveform into PS
        self.ps_waveform = pulse_streamer_pb2.SequenceMessage(pulse=msgs, n_runs=0, initial=laser_on,
                                                              final=laser_and_uw_on, underflow=blank_pulse, start=1)

        self.current_loaded_asset = load_dict, 'waveform'
        return self.current_loaded_asset


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

        # Check is sequence already exists on PulseStreamer Hardware, and clears it.
        if self.pulse_streamer.hasSequence():
            self.reset()
        self.to_be_streamed = self.pulse_streamer.createSequence()
        return_dict = {}
        for seq_chnl_name, ps_seq in self.seq_master_dict[sequence_name].items():
            self.to_be_streamed.setDigital(self._channel_to_index(seq_chnl_name), ps_seq)
            return_dict[self._channel_to_index(seq_chnl_name) + 1] = seq_chnl_name

        pulse_sequence_raw = self.to_be_streamed.getData()

        pulse_sequence = []
        for pulse in pulse_sequence_raw:
            pulse_sequence.append(pulse_streamer_pb2.PulseMessage(ticks=pulse[0], digi=pulse[1], ao0=pulse[2], ao1=pulse[3]))

        blank_pulse = pulse_streamer_pb2.PulseMessage(ticks=0, digi=0, ao0=0, ao1=0)
        laser_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=self._convert_to_bitmask([self._laser_channel]), ao0=0,
                                                   ao1=0)
        laser_and_uw_channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        laser_and_uw_on = pulse_streamer_pb2.PulseMessage(ticks=0, digi=laser_and_uw_channels, ao0=0, ao1=0)
        self._sequence = pulse_streamer_pb2.SequenceMessage(pulse=pulse_sequence, n_runs=0, initial=laser_on,
                                                            final=laser_and_uw_on, underflow=blank_pulse, start=1)
        #self.current_loaded_asset = sequence_name, 'sequence'

        # update, confirmation and return
        if self.pulse_streamer.hasSequence():
            # self.to_be_streamed.plot()
            self.current_loaded_asset = return_dict, 'sequence'
            return self.current_loaded_asset
        else:
            self.log.error('No sequence with name "{0}" uploaded to PulseStreamer.\n'
                           '"load_sequence" call ignored.'.format(sequence_name))

        return sequence_name

#     sequencegeneratorlogic.sample_pulse_sequence('dummy_seq')
#     dummy_seq[(('dummy_ens_ch1', 'dummy_ens_ch2', 'dummy_ens_ch3', 'dummy_ens_ch4', 'dummy_ens_ch5', 'dummy_ens_ch6',
#                 'dummy_ens_ch7', 'dummy_ens_ch8'),
#                {'ensemble': 'dummy_ens', 'repetitions': 1, 'go_to': -1, 'event_jump_to': -1, 'event_trigger': 'OFF',
#                 'wait_for': 'OFF', 'flag_trigger': [], 'flag_high': []}), (('dummy_ens_ch1', 'dummy_ens_ch2',
#                                                                             'dummy_ens_ch3', 'dummy_ens_ch4',
#                                                                             'dummy_ens_ch5', 'dummy_ens_ch6',
#                                                                             'dummy_ens_ch7', 'dummy_ens_ch8'),
#                                                                            {'ensemble': 'dummy_ens', 'repetitions': 0,
#                                                                             'go_to': -1, 'event_jump_to': -1,
#                                                                             'event_trigger': 'OFF', 'wait_for': 'OFF',
#                                                                             'flag_trigger': [], 'flag_high': []})]


    def get_waveform_names(self):
        """ Retrieve the names of all uploaded waveforms on the device.

        @return list: List of all uploaded waveform name strings in the device workspace.
        """

        dict, format = self.current_loaded_asset
        if format == 'waveform':
            for key, val in dict.items():
            return list(dict.values())
        elif format == 'sequence':
            pass
            # print(self.current_loaded_asset) # for key, val in dict.items():
        else:
            return []


    def get_sequence_names(self):
        """ Retrieve the names of all uploaded sequence on the device.

        @return list: List of all uploaded sequence name strings in the device workspace.
        """

        #TODO: not entirely right...i guess...
        dict, format = self.current_loaded_asset
        # print(self.current_loaded_asset)
        if format == 'sequence':
            for key, val in dict.items():
                return [val]
        else:
            return []


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

    def _convert_digi_samples_to_quick_blocks(self, digi_samps_dict):
        """

        :param digi_samps_dict: Dictionary of format {'d_ch1': array([0,1,1,1,0,0,1,0,1,0,0,...,0,1]) , ...}
        :return:
        """
        for key in digi_samps_dict.keys():
            temp_dict = [int(b) for b in digi_samps_dict[key]]
            digi_samps_dict[key] = temp_dict

        bitwise_matrix = np.column_stack((digi_samps_dict['d_ch1'], digi_samps_dict['d_ch2'],
                                          digi_samps_dict['d_ch3'], digi_samps_dict['d_ch4'],
                                          digi_samps_dict['d_ch5'], digi_samps_dict['d_ch6'],
                                          digi_samps_dict['d_ch7'], digi_samps_dict['d_ch8']))
        return self._bitarray_to_ps_seq(np.packbits(bitwise_matrix, bitorder='little'))

    def _bitarray_to_ps_seq(self, array):
        twin_list = [(ind, array[ind-1], array[ind]) for ind in range(len(array[1:])) if array[ind-1] != array[ind]]
        seq = [(i, j) for i,j,k in twin_list[1:]]
        seq.append((len(array)-twin_list[-1][0] ,twin_list[-1][2])) # last element
        return seq

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
            return int(ch[-1]) - 1
        elif ch.startswith('a_ch'):
            return int(ch[-1]) - 1
        elif ch[-4:-1] == '_ch':
            return int(ch[-1]) - 1
        elif ch in self.ps_sequence_dict.keys():
            return int(ch[-1]) - 1
        else:
            self.log.error('Channel not given or ill-defined')






    # def _array_to_ps_sequence(self, array):
    #     t0 = time.clock()
    #     seq = []
    #     state = [array[0]]
    #     for ind, i in enumerate(array[1:]):
    #         # print(ind, len(array[1:]))
    #         if i == state[-1]:
    #             state.append(i)
    #         else:
    #             seq.append((len(state), state[0]))
    #             state = [i]
    #         if ind == len(array[1:]) - 1:
    #             seq.append((len(state), state[0]))
    #     print(time.clock() - t0)
    #     return seq
    #
    # def bitarray_to_ps_seq(self, array):
    #     t0 = time.clock()
    #     seq = []
    #     state = array[0]
    #     count = 1
    #     for ind, i in enumerate(array[1:]):
    #         if i == state:
    #             count += 1
    #         else:
    #             seq.append((count, state))
    #             state = i
    #             count = 1
    #         if ind == len(array[1:]) - 1:
    #             seq.append((count, state))
    #     print(time.clock() - t0)
    #     return seq
    #
    # def bitarray_to_ps_seq_v2(self, array):
    #     #t0 = time.clock()
    #     twin_list = [(ind, array[ind-1], array[ind]) for ind in range(len(array[1:])) if array[ind-1] != array[ind]]
    #     #print(twin_list)
    #     seq = [(i, j) for i,j,k in twin_list[1:]]
    #     seq.append((len(array)-twin_list[-1][0] ,twin_list[-1][2])) # last element
    #     #print(time.clock() - t0)
    #     return seq
    #




    # belongs in write_waveform:

    # msgs = []
    # chnl_states = []
    # chnl_states_new = []
    # chnl_on_list = []
    # duration = 0
    # for i in range(total_number_of_samples):
    #     duration = duration + 1
    #     diff = False
    #     for j in range(1, 9):
    #         key = 'd_ch' + str(j)
    #         if i == 0:
    #             chnl_states.append(digital_samples[key][i])
    #         elif digital_samples[key][i] != digital_samples[key][i-1]:
    #             diff = True
    #             chnl_states_new.append(digital_samples[key][i])
    #             if j == 8 and diff:
    #                 #update list of channelnumbers, where channel state is on/"True"
    #                 for num, state in enumerate(chnl_states):
    #                     if state:
    #                         chnl_on_list.append(num)
    #                 msgs.append(pulse_streamer_pb2.PulseMessage(ticks=duration,
    #                                                             digi=self._convert_to_bitmask(chnl_on_list),
    #                                                             ao0=0, ao1=0))
    #                 duration = 0
    #                 chnl_states = chnl_states_new
    #         if i == total_number_of_samples - 1:
    #             msgs.append(pulse_streamer_pb2.PulseMessage(ticks=duration,
    #                                                         digi=self._convert_to_bitmask(chnl_on_list),
    #                                                         ao0=0, ao1=0))





    #        pulse_waveform = []
    #        for pulse in pulse_waveform_raw:
    #            pulse_waveform.append(pulse_streamer_pb2.PulseMessage(ticks=pulse[0], digi=pulse[1], ao0=0, ao1=1))
    #
    #        self._waveform = pulse_streamer_pb2.SequenceMessage(pulse=pulse_waveform, n_runs=0, initial=laser_on,
    #                                                            final=laser_and_uw_on, underflow=blank_pulse, start=1)

    #       make list of loaded waveform in pulse streamer format to help load them into ps from write_sequence