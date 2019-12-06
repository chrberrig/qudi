# -*- coding: utf-8 -*-

"""
This file contains the Qudi Logic module base class.

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

from qtpy import QtCore
from collections import OrderedDict
from interface.microwave_interface import MicrowaveMode
from interface.microwave_interface import TriggerEdge
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt

from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex
from core.connector import Connector
from core.configoption import ConfigOption
from core.statusvariable import StatusVar


class T1Logic(GenericLogic):
    """This is the Logic class for ODMR."""

    # declare connectors
    laser = Connector(interface='SimpleLaserInterface')
    fastcounter = Connector(interface='FastCounterInterface')
    pulser = Connector(interface='PulserInterface')
    # mw_generator = Connector(interface='MicrowaveInterface')
    fitlogic = Connector(interface='FitLogic')
    savelogic = Connector(interface='SaveLogic')

    # config option
    # clock_frequency = StatusVar('clock_frequency', 200)
    default_repeat = 200
    init_duration = 1000
    counter_delay = init_duration/2
    interleave_duration = 1000 # tau
    collection_duration = 200
    intersequence_delay = 10000
    # the above are all in units of ns.


    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)
        self.threadlock = Mutex()

    def on_activate(self):
        """
        Initialisation performed during activation of the module.
        """
        # Get connectors
        self.laser_hw = self.laser()
        self.fastcounter_hw = self.fastcounter()
        self.pulser_hw = self.pulser()
        # self.mw_gen_hw = self.mw_generator()
        self._fit_logic = self.fitlogic()
        self._save_logic = self.savelogic()

        self._laser_channel = self.pulser_hw._laser_channel
        self._trigger_channel = self.fastcounter_hw._channel_detect
        self._counter_channel = self.fastcounter_hw._channel_apd

        self.set_up_laser()


    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        # Stop measurement if it is still running
        if self.module_state() == 'locked':
            self.stop_measurement()
        timeout = 30.0
        start_time = time.time()
        while self.module_state() == 'locked':
            time.sleep(0.5)
            timeout -= (time.time() - start_time)
            if timeout <= 0.0:
                self.log.error('Failed to properly deactivate T1 measurement logic.' 
                                    'T1 measurement is still running but can not be stopped after 30 sec.')
                break
        # Switch off hardware
        self.laser_hw.off()

    def set_up_laser(self):
        self.laser_hw.on()

    # def set_up_pulser(self):
    #     pass

    def set_up_counter(self, bin_width_s=None, record_length_s=None, number_of_gates=0):
        """
        setting up fastcounter hardware as counter for photons.
        :return: bin_width_s, record_length_s, number_of_gates
        """
        #return self.fastcounter_hw.configure(self, bin_width_s=self.collection_duration, record_length_s=self.collection_duration, number_of_gates=0)

        if bin_width_s == None:
            bin_width_s = self.collection_duration
        if record_length_s == None:
            record_length_s = self.collection_duration

        return self.fastcounter_hw.configure(self, bin_width_s=bin_width_s, record_length_s=record_length_s, number_of_gates=number_of_gates)

    def set_up_mw_gen(self):
        pass


    def start_measurement(self):

        self.pulser_hw.pulser_on()
        return None


    def stop_measurement(self):
        self.pulser_hw.pulser_off()
        return None


    def generate_measurement_sequences(self, parameters, time_sweep, repeat=0, sequence_name='T1_measurement_sequence'):
        """
        Genarates pulse blocks and combine them to pulse sequence
        :param parameters:  tuple containing:
                            (initializasion_duration_in_ns, counter_delay_in_ns, interleave_duration or tau_in_ns, collection_duration_in_ns, intersequence_delay_in_ns) in units of ns
        :param time_sweep:  tuple containing:
                            (increment_size, number_of_increments)
        :param repeat:      Number of repeats for each choice of interleave-time
        :return str:        sequence name of created sequence.
        """
        if repeat == 0:
            repeat = self.default_repeat

        # defaultparameters = (self.init_duration, self.counter_delay, self.interleave_duration, self.collection_duration, self.intersequence_delay)
        # parameters = tuple(map(lambda x, y: y if y is not None else x, defaultparameters, parameters))
        if self.counter_delay > self.init_duration:
            self.log.warn('APD is not activated during initializasion')
            self.counter_delay = 0

        # crate sequence
        # TODO: This is pulsestreamer specific. This should of course be generalized to all pulsers... This might be doable using the pulser gui, looking for a specific file etc...
        # TODO: What is pulsestreamer specific, is that only the arguments necessay for the pulsestreamer_v2 is used in the write/load_waveform/sequence... this needs to be generalized.
        parameter_list = []
        timelist = self.make_time_list(parameters[2], time_sweep[0], time_sweep[1])
        for interleave in timelist:
            seq_after_init = [(interleave, 0), (parameters[3], 1), (parameters[4], 0)]
            _sequence_laser = [(parameters[0], 1)] + seq_after_init #)*repeat
            _sequence_counter = [(parameters[1], 0), (parameters[0] - parameters[1], 1)] + seq_after_init #)*repeat
            waveform_name = 'T1_measuerment_sequence_interleave_{0}ns'.format(interleave)

            # tab this section out if this should be done in waveform format "<= waveform" and "=> sequence"
            self.measurement_sequence_laser = self._seq_to_array(_sequence_laser)
            self.measurement_sequence_counter = self._seq_to_array(_sequence_counter)
            if len(self.measurement_sequence_laser) == len(self.measurement_sequence_laser):
                measurement_duration = len(self.measurement_sequence_laser)
            else:
                self.log.error('Error occurred in sequence writing process. Process terminated.')
            empty_array = np.asarray([0]*measurement_duration)
            # create waveform and write it.
            digi_samples = {}
            channel_seqs = []
            for chnl, val in self.pulser_hw.get_active_channels().items():
                channel_seqs.append(waveform_name + chnl[-4:])
                ch_index = self.pulser_hw._channel_to_index(chnl) - 1
                if ch_index == self._laser_channel:
                    digi_samples[chnl] = self.measurement_sequence_laser
                elif ch_index == self._counter_channel:
                    digi_samples[chnl] = self.measurement_sequence_counter
                else:
                    digi_samples[chnl] = empty_array
            self.pulser_hw.write_waveform(waveform_name, {}, digi_samples, True, True, measurement_duration)
            # make tuple format for loading into sequence.
            channel_seqs = tuple(channel_seqs)
            seq_parameters_dict = {}
            seq_parameters_dict['ensemble'] = waveform_name
            seq_parameters_dict['repetitions'] = repeat

            parameter_list.append((channel_seqs, seq_parameters_dict))
        # make sequence consisting of all waveforms and write it.
        self.pulser_hw.write_sequence(sequence_name, parameter_list)

        return sequence_name, parameter_list

    def load_sequence_to_pulser(self, sequence_name, parameter_list):
        """
        Load sequence into pulser.
        :param str sequence_name: name of already written sequence to load into pulsestreamer
        :return int err_code: error code; 0:OK, -1:err.
        """
        self.pulser_hw.load_sequence(sequence_name)
        return 0
        # for chnl, val in self.get_active_channels():
        #     ch_index = self.pulser_hw._channel_to_index(chnl)-1
        #     if ch_index == self._laser_channel:
        #         self.pulser_hw.write_waveform()
        #     elif ch_index == self._laser_channel:
        #         self.pulser_hw.write_waveform()
        #     else:
        #
        # pass

    def perform_measurement_routine(self, parameters, time_sweep, repeat=0, sequence_name='T1_measuerment_sequence'):

        if repeat == None:
            repeat = self.default_repeat

        seq_name, parameter_list = self.generate_measurement_sequences(self, parameters, time_sweep, repeat=repeat, sequence_name='T1_measurement_sequence')
        self.load_sequence_to_pulser(self, seq_name, parameter_list)




    # def convert_to_sequence_generator_logic_format(self, measurement_duration):
    #     digi_samples = {}
    #     channel_seqs = []
    #     for chnl, val in self.pulser_hw.get_active_channels():
    #         channel_seqs.append(waveform_name + '_' + chnl)
    #         ch_index = self.pulser_hw._channel_to_index(chnl) - 1
    #         if ch_index == self._laser_channel:
    #             digi_samples[chnl] = self.measurement_sequence_laser
    #         elif ch_index == self._laser_channel:
    #             digi_samples[chnl] = self.measurement_sequence_counter
    #         else:
    #             digi_samples[chnl] = empty_array
    #     self.pulser_hw.write_waveform(waveform_name, {}, digi_samples, True, True, measurement_duration)
    #
    #     # make tuple format for loading into sequence.
    #     channel_seqs = tuple(channel_seqs)
    #     seq_parameters_dict = {}
    #     seq_parameters_dict['ensemble'] = waveform_name
    #     seq_parameters_dict['repetitions'] = repeat
    #
    #     return channel_seqs, seq_parameters_dict

    def make_time_list(self, time_min, increment, num_incr):
        """
        Makes list consisting of timeintervals, from minimum time interval, increment size of time interval, and number of increments.
        :param int time_min:    minimum time interval
        :param int increment:   increment size of time interval
        :param int num_incr:    number of increments
        :return list:
        """
        return [time_min + i*increment for i in range(num_incr+1)]

    def _seq_to_array(self, seq):
        """
        Transform np array into pulse streamer sequence format

        @param list ps_seq: list consisting of tuples:
                            first element in tuple is duration in nsec
                            second element in tuple being the (digital) state of the channel

        @return array:  nparray containing the state of the digital channel of pullse streamer pr ns.
        """
        ret_list = []
        for duration, digi_state in seq:
            ret_list = ret_list + [digi_state] * duration

        return np.asarray(ret_list)

