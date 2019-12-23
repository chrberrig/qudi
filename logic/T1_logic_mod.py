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
    """This is the Logic class for T1 measurements."""

    # declare connectors
    laser = Connector(interface='SimpleLaserInterface')
    fastcounter = Connector(interface='FastCounterInterface')
    pulser = Connector(interface='PulserInterface')
    # mw_generator = Connector(interface='MicrowaveInterface')
    fitlogic = Connector(interface='FitLogic')
    savelogic = Connector(interface='SaveLogic')

    # config option
    # clock_frequency = StatusVar('clock_frequency', 200)
    default_repeat = 9

    size_increments = 0
    num_increments = 0

    init_duration = 2#1000
    counter_delay = init_duration
    interleave_duration = 1000#500 # tau
    collection_duration = 200
    intersequence_delay = 2#1000
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

        # collects relevant channels from hardware.
        #TODO: this is at the moment rather hardware specific. Maybe implement this as a configOption...
        self._laser_channel = self.pulser_hw._laser_channel
        self._trigger_channel = self.fastcounter_hw._channel_detect
        self._counter_channel = self.fastcounter_hw._channel_apd
        # self._mw_gen_channel =

        self.log.info('TimeTagger (fast counter) configured to use  channel {0} as counter input channel'
                      .format(self.fastcounter_hw._channel_apd))
        self.log.info('TimeTagger (fast counter) configured to use  channel {0} as gating channel'
                      .format(self.fastcounter_hw._channel_detect))
        self.log.info('PulseStreamer configured to use  channel {0} as laser channel'
                      .format(self._laser_channel))
        self.log.info('PulseStreamer configured to use  channel {0} as trigger channel fr TT gating'
                      .format(self._trigger_channel))


        print('self._laser_channel: ' + str(self._laser_channel))
        print('self._trigger_channel: ' + str(self._trigger_channel))
        print('self._counter_channel: ' + str(self._counter_channel))
        print("waiting_time: " + str(sum((self.init_duration, self.interleave_duration, self.collection_duration, self.intersequence_delay)) * self.default_repeat) + ' ns')


        self.set_up_laser()

        self.set_up_counter()
        self.fastcounter_hw.test_signal(self._counter_channel, True) #this line is to be uncommented only during testing !!!
        self.fastcounter_hw.start_measure()
        time.sleep(1)
        self.pulser_hw.set_constatnt_state([self._trigger_channel])
        time.sleep(sum((self.init_duration, self.interleave_duration, self.collection_duration, self.intersequence_delay)) * self.default_repeat * 1e-9)
        print(self.fastcounter_hw.get_data_trace())
        print(sum(self.fastcounter_hw.get_data_trace()[0]))
        self.fastcounter_hw.stop_measure()


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


    def set_up_counter(self, bin_width_s=0, record_length_s=0, number_of_gates=0):
        """
        setting up fastcounter hardware as counter for photons.
        :return: bin_width_s, record_length_s, number_of_gates
        """
        #return self.fastcounter_hw.configure(self, bin_width_s=self.collection_duration, record_length_s=self.collection_duration, number_of_gates=0)
        if bin_width_s == 0:
            bin_width_s = self.collection_duration*1e-10
        if record_length_s == 0:
            record_length_s = self.collection_duration*1e-9
        self.fastcounter_hw.configure(bin_width_s, record_length_s, number_of_gates)
        return 0

    def set_up_mw_gen(self):
        pass


    def generate_T1_waveform(self, wf_name, parameters):
        """
        Genarates pulse waveforms for T1 measuremens wrt. the given parameters.
        :param wf_name:  str defining the waveform name.
        :param parameters:  tuple containing:
                            (initializasion_duration_in_ns, counter_delay_in_ns, interleave_duration or tau_in_ns, collection_duration_in_ns, intersequence_delay_in_ns) in units of ns
        :return tuple:      tuple containing: waveform_name, tuple(channel_seqs_names), dict(digi_samples)
        """
        time.clock()
        seq_after_init = [(parameters[2], 0), (parameters[3], 1)] #, (parameters[4], 0)]
        _sequence_laser = [(parameters[4], 0), (parameters[0], 1)] + seq_after_init
        _sequence_counter = [((parameters[4] + parameters[1]), 0), (parameters[0] - parameters[1], 1)] + seq_after_init
        # _sequence_mw_gen = [(parameters[1], 0), (parameters[0] - parameters[1], 1)] + seq_after_init
        waveform_name = wf_name + '_interleave_{0}ns'.format(parameters[2])

        # tab this section out if this should be done in waveform format "<= waveform" and "=> sequence"
        self.measurement_sequence_laser = self._seq_to_array(_sequence_laser)
        self.measurement_sequence_counter = self._seq_to_array(_sequence_counter)
        # self.measurement_sequence_mw_gen = self._seq_to_array(_sequence_counter)
        if len(self.measurement_sequence_laser) == len(self.measurement_sequence_laser):
            measurement_duration = len(self.measurement_sequence_laser)
        else:
            self.log.error('Error occurred in sequence writing process. Process terminated.')
        # create waveform and write it.
        empty_array = np.array([0] * measurement_duration)
        digi_samples = {}
        channel_seqs = []
        for chnl, val in self.pulser_hw.get_active_channels().items():
            channel_seqs.append(waveform_name + chnl[-4:])
            ch_index = self._channel_to_index(chnl)
            if ch_index == self._laser_channel:
                digi_samples[chnl] = self.measurement_sequence_laser
            elif ch_index == self._trigger_channel:
                digi_samples[chnl] = self.measurement_sequence_counter
            # elif ch_index == self._mw_gen_channel:
            #     digi_samples[chnl] = self.measurement_sequence_mw_gen
            else:
                digi_samples[chnl] = empty_array
        self.pulser_hw.write_waveform(waveform_name, {}, digi_samples, True, True, measurement_duration)
        return waveform_name, tuple(channel_seqs), digi_samples


    def generate_T1_measurement_sequences(self, parameters, param_sweep, repeat=0, sequence_name='T1_measurement_sequence'):
        """
        Genarates pulse blocks and combine them to pulse sequence
        :param parameters:  tuple containing:
                            (initializasion_duration_in_ns, counter_delay_in_ns, interleave_duration or tau_in_ns, collection_duration_in_ns, intersequence_delay_in_ns) in units of ns
        :param param_sweep:  tuple containing:
                            (increment_size, number_of_increments)
        :param repeat:      Number of repeats for each choice of interleave-time
        :return str:        sequence name of created sequence.
        """
        # if repeat == 0:
        #     repeat = self.default_repeat
        if self.counter_delay > self.init_duration:
            self.log.warn('APD is not activated during initializasion')
            self.counter_delay = 0
        # crate sequence
        parameter_list = []
        timelist = self.make_param_sweep_list(parameters[2], param_sweep[0], param_sweep[1])
        for interleave in timelist:
            loopparameters = (parameters[0], parameters[1], interleave, parameters[3], parameters[4])
            waveform_name, channel_seqs, digi_samples = self.generate_T1_waveform(sequence_name, loopparameters)
            # print(digi_samples)
            seq_parameters_dict = {}
            seq_parameters_dict['ensemble'] = waveform_name
            seq_parameters_dict['repetitions'] = repeat
            parameter_list.append((channel_seqs, seq_parameters_dict))
        # make sequence consisting of all waveforms and write it.
        self.pulser_hw.write_sequence(sequence_name, parameter_list)
        return sequence_name, parameter_list


    def perform_T1_measurement_routine(self, parameters=None, param_sweep=None, repeat=0, sequence_name='T1_measuerment_sequence'):
        """
        Genarates pulse blocks and combine them to pulse sequence
        :param parameters:  tuple containing:
                            (initializasion_duration_in_ns, counter_delay_in_ns, interleave_duration or tau_in_ns, collection_duration_in_ns, intersequence_delay_in_ns) in units of ns
        :param param_sweep:  tuple containing:
                            (increment_size, number_of_increments)
        :param repeat:      Number of repeats for each choice of interleave-time
        :return str:        sequence name of created sequence.
        """
        if repeat == 0:
            repeat = self.default_repeat

        if param_sweep == None:
            param_sweep = (self.size_increments, self.num_increments)

        if parameters == None:
            parameters = (self.init_duration, self.counter_delay, self.interleave_duration, self.collection_duration,
                          self.intersequence_delay)
        t0 = time.time()
        self.fastcounter_hw.configure(parameters[3] * 1e-10, parameters[3] * 1e-9, repeat)
        data_dict = {}
        seq_name, parameter_list = self.generate_T1_measurement_sequences(parameters, param_sweep, repeat, sequence_name)
        # print(str(time.time() - t0))
        for channel_seqs, seq_parameters_dict in parameter_list:
            # print(seq_parameters_dict['ensemble'])
            # print(parameter_list)
            # self.load_waveform_to_pulser(seq_parameters_dict['ensemble'], parameter_list)
            self.load_sequence_to_pulser(seq_name)
            # print(str(time.time() - t0))
            # self.pulser_hw.to_be_streamed.plot()
            data_dict[seq_parameters_dict['ensemble']] = []
            print(str(time.time() - t0))
            self.fastcounter_hw.start_measure()
            self.pulser_hw.pulser_on()
            print(str(time.time() - t0))
            time.sleep(sum((parameters[0], parameters[2], parameters[3], parameters[4]))*repeat*1e-9)
            print(str(time.time() - t0))
            data_dict[seq_parameters_dict['ensemble']].append(self.fastcounter_hw.get_data_trace()[0])
            # print(str(time.time() - t0))
            self.pulser_hw.pulser_off()
            self.fastcounter_hw.stop_measure()
            # data_dict[seq_parameters_dict['ensemble']] = sum(data_dict[seq_parameters_dict['ensemble']][0][0]) #/(repeat+1)
            # print(str(time.time() - t0))
        return data_dict
        # why do we only get 0 counts?


# ============== internal funct.s ===============

    def make_param_sweep_list(self, time_min, increment, num_incr):
        """
        Makes list consisting of timeintervals, from minimum time interval, increment size of time interval, and number of increments.
        :param int time_min:    minimum time interval
        :param int increment:   increment size of time interval
        :param int num_incr:    number of increments
        :return list:
        """
        return [time_min + i*increment for i in range(num_incr+1)]


    def generate_waveform_load_dict(self, waveform_name, parameters):
        """
        generates load_dict for wavefom in accordence w. the parameters
        :param parameters:  str defining the waveform name.
        :param parameters:  tuple containing:
                            (initializasion_duration_in_ns, counter_delay_in_ns, interleave_duration or tau_in_ns, collection_duration_in_ns, intersequence_delay_in_ns) in units of ns
        :return:            Dictionaries containing as keys the generic channel indices and as values the corresponding waveform channel-names.
        """
        load_dict = {}
        if type(parameters) == tuple:
            waveform_name, channel_seqs, digi_samples = self.generate_waveform(waveform_name, parameters)
            for wf_chnl_name in channel_seqs:
                load_dict[int(wf_chnl_name[-1])] = wf_chnl_name
        elif type(parameters) == list:
            for channel_seqs, seq_parameters_dict in parameters:
                if waveform_name == seq_parameters_dict['ensemble']:
                    for ch in channel_seqs:
                        load_dict[self._channel_to_index(ch) + 1] = ch
        return load_dict


    def load_waveform_to_pulser(self, waveform_name, parameters):
        """
        Loads waveform into pulser.
        :param waveform_name str: name of (channel-collection of) waveform(s) to be loaded into pulser e.g. 'dummy_ens', not 'dummy_ens_ch1' etc.
        :param parameter_list: parameter list of the form generated in generate_measurement_sequences-function.
        :return: error code; 0: OK, -1: error
        """
        load_dict = self.generate_waveform_load_dict(waveform_name, parameters)
        return self.pulser_hw.load_waveform(load_dict)
        # return 0


    def load_sequence_to_pulser(self, seq_name):
        """
        Loads waveform into pulser.
        :param waveform_name str: name of (channel-collection of) waveform(s) to be loaded into pulser e.g. 'dummy_ens', not 'dummy_ens_ch1' etc.
        :param parameter_list: parameter list of the form generated in generate_measurement_sequences-function.
        :return: error code; 0: OK, -1: error
        """
        return self.pulser_hw.load_sequence(seq_name)
        # return 0


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
        else:
            self.log.error('Channel not given or ill-defined')

            # defaultparameters = (self.init_duration, self.counter_delay, self.interleave_duration, self.collection_duration, self.intersequence_delay)
            # parameters = tuple(map(lambda x, y: y if y is not None else x, defaultparameters, parameters))

            # def load_sequence_to_pulser(self, sequence_name):
            #     """
            #     Load sequence into pulser.
            #     :param str sequence_name: name of already written sequence to load into pulsestreamer
            #     :return int err_code: error code; 0:OK, -1:err.
            #     """
            #     self.pulser_hw.load_sequence(sequence_name)
            #     return 0

            # def load_waveform_from_seq_to_pulser(self, waveform_name, parameter_list):
            #     """
            #     Loads waveform into pulser.
            #     :param waveform_name str: name of (channel-collection of) waveform(s) to be loaded into pulser e.g. 'dummy_ens', not 'dummy_ens_ch1' etc.
            #     :param parameter_list: parameter list of the form generated in generate_measurement_sequences-function.
            #     :return: error code; 0: OK, -1: error
            #     """
            #     load_dict = {}
            #     for channel_seqs, seq_parameters_dict in parameter_list:
            #         if waveform_name == seq_parameters_dict['ensemble']:
            #             for ch in channel_seqs:
            #                 # print(self._channel_to_index(ch))
            #                 # print(ch)
            #                 load_dict[self._channel_to_index(ch)+1] = ch
            #
            #     self.pulser_hw.load_waveform(load_dict)
            #     return 0