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


class PopulationCalibrationLogic(GenericLogic):
    """This is the Logic class for T1 measurements."""

    # declare connectors
    laser = Connector(interface='SimpleLaserInterface')
    fastcounter = Connector(interface='FastCounterInterface')
    pulser = Connector(interface='PulserInterface')
    # mw_generator = Connector(interface='MicrowaveInterface')
    fitlogic = Connector(interface='FitLogic')
    savelogic = Connector(interface='SaveLogic')

    default_repeat = 999
    time_resolution = 100

    # Sweep parameters:
    size_increments = 1
    num_increments = 10

    # Sequence parameters:
    laser_delay = 5
    laser_ontime = 10
    laser_power = 5
    intersequence_delay = 1000
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

        self.set_up_laser()
        self.set_up_counter()

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
        return self.laser_hw.on()

    def set_up_counter(self, bin_width_s=None, record_length_s=None, number_of_gates=0):
        """
        setting up fastcounter hardware as counter for photons.
        :return: bin_width_s, record_length_s, number_of_gates
        """
        if bin_width_s == None:
            bin_width_s = self.laser_ontime
        if record_length_s == None:
            record_length_s = self.laser_ontime
        return self.fastcounter_hw.configure(bin_width_s, record_length_s, number_of_gates)

    def set_up_mw_gen(self):
        pass

    def generate_population_waveform(self, wf_name, parameters):
        """
        Genarates pulse waveforms for T1 measuremens wrt. the given parameters.
        :param wf_name:     str defining the waveform name.
        :param parameters:  tuple containing:
                            (laser_delay, laser_ontime_in_ns, laser_power, intersequence_delay_in_ns) in units of ns
        :return tuple:      tuple containing: waveform_name, tuple(channel_seqs_names), dict(digi_samples)
        """
        seq_common = [(parameters[0], 0), (parameters[1], 1), (parameters[3], 0)]
        _sequence_laser = seq_common
        _sequence_counter = seq_common
        # _sequence_mw_gen = [(parameters[1], 0), (parameters[0] - parameters[1], 1)] + seq_after_init
        waveform_name = wf_name + '_lasercurrent_{0}%'.format(parameters[2])

        # tab this section out if this should be done in waveform format "<= waveform" and "=> sequence"
        self.measurement_sequence_laser = self._seq_to_array(_sequence_laser)
        self.measurement_sequence_counter = self._seq_to_array(_sequence_counter)
        # self.measurement_sequence_mw_gen = self._seq_to_array(_sequence_counter)
        measurement_duration = len(self.measurement_sequence_laser)
        if not len(self.measurement_sequence_laser) == len(self.measurement_sequence_laser):
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

    def generate_population_measurement_sequences(self, parameters, param_sweep, repeat=0, sequence_name='population_measurement_sequence'):
        """
        Genarates pulse blocks and combine them to pulse sequence
        :param parameters:  tuple containing:
                            (laser_delay, laser_ontime_in_ns, laser_power, intersequence_delay_in_ns) in units of ns
        :param param_sweep:  tuple containing:
                            (increment_size, number_of_increments)
        :param repeat:      Number of repeats for each choice of interleave-time
        :return str:        sequence name of created sequence.
        """
        if param_sweep[1] > 100 or param_sweep[0] > 100:
            self.log.error('Error occurred in sequence writing process. Process terminated.'
                           'Laser sweep parameters are ill defined.')
        # crate sequence
        parameter_list = []
        laserlist = self.make_param_sweep_list(parameters[2], param_sweep[0], param_sweep[1])
        for laserpower in laserlist:
            loopparameters = (parameters[0], parameters[1], laserpower, parameters[3])
            waveform_name, channel_seqs, digi_samples = self.generate_population_waveform(sequence_name, loopparameters)
            seq_parameters_dict = {}
            seq_parameters_dict['ensemble'] = waveform_name
            seq_parameters_dict['repetitions'] = repeat
            seq_parameters_dict['laserpower'] = laserpower
            parameter_list.append((channel_seqs, seq_parameters_dict))
        # make sequence consisting of all waveforms and write it.
        self.pulser_hw.write_sequence(sequence_name, parameter_list)
        return sequence_name, parameter_list

    def perform_population_measurement_routine(self, parameters=None, param_sweep=None, repeat=0, sequence_name='T1_measuerment_sequence'):
        """
        Genarates pulse blocks and combine them to pulse sequence
        :param parameters:  tuple containing:
                            (laser_delay, laser_ontime_in_ns, laser_power, intersequence_delay_in_ns) in units of ns
        :param param_sweep: tuple containing:
                            (increment_size, number_of_increments)
        :param repeat:      Number of repeats for each choice of interleave-time
        :return str:        sequence name of created sequence.
        """
        if repeat == 0:
            repeat = self.default_repeat

        if param_sweep == None:
            param_sweep = (self.size_increments, self.num_increments)

        if parameters == None:
            parameters = (self.laser_delay, self.laser_ontime, self.laser_power, self.intersequence_delay)

        # laserlist = self.make_param_sweep_list(parameters[2], param_sweep[0], param_sweep[1])

        self.fastcounter_hw.configure(parameters[1]/self.time_resolution * 1e-9, parameters[1] * 1e-9, repeat)
        data_dict = {}
        seq_name, parameter_list = self.generate_population_measurement_sequences(parameters, param_sweep, repeat, sequence_name)
        for channel_seqs, seq_parameters_dict in parameter_list:
            self.load_sequence_to_pulser(seq_name)
            # self.pulser_hw.plot_loaded_asset()
            data_dict[seq_parameters_dict['ensemble']] = []
            # input("Press Enter to continue...")
            self.laser_hw.set_current(float(seq_parameters_dict['laserpower']))
            self.fastcounter_hw.start_measure()
            self.pulser_hw.pulser_on()
            # for rep in range(repeat+1):
            #     self.pulser_hw.pulser_on()
            #     time.sleep(parameters[3]*1e-9)
            #     data_dict[seq_parameters_dict['ensemble']].append(self.fastcounter_hw.get_data_trace()[0])
            time.sleep(1.5 + sum([parameters[0], parameters[1], parameters[3]])*1e-9*repeat)
            data_dict[seq_parameters_dict['ensemble']].append(self.fastcounter_hw.get_data_trace()[0])
            # data_dict[seq_parameters_dict['ensemble']] = sum(data_dict[seq_parameters_dict['ensemble']][0][0]) #/(repeat+1)
            self.pulser_hw.pulser_off()
            # self.pulser_hw.reset()
            self.fastcounter_hw.stop_measure()
            data_dict[seq_parameters_dict['ensemble']] = np.sum(data_dict[seq_parameters_dict['ensemble']][0], axis=0)#, keepdims=True)
        return data_dict

# ============== internal funct.s ===============

    def make_param_sweep_list(self, param_min, increment, num_incr):
        """
        Makes list consisting of timeintervals, from minimum time interval, increment size of time interval, and number of increments.
        :param int time_min:    minimum time interval
        :param int increment:   increment size of time interval
        :param int num_incr:    number of increments
        :return list:
        """
        return [param_min + i*increment for i in range(num_incr+1)]

    def generate_waveform_load_dict(self, waveform_name, parameters):
        """
        generates load_dict for wavefom in accordence w. the parameters
        :param parameters:  str defining the waveform name.
        :param parameters:  tuple containing parameters characterizing the waveform eg.:
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
