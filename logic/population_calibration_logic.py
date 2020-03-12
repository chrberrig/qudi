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
    slowcounter = Connector(interface='SlowCounterInterface')
    pulser = Connector(interface='PulserInterface')
    scanner = Connector(interface='ConfocalScannerInterface')
    # mw_generator = Connector(interface='MicrowaveInterface')
    fitlogic = Connector(interface='FitLogic')
    savelogic = Connector(interface='SaveLogic')
    optimizerlogic = Connector(interface='OptimizerLogic')

# ====== Global parameters ======
    # === Signals ===
    sigRunFinished = QtCore.Signal()
    sigNewParameter = QtCore.Signal()
    sigMeasurementRoutineFinished = QtCore.Signal()
    sigOptimizasionRoundDone = QtCore.Signal()

    # === channels ===
    _laser_channel = 0
    _trigger_channel = 0
    _counter_channel = 0
    # _mw_gen_channel =

    # === Internal parameters ===
    measurement_routine = "pop_cal"
    wf_name = ""
    default_repeat = int(1e4)
    time_resolution = 10
    num_optim = 2
    num_attempts = 3

    # === Initialization constants ===
    ref_count = 1
    optim_round = 0
    optim_attempt = 0

    reps_pr_optim = 0
    sleep_time = 0

    # === Sweep parameters ===
    init_routine = True
    size_increments = 0 #1
    num_increments = 0 #5
    param_sweep_list = list()

    # === Data ===
    data_dict = {}
    data_dict_opt = {}
    position_array = np.zeros((num_optim + 1, 4))

# ====== Measurement Routine specific parameters ======

    # === Sequence parameters for pop_cal ===
    laser_delay = 4*int(1e6)
    laser_ontime = 2*int(1e3)
    laser_power = 5
    intersequence_delay = 4*int(1e6)
    # the above are all in units of ns.

    # waveforms for single channels:
    # measurement_wf_laser = # self._seq_to_array([])
    # measurement_wf_counter
    # measurement_wf_mw_gen


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
        self.slowcounter_hw = self.slowcounter()
        self.pulser_hw = self.pulser()
        self.scanner_hw = self.scanner()
        # self.mw_gen_hw = self.mw_generator()
        self._fit_logic = self.fitlogic()
        self._save_logic = self.savelogic()
        self._optimizer_logic = self.optimizerlogic()

        self._optimizer_logic.sigRefocusFinished.connect(self.optimize_pos, QtCore.Qt.QueuedConnection)
        self.sigNewParameter.connect(self.perform_measurement_routine, QtCore.Qt.QueuedConnection)
        # self.sigMeasurementRoutineFinished.connect(self.save_data, QtCore.Qt.QueuedConnection)
        self.sigRunFinished.connect(self.optimize_pos, QtCore.Qt.QueuedConnection)
        self.sigOptimizasionRoundDone.connect(self.single_run, QtCore.Qt.QueuedConnection)

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

    def generate_population_waveform(self, wf_name): #, parameters):
        """
        Genarates pulse waveforms for T1 measuremens wrt. the given parameters.
        :param wf_name:     str defining the waveform name.
        :param parameters:  tuple containing:
                            (laser_delay, laser_ontime_in_ns, laser_power, intersequence_delay_in_ns) in units of ns
        :return tuple:      tuple containing: waveform_name, tuple(channel_seqs_names), dict(digi_samples)
        """
        # print('Entering generate_population_waveform: ' + str(time.clock()))
        seq_common = [(self.laser_delay, 0), (self.laser_ontime, 1), (self.intersequence_delay, 0)]
        _wf_laser = seq_common
        _wf_counter = seq_common
        _wf_mw_gen = [(self.laser_delay + self.laser_ontime + self.intersequence_delay, 0)] # [(parameters[1], 0), (parameters[0] - parameters[1], 1)] #+ seq_after_init ?
        waveform_name = wf_name + '_lasercurrent_{0}%'.format(self.laser_power)
        self.measurement_wf_laser = self._seq_to_array(_wf_laser)
        self.measurement_wf_counter = self._seq_to_array(_wf_counter)
        self.measurement_wf_mw_gen = self._seq_to_array(_wf_counter)
        return self.write_to_pulser(waveform_name)

    def perform_measurement_routine(self): #, measurement_routine):
        """
        Genarates all pulse waveforms according to sweep of parameters in measurement routine.
        :param parameters:  tuple containing:
                            (laser_delay, laser_ontime_in_ns, laser_power, intersequence_delay_in_ns) in units of ns
        :param param_sweep: tuple containing:
                            (increment_size, number_of_increments)
        :param repeat:      Number of repeats for each choice of interleave-time
        :return str:        sequence name of created sequence.
        """
        # print('Entering perform_population_measurement_routine: ' + str(time.clock()))

        # ================== Setting routine values, depending on global parameters. =====================
        # param_sweep = (self.size_increments, self.num_increments)
        # self.measurement_routine = measurement_routine
        self.reps_pr_optim = int((self.default_repeat + 1) / (self.num_optim + 1))
        num_bins = 1 + int(self.time_resolution)
        self.pulser_hw.set_num_runs(self.reps_pr_optim)
        self.sleep_time = 1 + sum([self.laser_delay, self.laser_ontime, self.intersequence_delay]) * 1e-9 * self.reps_pr_optim
        # position_array = np.zeros((self.num_optim + 1, 4))

        # ================== Decide measurement routine =====================
        if self.measurement_routine == 'pop_cal':
            # parameters = (self.laser_delay, self.laser_ontime, self.laser_power, self.intersequence_delay)
            if self.init_routine:
                self.param_sweep_list = self.make_param_sweep_list(self.laser_power,
                                                                   self.size_increments,
                                                                   self.num_increments)
                # ================== Sanity checks =====================
                if self.num_increments > 100 or self.size_increments > 100 or self.num_increments < 0 or self.size_increments < 0:
                    self.log.error('Error occurred in sequence writing process. Process terminated. '
                                   'Laser sweep parameters (power/current) are ill defined. Must be integer between 0 and 100')
                self.init_routine = False
            if not self.param_sweep_list:  # check if laserlist is empty:
                # self.sigMeasurementRoutineFinished.emit()
                # self.init_routine = True # Comment this line if the procedure should be repeatable...
                return self.save_data()
            else:
                # === set new laser power ===
                self.laser_power = self.param_sweep_list.pop(0)
                self.laser_hw.set_current(float(self.laser_power))
        # === generate and load waveform ===
            waveform_name, channel_wf_name, digi_samples = self.generate_population_waveform(self.measurement_routine)
            self.wf_name = waveform_name
            # self.pulser_hw.load_waveform(self.pulser_hw.get_waveform_names())  # loading waveform
            # === init data dicts + ref_count ===
            print('self.optim_round: ' + str(self.optim_round))
            # if self.optim_round == 0:
            self.data_dict[waveform_name] = []  # initialize data entries
            self.data_dict_opt[waveform_name] = np.zeros((self.num_optim + 1, num_bins))
            # --------------------------- reference count: --------------------------------
            self.slowcounter_hw.set_up_counter()
            self.ref_count = self.slowcounter_hw.get_counter()[0][0] + 1
            self.log.info('ref_count: ' + str(self.ref_count))
            self.optim_round = 0
            self.single_run()


        # update global variable according to next action:


    #   ============== internal funct.s ===============

    def single_run(self): #, sleep_time, reps):
        # === core functionality ===
        self.pulser_hw.load_waveform(self.pulser_hw.get_waveform_names())  # loading waveform
        self.fastcounter_hw.configure(self.laser_ontime / self.time_resolution * 1e-9, self.laser_ontime * 1e-9, self.reps_pr_optim)
        self.fastcounter_hw.start_measure()
        self.pulser_hw.pulser_on()
        time.sleep(self.sleep_time)
        # self.pulser_hw.reset()
        self.pulser_hw.pulser_off()
        self.fastcounter_hw.stop_measure()
        self.data_dict_opt[self.wf_name][self.optim_round] = np.sum(self.fastcounter_hw.get_data_trace()[0], axis=0)
        # Save partial data here?
        # === control flow: ===
        if self.optim_round < self.num_optim:
            self.log.info('optimizasion round: ' + str(self.optim_round))
            self.optim_attempt = 0
            self.sigRunFinished.emit()
            # self.optimize_pos()
        elif self.optim_round == self.num_optim:
            self.sigNewParameter.emit()
        return 1

    def optimize_pos(self):
        self._laser_on()
        self.slowcounter_hw.set_up_counter()
        temp_count = self.slowcounter_hw.get_counter()[0][0]
        # print('temp_count: ' + str(temp_count))
        self.log.info('optimizasion index: ' + str(temp_count / self.ref_count))
        # while temp_count <= 0.85*ref_count:
        if temp_count <= 0.90*self.ref_count and self.optim_attempt < self.num_attempts:
            self.log.info("Begin Optimization in optimizasion round {0}, attempt {1}:".format(self.optim_round, self.optim_attempt))
            self._optimizer_logic.start_refocus()
            self.log.info('optimizasion index: ' + str(temp_count/self.ref_count))
            self.optim_attempt += 1
            self.log.info('optimizasion attempt:' + str(self.optim_attempt))
        else:
            self.log.info("Optimization round {0} done.".format(self.optim_round))
            self.position_array[self.optim_round] = self.scanner_hw.get_scanner_position()
            self.log.info('scanner_pos : ' + str(self.scanner_hw.get_scanner_position()))
            self.optim_round += 1
            self.sigOptimizasionRoundDone.emit()
        return 1

    def save_data(self): # this is experiment dependent...
        if self.measurement_routine == 'pop_cal':
            waveform_name = self.wf_name
            self.data_dict[waveform_name] = np.sum(self.data_dict_opt[waveform_name], axis=0)
            self._save_logic.save_data({waveform_name: self.data_dict[waveform_name]},
                                       parameters={'laser_ontime': self.laser_ontime,
                                                   'laser_power': '0.1mW',  # self.laser_power,
                                                   'intersequence_delay': (self.intersequence_delay + self.laser_delay)},
                                       filetype="npz"
                                       )
            time.sleep(3)
            self._save_logic.save_data({waveform_name + '_pos': self.position_array},
                                       filetype="npz"
                                       )
        return self.data_dict

    def write_to_pulser(self, waveform_name):
        if not len(self.measurement_wf_laser) == len(self.measurement_wf_counter) and not len(self.measurement_wf_laser) == len(self.measurement_wf_mw_gen):
            self.log.error('Error occurred in waveform writing process. Uneven waveform lengths. Process terminated.')
        # create waveform and write it.
        measurement_duration = len(self.measurement_wf_laser)
        empty_array = np.zeros(measurement_duration)
        digi_samples = {}
        channel_wf_name = []
        for chnl, val in self.pulser_hw.get_active_channels().items():
            channel_wf_name.append(waveform_name + chnl[-4:])
            ch_index = self._channel_to_index(chnl)
            if ch_index == self._laser_channel:
                digi_samples[chnl] = self.measurement_wf_laser
            elif ch_index == self._trigger_channel:
                digi_samples[chnl] = self.measurement_wf_counter
            # elif ch_index == self._mw_gen_channel:
            #     digi_samples[chnl] = self.measurement_wf_mw_gen
            else:
                digi_samples[chnl] = empty_array
        self.pulser_hw.write_waveform(waveform_name, {}, digi_samples, True, True, measurement_duration)
        # print('end generate_population_waveform: ' + str(time.clock()))
        return waveform_name, tuple(channel_wf_name), digi_samples

    def _laser_on(self):
        """
        Sets pulser laser channel to constant on (pulse streamer specific.)
        :return:
        """
        self.pulser_hw.set_constant([self._laser_channel])

    def make_param_sweep_list(self, param_min, increment, num_incr):
        """
        Makes list consisting of timeintervals, from minimum time interval, increment size of time interval, and number of increments.
        :param int time_min:    minimum time interval
        :param int increment:   increment size of time interval
        :param int num_incr:    number of increments
        :return list:
        """
        return [param_min + i*increment for i in range(num_incr+1)]

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
                        where X corresponds to generic channel number
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

    # def _load_data(self, npz_file, inspect=False, file_entry=None):
    #     """
    #     Loads data from npz file (campatible w. savelogic)
    #     :param str npz_file: path to .npz file containing data to be loaded.
    #     # :param str file_entry: file-entry of data to be loaded
    #     :return: returns actual np.ndarrray containing data.
    #     """
    #
    #     if inspect == True:
    #         return np.load(npz_file).files
    #     else:
    #         if file_entry == None:
    #             data_load = np.load(npz_file)[np.load(npz_file).files[0]]
    #         elif file_entry in np.load(npz_file).files:
    #             data_load = np.load(npz_file)[file_entry]
    #         return data_load