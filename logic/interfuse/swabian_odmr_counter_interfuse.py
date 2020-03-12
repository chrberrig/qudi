# -*- coding: utf-8 -*-
"""
A hardware module for communicating with the fast counter FPGA.

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

import numpy as np
# import TimeTagger as tt
from core.module import Base
from core.configoption import ConfigOption
from core.connector import Connector
from interface.odmr_counter_interface import ODMRCounterInterface
# from interface.fast_counter_interface import FastCounterInterface
import os
import time


class SwabianODMRCounterInterfuse(Base, ODMRCounterInterface):
    """Inetrfuse of Timetagger and Pulstreamer from Swabian Instruments, to implement ODMR counter functionality """

    # from Dummy ODMR counter:
    # # connectors
    # fitlogic = Connector(interface='FitLogic')

    # config options
    _clock_frequency = ConfigOption('clock_frequency', 100, missing='warn')
    _switch_channel = ConfigOption('switch_channel', 2, missing='warn')
    #this is from dummy counter, and should not be used.
    # _number_of_channels = ConfigOption('number_of_channels', 1, missing='warn') # should this be set to 1 by default?

     # osc_period = 1


    # connectors
    fast_counter = Connector(interface='FastCounterInterface')
    pulser = Connector(interface='PulserInterface')
    laser = Connector(interface='SimpleLaserInterface')

    #from odmr_counter_dummy
    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self. osc_period = 1
        self._odmr_length = 100 # None
        # self._pulse_out_channel = ''
        self._lock_in_active = False
        self._oversampling = 10

    def on_activate(self):
        """ Connect and configure the access to the TimeTagger and pulsestreamer?.
        """
        # Instrument connectors
        self.fastcounter_hw = self.fast_counter()
        self.pulser_hw = self.pulser()
        self.laser_hw = self.laser()

        # ODMR variables:
        self._counter_channel = self.fastcounter_hw._channel_apd
        self._photon_source = self.pulser_hw._laser_channel
        self._clock_channel = self.fastcounter_hw._channel_detect
        self._odmr_trigger_channel = self.pulser_hw._uw_x_channel
        # self._switch_channel = 2

        # Instrument init.
        self.laser_hw.on()


    def get_constraints(self):
        """ Retrieve the hardware constrains from the Fast counting device.

        @return dict: dict with keys being the constraint names as string and
                      items are the definition for the constraints.

         The keys of the returned dictionary are the str name for the constraints
        (which are set in this method).

                    NO OTHER KEYS SHOULD BE INVENTED!

        If you are not sure about the meaning, look in other hardware files to
        get an impression. If still additional constraints are needed, then they
        have to be added to all files containing this interface.

        The items of the keys are again dictionaries which have the generic
        dictionary form:
            {'min': <value>,
             'max': <value>,
             'step': <value>,
             'unit': '<value>'}

        Only the key 'hardware_binwidth_list' differs, since they
        contain the list of possible binwidths.

        If the constraints cannot be set in the fast counting hardware then
        write just zero to each key of the generic dicts.
        Note that there is a difference between float input (0.0) and
        integer input (0), because some logic modules might rely on that
        distinction.

        ALL THE PRESENT KEYS OF THE CONSTRAINTS DICT MUST BE ASSIGNED!
        """

        constraints = dict()

        # the unit of those entries are seconds per bin. In order to get the
        # current binwidth in seonds use the get_binwidth method.
        constraints['hardware_binwidth_list'] = [1 / 1000e6]

        # TODO: think maybe about a software_binwidth_list, which will
        #      postprocess the obtained counts. These bins must be integer
        #      multiples of the current hardware_binwidth

        return constraints

    def on_deactivate(self):
        """ Deactivate the FPGA.
        """
        if self.module_state() == 'locked':
            pass

        self.log.debug('ODMR counter is shutting down.')

    def set_up_odmr_clock(self, clock_frequency=None, clock_channel=None):
        """ Configures the hardware clock of the NiDAQ card to give the timing.
            @param float clock_frequency: if defined, this sets the frequency of the clock
            @param str clock_channel: if defined, this is the physical channel of the clock
                @return int: error code (0:OK, -1:error)
        """
        if clock_frequency is not None:
            self._clock_frequency = float(clock_frequency)

        if clock_channel is not None:
            self._clock_channel = int(clock_channel)



        waveform_name, seq, measurement_duration = self._def_osc(self._clock_frequency)
        waveform_name, channel_wf_name, digi_samples = self.write_to_pulser(waveform_name, seq, measurement_duration)
        # self.pulser_hw.write_waveform(waveform_name, {}, digi_samples, True, True, measurement_duration) # TODO Writ helper funct. that makes waveform and print it to the desired chnl...
        self.pulser_hw.load_waveform(self.pulser_hw.get_waveform_names())
        self.log.info('SwabianODMRCounterInterfuse>set_up_odmr_clock')

        return 0

    def set_up_odmr(self, counter_channel=None, photon_source=None,
                    clock_channel=None, odmr_trigger_channel=None):
        """ Configures the actual counter with a given clock.
            @param str counter_channel: if defined, this is the physical channel of the counter
            @param str photon_source: if defined, this is the physical channel where the photons are to count from
            @param str clock_channel: if defined, this specifies the clock for the counter
            @param str odmr_trigger_channel: if defined, this specifies the trigger output for the microwave
                @return int: error code (0:OK, -1:error)
        """

        self.log.info('SwabianODMRCounterInterfuse>set_up_odmr')

        if self.module_state() == 'locked':
            self.log.error('A scan_line is already running, close this one first.')
            return -1

        if counter_channel is not None:
            self._counter_channel = counter_channel
        if photon_source is not None:
            self._photon_source = photon_source
        if clock_channel is not None:
            self._clock_channel = clock_channel
        if odmr_trigger_channel is not None:
            self._odmr_trigger_channel = odmr_trigger_channel

        # self.fastcounter_hw.configure(bin_width_s=1. / self._clock_frequency,
        #                             record_length_s=1. / self._clock_frequency,
        #                             number_of_gates=self._odmr_length)# 1)

        return 0

    def set_odmr_length(self, length=100):
        """ Set up the trigger sequence for the ODMR and the triggered microwave.
            @param int length: length of microwave sweep in pixel
                @return int: error code (0:OK, -1:error)
        """
        self._odmr_length = length

        self.fastcounter_hw.configure(bin_width_s=1. / self._clock_frequency,
                                      record_length_s=1. / self._clock_frequency,
                                      number_of_gates=self._odmr_length)
        self.pulser_hw.set_num_runs(self._odmr_length)

        # self.fastcounter_hw.configure_gated(number_of_gates=self._odmr_length)

        return 0

    def count_odmr(self, length=100):
        """ Sweeps the microwave and returns the counts on that sweep.
            @param int length: length of microwave sweep in pixel
                @return float[]: the photon counts per second
        """
        if self.module_state() == 'locked':
            self.log.error('A scan_line is already running, close this one first.')
            return -1

        self.module_state.lock()
        self._odmr_length = length

        self.set_odmr_length(self._odmr_length)

        ret = np.empty((len(self.get_odmr_channels()), self._odmr_length))

        self.pulser_hw.pulser_on()
        self.fastcounter_hw.start_measure()
        time.sleep(self._odmr_length * 1. / self._clock_frequency)
        self.fastcounter_hw.stop_measure()
        self.pulser_hw.pulser_off()

        for chnl_index in range(len(self.get_odmr_channels())):
            # print(self.fastcounter_hw.get_data_trace()[0][:, 0])# * self._clock_frequency )
            ret[chnl_index] = self.fastcounter_hw.get_data_trace()[0][:, 0] * self._clock_frequency
            # ret[chnl_index] = self.fastcounter_hw.get_data_trace()[0]

        self.module_state.unlock()

        return False, ret

    # for chnl_index in range(len(self.get_odmr_channels())):
    #     # line_data = np.empty((1, self._odmr_length))
    #     t0 = time.time()
    #     for i in range(self._odmr_length):
    #         t1 = time.time()
    #         time.sleep(1. / self._clock_frequency)
    #         # line_data[i] = self.fastcounter_hw.get_data_trace()[0]
    #         print(1. / self._clock_frequency)
    #         print(self.fastcounter_hw.get_data_trace()[0][0])
    #         ret[chnl_index][i] = self.fastcounter_hw.get_data_trace()[0][0]
    #         print(time.time() - t1)
    #     print(time.time() - t0)
    #     print(ret[chnl_index])


    def close_odmr(self):
        """ Close the odmr and clean up afterwards.
                @return int: error code (0:OK, -1:error)
        """

        self.log.info('SwabianODMRCounterInterfuse>close_odmr')
        if self.module_state() == 'locked':
            self.log.error('A scan_line is already running, close this one first.')
            return -1

        return 0

    def close_odmr_clock(self):
        """ Close the odmr and clean up afterwards.
                @return int: error code (0:OK, -1:error)
        """

        self.log.info('SwabianODMRCounterInterfuse>close_odmr_clock')
        if self.module_state() == 'locked':
            self.log.error('A scan_line is already running, close this one first.')
            return -1

        return 0

    def get_odmr_channels(self):
        """ Return a list of channel names.
                @return list(str): channels recorded during ODMR measurement
        """
        ch_list = []
        for entry in [self._counter_channel]: # , self._photon_source, self._clock_channel, self._odmr_trigger_channel]:
            if entry is not None:
                ch_list.append(entry)
        return ch_list

    @property
    def oversampling(self):
        return self._oversampling

    @oversampling.setter
    def oversampling(self, val):
        if not isinstance(val, (int, float)):
            self.log.error('oversampling has to be int of float.')
        else:
            self._oversampling = int(val)

    @property
    def lock_in_active(self):
        return self._lock_in_active

    @lock_in_active.setter
    def lock_in_active(self, val):
        if not isinstance(val, bool):
            self.log.error('lock_in_active has to be boolean.')
        else:
            self._lock_in_active = val
            if self._lock_in_active:
                self.log.warn('Lock-In is not implemented')

    # === internal_functions ===

    def write_to_pulser(self, waveform_name, single_osc_seq, measurement_duration):
        # create waveform and write it.
        # measurement_duration = len(self.measurement_wf_laser)
        empty_array = np.zeros(measurement_duration)
        digi_samples = {}
        channel_wf_name = []
        for chnl, val in self.pulser_hw.get_active_channels().items():
            channel_wf_name.append(waveform_name + chnl[-4:])
            ch_index = self._channel_to_index(chnl)
            if ch_index == self._clock_channel:
                digi_samples[chnl] = single_osc_seq
            elif ch_index == self._odmr_trigger_channel:
                digi_samples[chnl] = single_osc_seq
            # elif ch_index == self._switch_channel:
            #     digi_samples[chnl] = self._seq_to_array([(measurement_duration, 1)])
            # elif ch_index == self._photon_source:
            #     digi_samples[chnl] = self._seq_to_array([(measurement_duration, 1)])
            else:
                digi_samples[chnl] = empty_array
        self.pulser_hw.write_waveform(waveform_name, {}, digi_samples, True, True, measurement_duration)
        # print('end generate_population_waveform: ' + str(time.clock()))
        return waveform_name, tuple(channel_wf_name), digi_samples

    def _seq_to_array(self, seq):
        """
        Transform np array into pulse streamer sequence format

        @param list seq: list consisting of tuples:
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

    def _def_osc(self, freq=None):
        """ Define a single oscillation of a TTL waveform defined from a frequency

        @:param float freq: frequency defining the oscillation
        """
        if freq is not None:
            self._clock_frequency = float(freq)

        period = int(1e9 * 1./self._clock_frequency) # in ns
        ontime = int(period / 2)
        single_osc_seq = [(ontime, 1), (period - ontime, 0)]
        self.osc_period = period

        waveform_name = 'osc_period_{0}_ns'.format(period)

        return waveform_name, self._seq_to_array(single_osc_seq), period




