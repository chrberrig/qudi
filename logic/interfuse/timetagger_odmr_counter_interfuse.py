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
import TimeTagger as tt
from core.module import Base
from core.configoption import ConfigOption
from core.connector import Connector
from interface.odmr_counter_interface import ODMRCounterInterface
from interface.fast_counter_interface import FastCounterInterface
import os


class TimeTaggerODMRCounterInterfuse(Base, ODMRCounterInterface): # , FastCounterInterface):
    # _modclass = 'TimeTaggerFastCounter'
    # _modtype = 'hardware'

    #config options
    _channel_apd_0 = ConfigOption('timetagger_channel_apd_0', missing='error')
    _channel_apd_1 = ConfigOption('timetagger_channel_apd_1', missing='error')
    _channel_detect = ConfigOption('timetagger_channel_detect', missing='error') # what is channel detect used for exactly? External Trigger events.
#     _channel_sequence = ConfigOption('timetagger_channel_sequence', missing='error')# what is channel_sequence used for exactly? nothing according to fastcounter hardware module...
    _sum_channels = ConfigOption('timetagger_sum_channels', True, missing='warn')
    _clock_frequency = ConfigOption('clock_frequency', 100, missing='warn')

    # connectors
    fast_counter = Connector(interface='FastCounterInterface')
    # odmr_counter1 = Connector(interface='ODMRCounterInterface')

    def on_activate(self):
        """ Connect and configure the access to the TimeTagger.
        """
        # self._tagger = tt.createTimeTagger()
        # self._tagger.reset()
        #
        # # defining the APD-channel
        # if self._sum_channels == True:
        #     self._channel_combined = tt.Combiner(self._tagger, channels=[self._channel_apd_0, self._channel_apd_1])
        #     self._channel_apd = self._channel_combined.getChannel()
        # else:
        #     self._channel_apd = self._channel_apd_0
        #
        # self.log.info('TimeTagger (fast counter) configured to use  channel {0}'
        #               .format(self._channel_apd))
        #
        # self.statusvar = 0

        self.fastcount_hw = self.fast_counter()
        # self.odmrcount_hw = self.odmr_counter1()

        # ODMR variables:
        self._counter_channel = self._channel_apd
        self._photon_source = None
        self._clock_channel = None # self._channel_detect ?
        self._odmr_trigger_channel = None

        self._odmr_length = None


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
            self._tagger.stop()
        self._tagger.clear()
        self._tagger = None

    def set_up_odmr_clock(self, clock_frequency=None, clock_channel=None):
        """ Configures the hardware clock of the NiDAQ card to give the timing.
            @param float clock_frequency: if defined, this sets the frequency of the clock
            @param str clock_channel: if defined, this is the physical channel of the clock
                @return int: error code (0:OK, -1:error)
        """
        if clock_frequency is not None:
            self._clock_frequency = float(clock_frequency)

        self.log.info('ODMRCounterDummy>set_up_odmr_clock')

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

        self.log.info('TimeTaggerODMRCounter>set_up_odmr')

        if self.module_state() == 'locked':
            self.log.error('A scan_line is already running, close this one '
                           'first.')
            return -1

        if counter_channel is not None:
            self._counter_channel = counter_channel
        if photon_source is not None:
            self._photon_source = photon_source
        if clock_channel is not None:
            self._clock_channel = clock_channel
        if odmr_trigger_channel is not None:
            self._odmr_trigger_channel = odmr_trigger_channel

        self.fastcount_hw.configure(bin_width_s=0.01/self._clock_frequency,
                                    record_length_s=1/self._clock_frequency,
                                    number_of_gates=0)

        self.module_state.unlock()

        return 0

    def set_odmr_length(self, length=100):
        """ Set up the trigger sequence for the ODMR and the triggered microwave.
            @param int length: length of microwave sweep in pixel
                @return int: error code (0:OK, -1:error)
        """
        self._odmr_length = length
        return 0

    def count_odmr(self, length=100):
        """ Sweeps the microwave and returns the counts on that sweep.
            @param int length: length of microwave sweep in pixel
                @return float[]: the photon counts per second
        """
        if self.module_state() == 'locked':
            self.log.error('A scan_line is already running, close this one '
                           'first.')
            return -1



        self.module_state.lock()
        self._odmr_length = length

        ret = np.empty((1, length))


        self.module_state.unlock()

    def close_odmr(self):
        """ Close the odmr and clean up afterwards.
                @return int: error code (0:OK, -1:error)
        """
        pass

    def close_odmr_clock(self):
        """ Close the odmr and clean up afterwards.
                @return int: error code (0:OK, -1:error)
        """
        pass

    def get_odmr_channels(self):
        """ Return a list of channel names.
                @return list(str): channels recorded during ODMR measurement
        """
        ch_list = []
        for entry in [self._counter_channel, self._photon_source, self._clock_channel, self._odmr_trigger_channel]:
            if entry is not None:
                ch_list.append(entry)
        return ch_list

# =============================== Fastcounter functionality ===============================

    # def get_constraints(self):
    #     """ Retrieve the hardware constrains from the Fast counting device.
    #
    #     @return dict: dict with keys being the constraint names as string and
    #                   items are the definition for the constaints.
    #
    #      The keys of the returned dictionary are the str name for the constraints
    #     (which are set in this method).
    #
    #                 NO OTHER KEYS SHOULD BE INVENTED!
    #
    #     If you are not sure about the meaning, look in other hardware files to
    #     get an impression. If still additional constraints are needed, then they
    #     have to be added to all files containing this interface.
    #
    #     The items of the keys are again dictionaries which have the generic
    #     dictionary form:
    #         {'min': <value>,
    #          'max': <value>,
    #          'step': <value>,
    #          'unit': '<value>'}
    #
    #     Only the key 'hardware_binwidth_list' differs, since they
    #     contain the list of possible binwidths.
    #
    #     If the constraints cannot be set in the fast counting hardware then
    #     write just zero to each key of the generic dicts.
    #     Note that there is a difference between float input (0.0) and
    #     integer input (0), because some logic modules might rely on that
    #     distinction.
    #
    #     ALL THE PRESENT KEYS OF THE CONSTRAINTS DICT MUST BE ASSIGNED!
    #
    #     # Example for configuration with default values:
    #
    #     constraints = dict()
    #
    #     # the unit of those entries are seconds per bin. In order to get the
    #     # current binwidth in seonds use the get_binwidth method.
    #     constraints['hardware_binwidth_list'] = []
    #
    #     """
    #     self.fast_dev.get_constraints()
    #
    # def configure(self, bin_width_s, record_length_s, number_of_gates=0):
    #     """ Configuration of the fast counter.
    #
    #     @param float bin_width_s: Length of a single time bin in the time
    #                               trace histogram in seconds.
    #     @param float record_length_s: Total length of the timetrace/each
    #                                   single gate in seconds.
    #     @param int number_of_gates: optional, number of gates in the pulse
    #                                 sequence. Ignore for not gated counter.
    #
    #     @return tuple(binwidth_s, record_length_s, number_of_gates):
    #                 binwidth_s: float the actual set binwidth in seconds
    #                 gate_length_s: the actual record length in seconds
    #                 number_of_gates: the number of gated, which are accepted, None if not-gated
    #     """
    #     self.fast_dev.configure()
    #
    # @abstract_interface_method
    # def get_status(self):
    #     """ Receives the current status of the Fast Counter and outputs it as
    #         return value.
    #
    #     0 = unconfigured
    #     1 = idle
    #     2 = running
    #     3 = paused
    #   -1 = error state
    #     """
    #     pass
    #
    # @abstract_interface_method
    # def start_measure(self):
    #     """ Start the fast counter. """
    #     pass
    #
    # @abstract_interface_method
    # def stop_measure(self):
    #     """ Stop the fast counter. """
    #     pass
    #
    # @abstract_interface_method
    # def pause_measure(self):
    #     """ Pauses the current measurement.
    #
    #     Fast counter must be initially in the run state to make it pause.
    #     """
    #     pass
    #
    # @abstract_interface_method
    # def continue_measure(self):
    #     """ Continues the current measurement.
    #
    #     If fast counter is in pause state, then fast counter will be continued.
    #     """
    #     pass
    #
    # @abstract_interface_method
    # def is_gated(self):
    #     """ Check the gated counting possibility.
    #
    #     @return bool: Boolean value indicates if the fast counter is a gated
    #                   counter (TRUE) or not (FALSE).
    #     """
    #     pass
    #
    # @abstract_interface_method
    # def get_binwidth(self):
    #     """ Returns the width of a single timebin in the timetrace in seconds.
    #
    #     @return float: current length of a single bin in seconds (seconds/bin)
    #     """
    #     pass
    #
    # @abstract_interface_method
    # def get_data_trace(self):
    #     """ Polls the current timetrace data from the fast counter.
    #
    #     Return value is a numpy array (dtype = int64).
    #     The binning, specified by calling configure() in forehand, must be
    #     taken care of in this hardware class. A possible overflow of the
    #     histogram bins must be caught here and taken care of.
    #     If the counter is NOT GATED it will return a tuple (1D-numpy-array, info_dict) with
    #         returnarray[timebin_index]
    #     If the counter is GATED it will return a tuple (2D-numpy-array, info_dict) with
    #         returnarray[gate_index, timebin_index]
    #
    #     info_dict is a dictionary with keys :
    #         - 'elapsed_sweeps' : the elapsed number of sweeps
    #         - 'elapsed_time' : the elapsed time in seconds
    #
    #     If the hardware does not support these features, the values should be None
    #     """
    #     pass