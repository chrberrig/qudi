# -*- coding: utf-8 -*-
"""
This file contains the Qudi Dummy file for ODMRCounter.

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
import time

from core.module import Base
from core.connector import Connector
from core.configoption import ConfigOption
from interface.odmr_counter_interface import ODMRCounterInterface


class ODMRCounterDummy(Base, ODMRCounterInterface):
    """ Dummy hardware class to simulate the controls for a simple ODMR.

    Example config for copy-paste:

    odmr_counter_dummy:
        module.Class: 'odmr_counter_dummy.ODMRCounterDummy'
        clock_frequency: 100 # in Hz
        number_of_channels: 2
        fitlogic: 'fitlogic' # name of the fitlogic module, see default config

    """

    # connectors
    fitlogic = Connector(interface='FitLogic')

    # config options
    _clock_frequency = ConfigOption('clock_frequency', 100, missing='warn')
    _number_of_channels = ConfigOption('number_of_channels', 2, missing='warn')
    _gate_detection_channel = ConfigOption('gate_detection_channel', 1, missing='warn')

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self._scanner_counter_daq_task = None
        self._odmr_length = None
        self._pulse_out_channel = 'dummy'
        self._lock_in_active = False
        self._oversampling = 10

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """
        self._fit_logic = self.fitlogic()

        self._counter_channel = self._channel_apd
        self._photon_source = None
        self._clock_channel = None  # self._channel_detect ?
        self._odmr_trigger_channel = None
        self._channel_detect = None

        self._odmr_length = None

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        self.log.debug('ODMR counter is shutting down.')
        if self.module_state() == 'locked':
            self._tagger.stop()
        self._tagger.clear()
        self._tagger = None

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

        self.log.info('ODMRCounterDummy>set_up_odmr')

        if self.module_state() == 'locked' or self._scanner_counter_daq_task is not None:
            self.log.error('Another odmr is already running, close this one '
                    'first.')
            return -1

        self._number_of_gates = number_of_gates
        self._bin_width = bin_width_s * 1e9
        self._record_length = 1 + int(record_length_s / bin_width_s)
        self.statusvar = 1

        self.pulsed = tt.TimeDifferences(
            tagger=self._tagger,
            click_channel=self._channel_apd,
            start_channel=self._channel_detect,
            next_channel=self._channel_detect,
            sync_channel=tt.CHANNEL_UNUSED,
            binwidth=int(np.round(self._bin_width * 1000)),
            n_bins=int(self._record_length),
            n_histograms=number_of_gates)

        self.pulsed.stop()

        return 0

    def set_odmr_length(self, length=100):
        """ Sets up the trigger sequence for the ODMR and the triggered microwave.

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



        self.module_state.unlock()
        return False, ret


    def close_odmr(self):
        """ Closes the odmr and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """

        self.log.info('ODMRCounterDummy>close_odmr')



        return 0

    def close_odmr_clock(self):
        """ Closes the odmr and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """

        self.log.info('ODMRCounterDummy>close_odmr_clock')

        return 0

    def get_odmr_channels(self):
        """ Return a list of channel names.

        @return list(str): channels recorded during ODMR measurement
        """
        return ['ch{0:d}'.format(i) for i in range(1, self._number_of_channels + 1)]

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
