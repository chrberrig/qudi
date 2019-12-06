# -*- coding: utf-8 -*-
"""
Interfuse to do confocal scans with spectrometer data rather than APD count rates.

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

import time
import numpy as np

from core.module import Base, Connector
from core.configoption import ConfigOption
from interface.confocal_scanner_interface import ConfocalScannerInterface


class ConfocalScanner_PI_Swabian_Interfuse(Base, ConfocalScannerInterface):

    """This is the Interface class to define the controls for the simple
    microwave hardware.
    """
    # _modclass = 'confocalscannerinterface'
    # _modtype = 'hardware'

    # connectors
    fitlogic = Connector(interface='FitLogic')
    PI_E727_controller = Connector(interface='ConfocalScannerInterface')
    timetagger_counter = Connector(interface='SlowCounterInterface')
#    pulsestreamer = Connector(interface='PulserInterface')

    # config options
    #clock_frequency = ConfigOption('clock_frequency', 100, missing='warn')
    clock_frequency = 300
    #pixel_exposure_time = 1/float(clock_frequency) # exposure time for photon collection pr. pixel in sec.
                                            # Python time.sleep is reliable within ~ 0.005s...
                                            # not suitable for times below 0.1 sec
                                            # Dont use Python time.sleep, but timeres of timetagger itself,
                                            # reliable to ps resolution

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        # Internal parameters
        self._line_length = None
        self._scanner_counter_daq_task = None
        self._voltage_range = [-10., 10.]

        self._position_range = [[0, 100e-6], [0, 100e-6], [0, 10e-6], [0, 1e-6]]
        #self._position_range = [[0., 100.], [0., 100.], [0., 10.], [0., 1.]]
        self._current_position = [0., 0., 0., 0.]

        self._num_points = 500

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """

        self._fit_logic = self.fitlogic()
        self._scanner_hw = self.PI_E727_controller()
        self._counter_hw = self.timetagger_counter()
        # self._pulser_hw = self.pulsestreamer()

        # self._counter_hw.set_up_clock(self.clock_frequency, None)
        self.set_up_scanner_clock(self.clock_frequency)
        self._counter_hw.set_up_counter()
        # self._counter_hw.test_signal([self._counter_hw.get_counter_channels()], True)

    def on_deactivate(self):
        self.reset_hardware()

    def reset_hardware(self):
        """ Resets the hardware, so the connection is lost and other programs can access it.

        @return int: error code (0:OK, -1:error)
        """
        self.log.warning('Scanning Device will be reset.')
        return 0

    def get_position_range(self):
        """ Returns the physical range of the scanner.
        This is a direct pass-through to the scanner HW.

        @return float [4][2]: array of 4 ranges with an array containing lower and upper limit
        """
        return self._scanner_hw.get_position_range()

    def set_position_range(self, myrange=None):
        """ Sets the physical range of the scanner.
        This is a direct pass-through to the scanner HW

        @param float [4][2] myrange: array of 4 ranges with an array containing lower and upper limit

        @return int: error code (0:OK, -1:error)
        """
        if myrange is None:
            myrange = [[0,1],[0,1],[0,1],[0,1]]

        self._scanner_hw.set_position_range(myrange=myrange)

        return 0

    def set_voltage_range(self, myrange=None):
        """ Sets the voltage range of the NI Card.
        This is a direct pass-through to the scanner HW

        @param float [2] myrange: array containing lower and upper limit

        @return int: error code (0:OK, -1:error)
        """
        if myrange is None:
            myrange = [-10., 10.]

        self._scanner_hw.set_voltage_range(myrange=myrange)
        return 0

    def set_up_scanner_clock(self, clock_frequency=None, clock_channel = None):
        """ Configures the hardware clock of the NiDAQ card to give the timing.
        This is a direct pass-through to the scanner HW

        @param float clock_frequency: if defined, this sets the frequency of the clock
        @param string clock_channel: if defined, this is the physical channel of the clock

        @return int: error code (0:OK, -1:error)
        """
        self.clock_frequency = clock_frequency
        self._counter_hw.set_up_clock(self.clock_frequency, None)
        # print('interfuse freq = ' + str(self.clock_frequency))

        # self._pulser_hw.clock_channel = # clock_channel
        # if clock_frequency =! 167e6:
        #     self._pulser_hw.clock_frequency = 167e6 # clock_frequency determined by allowed external clock frequency of timetagger.


        #return self._scanner_hw.set_up_scanner_clock(clock_frequency=clock_frequency, clock_channel=clock_channel)
        return 0

    def set_up_scanner(self, counter_channel=None, photon_source=None, clock_channel=None, scanner_ao_channels=None):
        """ Configures the actual scanner with a given clock.

        TODO this is not technically required, because the spectrometer scanner does not need clock synchronisation.

        @param string counter_channel: if defined, this is the physical channel of the counter
        @param string photon_source: if defined, this is the physical channel where the photons are to count from
        @param string clock_channel: if defined, this specifies the clock for the counter
        @param string scanner_ao_channels: if defined, this specifies the analoque output channels

        @return int: error code (0:OK, -1:error)
        """

        # self._counter_hw.set_up_clock(1000)

        self.log.warning('ConfocalScannerInterfaceDummy>set_up_scanner')
        return 0

    def get_scanner_axes(self):
        """ Pass through scanner axes. """
        return self._scanner_hw.get_scanner_axes()

    def scanner_set_position(self, x = None, y = None, z = None, a = None):
        """Move stage to x, y, z, a (where a is the fourth voltage channel).
        This is a direct pass-through to the scanner HW

        @param float x: postion in x-direction (volts)
        @param float y: postion in y-direction (volts)
        @param float z: postion in z-direction (volts)
        @param float a: postion in a-direction (volts)

        @return int: error code (0:OK, -1:error)
        """

        self._scanner_hw.scanner_set_position(x=x, y=y, z=z, a=a)
        return 0

    def get_scanner_position(self):
        """ Get the current position of the scanner hardware.

        @return float[]: current position in (x, y, z, a).
        """

        return self._scanner_hw.get_scanner_position()

    def _set_up_line(self, length=100):
        """ Set the line length
        Nothing else to do here, because the line will be scanned using multiple scanner_set_position calls.

        @param int length: length of the line in pixel

        @return int: error code (0:OK, -1:error)
        """
        self._line_length = length
        return 0

    def get_scanner_count_channels(self):
        """ Counting channels in confocal: normal, negative and a ramp."""
        return self._scanner_hw.get_scanner_count_channels()

    def scan_line(self, line_path=None, pixel_clock=False):
        """ Scans a line and returns the counts on that line.

        @param float[][4] line_path: array of 4-part tuples defining the voltage points
        @param bool pixel_clock: whether we need to output a pixel clock for this line

        @return float[]: the photon counts per second
        """

        # if self._counter_hw.module_state() == 'locked':
        #     self.log.error('A scan_line is already running, close this one first.')
        #     return -1
        # self._counter_hw.module_state.lock()


        if not isinstance(line_path, (frozenset, list, set, tuple, np.ndarray, )):
            self.log.error('Given voltage list is no array type.')
            return np.array([[-1.]])

        if np.shape(line_path)[1] != self._line_length:
            self._set_up_line(np.shape(line_path)[1])

        count_data = np.zeros(self._line_length)
        # self._counter_hw.configure(1e-6, self._line_length*self.pixel_exposure_time, self._line_length)
        # self._counter_hw.test_signal([self._counter_hw._channel_apd], True)
        # self._counter_hw.pause_measure()
        #self._counter_hw.start_measure()
        # print(self._counter_hw.module_state())

        self._counter_hw.start_measure()
        for i in range(self._line_length):
            #t0 = time.clock()
            coords = line_path[:, i]
            # self._scanner_hw.scanner_set_position(x=coords[0], y=coords[1], z=coords[2], a=coords[3])
            self.scanner_set_position(x=coords[0], y=coords[1], z=coords[2], a=coords[3])
            # print(time.clock() - t0)
            # time.sleep(4/self.clock_frequency)
            # print(time.clock() - t0)
            #using the slow counter: (requires only timetagger)
            # count_data[i] = self._counter_hw.measure_for_gate()
            count_data[i] = self._counter_hw.get_counter()[0][0]/self.clock_frequency
            # self._counter_hw.stop_measure()

            # updating_current position.
            self._current_position = list(line_path[:, -1])
            #print(time.clock() - t0)
            #print(count_data[i])
        self._counter_hw.stop_measure()

        return np.array([[count_data[i]] for i in range(self._line_length)])

    def close_scanner(self):
        """ Closes the scanner and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """

        self._scanner_hw.close_scanner()

        return 0

    def close_scanner_clock(self):
        """ Closes the clock and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """

        self._scanner_hw.close_scanner_clock()

        return 0
