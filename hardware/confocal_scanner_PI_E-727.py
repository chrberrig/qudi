# -*- coding: utf-8 -*-
"""
This file contains the Qudi module for the PI E-727 controller module's confocal scanner functionality.

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
from pipython import GCSDevice, pitools

from core.module import Base, Connector, ConfigOption
from interface.confocal_scanner_interface import ConfocalScannerInterface

class ConfocalScannerPI_E727(Base, ConfocalScannerInterface):

    """ Confocal scanner for PI E727 controller.
    """

    _modclass = 'ConfocalScannerPI_E727'
    _modtype = 'hardware'

    # connectors
    fitlogic = Connector(interface='FitLogic')

    # config
    # clock_frequency = ConfigOption('clock_frequency', 100, missing='warn') # I have no idea if this is reqired for E727 or not...
    # The clock freq. sets the "resolution" for the picture along the scanned line...
    #   ... it dosn't... maybe it is used for sync.-ing?

    E727_USBserial = ConfigOption('E727_USBserial', '0119019672', missing='warn')

    CONTROLLERNAME = 'E-727'
    STAGES = None # stage model not needed, but stage used is 'P-733.3CD'
    REFMODES = None

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        # Internal parameters
        self._line_length = None
        self._voltage_range = [-10., 10.]

        self._position_range = [[0, 100e-6], [0, 100e-6], [0, 10e-6], [0, 1e-6]]
        # self._position_range = [[self.get_position_range()[0][0], self.get_position_range()[0][1]],
        #                         [self.get_position_range()[1][0], self.get_position_range()[1][1]],
        #                         [self.get_position_range()[2][0], self.get_position_range()[2][1]], [0, 1e-6]]
        self._current_position = [0, 0, 0, 0][0:len(self.get_scanner_axes())]
        self._num_points = 500

    def on_activate(self):
        """ Initialisation performed during activation of the module.
            Connect to a PIPython device, using GCSDevice as context manager with "with".
            Different options for establishing connection to E-727 controller module, are given below.
        """
        self.e727_controller = GCSDevice(self.CONTROLLERNAME) #.ConnectUSB(serialnum=self.E727_USBserial)
        # pidevice.ConnectTCPIP(ipaddress='192.168.178.42')
        self.e727_controller.ConnectUSB(serialnum=self.E727_USBserial)
        # pidevice.ConnectRS232(comport=1, baudrate=115200)
        pitools.startup(self.e727_controller, stages=self.STAGES, refmodes=self.REFMODES)

    def on_deactivate(self):
        """ Deactivate properly the confocal scanner dummy.
        """
        self.reset_hardware()

    def reset_hardware(self):
        """ Resets the hardware, so the connection is lost and other programs
            can access it.

        @return int: error code (0:OK, -1:error)
        """
        if self.e727_controller.IsConnected() == True:
            self.e727_controller.close()

        if self.e727_controller.IsConnected() == False:
            self.log.warning('Scanning Device will be reset.')
            return 0
        else:
            return -1

    def get_position_range(self):
        """ Returns the physical range of the scanner.

        @return float [4][2]: array of 4 ranges with an array containing lower
                              and upper limit
        """
        rangemin = self.e727_controller.qTMN()
        rangemax = self.e727_controller.qTMX()
        return [[round(rangemin['1']*1e-6, 9), round(rangemax['1']*1e-6, 9)],
                [round(rangemin['2']*1e-6, 9), round(rangemax['2']*1e-6, 9)],
                [round(rangemin['3']*1e-6, 9), round(rangemax['3']*1e-6, 9)], [0, 1e-6]]
        # return [[rangemin['1']* 1e-6, rangemax['1']*1e-6],
        #         [rangemin['2']* 1e-6, rangemax['2']*1e-6],
        #         [rangemin['3']* 1e-6, rangemax['3']*1e-6], [0, 1e-6]]
        #return self._position_range

    def set_position_range(self, myrange=None):
        """ Sets the physical range of the scanner.

        @param float [4][2] myrange: array of 4 ranges with an array containing
                                     lower and upper limit

        @return int: error code (0:OK, -1:error)
        """
        if myrange is None:
            myrange = [[0, 100e-6], [0, 100e-6], [0, 10e-6], [0, 1e-6]]

        if not isinstance(myrange, (frozenset, list, set, tuple, np.ndarray, )):
            self.log.error('Given range is no array type.')
            return -1

        if len(myrange) != 4:
            self.log.error('Given range should have dimension 4, but has '
                    '{0:d} instead.'.format(len(myrange)))
            return -1

        for pos in myrange:
            if len(pos) != 2:
                self.log.error('Given range limit {1:d} should have '
                        'dimension 2, but has {0:d} instead.'.format(
                            len(pos),
                            pos))
                return -1
            if pos[0]>pos[1]:
                self.log.error('Given range limit {0:d} has the wrong '
                        'order.'.format(pos))
                return -1

        self._position_range = myrange

        return 0

    def set_voltage_range(self, myrange=None):
        """ Sets the voltage range of the E727 controller.

        @param float [2] myrange: array containing lower and upper limit

        @return int: error code (0:OK, -1:error)
        """
        if myrange is None:
            myrange = [-10., 10.]

        if not isinstance(myrange, (frozenset, list, set, tuple, np.ndarray, )):
            self.log.error('Given range is no array type.')
            return -1

        if len(myrange) != 2:
            self.log.error('Given range should have dimension 2, but has '
                    '{0:d} instead.'.format(len(myrange)))
            return -1

        if myrange[0]>myrange[1]:
            self.log.error('Given range limit {0:d} has the wrong '
                    'order.'.format(myrange))
            return -1

        if self.module_state() == 'locked':
            self.log.error('A Scanner is already running, close this one '
                    'first.')
            return -1

        self._voltage_range = myrange

        return 0

    def get_scanner_axes(self):
        """ Dummy scanner is always 3D cartesian.
        """

        return ['x', 'y', 'z', 'a']

    def get_scanner_count_channels(self):
        """ Counting channels in confocal: normal, negative and a ramp."""
        return ['Norm']

    def set_up_scanner_clock(self, clock_frequency=None, clock_channel=None):
        """ Configures the hardware clock of the NiDAQ card to give the timing.

        @param float clock_frequency: if defined, this sets the frequency of the
                                      clock
        @param str clock_channel: if defined, this is the physical channel of
                                  the clock

        @return int: error code (0:OK, -1:error)
        """

        if clock_frequency is not None:
            self._clock_frequency = float(clock_frequency)

        self.log.debug('ConfocalScannerDummy>set_up_scanner_clock')
        # time.sleep(0.2)
        return 0


    def set_up_scanner(self, counter_channels=None, sources=None,
                       clock_channel=None, scanner_ao_channels=None):
        """ Configures the actual scanner with a given clock.

        @param str counter_channel: if defined, this is the physical channel of
                                    the counter
        @param str photon_source: if defined, this is the physical channel where
                                  the photons are to count from
        @param str clock_channel: if defined, this specifies the clock for the
                                  counter
        @param str scanner_ao_channels: if defined, this specifies the analoque
                                        output channels

        @return int: error code (0:OK, -1:error)
        """

        self.log.debug('ConfocalScannerDummy>set_up_scanner')
        # time.sleep(0.2)
        return 0


    def scanner_set_position(self, x=None, y=None, z=None, a=None):
        """Move stage to x, y, z, a (where a is the fourth voltage channel).

        @param float x: postion in x-direction (volts)
        @param float y: postion in y-direction (volts)
        @param float z: postion in z-direction (volts)
        @param float a: postion in a-direction (volts)

        @return int: error code (0:OK, -1:error)
        """

        if self.module_state() == 'locked':
            self.log.error('A Scanner is already running, close this one first.')
            return -1
        # coord_list = [round(float(x),8), round(float(y),8), round(float(z),8), round(float(a),8)]
        # coord_list = [x, y, z, a]
        coord_list = [x*1e6, y*1e6, z*1e6, a*1e6]
        # print("coord_list:" + str(coord_list))
        # print(self._volt_to_position(coord_list))
        t0 = time.clock()
        for axis, target in zip(self.e727_controller.axes, coord_list[:-1]):
            self.e727_controller.MOV(axis, target)
        pitools.waitontarget(self.e727_controller, axes=axis)
        self._current_position = [x, y, z, a][0:len(self.get_scanner_axes())]
        wait_time = time.clock() - t0
        time.sleep(8*wait_time)
        # for axis, target in zip(self.e727_controller.axes, self._volt_to_position(coord_list)[:-1]):
        #     self.e727_controller.MOV(axis, target*1e6)
        #pitools.waitontarget(self.e727_controller, axes=axis)
        #time.sleep(0.05)
        #self._current_position = [x, y, z, a][0:len(self.get_scanner_axes())]
        # print("current_pos:" + str(self._current_position))
        # print("get_scanner_pos" + str(self.get_scanner_position()))
        return 0

    def get_scanner_position(self):
        """ Get the current position of the scanner hardware.

        @return float[]: current position in (x, y, z, a).
        """

        curpos = self.e727_controller.qPOS()
        # print(curpos)
        self._current_position = [curpos['1'], curpos['2'], curpos['3'], self._current_position[-1]]

        return self._current_position[0:len(self.get_scanner_axes())]

    def _set_up_line(self, length=100):
        """ Sets up the analoque output for scanning a line.

        @param int length: length of the line in pixel

        @return int: error code (0:OK, -1:error)
        """

        self._line_length = length
        self.log.debug('ConfocalScannerPI_E-727>set_up_line')
        return 0

    def scan_line(self, line_path=None, pixel_clock=False):
        """ Scans a line and returns the counts on that line.

        @param float[][4] line_path: array of 4-part tuples defining the voltage points
        @param bool pixel_clock: whether we need to output a pixel clock for this line

        @return float[]: the photon counts per second
        """

        if not isinstance(line_path, (frozenset, list, set, tuple, np.ndarray, )):
            self.log.error('Given voltage list is no array type.')
            return np.array([[-1.]])

        if np.shape(line_path)[1] != self._line_length:
            self._set_up_line(np.shape(line_path)[1])

        self._current_position = list(line_path[:, -1])
        # time.sleep(0.001)
        return np.array([[i] for i in range(self._line_length)])#.transpose()


    def close_scanner(self):
        """ Closes the scanner and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """

        self.log.debug('ConfocalScannerDummy>close_scanner')
        return 0

    def close_scanner_clock(self, power=0):
        """ Closes the clock and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """

        self.log.debug('ConfocalScannerDummy>close_scanner_clock')
        return 0


    def _volt_to_position(self, volts=None):
        """ Converts a set of position pixels to actual voltages.

        @param float[][n]: array of n-part tuples of corresponding voltages

        @return float[] positions: array of n-part tuples defining the pixels

        The positions is typically a matrix like
            [[x_values], [y_values], [z_values], [a_values]]
            but x, xy, xyz and xyza are allowed formats.
        """

        if not isinstance(volts, (frozenset, list, set, tuple, np.ndarray, )):
            self.log.error('Given voltage list is no array type.')
            return np.array([np.NaN])

        poslist = []
        for i, volt in enumerate(volts):
            poslist.append(
                round((self._position_range[i][1] - self._position_range[i][0])
                / (self._voltage_range[1] - self._voltage_range[0])
                * (volt - self._voltage_range[0])
                + self._position_range[i][0], 9)
            )
        positions = np.vstack(poslist)

        for i, pos in enumerate(positions):
            if pos.min() < self._position_range[i][0] or pos.max() > self._position_range[i][1]:
                self.log.error(
                    'Positions ({0}, {1}) exceed the limit, the positions have to '
                    'be adjusted to stay in the given range.'.format(pos.min(), pos.max()))
                return np.array([np.NaN])

        positions = [i[0] for i in positions]
        # print(positions)

        return positions
