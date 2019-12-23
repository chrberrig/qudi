# -*- coding: utf-8 -*-
"""
This module acts like a laser.

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
from interface.simple_laser_interface import SimpleLaserInterface
from interface.simple_laser_interface import LaserState
from interface.simple_laser_interface import ShutterState
from interface.simple_laser_interface import ControlMode

import serial


class DLnsec(Base, SimpleLaserInterface):
    """
    Swabian instrumants / LABS electronics DLnsec laser

    Config ex:

    DLnSec_laser:
        module.Class: 'laser.DLnsec_laser.DLnsec'
        com_port: 'COM3'
    """

    # _modclass = 'DLnsec_laser'
    # _modtype = 'hardware'

    _com_port = ConfigOption('com_port', missing='error')

    # eol = '\r'
    # _model_name = 'UNKNOWN'

    def __init__(self, **kwargs):
        """ """
        super().__init__(**kwargs)
        self.lstate = LaserState.OFF
        self.shutter = ShutterState.NOSHUTTER
        self.mode = ControlMode.MIXED
        self.current_setpoint = 0
        self.power_setpoint = 0

# --- Activation and connection ---

    def on_activate(self):
        """ Activate module.
        """
        self.ser = serial.Serial(self._com_port, timeout=1)
        self.ser.write_timeout = 1
        self.ser.read_timeout = 1
        connected = self.connect_laser()

        if not connected:
            self.log.error('Laser does not seem to be connected.')
            return -1
        else:
            self._model_name = self.read_serial('*IDN')
            self.log.info('Laser {0} now connected via. serial {1}'.format(self._model_name, self._com_port))
            self.get_current()
            return 0

    def on_deactivate(self):
        """ Deactivate module.
        """
        self.disconnect_laser()

    def connect_laser(self):
        """ Connect to Instrument.
        @return bool: connection success
        """
        response = self.read_serial('*IDN')
        if 'DLNS' in response.upper():
            return True
        else:
            return False

    def disconnect_laser(self):
        """ Close the connection to the instrument.
        """
        self.off()
        self.ser.close()

    def write_serial(self, command):
        """ Write command to serialport of instrument.
        """
        assert type(command) == str
        self.ser.write(command.encode() + b'\n')

    def read_serial(self, command):
        """ Write command to serialport of instrument and get response.
        """
        assert type(command) == str
        self.write_serial(command)
        response = self.ser.readline().strip()
        self.ser.read()  # fetch the additional '\r' char
        response = response.decode('ascii', 'ignore')
        return response

# --- Power Settings ---

    def get_power_range(self):
        """ Return optical power range
            @return (float, float): power range
                value taken from spec. table
        """
        return (0, 0.040)

    def get_power(self):
        """ Return laser power
            @return float: Laser power in watts
                note that setpoint and actual val. is identically defined
        """
        if self.lstate == LaserState.ON:
            self.current_setpoint = float(self.read_serial('PWR?'))
            self.power_setpoint = self.current_to_power(self.current_setpoint)
            return self.power_setpoint
        else:
            return float(0)

    def get_power_setpoint(self):
        """ Return optical power setpoint.
            @return float: power setpoint in watts
                note that setpoint and actual val. is identically defined
        """
        # return int(self.read_serial('PWR?'))*self.get_power_range()[1]
        return self.power_setpoint

    def set_power(self, power):
        """ Set power setpoint.
            @param float power: power setpoint
            @return float: actual new power setpoint
        """
        assert type(power) == float and self.get_power_range()[0] <= power <= self.get_power_range()[1]
        current = self.power_to_current(power)
        self.set_current(current)
        self.power_setpoint = self.current_to_power(current)
        # self.get_extra_info()
        return self.power_setpoint

# --- Current settings ---

    def get_current_unit(self):
        """ Get unit for laser current.
            @return str: unit
        """
        return '%'

    def get_current_range(self):
        """ Get laser current range.
            @return (float, float): laser current range

        """
        return (0, 100)

    def get_current(self):
        """ Get current laser current
            @return float: laser current in current current units
                note that setpoint and actual val. is identically defined
        """
        if self.lstate == LaserState.ON:
            self.current_setpoint = float(self.read_serial('PWR?'))
            self.power_setpoint = self.current_to_power(self.current_setpoint)
            return self.current_setpoint
        else:
            return float(0)

    def get_current_setpoint(self):
        """ Get laser curent setpoint
            @return float: laser current setpoint
                note that setpoint and actual val. is identically defined
        """
        return self.current_setpoint

    def set_current(self, current):
        """ Set laser current setpoint
            @prarm float current: desired laser current setpoint
            @return float: actual laser current setpoint
        """
        assert type(current) == float and self.get_current_range()[0] <= current <= self.get_current_range()[1]
        self.write_serial('PWR ' + str(int(current)))
        self.current_setpoint = float(self.read_serial('PWR?'))
        self.power_setpoint = self.current_to_power(self.power_setpoint)
        return self.current_setpoint

# --- Control modes and laser states ---

    def allowed_control_modes(self):
        """ Get supported control modes
            @return list(): list of supported ControlMode
        """
        # return [ControlMode.POWER, ControlMode.CURRENT]
        return ControlMode.MIXED

    def get_control_mode(self):
        """ Get the currently active control mode
            @return ControlMode: active control mode
        """
        return self.mode

    def set_control_mode(self, control_mode):
        """ Set the active control mode
            @param ControlMode control_mode: desired control mode
            @return ControlMode: actual active ControlMode
        """
        # self.mode = control_mode
        # return self.mode
        return ControlMode.MIXED

    def on(self):
        """ Turn on laser.
            @return LaserState: actual laser state
        """
        self.write_serial('*ON')
        self.write_serial('LASE')
        self.get_extra_info()
        self.lstate = LaserState.ON
        return self.lstate

    def off(self):
        """ Turn off laser.
            @return LaserState: actual laser state
        """
        self.write_serial('STOP')
        self.write_serial('*OFF')
        self.get_extra_info()
        self.lstate = LaserState.OFF
        return self.lstate

    def get_laser_state(self):
        """ Get laser state
            @return LaserState: actual laser state
        """
        return self.lstate

    def set_laser_state(self, status):
        """ Set desited laser state.
        @param LaserState status: desired laser state
        @return LaserState: actual laser state
        """
        # TODO: this is big. cannot be called without having LaserState,
        #       which is only defined in the simple laser interface.
        #       I think this should be a private method.
        if self.get_laser_state() != status:
            if status == LaserState.ON:
                self.on()
            elif status == LaserState.OFF:
                self.off()
        return self.get_laser_state()

    def get_shutter_state(self):
        """ Get laser shutter state
            @return ShutterState: actual laser shutter state
        """
        return self.shutter

    def set_shutter_state(self, state):
        """ Set laser shutter state.
            @param ShutterState state: desired laser shutter state
            @return ShutterState: actual laser shutter state
        """
        # self.shutter = state
        return self.shutter

# --- Temperature settings ---

    def get_temperatures(self):
        """ Get all available temperatures.
            @return dict: dict of temperature name and value in degrees Celsius
        """
        return {}

    def set_temperatures(self, temps):
        """ Set temperatures for lasers with tunable temperatures.
            @return {}: empty dict, dummy not a tunable laser
        """
        return {}

    def get_temperature_setpoints(self):
        """ Get temperature setpoints.
            @return dict: temperature setpoints for temperature tunable lasers
        """
        return {}

    def get_extra_info(self):
        """ Multiple lines of diagnostic information
            @return str: much laser, very useful
        """
        response = self.read_serial('*IDN')
        model, sernum = response.split('_')
        sernum = model[5:] + '_' + sernum

        extra = ('System Identification: '   +  self.read_serial('*IDN') + '%' + '\n'
                'System Model Name: '        + model[:4] + '\n'
                'System Serial Number: '     + sernum + '\n'
                'System Wavelength: '        + model[5:] + '\n'
                'System Power Rating: '      + self.read_serial('POW?') + '%' + '\n'
                'System Memory writes: '     + self.read_serial('NSAV?') + '%' + '\n'
                )
#         return model.upper(), sernum.upper()
        return extra

    def get_error(self):
        return self.read_serial('ERR?')

    def current_to_power(self, current):
        """
        :param current float: integer from 0 to 100, describing currnet in % of max current.
        :return float: power corresponding to current
        """
        power = float(self.get_power_range()[1] * current / self.get_current_range()[1])
        return power

    def power_to_current(self, power):
        """

        :param power float: float corresponding to power output in W
        :return float: integer from 0 to 100, describing currnet in % of max current.
        """
        current = float(int(self.get_current_range()[1] * power / self.get_power_range()[1]))
        return current