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

from core.module import Base, Connector, ConfigOption
from interface.fast_counter_interface import FastCounterInterface
from interface.pulser_interface import PulserInterface , PulserConstraints
# from interface.slow_counter_interface import SlowCounterInterface
# from interface.slow_counter_interface import SlowCounterConstraints
# from interface.slow_counter_interface import CountingMode
from interface.odmr_counter_interface import ODMRCounterInterface
# from interface.confocal_scanner_interface import ConfocalScannerInterface


class ConfocalScanner_PI_Swabian_Interfuse(Base, FastCounterInterface, PulserInterface):

    """This is the Interface class to define the controls for the simple
    microwave hardware.
    """
    _modclass = 'Swabian I/O Interfuse'
    _modtype = 'hardware'

    # connectors
    fitlogic = Connector(interface='FitLogic')
    pulsestreamer = Connector(interface='PulserInterface')
    timetagger = Connector(interface='FastCounterInterface')

    # config options
    _clock_frequency = ConfigOption('clock_frequency', 100, missing='warn')


    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        # Internal parameters


    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """

        self._fit_logic = self.fitlogic()
        self._output_hw = self.pulsestreamer()
        self._input_hw = self.timetagger()

        self._input_hw.on_activate()
        self._output_hw.on_activate()
        # configure(self, bin_width_s, record_length_s, number_of_gates=0)
        self._input_hw.configure(self._output_hw._channel_detect)

    def on_deactivate(self):

        self._input_hw.on_deactivate()
        self._output_hw.on_deactivate()

# === timetagger / input functionality ===
    def get_constraints(self):
        """ Retrieve the hardware constrains from the Fast counting device.

        @return dict: dict with keys being the constraint names as string and
                      items are the definition for the constaints.
        """

        self._input_hw.get_constraints()

    def configure(self, bin_width_s, record_length_s, number_of_gates=0):
        """ Configuration of the fast counter.

        @param float bin_width_s: Length of a single time bin in the time
                                  trace histogram in seconds.
        @param float record_length_s: Total length of the timetrace/each
                                      single gate in seconds.
        @param int number_of_gates: optional, number of gates in the pulse
                                    sequence. Ignore for not gated counter.

        @return tuple(binwidth_s, record_length_s, number_of_gates):
                    binwidth_s: float the actual set binwidth in seconds
                    gate_length_s: the actual record length in seconds
                    number_of_gates: the number of gated, which are accepted, None if not-gated
        """
        self._input_hw.configure(bin_width_s, record_length_s, number_of_gates)

    def get_status(self):
        """ Receives the current status of the Fast Counter and outputs it as
            return value.

        0 = unconfigured
        1 = idle
        2 = running
        3 = paused
        -1 = error state
        """
        self._input_hw.get_status()

    def start_measure(self):
        """ Start the fast counter. """
        self._input_hw.start_measure()

    def stop_measure(self):
        """ Stop the fast counter. """
        self._input_hw.stop_measure()

    def pause_measure(self):
        """ Pauses the current measurement.

        Fast counter must be initially in the run state to make it pause.
        """
        self._input_hw.pause_measure()

    def continue_measure(self):
        """ Continues the current measurement.

        If fast counter is in pause state, then fast counter will be continued.
        """
        self._input_hw.continue_measure()

    def is_gated(self):
        """ Check the gated counting possibility.

        @return bool: Boolean value indicates if the fast counter is a gated
                      counter (TRUE) or not (FALSE).
        """
        self._input_hw.is_gated()

    def get_binwidth(self):
        """ Returns the width of a single timebin in the timetrace in seconds.

        @return float: current length of a single bin in seconds (seconds/bin)
        """
        self._input_hw.get_binwidth()


    def get_data_trace(self):
        """ Polls the current timetrace data from the fast counter.

        Return value is a numpy array (dtype = int64).
        The binning, specified by calling configure() in forehand, must be
        taken care of in this hardware class. A possible overflow of the
        histogram bins must be caught here and taken care of.
        If the counter is NOT GATED it will return a 1D-numpy-array with
            returnarray[timebin_index]
        If the counter is GATED it will return a 2D-numpy-array with
            returnarray[gate_index, timebin_index]
        """
        self._input_hw.get_data_trace()


# === PulseStreamer / output ===

    def get_constraints(self):
        """
        Retrieve the hardware constrains from the Pulsing device.

        @return constraints object: object with pulser constraints as attributes.
        """
        self._output_hw.get_constraints()

    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.pulser_on()

    def pulser_off(self):
        """ Switches the pulsing device off.

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.pulser_off()

    def upload_asset(self, asset_name=None):
        """ Upload an already hardware conform file to the device mass memory.
            Also loads these files into the device workspace if present.
            Does NOT load waveforms/sequences/patterns into channels.

        @param asset_name: string, name of the ensemble/sequence to be uploaded

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.upload_asset(asset_name)

    def load_asset(self, asset_name, load_dict=None):
        """ Loads a sequence or waveform to the specified channel of the pulsing device.
        For devices that have a workspace (i.e. AWG) this will load the asset from the device
        workspace into the channel.
        For a device without mass memory this will transfer the waveform/sequence/pattern data
        directly to the device so that it is ready to play.

        @param str asset_name: The name of the asset to be loaded

        @param dict load_dict:  a dictionary with keys being one of the available channel numbers
                                and items being the name of the already sampled waveform/sequence
                                files.
                                Examples:   {1: rabi_Ch1, 2: rabi_Ch2}
                                            {1: rabi_Ch2, 2: rabi_Ch1}
                                This parameter is optional. If none is given then the channel
                                association is invoked from the file name, i.e. the appendix
                                (_ch1, _ch2 etc.)

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.load_asset(asset_name, load_dict)

    def get_loaded_asset(self):
        """ Retrieve the currently loaded asset name of the device.

        @return str: Name of the current asset ready to play. (no filename)
        """
        self._output_hw.get_loaded_asset()

    def clear_all(self):
        """ Clears all loaded waveforms from the pulse generators RAM/workspace.

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.clear_all()

    def get_status(self):
        """ Retrieves the status of the pulsing hardware

        @return (int, dict): tuple with an interger value of the current status and a corresponding
                             dictionary containing status description for all the possible status
                             variables of the pulse generator hardware.
        """
        self._output_hw.get_status()

    def get_sample_rate(self):
        """ Get the sample rate of the pulse generator hardware

        @return float: The current sample rate of the device (in Hz)
        """
        self._output_hw.get_sample_rate()

    def set_sample_rate(self, sample_rate):
        """ Set the sample rate of the pulse generator hardware.

        @param float sample_rate: The sampling rate to be set (in Hz)

        @return float: the sample rate returned from the device (in Hz).
        """
        self._output_hw.set_sample_rate(sample_rate)

    def get_analog_level(self, amplitude=None, offset=None):
        """ Retrieve the analog amplitude and offset of the provided channels.

        @param list amplitude: optional, if the amplitude value (in Volt peak to peak, i.e. the
                               full amplitude) of a specific channel is desired.
        @param list offset: optional, if the offset value (in Volt) of a specific channel is
                            desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel descriptor string
                               (i.e. 'a_ch1') and items being the values for those channels.
                               Amplitude is always denoted in Volt-peak-to-peak and Offset in volts.
        """
        self._output_hw.get_analog_level(amplitude, offset)

    def set_analog_level(self, amplitude=None, offset=None):
        """ Set amplitude and/or offset value of the provided analog channel(s).

        @param dict amplitude: dictionary, with key being the channel descriptor string
                               (i.e. 'a_ch1', 'a_ch2') and items being the amplitude values
                               (in Volt peak to peak, i.e. the full amplitude) for the desired
                               channel.
        @param dict offset: dictionary, with key being the channel descriptor string
                            (i.e. 'a_ch1', 'a_ch2') and items being the offset values
                            (in absolute volt) for the desired channel.

        @return (dict, dict): tuple of two dicts with the actual set values for amplitude and
                              offset for ALL channels.
        """
        self._output_hw.set_analog_level(amplitude, offset)

    def get_digital_level(self, low=None, high=None):
        """ Retrieve the digital low and high level of the provided/all channels.

        @param list low: optional, if the low value (in Volt) of a specific channel is desired.
        @param list high: optional, if the high value (in Volt) of a specific channel is desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel descriptor strings
                               (i.e. 'd_ch1', 'd_ch2') and items being the values for those
                               channels. Both low and high value of a channel is denoted in volts.
        """
        self._output_hw.get_digital_level(low, high)

    def set_digital_level(self, low=None, high=None):
        """ Set low and/or high value of the provided digital channel.

        @param dict low: dictionary, with key being the channel descriptor string
                         (i.e. 'd_ch1', 'd_ch2') and items being the low values (in volt) for the
                         desired channel.
        @param dict high: dictionary, with key being the channel descriptor string
                          (i.e. 'd_ch1', 'd_ch2') and items being the high values (in volt) for the
                          desired channel.

        @return (dict, dict): tuple of two dicts where first dict denotes the current low value and
                              the second dict the high value for ALL digital channels.
                              Keys are the channel descriptor strings (i.e. 'd_ch1', 'd_ch2')
        """
        self._output_hw.set_digital_level(low, high)

    def get_active_channels(self, ch=None):
        """ Get the active channels of the pulse generator hardware.

        @param list ch: optional, if specific analog or digital channels are needed to be asked
                        without obtaining all the channels.

        @return dict:  where keys denoting the channel string and items boolean expressions whether
                       channel are active or not.
        """
        self._output_hw.get_active_channels(ch)

    def set_active_channels(self, ch=None):
        """ Set the active channels for the pulse generator hardware.

        @param dict ch: dictionary with keys being the analog or digital string generic names for
                        the channels (i.e. 'd_ch1', 'a_ch2') with items being a boolean value.
                        True: Activate channel, False: Deactivate channel

        @return dict: with the actual set values for ALL active analog and digital channels
        """
        self._output_hw.set_active_channels(ch)

    def get_uploaded_asset_names(self):
        """ Retrieve the names of all uploaded assets on the device.

        @return list: List of all uploaded asset name strings in the current device directory.
                      This is no list of the file names.

        Unused for pulse generators without sequence storage capability (PulseBlaster, FPGA).
        """
        self._output_hw.get_uploaded_asset_names()

    def get_saved_asset_names(self):
        """ Retrieve the names of all sampled and saved assets on the host PC. This is no list of
            the file names.

        @return list: List of all saved asset name strings in the current
                      directory of the host PC.
        """
        self._output_hw.get_saved_asset_names()

    def delete_asset(self, asset_name):
        """ Delete all files associated with an asset with the passed asset_name from the device
            memory (mass storage as well as i.e. awg workspace/channels).

        @param str asset_name: The name of the asset to be deleted
                               Optionally a list of asset names can be passed.

        @return list: a list with strings of the files which were deleted.

        Unused for pulse generators without sequence storage capability (PulseBlaster, FPGA).
        """
        self._output_hw.delete_asset(asset_name)

    def set_asset_dir_on_device(self, dir_path):
        """ Change the directory where the assets are stored on the device.

        @param str dir_path: The target directory

        @return int: error code (0:OK, -1:error)

        Unused for pulse generators without changeable file structure (PulseBlaster, FPGA).
        """
        self._output_hw.set_asset_dir_on_device(dir_path)

    def get_asset_dir_on_device(self):
        """ Ask for the directory where the hardware conform files are stored on the device.

        @return str: The current file directory

        Unused for pulse generators without changeable file structure (i.e. PulseBlaster, FPGA).
        """
        self._output_hw.get_asset_dir_on_device()

    def get_interleave(self):
        """ Check whether Interleave is ON or OFF in AWG.

        @return bool: True: ON, False: OFF

        Will always return False for pulse generator hardware without interleave.
        """
        self._output_hw.get_interleave()

    def set_interleave(self, state=False):
        """ Turns the interleave of an AWG on or off.

        @param bool state: The state the interleave should be set to
                           (True: ON, False: OFF)

        @return bool: actual interleave status (True: ON, False: OFF)

        Note: After setting the interleave of the device, retrieve the
              interleave again and use that information for further processing.

        Unused for pulse generator hardware other than an AWG.
        """
        self._output_hw.set_interleave(state)

    def tell(self, command):
        """ Sends a command string to the device.

        @param string command: string containing the command

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.tell(command)

    def ask(self, question):
        """ Asks the device a 'question' and receive and return an answer from it.

        @param string question: string containing the command

        @return string: the answer of the device to the 'question' in a string
        """
        self._output_hw.ask(question)

    def reset(self):
        """ Reset the device.

        @return int: error code (0:OK, -1:error)
        """
        self._output_hw.reset()

    def has_sequence_mode(self):
        """ Asks the pulse generator whether sequence mode exists.

        @return: bool, True for yes, False for no.
        """
        self._output_hw.has_sequence_mode()