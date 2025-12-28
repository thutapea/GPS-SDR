This is a C++ implementation of a GPS SDR Receiver. It's created with the help of "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach" by K Borre. You'll also need the FFTW library and Make. 

The acquisition script only works on 8-bit data, so using an SDR like the HackRF or RTL-SDR Dongle, wouldn't work. However, the change should be minimal.
The acquisition script was tested with a prerecorded bit stream, and the tracking script is a WIP.
