# LongHCPulse
## V. 1.3.2

Code for processing long pulse heat capacity data taken on a Quantum Design PPMS.
Please cite https://arxiv.org/pdf/1705.07129.pdf

The main file here is "PPMS_LongPulse.py". This contains the class LongHCPulse, which computes and plots heat capacity.
This is the file you should import into your python scripts.

"Minimal_Working_Example.ipynb" is an ipython notebook giving a simple introduction to the code and how to use the plotting functions.

"YbTiO_HeatCapacity.ipynb" is an ipython notebook showing how LongHCPulse was used in a more complicated way to process heat capacity data in DOI:10.1103/PhysRevLett.119.127201.

"DRPuck27.cal", "Yb2Ti2O7_longpulse.raw", and "YbTiO_MvsB_100mK.txt" are all data files used in YbtiO_HeatCapacity.ipynb.

A runnable version of YbTiO_HeatCapacity.ipynb can be found at:
[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/asche1/longhcpulse)
 [Note: MyBinder.org currently not working 6/30/17]

### Updates in version 1.3.2 (June 7, 2018)
- Fixed bug in savetraces function so that python closes the file properly.
- Added NonPPMS version of LongHCPulse which can process data without a calibration file, so that LongHCPulse can be used with data not taken on a Quantum Design PPMS.

### Updates in version 1.3.1 (October 16, 2017)
- Updated meshgrid so that heat pulses can be included in the binning (with the command useHeatPulses=True).
- Updated lineplotCombine so that heat pulses can be included in the final data (useHeatPulses=True), or so that heat pulses can be used instead of cooling pulses (onlyHeatPulses=True).

### Updates in version 1.3 (October 4, 2017)
- Re-wrote the combineTraces function so that it averages all overlapping regions of data, instead of concatenating, sorting, and then applying a moving average. The result is fewer jagged discontinuities when combining multiple heating and cooling pulses.

### Updates in version 1.2.1 (August 28, 2017)
- Added option in "lineplot" command to plot heating pulse data as well (using "plotHeatPulses=True")
- Added command to save processed data into a comma-separated text file: "savetraces". This allows further processing using external software.

### Updates in version 1.2 (June 30, 2017): 
- Updated LongHCPulse.py to be compatible with Python 3.6 as well as Python 2.7
- Fixed bug so that code runs properly on Windows
- Added minimal working example that includes an explanation of the plotting functions
