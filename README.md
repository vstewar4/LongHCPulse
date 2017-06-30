# LongHCPulse
## V. 1.2

Code for processing long pulse heat capacity data taken on a Quantum Design PPMS
Please cite https://arxiv.org/pdf/1705.07129.pdf

The main file here is "PPMS_LongPulse.py". This contains the class LongHCPulse, which computes and plots heat capacity.
This is the file you should import into your python scripts.

YbTiO_HeatCapacity.ipynb is an ipython notebook that gives an example of LongHCPulse being used to process heat capacity data
in https://arxiv.org/pdf/1703.06904.pdf.

"DRPuck27.cal", "Yb2Ti2O7_longpulse.raw", and "YbTiO_MvsB_100mK.txt" are all data files used in YbtiO_HeatCapacity.ipynb.

A runnable version of YbTiO_HeatCapacity.ipynb can be found at:
[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/asche1/longhcpulse)
[Note: MyBinder.org currently not working 6/30/17]

## V. 1.2 updates (June 30, 2017): 
- Updated LongHCPulse.py to be compatible with Python 3.6 as well as Pythin 2.7
- Fixed bug so that code runs properly on Windows
