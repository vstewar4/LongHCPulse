# Code to extract heat capacity from a long pulse in PPMS measurement
# Allen Scheie
# August, 2016

import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
import copy
import pickle
cc = matplotlib.colors.ColorConverter()


class LongHCPulse:
	"""For importing long pulse heat capacity data"""
	def __init__(self,datafile,calfile=None,sampmass=None,molarmass=1,scaleshortpulse=1,
				AdiabaticCriterion=0.1):
		# If only rawfile is specified, it assumes that you're feeding it a pickle object
		print('**************** LongHCPulse v 1.3.3 *****************\n'+\
			' please cite   https://doi.org/10.1007/s10909-018-2042-9\n'+\
			'******************************************************')


		if ((calfile == None and sampmass == None) and datafile.endswith('.pickle')):
			with open(datafile, 'rb') as pf:
				pickdict = pickle.load(pf)
			self._importPickle(**pickdict)

		else: self._importData(datafile,calfile,sampmass,molarmass,scaleshortpulse,
								AdiabaticCriterion)



	def _importPickle(self, **entries):
		'''file must have been saved with savefile function'''
		self.__dict__.update(entries)



	def _importData(self,rawfile,calfile,sampmass,molarmass,scaleshortpulse, 
					AdiabaticCriterion):
		# Import thermal conductivity and thermometer resistivity values from calibration file
		print(" - Importing data...")
		self.importCalibration(calfile)

		# Get the number of points per heating/cooling curve, and the number of curves in the file
		numpts, numcurves = self._countpoints(rawfile)

		DataHeader = False
		DataLines = False	#index for keeping track of when data begins in file
		i = 0 		# Index for keeping track of data points
		j = 0-1 	# Index for keeping track of heating curves

		self.rawdata = np.empty((numcurves,3,numpts,2))*np.nan  # Array which holds all raw data in file
			# First index: pulse number
			# Second index: data type (Time, Sample Temperature, Heater Power)
			# Third index: actual data points
			# Fourth index: 0 = heating curve, 1 = cooling curve
		self.Tb = np.empty((numcurves))*np.nan	# Array of system temperatures
		self.Tsamp = np.empty((numcurves))*np.nan	# Array of average sample temperatures
		self.Tstart = np.empty((numcurves))*np.nan
		self.Bfield = np.empty((numcurves))*np.nan	# Array for holding magnetic field
		self.ShortPulse = np.zeros((numcurves))	# gives HC if pulse is a short pulse, zero if long pulse
		self.AddendaHC = np.zeros((2,numcurves))	# Subtracted from HC at the very end.
		
		self.ThermCondWire = np.zeros((2,numcurves))	# Used to compare 
		##self.ThermCondWire = []

		# Import data
		for line in open(rawfile):

			if DataLines==True:
				d = line.strip('\n').split(",")
				if d[1] == '':  #only accept lines for which there are no comments
					if d[4] == '0':		
						kk = 1 # If heater power is zero, it's a cooling curve.
						self.rawdata[j,0,i,kk] = d[0] #time (s)
						#self.rawdata[j,1,i,kk] = d[3] #Temp (K)  BAD!! UNCORRECTED!!
						self.rawdata[j,2,i,kk] = d[4] #Heater Power (W)
						self.rawdata[j,1,i,kk] = self._resisToTemp(float(d[2]), self.Bfield[j])
						i+=1
					else:	
						kk = 0  # heating curve.
						self.rawdata[j,0,ii,kk] = d[0] #time (s)
						#self.rawdata[j,1,ii,kk] = d[3] #Temp (K) BAD!! UNCORRECTED!!
						self.rawdata[j,2,ii,kk] = d[4] #Heater Power (W)
						self.rawdata[j,1,ii,kk] = self._resisToTemp(float(d[2]), self.Bfield[j])

						#  Attempt to correct heating pulses for improper power values. Didn't work.
						# Assumes that voltage is measured across heater 
						startT = self.rawdata[j,1,0,0]
						#self.rawdata[j,2,ii,kk] *= self._HeaterRes(startT) /\
						#		 self._HeaterRes(self.rawdata[j,1,ii,kk])
						ii +=1

			# find the information needed to compute heat capacity
			if DataLines==False:
				if "SystemTemp=" in line:
					self.Tb[j] = float(line.strip('\n').split(',')[1][11:])
					# Note that the log files typically show this to be constant over each pulse
				#if "StableStartTemperature=" in line:
				#	self.Tb[j] = float(line.split(',')[1][len("StableStartTemperature="):-2])


				if "TempMinMidMax=" in line:
					self.Tstart[j] = float(line.strip('\n').split(',')[1][len("TempMinMidMax="):].split('|')[0])


				if "Field=" in line:
					self.Bfield[j] = round(float(line.strip('\n').split(',')[1][6:]),0)

				if "SampHC=" in line:
					SampHC = float(line.strip('\n').split(',')[1][7:])

				# determine if it's long pulse or short pulse
				if "SampleTemp=" in line:
					self.Tsamp[j] = float(line.strip('\n').split(',')[1][11:])
				if "TempRise=" in line:
					TempRise = float(line.strip('\n').split(',')[1][9:])
					Tratio = TempRise / self.Tsamp[j]
				# CRITERION FOR ADIABATIC PULSE: Trise/Tsamp < AdiabaticCriterion
				# If this is true, then save the sample heat capacity. If not, value = 0
				if ("SampHC=" in line and Tratio < AdiabaticCriterion):
					self.ShortPulse[j] = float(line.strip('\n').split(',')[1][7:])

				if "AddendaHC=" in line:
					self.AddendaHC[1,j] = float(line.strip('\n').split(',')[1][10:])	# microJ/K
					self.AddendaHC[0,j] = self.Tsamp[j]
					# Convert to J/K from microJ/K
					self.AddendaHC[1,j] *= 1e-6

				if "ThermCondWire=" in line: #and Tratio < AdiabaticCriterion):
					self.ThermCondWire[1,j] = float(line.strip('\n').split(',')[1][len("ThermCondWire="):])
					self.ThermCondWire[0,j] = self.Tsamp[j]
					#self.ThermCondWire.append( [ self.Tsamp[j],
					#	float(line.split(',')[1][len("ThermCondWire="):-2]) ])

			# These two statements check for data lines
			if "END:PULSE:PARAMS" in line:	
				DataLines=True
				# prepare to take data
				i=0
				ii=0
			if "BEGIN:PULSE:PARAMS" in line:
				j+=1
				DataLines=False
				# Print progress
				sys.stdout.write('\r %d%% ' % (100*j/numcurves))
				sys.stdout.flush() # important for printing progress

		print("\r 100%")


		# Smooth data
		self._smooth(n=5)
		self.sampmass = sampmass
		self.molarmass = molarmass

		# Set all "zero-field values" to zero
		self.Bfield[np.where(np.around(self.Bfield,0) == 0)] = 0
		# Round all Bfield values to nearest 10 Oe
		self.Bfield = np.around(self.Bfield,-1)

		# Scale short pulse by user-provided factor (depends on how the measurement
		# was taken).
		self.ShortPulse *=scaleshortpulse
		self.shortPulseAverage()

		# Set scaled Bfield to be Bfield (this is changed if a demagnetization factor
		# is applied.)
		self.ScaledBfield = copy.deepcopy(self.Bfield)

		#Define label array for later
		self.labels = []
		self.shortpulselabels = []
		self.Bflag = []
		self.shortpulseBflag = []
		self.Blabels = []



	def _countpoints(self, infile):
		numberpoints = 0
		numberpulses = 0
		for line in open(infile):
			if "NBinsOn" in line:
				pointson = int(line.strip('\n').split(",")[1][8:])
				numberpulses += 1
			if "NBinsOff" in line:
				pointsoff = int(line.strip('\n').split(",")[1][9:])
				if (pointson) > numberpoints:
					numberpoints = pointson
		return numberpoints, numberpulses


	def _smooth(self, n):
		'''Take moving average of heating and cooling pulses'''
		self.smoothedData = np.zeros((len(self.rawdata[:,1,:,:]),3,len(self.rawdata[0,1,:,:]),2))
		for i in range(len(self.rawdata[:,1,:,:])):
			self.smoothedData[i,1,:,0] = self._movingaverage(self.rawdata[i,1,:,0],n)
			self.smoothedData[i,1,:,1] = self._movingaverage(self.rawdata[i,1,:,1],n)
			self.smoothedData[i,0,:,:] = self.rawdata[i,0,:,:]

	def heatcapacity(self, smoothlevel=0, StaticOffset = 0.1):
		print(" - Computing Heat Capacity...")
		# Initialize arrays
		lenSD = len(self.smoothedData)
		self.HC = np.zeros((lenSD, len(self.smoothedData[0,0]),2))
		self.T = np.zeros((lenSD, len(self.smoothedData[0,0]),2))
		self.HC_uncertainty = np.zeros((lenSD, len(self.smoothedData[0,0]),2))

		# Sort the addenda heat capacity curve (so it can be interpolated)
		self.AddendaHC = self.AddendaHC[:,self.AddendaHC[0].argsort()]
		self.ThermCondWire = self.ThermCondWire[:,self.ThermCondWire[0].argsort()]

		SData = self.smoothedData #create convenient shorthand for smoothed data

		# Loop through all curves and compute heat capacity
		for i in range(lenSD):
			maxT = self.rawdata[i,1,-1,0]  #maximum temperature reached in heating pulse
			for k in range(2):
				for j in range(1,len(self.smoothedData[i,0,:,0])-2):
					
					Ts = self.smoothedData[i,1,j,k]
					Ph = self.rawdata[i,2,j,k]
					# Let thermal conductivity be the integral of conductivity w.r.t. temperature
					# plus a static offset parameter (as defined by eq. 4.7b in the PPMS manual)

					KwIntegral = self._WireCondIntegral(self.Tb[i], Ts)
					KwdT = KwIntegral  + self._WireCond(self.Tb[i])*(Ts- self.Tb[i])*StaticOffset
					
					# compute dT/dt using 2nd order central finite difference
					dTdt = (8.*(SData[i,1,j+1,k]- SData[i,1,j-1,k]) - (SData[i,1,j+2,k]- SData[i,1,j-2,k]) )/\
						(12.*(SData[i,0,j,k]-SData[i,0,j-1,k]))
					# compute dT/dt using 1st order central finite difference
					#dTdt = (SData[i,1,j+1,k]- SData[i,1,j-1,k])/\
					#	(2*(SData[i,0,j,k]-SData[i,0,j-1,k]))

					# compute heat capacity
					#################################################
					SData[i,2,j,k] = (-KwdT + Ph)/dTdt
					#################################################

					# subtract addenda
					SData[i,2,j,k] -= self._addenda(Ts)

					###############################
					# Compute uncertainty
					###############################
					deltaTb = 0.0001		
					deltaTs = 0.00003
					if k == 0 : deltaP = 0.1e-12	# typically negligible
					else: deltaP = 0
					# We approximate deltaKw from the spread of any short pulse data that exists.
					try: deltaKw = np.interp(Ts, self.avgThermCondWire[0], self.avgThermCondWire[2])
					except AttributeError: deltaKw = 4e-10   # <-- Dominant term
					deltaS = StaticOffset*0.1
					#######################################
					self.HC_uncertainty[i,j,k] = 1/dTdt * np.sqrt(
							((1-StaticOffset)*self._WireCond(self.Tb[i])*  deltaTb)**2  +\
							((1+StaticOffset)*(Ts - self.Tb[i])*  deltaKw)**2  +\
							deltaP**2  +\
							(self._WireCond(self.Tb[i])*(Ts - self.Tb[i])*deltaS)**2 +\
							((self._WireCond(Ts)+StaticOffset*self._WireCond(self.Tb[i]))**2 +\
									(65*SData[i,2,j,k]**2)/
									(72*(SData[i,0,j,k]-SData[i,0,j-1,k])**2)) * deltaTs**2
						)
					####################################################

					#Eliminate values too close to max T or min T
					if (((Ts - self.Tb[i])/self.Tb[i] < 0.1) or ((Ts - self.Tb[i]) < 0.025)):  
						SData[i,2,j,k] = np.nan

					if k == 1: #Eliminate values too close to max T on heating
						if (maxT  - Ts)/(maxT-self.Tb[i]) < 0.15:  
							SData[i,2,j,k] = np.nan
					elif k == 0:  #Eliminate values to close to min T on cooling
						if (maxT  - Ts)/(maxT-self.Tb[i]) < 0.16:  
							SData[i,2,j,k] = np.nan

				# first two and last two data points can't have value: no n-1 term
				SData[i,2,0,k] = np.nan  
				SData[i,2,1,k] = np.nan  
				SData[i,2,-1,k] = np.nan  
				SData[i,2,-2,k] = np.nan
				SData[i,2,-3,k] = np.nan
				# eliminate the first few points from a heating pulse (never reliable)
				if k ==0: SData[i,2,2:5,k] *= np.nan 

				# Apply this same elimination to uncertainty
				self.HC_uncertainty[i,np.argwhere(np.isnan(SData[i,2,:,k])),k] *= np.nan

				# Take moving average of data
				self.T[i,:,k] = SData[i,1,:,k]
				self.HC[i,:,k] =SData[i,2,:,k]
				#self.HC[i,:,k] = self._movingaverage( SData[i,2,:,k], smoothlevel)
			#Print progress
			sys.stdout.write('\r %d%% ' % (100*i/lenSD))
			sys.stdout.flush() # important

		# Convert to J/K/mol from J/K
		self.HC *= self.molarmass/(self.sampmass*1e-3)
		self.HC_uncertainty *= self.molarmass/(self.sampmass*1e-3) 

		print("\r 100%")

	def shortPulseAverage(self):
		"""Combine short pulse data points taken at the same temperatures. Threshold
		is defined in terms of percent change. Adjust this as necessary."""
		threshold = 0.01

		avgSpHc = [[]]
		avgSpT = [[]]
		avgTCW = []
		Bindex = 0
		# Initialize the first values:
		firstnonzero = next((i for i, x in enumerate(self.ShortPulse) if x), None)
		if firstnonzero is None:
			self.avgSpHc = []
			self.avgSpT = [[]]
			self.avgSpB = [[]]
			self.ScaledavgSpB = [[]]
			return 0
		Bprev = self.Bfield[firstnonzero]
		avgSpB = [Bprev]
		Tprev = round(self.Tsamp[firstnonzero],3)
		Tavg = [self.Tsamp[firstnonzero]]
		Cavg = [self.ShortPulse[firstnonzero]]
		KWavg = [self.ThermCondWire[1,firstnonzero]]
		for i, sp in enumerate(self.ShortPulse):
			if sp != 0:
				if self.Bfield[i] not in avgSpB:
					avgSpHc[Bindex].append(np.mean(Cavg))
					avgSpT[Bindex].append(np.mean(Tavg))
					avgTCW.append([np.mean(Tavg), np.mean(KWavg), np.std(KWavg), len(Tavg)])
					Tavg = [self.Tsamp[i]]
					Cavg = [sp]
					KWavg = [self.ThermCondWire[1,i]]
					# start new data set
					Tprev = np.mean(Tavg)
					avgSpB.append(self.Bfield[i])
					avgSpHc.append([])
					avgSpT.append([])
					Bindex = np.where(avgSpB == self.Bfield[i])[0][0]
				else:
					if np.abs((self.Tsamp[i] - Tprev)/Tprev) < threshold:
						Tavg.append(self.Tsamp[i])
						Cavg.append(sp)
						KWavg.append(self.ThermCondWire[1,i])
					else:
						avgSpHc[Bindex].append(np.mean(Cavg))
						avgSpT[Bindex].append(np.mean(Tavg))
						avgTCW.append([np.mean(Tavg), np.mean(KWavg), np.std(KWavg), len(Tavg)])
						Tavg = [self.Tsamp[i]]
						Cavg = [sp]
						KWavg = [self.ThermCondWire[1,i]]
						Tprev = np.mean(Tavg)
				#Bprev = self.Bfield[i]
				#Tprev = round(self.Tsamp[i],3)
		avgSpHc[Bindex].append(np.mean(Cavg))
		avgSpT[Bindex].append(np.mean(Tavg))
		avgTCW.append([np.mean(Tavg), np.mean(KWavg), np.std(KWavg), len(Tavg)])

		self.avgSpHc = np.array(avgSpHc)
		self.avgSpT = avgSpT
		self.avgSpB = avgSpB
		self.ScaledavgSpB = copy.deepcopy(self.avgSpB)
		self.avgThermCondWire = np.array(avgTCW)

		#sort the thermal conductivity
		Tarrinds = self.avgThermCondWire[:,0].argsort()
		self.avgThermCondWire = self.avgThermCondWire[Tarrinds].T
		#self.avgThermCondWire[1] = self._movingaverage(self.avgThermCondWire[1], 2)
		#self.avgThermCondWire[2] = self._movingaverage(self.avgThermCondWire[2], 2)


	def importCalibration(self,calfile):
		"""Imports thermal conductivity, heater resistivity,
		and thermometer resistance from calibration file"""
		lines = open(calfile,'r').readlines()
		for ll in range(len(lines)):
			line = lines[ll]

			if "[Temp_Cond]" in line:
				self.Kw = np.array(self._recordCalData(lines, ll)).T

			elif "[Temp_HtrRes]" in line:
				self.HtrRes = np.array(self._recordCalData(lines, ll)).T

			# Import thermometer resistance at various fields (rest of function)
			elif "[CalibrationFields]" in line:
				self.CalFields = {'f0':0.0}
				FieldNames = ['f0']  # there's always zero field
				numfields = int(lines[ll+1].split('=')[1])
				for i in range(numfields):
					cf = lines[ll+2+i].split('=')
					self.CalFields[cf[0]] = float(cf[1])
					FieldNames.append(cf[0])
				ThRes = {fn: [] for fn in FieldNames}

			elif "[Temp_ThRes" in line:
				if any([x in line for x in FieldNames[1:]]):
					CalField = [x for x in FieldNames if x in line][0]
					newdata = self._recordCalData(lines, ll)
					ThRes[CalField].extend(newdata) 

				# Import zero field data
				elif "f" not in line:
					newdata = self._recordCalData(lines, ll)
					ThRes['f0'].extend(newdata) 

		for f in FieldNames:  # convert to sorted numpy array
			ThRes[f] = np.array(ThRes[f])
			ThRes[f] = ThRes[f][ThRes[f][:,0].argsort()].T

		# Combine close repeated values in thermometer resistance
		self.AvgThRes = {fn: [[],[]] for fn in FieldNames}
		AvgThTemp = {fn: [[],[]] for fn in FieldNames}
		threshold = 0.005 #(K)
		for f in FieldNames:
			Tavg = []
			Ravg = []
			Tprev = ThRes[f][0,0]
			for i, Tem in enumerate(ThRes[f][0]):
				if (np.abs(Tem-Tprev) <= threshold) | (np.abs(Tem-Tprev)/Tem <= threshold):
					Tavg.append(Tem)
					Ravg.append(ThRes[f][1,i])
				else:
					self.AvgThRes[f][0].append(np.mean(Tavg))
					self.AvgThRes[f][1].append(np.mean(Ravg))
					Tavg = [Tem]
					Ravg = [ThRes[f][1,i]]
				Tprev = Tem

			# Add a few extra points at the end so the cubic spline doesn't get messed up
			# atr = self.AvgThRes[f]
			# beginslope = (atr[1][0] - atr[1][1])/\
			# 			(atr[0][0] - atr[0][1])
			# xbeg1 = atr[0][0] - (atr[0][1] - atr[0][0])
			# xbeg2 = (atr[0][0] + xbeg1)/2.
			# xbeg3 = (xbeg1 + xbeg2)/2.
			# xbeg4 = (xbeg1 + xbeg3)/2.
			# for x in [xbeg2, xbeg3, xbeg4, xbeg1]:
			# 	ynew = beginslope*(x-atr[0][0]) + atr[1][0]
			# 	self.AvgThRes[f][0].insert(0,x)
			# 	self.AvgThRes[f][1].insert(0,ynew)	

			# endslope = (atr[1][-1] - atr[1][-2])/\
			# 			(atr[0][-1] - atr[0][-2])
			# xend1 = atr[0][-1] + (atr[0][-1] - atr[0][-2])
			# xend2 = (atr[0][-1] + xend1)/2.
			# xend3 = (xend2 + xend1)/2.
			# xend4 = (xend2 + xend3)/2.
			# for x in [xend4, xend2, xend3, xend1]:
			# 	ynew = endslope*(x-atr[0][-1]) + atr[1][-1]
			# 	self.AvgThRes[f][0].append(x)
			# 	self.AvgThRes[f][1].append(ynew)

		# Prepare functions for interpolation
		self.ResTempFunc = {}
		for f in FieldNames:
			self.AvgThRes[f] = np.array(self.AvgThRes[f])   # Used for "_tempToResis"
			# Smooth data (unhelpful)
			#self.AvgThRes[f][0] = self._movingaverage( self.AvgThRes[f][0], 3)
			AvgThTemp[f] = self.AvgThRes[f][:,self.AvgThRes[f][1].argsort()]  # Used for "_resisToTemp"
			# self.ResTempFunc[f] = interp1d(AvgThTemp[f][1], AvgThTemp[f][0], kind='cubic',
			# 	bounds_error = False)
			self.ResTempFunc[f] = CubicSpline(AvgThTemp[f][1], AvgThTemp[f][0], 
				bc_type = 'natural',extrapolate=True)

		# plt.figure()
		# plt.plot(AvgThTemp[FieldNames[0]][1], 1/AvgThTemp[FieldNames[0]][0], '.')
		# xdat = np.linspace(AvgThTemp[FieldNames[0]][1][0],AvgThTemp[FieldNames[0]][1][-1], 1000)
		# plt.plot(xdat, 1/self.ResTempFunc[FieldNames[0]](xdat))

	def _recordCalData(self,datalines,start):
		"""Used with the importCalibration function"""
		count = int(datalines[start + 6].split('=')[1])
		data = []
		for i in range(start+7, start+7+count):
			data.append([float(d) for d in datalines[i].strip('\n').strip(',').split(',')])
		return data

	def _movingaverage(self, datay, n):
		"""Computes moving average of the y data set"""
		if not (n & 0x1):
			n+=1 	# Force the average number to be odd (to keep x values accurate)
		newdatay = np.convolve(datay, np.ones((n,))/n, mode='same')
		for i in range(int((n-1)/2)):
			newdatay[i] = np.average(datay[:(2*i)+1])
			newdatay[-i-1] = np.average(datay[-(2*i)-1:])
		return newdatay

	def _movingaverage_xspacing(self, datax, datay, smooth, xSpacing):
		"""Computes moving average of the y data set, but limited only to
		data points under xSpacing distance."""
		length = len(datay)
		newy  = np.zeros(length)
		for i in range(length):
			if (i-smooth) < 0:
				llim = 0
				ulim = 2*i
			elif (i+smooth) >= length:
				llim = 2*i - length
				ulim = length-1
			else:
				llim = i-smooth
				ulim = i+smooth
			indexrange = np.arange(llim, ulim+1)
			xspacerange = np.where((datax < datax[i]+xSpacing) & (datax > datax[i]-xSpacing))[0]

			indices = np.intersect1d(indexrange, xspacerange)
			#TODO: Force range to be symmetric about i

			lowind, highind = indices[0], indices[-1]
			if i-lowind < highind-i:
				highind = i + (i-lowind)
			elif i-lowind > highind-i:
				lowind = i - (highind-i)

			#print i, indexrange, xspacerange, lowind, highind
			#print np.array(datay)[indices]

			newy[i]= np.average(datay[lowind:highind+1])
		return newy





	def _WireCond(self,temp):
		"""Returns the wire conductivity for a given temperature
		(must run 'importCalibration' first)"""
		tcw= self.Kw
		return np.interp(temp, tcw[0], tcw[1])

		# #atcw = self.avgThermCondWire	
		# if np.logical_and(temp < tcw[0,-1], temp>tcw[0,0]):
		# 	return np.interp(temp, tcw[0], tcw[1])
		# #If temp is not in range of avgThermCondWire, extrapolate last two points.
		# elif temp < atcw[0,0]:
		# 	x1, y1, x2, y2 = tcw[0,0], tcw[1,0], tcw[0,1], tcw[1,1]
		# 	AA = (y2-y1)/(x2-x1)
		# 	BB = y1-AA*x1
		# 	return AA*temp + BB
		# elif temp > atcw[0,-1]:
		# 	x1, y1, x2, y2 = tcw[0,-2], tcw[1,-2], tcw[0,-1], tcw[1,-1]
		# 	AA = (y2-y1)/(x2-x1)
		# 	BB = y1-AA*x1
		# 	return AA*temp + BB
		# #return np.interp(temp, self.ThermCondWire[0], self.ThermCondWire[1])

	def _WireCondIntegral(self, temp1, temp2):
		"""Used for computing heat capacity"""
		tcw = self.Kw

		startK = self._WireCond(temp1)
		endK = self._WireCond(temp2)
		intermediate = np.where(np.logical_and(tcw[0]>temp1, tcw[0]<temp2))[0]
		if len(intermediate) == 0:
			Kint = 0.5*(endK + startK)*(temp2-temp1)
		else:
			interT = tcw[0][intermediate]
			interK = tcw[1][intermediate]
			# Integrate with the trapezoid rule
			Kint = 0.5*(startK + interK[0])*(interT[0]-temp1)
			for i in range(len(intermediate)-1):
				Kint += 0.5*(interK[i] + interK[i+1])*(interT[i+1]-interT[i])
			Kint += 0.5*(endK + interK[-1])*(temp2-interT[-1])
		return Kint


	def _HeaterRes(self,temp):
		"""Returns the heater resistance for a given temperature
		(must run 'importCalibration' first)"""
		return np.interp(temp, self.HtrRes[0], self.HtrRes[1])

	def _tempToResis(self, Temp, Field):
		"""Returns the thermometer resistance for a given temperature"""
		## Interpolate temperatures
		FieldArray = np.array(list(self.CalFields.values()))
		IntrpTRes = np.zeros(len(FieldArray))
		for i,f in enumerate(FieldArray):
			IntrpTRes[i] = np.interp(Temp, self.AvgThRes['f'+str(i)][0], self.AvgThRes['f'+str(i)][1])
		## Interpolate resistance between two closest fields
		return np.interp(Field, FieldArray, IntrpTRes)

	def _resisToTemp(self, Resis, Field):
		"""Used to compute temperature from thermometer resistance.
		For long pulses, the PPMS calculates temperature incorrectly."""
		## Interpolate temperatures
		FieldArray = np.array(list(self.CalFields.values()))
		IntrpTTemp = np.zeros(len(FieldArray))
		for i,f in enumerate(FieldArray):
			IntrpTTemp[i] = self.ResTempFunc['f'+str(i)](Resis)
		## Interpolate resistance between two closest fields
		return np.interp(Field, FieldArray, IntrpTTemp)

	def _addenda(self,temp):
		"""Returns the addenda Hc (in J/K) for a given temperature
		(also works for an array of temperatures)"""
		return np.interp(temp, self.AddendaHC[0], self.AddendaHC[1])


	def _combineTraces(self,smooth, FieldBinSize=10, useHeatPulses=False, onlyHeatPulses=False):
		"""Combines all heat capacity traces into a single line.
		Used for lineplotCombine and Entropy calculations.
		FieldBinSize in in Oe, default setting is to only use cooling pulses."""

		Barray = np.sort(list(set(self.Bfield)))
		self.CombB = Barray
		ScaledBarray = np.zeros(len(Barray))
		for ii in range(len(Barray)):
			B = Barray[ii]
			Bindex = np.where(self.Bfield==B)[0][0]
			ScaledBarray[ii] = self.ScaledBfield[Bindex]
		#ScaledBarray = np.sort(list(set(np.around(self.ScaledBfield,1))))
		self.ScaledCombB = ScaledBarray
		if len(ScaledBarray) != len(Barray):
			raise IndexError("ScaledBarray != Barray")


		HCtoComb = []
		TtoComb = []
		self.CombHC = []
		self.CombT = []

		#plt.figure()

		# Find where the short pulses are
		LongPindices = np.where(self.ShortPulse==0)[0]
		# loop through fields and combine traces
		for jj in range(len(Barray)):
			B = Barray[jj]

			# find all indices where the magnetic field is of the specified value
			Bindices = np.where(np.abs(self.Bfield - B)<FieldBinSize)[0]
			# Take intersection of long pulses and field to eliminate short pulse data
			Bindices = np.intersect1d(Bindices, LongPindices)
			if len(Bindices) == 0:		
				continue	# Skip if no data exists for this field

			combinedHc = []
			combinedT = []
			# Average the points which overlap with another curve
			for bi in Bindices:
				if onlyHeatPulses: # Combine all heating pulses into a single trace
					nonnan = np.where(~np.isnan(np.array(self.HC[bi][:,0])))
					overlapdataHC = self.HC[bi][nonnan,0].flatten()
					for bj in Bindices[np.where(Bindices != bi)]:
						overlapdataHC = np.vstack((overlapdataHC,
							np.interp(self.T[bi][nonnan,0].flatten(), 
							self.T[bj][:,0], self.HC[bj][:,0],
							right = np.nan, left = np.nan)))
					combinedHc.extend(np.nanmean(np.array(overlapdataHC), axis=0))
					combinedT.extend(self.T[bi][nonnan,0].flatten())

				else:
					#Combine all cooling pulses into a single trace
					nonnan = np.where(~np.isnan(np.array(self.HC[bi][:,1])))
					overlapdataHC = self.HC[bi][nonnan,1].flatten()
					for bj in Bindices[np.where(Bindices != bi)]:
						overlapdataHC = np.vstack((overlapdataHC,
							np.interp(self.T[bi][nonnan,1].flatten(), 
							self.T[bj][:,1][::-1], self.HC[bj][:,1][::-1],
							right = np.nan, left = np.nan)))

						if useHeatPulses: # include heating pulses
							overlapdataHC = np.vstack((overlapdataHC,
								np.interp(self.T[bi][nonnan,1].flatten(), 
								self.T[bj][:,0], self.HC[bj][:,0],
								right = np.nan, left = np.nan)))

					if useHeatPulses:
						# Compute overlap with same index heating pulse
						overlapdataHC = np.vstack((overlapdataHC,
							np.interp(self.T[bi][nonnan,1].flatten(), 
							self.T[bi][:,0], self.HC[bi][:,0],
							right = np.nan, left = np.nan)))

					# Concatenate data to array, to be sorted later
					combinedHc.extend(np.nanmean(np.array(overlapdataHC), axis=0))
					combinedT.extend(self.T[bi][nonnan,1].flatten())

					if useHeatPulses:
						# Compute overlap with same index heating pulse
						overlapdataHC = np.vstack((overlapdataHC,
							np.interp(self.T[bi][nonnan,1].flatten(), 
							self.T[bi][:,0], self.HC[bi][:,0],
							right = np.nan, left = np.nan)))
						# Now repeat, but computing overlap with heating pulse
						nonnan = np.where(~np.isnan(np.array(self.HC[bi][:,0])))
						overlapdataHC = self.HC[bi][nonnan,0].flatten()

						for bj in Bindices[np.where(Bindices != bi)]:
							overlapdataHC = np.vstack((overlapdataHC,
								np.interp(self.T[bi][nonnan,0].flatten(), 
								self.T[bj][:,1][::-1], self.HC[bj][:,1][::-1],
								right = np.nan, left = np.nan)))
							overlapdataHC = np.vstack((overlapdataHC,
								np.interp(self.T[bi][nonnan,0].flatten(), 
								self.T[bj][:,0], self.HC[bj][:,0],
								right = np.nan, left = np.nan)))

						# Concatenate data to array, to be sorted later
						combinedHc.extend(np.nanmean(np.array(overlapdataHC), axis=0))
						combinedT.extend(self.T[bi][nonnan,0].flatten())

			combinedHc = np.array(combinedHc)
			combinedT = np.array(combinedT)

			# Sort data by temperature
			Tarrinds = combinedT.argsort()
			combinedHc = combinedHc[Tarrinds]
			combinedT = combinedT[Tarrinds]

			AvgCombHc= []
			AvgCombT = []
			# Combine close repeated values
			threshold = 0.002 #(K)
			Tprev = combinedT[0]
			HCprev = combinedHc[0]
			i=1
			while i < len(combinedT):
				Tem = combinedT[i]
				HCap = combinedHc[i]
				if np.abs(Tem-Tprev) <= threshold:
					AvgCombT.append(0.5*(Tem + Tprev))
					AvgCombHc.append(0.5*(HCap + HCprev))
					HCprev = combinedHc[i]
					Tprev = combinedT[i]
					i+=1
				else:
					AvgCombT.append(Tprev)
					AvgCombHc.append(HCprev)
					HCprev = combinedHc[i]
					Tprev = combinedT[i]
				i+=1

			# Smooth data
			AvgCombHc = self._movingaverage_xspacing(AvgCombT,AvgCombHc, smooth, 0.01)
			# Append to object list
			self.CombHC.append(AvgCombHc)
			self.CombT.append(AvgCombT)


	def _combineTracesOld(self,smooth, FieldBinSize=10):
		"""Combines all heat capacity traces into a single line.
		Used for lineplotCombine and Entropy calculations.
		FieldBinSize in in Oe"""

		Barray = np.sort(list(set(self.Bfield)))
		self.CombB = Barray
		ScaledBarray = np.zeros(len(Barray))
		for ii in range(len(Barray)):
			B = Barray[ii]
			Bindex = np.where(self.Bfield==B)[0][0]
			ScaledBarray[ii] = self.ScaledBfield[Bindex]
		#ScaledBarray = np.sort(list(set(np.around(self.ScaledBfield,1))))
		self.ScaledCombB = ScaledBarray
		if len(ScaledBarray) != len(Barray):
			raise IndexError("ScaledBarray != Barray")

		self.CombHC = []
		self.CombT = []

		# Find where the short pulses are
		LongPindices = np.where(self.ShortPulse==0)[0]
		# loop through fields and combine traces
		for jj in range(len(Barray)):
			B = Barray[jj]

			# find all indices where the magnetic field is of the specified value
			Bindices = np.where(np.abs(self.Bfield - B)<FieldBinSize)[0]
			# Take intersection of long pulses and field to eliminate short pulse data
			Bindices = np.intersect1d(Bindices, LongPindices)
			if len(Bindices) == 0:		
				continue	# Skip if no data exists for this field

			# Concatenate data
			combinedHc = self.HC[Bindices[0]][:,1]
			combinedT = self.T[Bindices[0]][:,1]
			for bb in Bindices[1:]:
				combinedHc = np.hstack((combinedHc,self.HC[bb][:,1]))
				combinedT = np.hstack((combinedT,self.T[bb][:,1]))

			# Eliminate nan values
			nonnan = np.where(~np.isnan(combinedHc))
			combinedHc = combinedHc[nonnan]
			combinedT = combinedT[nonnan]

			# Sort data by temperature
			Tarrinds = combinedT.argsort()
			combinedHc = combinedHc[Tarrinds]
			combinedT = combinedT[Tarrinds]

			AvgCombHc= []
			AvgCombT = []
			# Combine close repeated values
			threshold = 0.005 #(K)
			Tprev = combinedT[0]
			HCprev = combinedHc[0]
			i=1
			while i < len(combinedT):
				Tem = combinedT[i]
				HCap = combinedHc[i]
				if np.abs(Tem-Tprev) <= threshold:
					AvgCombT.append(0.5*(Tem + Tprev))
					AvgCombHc.append(0.5*(HCap + HCprev))
					HCprev = combinedHc[i]
					Tprev = combinedT[i]
					i+=1
				else:
					AvgCombT.append(Tprev)
					AvgCombHc.append(HCprev)
					HCprev = combinedHc[i]
					Tprev = combinedT[i]
				i+=1

			# Smooth data
			AvgCombHc = self._movingaverage_xspacing(AvgCombT,AvgCombHc, smooth, 0.01)
			# Append to object list
			self.CombHC.append(AvgCombHc)
			self.CombT.append(AvgCombT)






	# Functions called by user:

	def plotHC(self,axes,index,heatingcolor=None,coolingcolor=None,shortpulsecolor=None,
			Blabels=False, demag=True, marker = 'o', PlotUncertainty=False, **kwargs):
		"""If the pulse is a long pulse, plot it. If it is a short pulse,
		plot a single point. If you don't want to plot a curve, leave
		the color value as blank."""
		ShortHC = self.ShortPulse[index]
		if ShortHC == 0:
			if heatingcolor is not None:
				if 'color' in kwargs:
					axes.plot(self.T[index][:,0], self.HC[index][:,0],**kwargs)
				else:
					axes.plot(self.T[index][:,0], self.HC[index][:,0],color=heatingcolor, **kwargs)
				if PlotUncertainty == True:
					axes.fill_between(self.T[index][:,0], 
						self.HC[index][:,0]+ self.HC_uncertainty[index][:,0],
						self.HC[index][:,0]- self.HC_uncertainty[index][:,0],
						facecolor=[(2*i+0.5)/3. for i in cc.to_rgba(heatingcolor, alpha=0.5)], 
						edgecolor='none', interpolate=True)
			if coolingcolor is not None:
				if 'color' in kwargs:
					axes.plot(self.T[index][:,1], self.HC[index][:,1], **kwargs)
				else:
					axes.plot(self.T[index][:,1], self.HC[index][:,1],color=coolingcolor, **kwargs)
				if PlotUncertainty == True:
					axes.fill_between(self.T[index][:,1], 
						self.HC[index][:,1]+ self.HC_uncertainty[index][:,1],
						self.HC[index][:,1]- self.HC_uncertainty[index][:,1],
						facecolor=[(2*i+0.5)/3. for i in cc.to_rgba(coolingcolor, alpha=0.5)], 
						edgecolor='none', alpha=0.3, interpolate=True)
		else: 
			if shortpulsecolor is not None:
				axes.plot(self.Tsamp[index],self.ShortPulse[index], marker=marker,
					markeredgecolor=shortpulsecolor, markerfacecolor='none')

		# Create labels for a data legend
		if Blabels == False:
			heating_curve = mlines.Line2D([], [], color=heatingcolor, label='Heating Curve', **kwargs)
			cooling_curve = mlines.Line2D([], [], color=coolingcolor, label='Cooling Curve', **kwargs)
			short_pulses = mlines.Line2D([], [], color=shortpulsecolor, marker=marker,
						markeredgecolor=shortpulsecolor, markerfacecolor='none', 
						label='Adiabatic Pulses', linestyle="None")

			self.labels = []
			if heatingcolor is not None:
				self.labels.append(heating_curve)
			if coolingcolor is not None:
				self.labels.append(cooling_curve)
			if (shortpulsecolor is not None) and (ShortHC != 0):
				self.labels.append(short_pulses)

		if Blabels == True:
			if demag == True:
				Bvalue = self.ScaledBfield[index]
			else:
				Bvalue = self.Bfield[index]
			if len(self.labels) == 0: self.Bflag = [] # This means labels has been reset
			if Bvalue not in self.Bflag:
				# Right now, it uses the cooling curve color for the label.

				labl = str(abs(Bvalue)/10000)+' T'
				#Blab = mlines.Line2D([], [], color=coolingcolor, label=labl)

				self.Bflag.append(Bvalue)

				if (coolingcolor is not None) and (heatingcolor is not None):
					self.labels.append(mlines.Line2D([], [], color=coolingcolor, 
						label=labl + ' (cooling)'))
					self.labels.append(mlines.Line2D([], [], color=heatingcolor, 
						label=labl+ ' (heating)'))
					self.Bflag.append(Bvalue)
				elif coolingcolor is not None:
					self.labels.append(mlines.Line2D([], [], color=coolingcolor, label=labl))
				elif heatingcolor is not None:
					self.labels.append(mlines.Line2D([], [], color=heatingcolor, label=labl))
				else:
					self.labels.append(mlines.Line2D([], [], color=shortpulsecolor, label=labl))


				#Sort Bflags, and make into array of labels
				self.labels = [x for (y,x) in sorted(zip(self.Bflag,self.labels), 
					key=lambda pair: pair[0])]
				self.Bflag = sorted(self.Bflag)

			if (Bvalue not in self.shortpulseBflag and ShortHC != 0):
				labl = str(Bvalue/10000)+' T'
				self.shortpulseBflag.append(Bvalue)

				if shortpulsecolor is not None:
					self.shortpulselabels.append(mlines.Line2D([], [], color=shortpulsecolor, 
						marker='o', linestyle="None", markeredgecolor=shortpulsecolor, 
						markerfacecolor='none', label=labl))

				self.shortpulselabels = [x for (y,x) in sorted(zip(self.shortpulseBflag,
					self.shortpulselabels), 
					key=lambda pair: pair[0])]
				self.shortpulseBflag = sorted(self.shortpulseBflag)



	def plotHCT(self,axes,index,heatingcolor=None,coolingcolor=None,shortpulsecolor=None,
			Blabels=False, demag=True, marker = 'o', PlotUncertainty=False, **kwargs):
		"""If the pulse is a long pulse, plot it. If it is a short pulse,
		plot a single point. If you don't want to plot a curve, leave
		the color value as blank."""
		ShortHC = self.ShortPulse[index]
		if ShortHC == 0:
			if heatingcolor is not None:
				if 'color' in kwargs:
					axes.plot(self.T[index][:,0], self.HC[index][:,0]/self.T[index][:,0], **kwargs)
				else:
					axes.plot(self.T[index][:,0], self.HC[index][:,0]/self.T[index][:,0],color=heatingcolor, **kwargs)
				if PlotUncertainty == True:
					axes.fill_between(self.T[index][:,0], 
						(self.HC[index][:,0]+ self.HC_uncertainty[index][:,0])/self.T[index][:,0],
						(self.HC[index][:,0]- self.HC_uncertainty[index][:,0])/self.T[index][:,0],
						facecolor=[(2*i+0.5)/3. for i in cc.to_rgba(heatingcolor, alpha=0.5)], 
						edgecolor='none', interpolate=True)
			if coolingcolor is not None:
				if 'color' in kwargs:
					axes.plot(self.T[index][:,1], self.HC[index][:,1]/self.T[index][:,1], **kwargs)
				else:
					axes.plot(self.T[index][:,1], self.HC[index][:,1]/self.T[index][:,1],color=coolingcolor, **kwargs)
				if PlotUncertainty == True:
					axes.fill_between(self.T[index][:,1], 
						(self.HC[index][:,1]+ self.HC_uncertainty[index][:,1])/self.T[index][:,1],
						(self.HC[index][:,1]- self.HC_uncertainty[index][:,1])/self.T[index][:,1],
						facecolor=[(2*i+0.5)/3. for i in cc.to_rgba(coolingcolor, alpha=0.5)], 
						edgecolor='none', alpha=0.3, interpolate=True)
		else: 
			if shortpulsecolor is not None:
				axes.plot(self.Tsamp[index],self.ShortPulse[index]/self.Tsamp[index], marker=marker,
					markeredgecolor=shortpulsecolor, markerfacecolor='none')

		# Create labels for a data legend
		if Blabels == False:
			heating_curve = mlines.Line2D([], [], color=heatingcolor, label='Heating Curve', **kwargs)
			cooling_curve = mlines.Line2D([], [], color=coolingcolor, label='Cooling Curve', **kwargs)
			short_pulses = mlines.Line2D([], [], color=shortpulsecolor, marker=marker,
						markeredgecolor=shortpulsecolor, markerfacecolor='none', 
						label='Adiabatic Pulses', linestyle="None")

			self.labels = []
			if heatingcolor is not None:
				self.labels.append(heating_curve)
			if coolingcolor is not None:
				self.labels.append(cooling_curve)
			if (shortpulsecolor is not None) and (ShortHC != 0):
				self.labels.append(short_pulses)

		if Blabels == True:
			if demag == True:
				Bvalue = self.ScaledBfield[index]
			else:
				Bvalue = self.Bfield[index]
			if len(self.labels) == 0: self.Bflag = [] # This means labels has been reset
			if Bvalue not in self.Bflag:
				# Right now, it uses the cooling curve color for the label.

				labl = str(abs(Bvalue)/10000)+' T'
				#Blab = mlines.Line2D([], [], color=coolingcolor, label=labl)

				self.Bflag.append(Bvalue)

				if (coolingcolor is not None) and (heatingcolor is not None):
					self.labels.append(mlines.Line2D([], [], color=coolingcolor, 
						label=labl + ' (cooling)'))
					self.labels.append(mlines.Line2D([], [], color=heatingcolor, 
						label=labl+ ' (heating)'))
					self.Bflag.append(Bvalue)
				elif coolingcolor is not None:
					self.labels.append(mlines.Line2D([], [], color=coolingcolor, label=labl))
				elif heatingcolor is not None:
					self.labels.append(mlines.Line2D([], [], color=heatingcolor, label=labl))
				else:
					self.labels.append(mlines.Line2D([], [], color=shortpulsecolor, label=labl))


				#Sort Bflags, and make into array of labels
				self.labels = [x for (y,x) in sorted(zip(self.Bflag,self.labels), 
					key=lambda pair: pair[0])]
				self.Bflag = sorted(self.Bflag)

			if (Bvalue not in self.shortpulseBflag and ShortHC != 0):
				labl = str(Bvalue/10000)+' T'
				self.shortpulseBflag.append(Bvalue)

				if shortpulsecolor is not None:
					self.shortpulselabels.append(mlines.Line2D([], [], color=shortpulsecolor, 
						marker='o', linestyle="None", markeredgecolor=shortpulsecolor, 
						markerfacecolor='none', label=labl))

				self.shortpulselabels = [x for (y,x) in sorted(zip(self.shortpulseBflag,
					self.shortpulselabels), 
					key=lambda pair: pair[0])]
				self.shortpulseBflag = sorted(self.shortpulseBflag)



	def lineplot(self,axes,Barray, plotHeatPulses=False, markers = ['s','^','o','x'], **kwargs):
		""" Plots all the traces of long-pulse heat capacity and points of short-pulse
		Currently uses 'gist_rainbow' colormap. If you don't like it, change it."""
		if (Barray == 'All' or Barray == 'all'):
			Barray = np.sort(list(set(self.Bfield)))
		else:
			Barray = np.array(Barray)

		# determine the color map
		#colormap = plt.cm.hsv((Barray*1.0)/Barray[-1]) * 0.6   #based on field

		colormap = plt.cm.hsv(np.arange(len(Barray))*1.0/len(Barray))* 0.75   #based on index
		colors = dict(zip(Barray, colormap))

		for jj in range(len(self.Bfield)):
			B = self.Bfield[jj]
			for b in Barray:
				if  B == b:
					self.plotHC(axes=axes,index=jj,coolingcolor=colors[B],Blabels=True, **kwargs)
					if plotHeatPulses == True:
						heatpulsecolor = colors[B]*0.6
						heatpulsecolor[-1] = 0.9
						self.plotHC(axes=axes,index=jj,heatingcolor=heatpulsecolor,Blabels=True, **kwargs)
		
		# Plot short pulse data 
		if np.count_nonzero(self.ShortPulse) == 0: return
		for jj, b in enumerate(Barray):
			for i in range(len(self.avgSpB)):
				spB = self.avgSpB[i]
				if spB == b:
					axes.plot(self.avgSpT[i], self.avgSpHc[i],color=colors[spB], 
						marker=markers[i%len(markers)],
						markeredgecolor=colors[spB], markerfacecolor='none', 
						label='Adiabatic Pulses', linestyle="None")
					try: 
						if demag == True:
							labl = str(abs(self.ScaledavgSpB[i])/10000.)+' T'
						else: labl = str(round(self.avgSpB[i],1)/10000.)+' T'
					except NameError: labl = str(round(self.avgSpB[i],1)/10000)+' T'
					self.shortpulselabels.append(mlines.Line2D([], [], color=colors[spB], 
							marker=markers[i%len(markers)], linestyle="None", markeredgecolor=colors[spB], 
							markerfacecolor='none', label=labl))


	def lineplotCombine(self,axes,Barray,smooth, demag=True, plotShortPulse=True, 
		markers = ['s','^','o','x'], FieldBinSize = 10, onlyHeatPulses=False,
		useHeatPulses=False, **kwargs):
		"""Combines all the heat capacity traces in a given field so that 
		there is only one line plotted"""
		self.labels = []

		if Barray in ['All', 'all']:
			Barray = np.sort(list(set(self.Bfield)))
		else:
			Barray = np.array(Barray)

		# determine the color map
		colormap = plt.cm.hsv(np.arange(len(Barray))*1.0/len(Barray))* 0.75   #based on index
		colors = dict(zip(Barray, colormap))

		# Combine traces into single line
		self._combineTraces(smooth, FieldBinSize, useHeatPulses, onlyHeatPulses)

		# plot the long pulse data
		for jj in range(len(Barray)):
			B = Barray[jj]
			if B == 0: B=0.0   # get rid of negative sign which may be there
			# find all indices where the magnetic field is of the specified value
			try:
				Bindex = np.where(self.CombB==B)[0][0]
			except IndexError:
				continue
			# Plot
			if demag == True:
				fieldval = round(self.ScaledCombB[Bindex]/10000 , 3)
				if fieldval == 0: fieldval = 0.0
				labl = str(fieldval)+' T'
			else:
				labl = str(Barray[jj]/10000.)+' T'

			if ('color' in kwargs) or ('c' in kwargs) :
				axes.plot(self.CombT[Bindex], self.CombHC[Bindex], **kwargs)
			else:
				axes.plot(self.CombT[Bindex], self.CombHC[Bindex], color=colors[B], **kwargs)
			self.labels.append(mlines.Line2D([], [], color=colors[B], label=labl))

		#Plot short pulse data
		if np.count_nonzero(self.ShortPulse) == 0: return
		if plotShortPulse == True:
			edgewidth = 0.8
			for jj, b in enumerate(Barray):
				for i in range(len(self.avgSpB)):
					spB = self.avgSpB[i]
					if spB == b:
						axes.plot(self.avgSpT[i], self.avgSpHc[i],color=colors[spB], 
							marker=markers[i%len(markers)], markeredgewidth = edgewidth,
							markeredgecolor=colors[spB], markerfacecolor='none', 
							label='Adiabatic Pulses', linestyle="None")
						if demag == True:
							labl = str(round(abs(self.ScaledavgSpB[i])/10000 , 3))+' T'
						else: labl = str(round(self.avgSpB[i],1)/10000)+' T'
						self.shortpulselabels.append(mlines.Line2D([], [], color=colors[spB], 
								marker=markers[i%len(markers)], linestyle="None", 
								markeredgewidth = edgewidth, markeredgecolor=colors[spB], 
								markerfacecolor='none', label=labl))


	def _computeEntropy(self,smooth):
		'''computes entropy for all field values'''
		Barray = np.sort(list(set(self.Bfield)))

		#Combine traces into single line (if not done already)
		try:
			self.CombHC
		except AttributeError:
			print(" combining traces...")
			self._combineTraces(smooth)

		self.Entropy = np.zeros_like(self.CombT)
		# loop through magnetic fields
		for jj in range(len(Barray)):
			B = Barray[jj]
			if B == 0: B=0.0   # get rid of negative sign which may be there
			# find all indices where the magnetic field is of the specified value
			try:
				Bindex = np.where(self.CombB==B)[0][0]
			except IndexError:		
				continue	# Skip if no data exists for this field

			# Compute Entropy

			T = self.CombT[Bindex]
			C = self.CombHC[Bindex]

			entropy = np.zeros(C.size)
			en = (T[0])*0.5*(C[0]/T[0])
			entropy[0] = en
			
			for i in range(0,C.size-1):
				ds = (T[i+1]-T[i])*0.5*(C[i]/T[i]+C[i+1]/T[i+1])
				en = en + ds
				entropy[i+1] = en

			self.Entropy[Bindex]=entropy

	def plotEntropy(self, axes,Barray,smooth):
		"""Plots entropy vs. T for various magnetic fields"""
		self.entropylabels = []
		if (Barray == 'All' or Barray == 'all'):
			Barray = np.sort(list(set(self.Bfield)))
		else:
			Barray = np.array(Barray)

		# determine the color map
		colormap = plt.cm.hsv(np.arange(len(Barray))*1.0/len(Barray))* 0.85   #based on index
		colors = dict(zip(Barray, colormap))

		#Compute entropy (if not done already)
		try:
			self.Entropy
		except AttributeError:
			print(" computing entropy...")
			self._computeEntropy(smooth)

		# loop through magnetic fields
		for jj in range(len(Barray)):
			B = Barray[jj]
			if B == 0: B=0.0   # get rid of negative sign which may be there
			# find all indices where the magnetic field is of the specified value
			try:
				Bindex = np.where(self.CombB==B)[0][0]
			except IndexError:		
				continue	# Skip if no data exists for this field

			# Plot
			fieldval = round(self.ScaledCombB[Bindex]/10000 , 3)
			if fieldval == 0: fieldval = 0.0
			labl = str(fieldval)+' T'
			axes.plot(self.CombT[Bindex], self.Entropy[Bindex], color=colors[B])
			self.entropylabels.append(mlines.Line2D([], [], color=colors[B], label=labl))
			


	def meshgrid(self,Tarray,Barray, useHeatPulses=False):
		"""Tarray is the array of x values to be binned to
		Barray is the array of magnetic fields to loop through
		Set Barray to 'all' if you want to plot all the fields."""

		if (Barray == 'All' or Barray == 'all'):
			Barray = np.sort(list(set(self.Bfield)))
			# Create array of scaledB:
			ScaledBarray = np.zeros(len(Barray))
			for ii in range(len(Barray)):
				B = Barray[ii]
				Bindex = np.where(self.Bfield==B)[0][0]
				ScaledBarray[ii] = np.around(self.ScaledBfield[Bindex],1)


		Intensity = np.empty((len(Tarray)-1,len(Barray)))*np.nan
		for ii in range(len(Barray)):
			B = Barray[ii]
			# find all indices where the magnetic field is of the specified value
			Bindices = np.where(self.Bfield==B)[0]
			if len(Bindices) == 0:		
				continue	# Skip if no data exists for this field
			# Concatenate data into single array
			binnedhc = self.binAndInterpolate(Tarray, self.T[Bindices[0]][:,1], 
					self.HC[Bindices[0]][:,1])

			for bb in Bindices:
				if self.ShortPulse[bb] != 0:
					continue	# skip this value if it is a short pulse
				newbinnd = self.binAndInterpolate(Tarray, self.T[bb][:,1] , self.HC[bb][:,1])
				binnedhc = np.vstack((binnedhc,newbinnd))

				if useHeatPulses:
					newbinnd = self.binAndInterpolate(Tarray, self.T[bb][:,0] , self.HC[bb][:,0])
					binnedhc = np.vstack((binnedhc,newbinnd))

			# Take average along given axis
			AvgBinnedHc = np.nanmean(binnedhc, axis=0)

			# Bin into x array
			Intensity[:,ii] = AvgBinnedHc

		#define edges of magnetic field bins
		d = np.diff(Barray)/2.
		Bedges = np.hstack([Barray[0]-d[0],Barray[0:-1]+d,Barray[-1]+d[-1]])

		d = np.diff(ScaledBarray)/2.
		Bedges = np.hstack([ScaledBarray[0]-d[0],ScaledBarray[0:-1]+d,ScaledBarray[-1]+d[-1]])

		# Mask all the zero elements
		Intensity = np.ma.masked_where(np.isnan(Intensity), Intensity)

		return Intensity.T, Bedges   #transpose


	def binAndInterpolate(self, x_new, x_old, y_old):
		###BIN
		edges = np.array(x_new)
		centers = edges[:-1] + (edges[1:]-edges[:-1])/2.

		# check that values in x_old fall into the range of edges
		# discard any that are outside the range
		if x_old[0] > x_old[-1]:  # reverse order so that low temperature comes first
			x_old = x_old[::-1]
			y_old = y_old[::-1]
		if x_old[-1] > edges[-1]:  # check if x_old goes higher than edges
			high_idx = np.nonzero(x_old > edges[-1])[0][0]
			x_old = x_old[:high_idx]
			y_old = y_old[0:high_idx]
			if len(x_old) == 0:    # return a masked array if no data falls in range
				y_new = np.zeros(np.size(x_new)-1)
				return np.ma.masked_where(y_new==0 , y_new)
		if x_old[0] < edges[0]:
			low_idx = np.nonzero(x_old <= edges[0])[0][-1]
			x_old = x_old[low_idx:]
			y_old = y_old[low_idx:]
			if len(x_old) == 0:
				y_new = np.zeros(np.size(x_new)-1)
				return np.ma.masked_where(y_new==0 , y_new)

		bin_idx = np.digitize(x_old,edges)
		bin_count, b = np.histogram(x_old, edges)

		y_new = np.zeros(np.size(x_new)-1)

		mask_idx = bin_count < 1. 
		mask_ct = np.ma.array(bin_count,mask = mask_idx)

		for ii, idx in enumerate(bin_idx):
			#print ii, idx, x_new[idx-1], x_old[ii]
			y_new[idx-1] += y_old[ii]

		# mask zeros then divide by bin_counts
		y_new = np.ma.masked_where(y_new==0 , y_new)/mask_ct

		###INTERPOLATE
		# remove masked elements of array
		unmaskedindices = np.where(np.logical_not(y_new.mask))[0]
		y_unmasked = np.ma.compressed(y_new)
		x_unmasked = centers[unmaskedindices]

		if len(y_unmasked) ==0:
			return np.zeros(len(centers))*np.nan
		#print y_unmasked, x_unmasked

		#interpolate
		y_interp = np.interp(centers, x_unmasked, y_unmasked, left =np.nan, right =np.nan)
		return y_interp


	def scale(self, factor):
		"""scale the heat capacity data by a constant.
		Useful for going between per FU and per ion"""
		self.HC = self.HC/factor
		self.HC_uncertainty = self.HC_uncertainty/factor
		self.ShortPulse = self.ShortPulse/factor
		for i, lst in enumerate(self.avgSpHc):
			for j, val in enumerate(lst):
				self.avgSpHc[i][j] = val/factor

	def scaleDemagFactor(self, demagFac, magvsfield):
		"""Scale the magnetic field values by a demagnetization factor.
		magvsfield should be a numpy array giving M as a function of INTERNAL field. 
		Note that this correction must be done numerically:
			H_int = H_0 - mu0*(D)*M(H_int)
		We use the bisection method here."""

		mu0 = 4*np.pi*1e-7
		def H0DM(field):
			"""Computes H + mu0*(D)*M(H)"""
			return field + mu0*demagFac*np.interp(field, magvsfield[0], magvsfield[1])*10000 #mu_B / F.U.

		print("Scaling demagnetization...")
		self.ScaledBfield = np.zeros(len(self.Bfield))
		self.ScaledavgSpB = np.zeros(len(self.avgSpB))

		for i, H0 in enumerate(self.Bfield):
			if H0 == 0.0:
				self.ScaledBfield[i] = 0.0
				continue
			# First, obtain two H values that bound the problem.
			H1 = H0*1.0
			mag1 = H0DM(H0) 
			if mag1 < H0:
				while mag1 < H0:
					H1 += 0.1*H0
					mag1 = H0DM(H1)
				H2 = H1*1.0
				H1 = H2 - 0.1*H0
			elif mag1 > H0:
				while mag1 > H0:
					H1 -= 0.1*H0
					mag1 = H0DM(H1)   #mu_B / F.U.
				H2 = H1 + 0.1*H0
			else:
				H2 = H0
			#print H0DM(H1)-H0, H0DM(H2)-H0

			# Next, use the bisection method to determine the value
			Hint = H1
			Hintold = H0*1.0
			while np.abs((Hint - Hintold)/H0) > 0.00001:
				Hintold = Hint*1.0
				Hint = 0.5*(H1+H2)
				mag1 = H0DM(H1)
				mag2 = H0DM(H2)
				magbi = H0DM(Hint)
				if magbi < H0:
					H1 = Hint
				else:
					H2 = Hint
			self.ScaledBfield[i] = Hint

		# Do the same thing, but for the short pulse field values
		for i, H0 in enumerate(self.avgSpB):
			if H0 == 0.0:
				self.ScaledBfield[i] = 0.0
				continue
			# First, obtain two H values that bound the problem.
			H1 = H0*1.0
			mag1 = H0DM(H0) 
			if mag1 < H0:
				while mag1 < H0:
					H1 += 0.1*H0
					mag1 = H0DM(H1)
				H2 = H1*1.0
				H1 = H2 - 0.1*H0
			elif mag1 > H0:
				while mag1 > H0:
					H1 -= 0.1*H0
					mag1 = H0DM(H1)   #mu_B / F.U.
				H2 = H1 + 0.1*H0
			else:
				H2 = H0
			#print H0DM(H1)-H0, H0DM(H2)-H0

			# Next, use the bisection method to determine the value
			Hint = H1
			Hintold = H0*1.0
			while np.abs((Hint - Hintold)/H0) > 0.00001:
				Hintold = Hint*1.0
				Hint = 0.5*(H1+H2)
				mag1 = H0DM(H1)
				mag2 = H0DM(H2)
				magbi = H0DM(Hint)
				if magbi < H0:
					H1 = Hint
				else:
					H2 = Hint
			self.ScaledavgSpB[i] = Hint

		try:
			# Do the same thing, but for the combined traces (if they exist)
			for i, H0 in enumerate(self.CombB):
				if H0 == 0.0:
					self.ScaledCombB[i] = 0.0
					continue
				# First, obtain two H values that bound the problem.
				H1 = H0*1.0
				mag1 = H0DM(H0) 
				if mag1 < H0:
					while mag1 < H0:
						H1 += 0.1*H0
						mag1 = H0DM(H1)
					H2 = H1*1.0
					H1 = H2 - 0.1*H0
				elif mag1 > H0:
					while mag1 > H0:
						H1 -= 0.1*H0
						mag1 = H0DM(H1)   #mu_B / F.U.
					H2 = H1 + 0.1*H0
				else:
					H2 = H0
				#print H0DM(H1)-H0, H0DM(H2)-H0

				# Next, use the bisection method to determine the value
				Hint = H1
				Hintold = H0*1.0
				while np.abs((Hint - Hintold)/H0) > 0.00001:
					Hintold = Hint*1.0
					Hint = 0.5*(H1+H2)
					mag1 = H0DM(H1)
					mag2 = H0DM(H2)
					magbi = H0DM(Hint)
					if magbi < H0:
						H1 = Hint
					else:
						H2 = Hint
				self.ScaledCombB[i] = Hint
		except AttributeError:
			pass


	def savetrace(self, index, outfile):
		heatarray = np.vstack((self.T[index][:,0],self.HC[index][:,0])).T
		np.savetxt(outfile+'_heating.txt',heatarray, fmt='%.9f', header = 'Temp(K)\tHC(J/[K mol-ion])', 
			delimiter=', ')
		coolarray = np.vstack((self.T[index][:,1],self.HC[index][:,1])).T
		np.savetxt(outfile+'_cooling.txt',coolarray, fmt='%.9f', header = 'Temp(K)\tHC(J/[K mol-ion])', 
			delimiter=', ')
		rawarray = np.vstack((np.hstack((self.rawdata[index][0,:,0],self.rawdata[index][0,:,1])),
			np.hstack((self.rawdata[index][1,:,0],self.rawdata[index][1,:,1])) ))
		np.savetxt(outfile+'_raw-pulse.txt',rawarray.T, fmt='%.6f', header = 'time(s)\tTemp(K)', 
			delimiter=', ')

	def savetraces(self, outfile, Barray = 'all'):
		if (Barray == 'All' or Barray == 'all'):
			Barray = np.sort(list(set(self.Bfield)))
		else:
			Barray = np.array(Barray)

		f = open(outfile, 'w')
		f.write('# Heat Capacity data from LongHCPulse\n')
		f.write('# Temp_heating (K),\t C_heating,\t Temp_cooling (K),\t C_cooling\n')
		for jj in range(len(self.Bfield)):
			B = self.Bfield[jj]
			for b in Barray:
				if  (B == b) and (self.ShortPulse[jj] == 0):
					f.write('\n# B='+str(B)+'K, curve number = '+str(jj)+'\n')
					for ii in range(len(self.T[jj][:,0])):
						f.write(str(self.T[jj][ii,0]) +',\t'+ str(self.HC[jj][ii,0])+',\t'+
							str(self.T[jj][ii,1]) +',\t'+ str(self.HC[jj][ii,1])+'\n')
		f.close()


	def saveData(self, outfile):
		#Combine traces into single line (if not done already)
		try:
			self.CombHC
		except AttributeError:
			print(" combining traces...")
			self._combineTraces(1)

		#so we don't save anything that hasn't been created
		members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]
		membersToBeSaved = ['rawdata','smoothedData','Bfield', 'Blabels', 'CombB', 'CombHC', 'CombT', 'HC', 
							'ScaledBfield', 'ThermCondWire', 'avgThermCondWire', 'Kw', 'AddendaHC',
							'ScaledCombB', 'ScaledavgSpB', 'ShortPulse', 'T', 'Tsamp', 'avgSpB', 'avgSpHc', 
							'avgSpT', 'molarmass', 'sampmass', 'shortpulseBflag', 'shortpulselabels', 
							'entropylabels', 'Entropy','labels']
		membersToBeSaved = list(set(members).intersection(membersToBeSaved))  

		# Make list of data to be saved
		dataList =eval('[self.'+ ', self.'.join(membersToBeSaved)+']')
		dataToSave = dict(zip(membersToBeSaved, dataList))

		# Edit outfile so it has the extension .pickle
		try:
			periodindex = outfile.index('.')
			outfile = outfile[:periodindex] + '.pickle'
		except ValueError:
			outfile = outfile + '.pickle'

		# Save data with pickle
		with open(outfile, 'w') as f:
			pickle.dump(dataToSave, f)

