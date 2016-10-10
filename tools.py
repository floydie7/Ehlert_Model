import numpy as np
import os
from astropy.io import fits

def ABToFlux(mag):
	'''Converts from AB magnitude to flux in uJy.
	
	INPUT
	------------------------------------------------------
	
	mag: float
		AB magnitude
	
	OUTPUT
	------------------------------------------------------
	
	flux: float
		Flux in uJy
	'''
	flux = 10**(0.4*(23.9-mag))
	return flux

def ABToVega(AB,band):
	'''Converts from an AB magnitude to a Vega magnitude in one of the \
supported filters.
	
	INPUT
	------------------------------------------------------
	
	AB: float
		Magnitude in AB
	band: str, kwarg
		The filter in which to convert. See 'vegaOffsets' for a list of \
		currently supported filters
	
	OUTPUT
	------------------------------------------------------
	
	Vega: float
		Vega magnitude in the given band
	'''
	Vega = AB - vegaOffsets(band)
	return Vega

def angDiaDist(z,H0 = 70.0,Om = 0.3,Ol = 0.7):
	'''Derivative of the comoving distance function that calculates \
the angular diameter scale at a given redshift in kpc/arcsec.
	
	INPUT
	------------------------------------------------------
	
	z: float
		The redshift in question
	H0: float
		The Hubble constant in units of km/s/Mpc. Default is 70.0
	Om: float
		Omega_matter. Default is 0.3
	Ol: float
		Omega_lambda. Default is 0.7
	
	OUTPUT
	------------------------------------------------------
	
	DA: float
		The angluar scale in kpc/" at the given redshift
	'''
	DAM = comovingDistance(z,H0,Om,Ol)/(1+z)
	DA = DAM/206.265 #Converting from Mpc/radian to kpc/"
	return DA
	
def area(covmap,mincov):
	'''Takes a coverage map and calculaes the total area \
(in square arcmins) that has been covered above a given level.

	INPUT
	------------------------------------------------------
	
	covmap: str
		The filename of the coverage map.
	nimcov: float
		The minimum number of exposures needed to be counted in \
		the area
	
	OUTPUT
	------------------------------------------------------
	
	ar: float
		The area in squae arcminutes
	'''
	cov = fits.getdata(covmap)
	goodnum = 0
	goodxpix = []
	goodypix = []
	for n in range(len(cov)):
		for m in range(len(cov[n])):
			if cov[n][m] >= mincov:
				goodnum += 1
				goodxpix.append(n)
				goodypix.append(m)
	try:
		xpix = fits.getval(covmap,'PXSCAL1')
		ypix = fits.getval(covmap,'PXSCAL2')
	except:
		KeyError
		#This is because a few fits headers haven't had 'PXSCAL' and
		#I don't want it to throw up an error after doing the hard bit
		try:
			xpix = fits.getval(covmap,'CDELT1')*3600
			ypix = fits.getval(covmap,'CDELT2')*3600
		except:
			KeyError
			return goodnum,'The number of pixels is {0}. The scale could not be determined.'.format(goodnum)
	ar = goodnum*xpix*ypix*-1/3600
	return ar

def comovingDistance(z,H0 = 70.0,Om = 0.3,Ol = 0.7):
	'''Calculates the comoving radial distance at a given redshift in Mpc.
	
	INPUT
	------------------------------------------------------
	
	z: float
		The redshift in question
	H0: float
		The Hubble constant in units of km/s/Mpc. Default is 70.0
	Om: float
		Omega_matter. Default is 0.3
	Ol: float
		Omega_lambda. Default is 0.7
	
	OUTPUT
	------------------------------------------------------
	
	DC: float
		The comoving distance in Mpc at the given redshift
	'''
	c = 300000 #km/s
	intpart = 0
	for zp in np.arange(0.0005,z,0.001):
		intpart += (np.sqrt(Om*(1+zp)**3 + Ol)**-1)*0.001
	DC = c*intpart/(H0)
	return DC

def coordList(coord,sep):
	'''Separates a sexagesimal coordinate string and returns \
a list of three floats. How the segments of the coordinate are \
separated must be specified.

	INPUT
	------------------------------------------------------
	
	coord: str
		Simple sexagesimal coordinate, either RA or Dec. \
		EG: '13:45:02' or '+43 12 05.2'
	sep: str
		Whatever string is separating the components of the \
		coordinate. EG: ':' or ' '
	
	OUTPUT
	------------------------------------------------------
	
	Coordinate: list
		A list of three floats: hours/degrees, minutes, seconds.
	'''
	hd = np.float128(coord.split(sep)[0])
	mm = np.float128(coord.split(sep)[1])
	ss = np.float128(coord.split(sep)[2])
	Coordinate = [hd,mm,ss]
	return Coordinate

def decToSex(RA, Dec):
	'''Takes RA and dec in decimal and converts to sexagesimal. Returns \
a tuple of two strings.

	INPUT
	------------------------------------------------------
	
	RA: float or int
		The decimal right ascension. Automatically converted into a \
		numpy float128.
	Dec: float or int
		The decimal right ascension. Automatically converted into a \
		numpy float128.
	
	OUTPUT
	------------------------------------------------------
	
	RightAsc: str
		The sexagesimal right ascension
	Declination: str
		The sexagesimal declination
	'''
	RAH = np.floor(np.float128(RA)*24/360)
	RAM = np.floor((np.float128(RA)*24/360 - RAH)*60)
	RASec = ((np.float128(RA)*24/360 - RAH)*60 - RAM)*60
	RAS = round(RASec,2)
	#Maths to convert RA into HMS. Note the variables are left as floats
	#until they are displayed.
	if RAH < 10:
		HrStr = '0'+str(int(RAH))
	else:
		HrStr = str(int(RAH))
	if RAM < 10:
		RMStr = '0'+str(int(RAM))
	else:
		RMStr = str(int(RAM))
	if RAS < 10:
		RSStr = '0'+str(RAS)
	else:
		RSStr = str(RAS)
	#All the if statements are just to add leading zeros where
	#aesthetically appropriate
	RightAsc = HrStr + ':' + RMStr + ':' + RSStr
	DD = np.floor(np.abs(np.float128(Dec)))
	DM = np.floor((np.abs(np.float128(Dec)) - DD)*60)
	DSec = ((np.abs(Dec) - DD)*60 - DM)*60
	DS = round(DSec,1)
	#Maths to convert to DMS. Variables left as floats as above, but also
	#note that the absolute value of the declinaion is used to ensure
	#correct maths. Negative values are returned at the display stage.
	if Dec > 0:
		if DD >= 10:
			DegStr = '+' + str(int(DD))
		else:
			DegStr = '+0' + str(int(DD))
	elif Dec < 0:
		if DD >= 10:
			DegStr = '-' + str(int(DD))
		else:
			DegStr = '-0' + str(int(DD))
	else:
		Declination = '+00:00:00.0'
	#In addition to getting the leading zero right, this also adds +/-
	if DM < 10:
		DMStr = '0'+str(int(DM))
	else:
		DMStr = str(int(DM))
	if DS < 10:
		DSStr = '0' + str(DS)
	else:
		DSStr = str(DS)
	#More hoops through which to jump for neat display
	Declination = DegStr + ':' + DMStr + ':' + DSStr
	return RightAsc,Declination

def distance(coord1,coord2):
	'''Calculates the true distance (in arcminutes) between two specified \
coordinates. Possibly vulnerable to floating point precision difficulties \
at very small distances.

	INPUT
	------------------------------------------------------
	
	coord1: list or tuple of two elements
		This can be a list, tuple or np array, but needs to have the RA \
		and Dec both in decimal form with the RA first.
	coord2: list or tuple of two elements
		Same as coord1
	
	OUTPUT
	------------------------------------------------------
	
	dmin: float
		Angular distance, in arcminutes
	'''
	ra1 = coord1[0]*np.pi/180
	ra2 = coord2[0]*np.pi/180
	dec1 = coord1[1]*np.pi/180
	dec2 = coord2[1]*np.pi/180
	#Putting the coordinates into radians
	sins = np.sin(dec1)*np.sin(dec2)
	coss = np.cos(dec1)*np.cos(dec2)
	dellon = np.abs(ra2-ra1)
	#Splitting up some of the terms in the equation
	drad = np.arccos(sins+coss*np.cos(dellon))
	#Distace in radians
	ddeg = drad*180/np.pi
	dmin = ddeg*60   
	#Converting from radiants to degrees then to arcmin
	return dmin

def fluxToAB(flux):
	'''Takes flux in uJy and converts to AB magnitude.

	INPUT
	------------------------------------------------------
	
	flux: float
		The flux in uJy
	
	OUTPUT
	------------------------------------------------------
	
	AB: float
		AB magnitude
	'''
	AB = -2.5*np.log10(flux) + 23.9
	return AB
    
def fluxToVega(flux,band='ch1'):
	'''Converts flux in uJy to Vega magnitude in Irac channel one or two. \
Only those two cahnnels are supported because I don't have the offsets to \
hand for the others.

	INPUT
	------------------------------------------------------
	
	flux: float
		The flux in uJy
	band: str, kwarg
		The filter in which to convert. See 'vegaOffsets' for a list of \
		currently supported filters
	
	OUTPUT
	------------------------------------------------------
	
	Vega: float
		Vega magnitude in the appropriate channel
	'''
	Vega = fluxToAB(flux) - vegaOffsets(band)
	return Vega

def lookbackTime(z,H0 = 70.0,Om = 0.3,Ol = 0.7):
	'''Calculates the lookback time to a given redshift.

	INPUT
	------------------------------------------------------
	
	z: float
		redshift
	H0: float
		The Hubble constant in units of km/s/Mpc. Default is 70.0
	Om: float
		Omega_matter. Default is 0.3
	Ol: float
		Omega_lambda. Default is 0.7
	
	OUTPUT
	------------------------------------------------------
	
	tL: float
		Lookback time
	'''
	H0t = H0*(3.16*10**7)*(10**9)/(3.09*10**19) #km/s/Mpc -> 1/Gyr
	tH = 1/H0t
	intpart = 0
	for zp in np.arange(0.0005,z,0.001):
		intpart += ((np.sqrt(Om*(1+zp)**3 + Ol)*(1+zp))**-1)*0.001
	tL = intpart*tH
	return tL

def luminosityDistance(z,H0 = 70.0,Om = 0.3,Ol = 0.7):
	'''A derivative of the comoving distance function that \
returns luminostiy distance in Mpc.

	INPUT
	------------------------------------------------------
	
	z: float
		redshift
	H0: float
		The Hubble constant in units of km/s/Mpc. Default is 70.0
	Om: float
		Omega_matter. Default is 0.3
	Ol: float
		Omega_lambda. Default is 0.7
	
	OUTPUT
	------------------------------------------------------
	
	DL: float
		Luminosity distance in Mpc
	'''
	DL = comovingDistance(z,H0,Om,Ol)*(1+z)
	return DL

def m500Tor500(m500,z,H0 = 70.0,Om = 0.3,Ol = 0.7):
	'''Finds r500 in Mpc for an m500 given in solar masses. Largely taken \
from Chris Greer's code.

	INPUT
	------------------------------------------------------
	
	m500: float
		The m500 in question, in solar mass units	
	z: float
		redshift
	H0: float
		The Hubble constant in units of km/s/Mpc. Default is 70.0
	Om: float
		Omega_matter. Default is 0.3
	Ol: float
		Omega_lambda. Default is 0.7
	
	OUTPUT
	------------------------------------------------------
	
	r500: float
		The equivalent r500, in Mpc
	'''
	Mpc = 3.08568e+24 #cm
	G = 6.67259e-8 #cm^3 g^-1 s^-2
	Msun = 1.9891e+33 #grams
	H0 = H0*1000.0*100.0/Mpc #km/s/Mpc -> s^-1
	Ez = np.sqrt((Om*(1.0+z)**(3.0))+Ol)
	Hz = H0*Ez
	rc = (3.0*Hz**2)/(8.0*np.pi*G)
	M = m500*Msun #Msun -> grams
	R = ((3.0 * M / (4.0 * np.pi * rc * 500.0))**(1.0/3.0)) #cm
	r500 = R/Mpc #cm -> Mpc
	return r500

def ptSlope(coord1,coord2):
	'''Takes (x1,y1) & (x2,y2) and returns the equation for a line \
that connects the points.

	INPUT
	------------------------------------------------------
	
	coord1: two-element list/tuple/array
		The x and y coordinates of the first point, each a float or an int
	coord2: two-element list/tuple/array
		The x and y coordinates of the second point, each a float or an int
	
	OUTPUT
	------------------------------------------------------
	
	eq: str
		The equation written in y = mx + b form
	'''
	x1 = float(coord1[0])
	y1 = float(coord1[1])
	x2 = float(coord2[0])
	y2 = float(coord2[1])
	m = (y2 - y1)/(x2 - x1)
	b = m*-1*x1 + y1
	eq = 'y = {0:2.2f}x + {1:2.2f}'.format(m,b)
	return eq

def r500Tom500(r500,z,H0 = 70.0,Om = 0.3,Ol = 0.7):
	'''Finds m500 in solar masses for an r500 given in Mpc. Largely taken \
from Chris Greer's code.

	INPUT
	------------------------------------------------------
	
	r500: float
		The r500 in question, in Mpc
	z: float
		redshift
	H0: float
		The Hubble constant in units of km/s/Mpc. Default is 70.0
	Om: float
		Omega_matter. Default is 0.3
	Ol: float
		Omega_lambda. Default is 0.7
	
	OUTPUT
	------------------------------------------------------
	
	m500: float
		The equivalent m500, in solar masses
	'''
	Mpc = 3.08568e+24 #cm
	G = 6.67259e-8 #cm^3 g^-1 s^-2
	Msun = 1.9891e+33 #grams
	H0 = H0*1000.0*100.0/Mpc #km/s/Mpc -> s^-1
	Ez = np.sqrt((Om*(1.0+z)**(3.0))+Ol)
	Hz = H0*Ez
	rc = (3.0*Hz**2)/(8.0*np.pi*G)
	R = r500*Mpc #Mpc -> cm
	M = (4.0/3.0)*np.pi*rc*500*R**3 #grams
	m500 = M/Msun #grams -> Msol
	return m500

def recArea(corner1,corner2):
	'''Takes two corners and returns the area (in sq degrees) of a \
rectangle defined by those ra and dec lines. That is, corners of \
(145,40) and (155,50) will give the area of a patch bounded by \
dec = 40d in the south, dec = 50d in the north, RA = 155d in the \
east and RA = 145d in the west.

	INPUT
	------------------------------------------------------
	
	corner1: two-element list/tuple/array
		The RA and Dec of the first corner, in decimal form
	corner2: two-element list/tuple/array
		The RA and Dec of the second corner, in decimal form
	
	OUTPUT
	------------------------------------------------------
	
	ard2: float
		The rectangular area in square degrees
	'''
	decN = np.max([corner1[1],corner2[1]])*np.pi/180
	decS = np.min([corner1[1],corner2[1]])*np.pi/180
	raE = np.max([corner1[0],corner2[0]])*np.pi/180
	raW = np.min([corner1[0],corner2[0]])*np.pi/180
	#Putting the RAs and Decs into radians
	arrad2 = (np.sin(decN) - np.sin(decS))*(raE - raW)
	#Equation for rectangular solid angle, in sq rads
	ard2 = arrad2*(180/np.pi)**2
	#Into sq degs
	return ard2

def sexToDec(RA, Dec):
	'''Takes sexagesimal RA and Dec strings and returns a tuple of RA and Dec \
decimal floats.

	INPUT
	------------------------------------------------------
	
	RA: str
		Simple RA coordinate string in sexagesimal.
	Dec: str
		Simple Dec string in sexagesimal
	
	OUTPUT
	------------------------------------------------------
	
	RightAsc: float
		The right ascension in degrees-decimal, rounded to five places
	Declination: float
		The declination in degrees-decimal, rounded as the ra
	'''
	RAdel = RA[2]
	Ddel = Dec[3]
	#Identifies how the coordinates are written, noting Dec has a +/-
	RAHMS = coordList(RA, RAdel)
	RightAsc = 360*(RAHMS[0] + RAHMS[1]/60 + RAHMS[2]/3600)/24
	DecDMS = coordList(Dec, Ddel)
	if Dec[0] == '+':
		Declination = DecDMS[0] + DecDMS[1]/60 + DecDMS[2]/3600
	elif Dec[0] == '-':
		Declination = DecDMS[0] - DecDMS[1]/60 - DecDMS[2]/3600
	#Easy maths, put /slightly/ different for positve and negative
	RightAsc = round(RightAsc,5)    
	Declination = round(Declination,5)
	#So as not to get a stupidly long output
	return RightAsc,Declination

def smallDist(coord1,coord2):
	'''Takes a pair of decimal coordinates and calculates the small \
angle approximation distance between them in arcmins.

	INPUT
	------------------------------------------------------
	
	coord1: list or tuple of two elements
		This can be a list, tuple or np array, but needs to have the RA \
		and Dec both in decimal form with the RA first.
	coord2: list or tuple of two elements
		Same as coord1
	
	OUTPUT
	------------------------------------------------------
	
	dm: float
		Distance in arcminutes
	'''
	if np.abs(coord1[1]) >= np.abs(coord2[1]):
		DecRad = coord1[1]*np.pi/180
	else:
		DecRad = coord2[1]*np.pi/180
	#Converts the higher dec into radians
	d = np.sqrt((coord2[0]-coord1[0])**2*(np.cos(DecRad))**2 + (coord2[1]-coord1[1])**2)
	#Small angle approximation for distance. 
	dm = d*60   
	return dm

def vegaOffsets(band):
	'''The offset between AB magnitude and Vega magnitude; ie: AB - Vega.
	
	INPUT
	------------------------------------------------------
	
	band: str, kwarg
		['ch1' | 'ch2' | 'b' | 'v' | 'r' | 'i' | 'j' | 'h' | 'ks' \
		| 'sloan_u' | 'sloan_g' | 'sloan_r' | 'sloan_i' | 'sloan_z' ]
	
	OUTPUT
	------------------------------------------------------
	
	vegoff: float
		The offset
	'''
	offset = {'ch1':2.79,'ch2':3.26,'b':-0.09,'v':0.02,'r':0.21,'i':0.45,\
	'j':0.91,'h':1.39,'ks':1.85,'sloan_u':0.91,'sloan_g':-0.08,\
	'sloan_r':0.16,'sloan_i':0.37,'sloan_z':0.54}
	vegoff = offset[band]
	return vegoff

def vegaToAB(Vega,band):
	'''Converts from a Vega magnitude to an AB magnitude in one of the \
supported filters.
	
	INPUT
	------------------------------------------------------
	
	Vega: float
		Magnitude in Vega
	band: str, kwarg
		The filter in which to convert. See 'vegaOffsets' for a list of \
		currently supported filters
	
	OUTPUT
	------------------------------------------------------
	
	AB: float
		AB magnitude in the given band
	'''
	AB = Vega + vegaOffsets(band)
	return AB

def vegaToFlux(mag,band='ch1'):
	'''Converts a Vega magnitude from a given filter into flux \
in uJy.

	INPUT
	------------------------------------------------------
	
	mag: float
		Magnitude in Vega
	band: str, kwarg, optional
		The filter in which to convert. See 'vegaOffsets' for a list of \
		currently supported filters
	
	OUTPUT
	------------------------------------------------------
	
	flux: float
		flux in uJy
	'''
	flux = 10**(0.4*(23.9 - vegaOffsets(band) - mag))
	return flux

def writeReg(catalogue,name,raind='ra',decind='dec'):
	'''Writes a list of coordinates from a catalogue to a file.

	INPUT
	------------------------------------------------------
	
	catalogue = numpy table or list of lists
		The catalogue with the coordinates. This is written for use \
		with a numpy table, but I think will also work with a list of \
		lists style table.
	name: str
		The name of the file to be written. No ending is needed, a \
		'.txt' is added automatically.
	raind: str or int, optional
		For a numpy table this is the column name of RAs. For a list \
		of lists it will be the numerical index. Defaults to 'ra'.
	decind: str or int, optional
		For a numpy table this is the column name of decs. For a list \
		of lists it will be the numerical index. Defaults to 'dec'
	
	OUTPUT
	------------------------------------------------------
	
	output: str
		A bit of text to confirm the file was written
	'''
	cwd = os.getcwd()
	filename = name + '.txt'
	reg = open(filename,'w')
	for obj in catalogue:
		coords = str(obj[raind]) + ',' + str(obj[decind]) + '\n'
		reg.write(coords)
	reg.close()
	os.chdir(cwd)
	output = 'Coordinates for ' + filename + ' has been written.'
	return output
    
