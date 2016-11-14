package gosolar

import (
	"math"
)

/////////////////////////////////////////////
//      Solar Position Algorithm (SPA)     //
//                   for                   //
//        Solar Radiation Application      //
//                                         //
//               May 12, 2003              //
//                                         //
//   Filename: SPA.C                       //
//                                         //
//   Afshin Michael Andreas                //
//   Afshin.Andreas@NREL.gov (303)384-6383 //
//                                         //
//   Measurement & Instrumentation Team    //
//   Solar Radiation Research Laboratory   //
//   National Renewable Energy Laboratory  //
//   1617 Cole Blvd, Golden, CO 80401      //
/////////////////////////////////////////////

/////////////////////////////////////////////
//   See the SPA.H header file for usage   //
//                                         //
//   This code is based on the NREL        //
//   technical report "Solar Position      //
//   Algorithm for Solar Radiation         //
//   Application" by I. Reda & A. Andreas  //
/////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//
//   NOTICE
//   Copyright (C) 2008-2011 Alliance for Sustainable Energy, LLC, All Rights Reserved
//
//The Solar Position Algorithm ("Software") is code in development prepared by employees of the
//Alliance for Sustainable Energy, LLC, (hereinafter the "Contractor"), under Contract No.
//DE-AC36-08GO28308 ("Contract") with the U.S. Department of Energy (the "DOE"). The United
//States Government has been granted for itself and others acting on its behalf a paid-up, non-
//exclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative
//works, and perform publicly and display publicly. Beginning five (5) years after the date
//permission to assert copyright is obtained from the DOE, and subject to any subsequent five
//(5) year renewals, the United States Government is granted for itself and others acting on
//its behalf a paid-up, non-exclusive, irrevocable, worldwide license in the Software to
//reproduce, prepare derivative works, distribute copies to the public, perform publicly and
//display publicly, and to permit others to do so. If the Contractor ceases to make this
//computer software available, it may be obtained from DOE's Office of Scientific and Technical
//Information's Energy Science and Technology Software Center (ESTSC) at P.O. Box 1020, Oak
//Ridge, TN 37831-1020. THIS SOFTWARE IS PROVIDED BY THE CONTRACTOR "AS IS" AND ANY EXPRESS OR
//IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CONTRACTOR OR THE
//U.S. GOVERNMENT BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//WHATSOEVER, INCLUDING BUT NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA OR PROFITS,
//WHICH MAY RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR OTHER TORTIOUS CLAIM THAT ARISES
//OUT OF OR IN CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
//
//The Software is being provided for internal, noncommercial purposes only and shall not be
//re-distributed. Please contact Jennifer Ramsey (Jennifer.Ramsey@nrel.gov) in the NREL
//Commercialization and Technology Transfer Office for information concerning a commercial
//license to use the Software.
//
//As a condition of using the Software in an application, the developer of the application
//agrees to reference the use of the Software and make this Notice readily accessible to any
//end-user in a Help|About screen or equivalent manner.
//
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// Revised 27-FEB-2004 Andreas
//         Added bounds check on inputs and return value for spa_calculate().
// Revised 10-MAY-2004 Andreas
//         Changed temperature bound check minimum from -273.15 to -273 degrees C.
// Revised 17-JUN-2004 Andreas
//         Corrected a problem that caused a bogus sunrise/set/transit on the equinox.
// Revised 18-JUN-2004 Andreas
//         Added a "function" input variable that allows the selecting of desired outputs.
// Revised 21-JUN-2004 Andreas
//         Added 3 new intermediate output values to SPA structure (srha, ssha, & sta).
// Revised 23-JUN-2004 Andreas
//         Enumerations for "function" were renamed and 2 were added.
//         Prevented bound checks on inputs that are not used (based on function).
// Revised 01-SEP-2004 Andreas
//         Changed a local variable from integer to double.
// Revised 12-JUL-2005 Andreas
//         Put a limit on the EOT calculation, so that the result is between -20 and 20.
// Revised 26-OCT-2005 Andreas
//         Set the atmos. refraction correction to zero, when sun is below horizon.
//         Made atmos_refract input a requirement for all "functions".
//         Changed atmos_refract bound check from +/- 10 to +/- 5 degrees.
// Revised 07-NOV-2006 Andreas
//         Corrected 3 earth periodic terms in the L_TERMS array.
//         Corrected 2 earth periodic terms in the R_TERMS array.
// Revised 10-NOV-2006 Andreas
//         Corrected a constant used to calculate topocentric sun declination.
//         Put a limit on observer hour angle, so result is between 0 and 360.
// Revised 13-NOV-2006 Andreas
//         Corrected calculation of topocentric sun declination.
//         Converted all floating point inputs in spa structure to doubles.
// Revised 27-FEB-2007 Andreas
//         Minor correction made as to when atmos. refraction correction is set to zero.
// Revised 21-JAN-2008 Andreas
//         Minor change to two variable declarations.
// Revised 12-JAN-2009 Andreas
//         Changed timezone bound check from +/-12 to +/-18 hours.
// Revised 14-JAN-2009 Andreas
//         Corrected a constant used to calculate ecliptic mean obliquity.
// Revised 01-APR-2013 Andreas
//		   Replace floor with new integer function for tech. report consistency, no affect on results.
//         Add "utility" function prototypes to header file for use with NREL's SAMPA.
//         Rename 4 "utility" function names (remove "sun") for clarity with NREL's SAMPA.
//		   Added delta_ut1 as required input, which the fractional second difference between UT and UTC.
//         Time must be input w/o delta_ut1 adjustment, instead of assuming adjustment was pre-applied.
// Revised 10-JUL-2014 Andreas
//         Change second in Spa_data structure from an integer to double to allow fractional second
// Revised 08-SEP-2014 Andreas
//         Corrected description of azm_rotation in header file
//         Limited azimuth180 to range of 0 to 360 deg (instead of -180 to 180) for tech report consistency
//         Changed all variables names from azimuth180 to azimuth_astro
//         Renamed 2 "utility" function names for consistency
///////////////////////////////////////////////////////////////////////////////////////////////
const (
	SPA_ZA     = iota //calculate zenith and azimuth
	SPA_ZA_INC        //calculate zenith, azimuth, and incidence
	SPA_ZA_RTS        //calculate zenith, azimuth, and sun rise/transit/set values
	SPA_ALL           //calculate all SPA output values
)

type Spa_data struct {
	//----------------------INPUT VALUES------------------------

	Year   int     // 4-digit year,      valid range: -2000 to 6000, error code: 1
	Month  int     // 2-digit month,         valid range: 1 to  12,  error code: 2
	Day    int     // 2-digit day,           valid range: 1 to  31,  error code: 3
	Hour   int     // Observer local hour,   valid range: 0 to  24,  error code: 4
	Minute int     // Observer local minute, valid range: 0 to  59,  error code: 5
	Second float64 // Observer local second, valid range: 0 to <60,  error code: 6

	Delta_ut1 float64 // Fractional second difference between UTC and UT which is used
	// to adjust UTC for earth's irregular rotation rate and is derived
	// from observation only and is reported in this bulletin:
	// http://maia.usno.navy.mil/ser7/ser7.dat,
	// where delta_ut1 = DUT1
	// valid range: -1 to 1 second (exclusive), error code 17

	Delta_t float64 // Difference between earth rotation time and terrestrial time
	// It is derived from observation only and is reported in this
	// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
	// where delta_t = 32.184 + (TAI-UTC) - DUT1
	// valid range: -8000 to 8000 seconds, error code: 7

	Timezone float64 // Observer time zone (negative west of Greenwich)
	// valid range: -18   to   18 hours,   error code: 8

	Longitude float64 // Observer longitude (negative west of Greenwich)
	// valid range: -180  to  180 degrees, error code: 9

	Latitude float64 // Observer latitude (negative south of equator)
	// valid range: -90   to   90 degrees, error code: 10

	Elevation float64 // Observer elevation [meters]
	// valid range: -6500000 or higher meters,    error code: 11

	Pressure float64 // Annual average local pressure [millibars]
	// valid range:    0 to 5000 millibars,       error code: 12

	Temperature float64 // Annual average local temperature [degrees Celsius]
	// valid range: -273 to 6000 degrees Celsius, error code 13

	Slope float64 // Surface slope (measured from the horizontal plane)
	// valid range: -360 to 360 degrees, error code: 14

	Azm_rotation float64 // Surface azimuth rotation (measured from south to projection of
	//     surface normal on horizontal plane, negative east)
	// valid range: -360 to 360 degrees, error code: 15

	Atmos_refract float64 // Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
	// valid range: -5   to   5 degrees, error code: 16

	Function int // Switch to choose functions for desired output (from enumeration)

	//-----------------Intermediate OUTPUT VALUES--------------------

	Jd float64 //Julian day
	jc float64 //Julian century

	jde float64 //Julian ephemeris day
	jce float64 //Julian ephemeris century
	jme float64 //Julian ephemeris millennium

	L float64 //earth heliocentric longitude [degrees]
	B float64 //earth heliocentric latitude [degrees]
	R float64 //earth radius vector [Astronomical Units, AU]

	theta float64 //geocentric longitude [degrees]
	beta  float64 //geocentric latitude [degrees]

	x0 float64 //mean elongation (moon-sun) [degrees]
	x1 float64 //mean anomaly (sun) [degrees]
	x2 float64 //mean anomaly (moon) [degrees]
	x3 float64 //argument latitude (moon) [degrees]
	x4 float64 //ascending longitude (moon) [degrees]

	Del_psi     float64 //nutation longitude [degrees]
	Del_epsilon float64 //nutation obliquity [degrees]
	epsilon0    float64 //ecliptic mean obliquity [arc seconds]
	Epsilon     float64 //ecliptic true obliquity  [degrees]

	del_tau float64 //aberration correction [degrees]
	lamda   float64 //apparent sun longitude [degrees]
	nu0     float64 //Greenwich mean sidereal time [degrees]
	nu      float64 //Greenwich sidereal time [degrees]

	alpha float64 //geocentric sun right ascension [degrees]
	delta float64 //geocentric sun declination [degrees]

	H           float64 //observer hour angle [degrees]
	xi          float64 //sun equatorial horizontal parallax [degrees]
	del_alpha   float64 //sun right ascension parallax [degrees]
	delta_prime float64 //topocentric sun declination [degrees]
	alpha_prime float64 //topocentric sun right ascension [degrees]
	h_prime     float64 //topocentric local hour angle [degrees]

	e0    float64 //topocentric elevation angle (uncorrected) [degrees]
	del_e float64 //atmospheric refraction correction [degrees]
	e     float64 //topocentric elevation angle (corrected) [degrees]

	eot  float64 //equation of time [minutes]
	srha float64 //sunrise hour angle [degrees]
	ssha float64 //sunset hour angle [degrees]
	sta  float64 //sun transit altitude [degrees]

	//---------------------Final OUTPUT VALUES------------------------

	Zenith        float64 //topocentric zenith angle [degrees]
	azimuth_astro float64 //topocentric azimuth angle (westward from south) [for astronomers]
	Azimuth       float64 //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
	Incidence     float64 //surface incidence angle [degrees]

	suntransit float64 //local sun transit time (or solar noon) [fractional hour]
	Sunrise    float64 //local sunrise time (+/- 30 seconds) [fractional hour]
	Sunset     float64 //local sunset time (+/- 30 seconds) [fractional hour]

}

const (
	PI         = 3.1415926535897932384626433832795028841971
	SUN_RADIUS = 0.26667

	L_COUNT = 6
	B_COUNT = 2
	R_COUNT = 5
	Y_COUNT = 63

	L_MAX_SUBCOUNT = 64
	B_MAX_SUBCOUNT = 5
	R_MAX_SUBCOUNT = 40
)

const (
	TERM_A = iota
	TERM_B
	TERM_C
	TERM_COUNT
)

const (
	TERM_X0 = iota
	TERM_X1
	TERM_X2
	TERM_X3
	TERM_X4
	TERM_X_COUNT
)

const (
	TERM_PSI_A = iota
	TERM_PSI_B
	TERM_EPS_C
	TERM_EPS_D
	TERM_PE_COUNT
)

const (
	JD_MINUS = iota
	JD_ZERO
	JD_PLUS
	JD_COUNT
)

const (
	SUN_TRANSIT = iota
	SUN_RISE
	SUN_SET
	SUN_COUNT
)

const TERM_Y_COUNT = TERM_X_COUNT

var l_subcount = [L_COUNT]int{64, 34, 20, 7, 3, 1}
var b_subcount = [B_COUNT]int{5, 2}
var r_subcount = [R_COUNT]int{40, 10, 6, 2, 1}

///////////////////////////////////////////////////
///  Earth Periodic Terms
///////////////////////////////////////////////////
// var L_TERMS = [L_COUNT][L_MAX_SUBCOUNT][TERM_COUNT]float64{{{175347046.0, 0, 0},
var L_TERMS = [][][]float64{{{175347046.0, 0, 0},
	{3341656.0, 4.6692568, 6283.07585},
	{34894.0, 4.6261, 12566.1517},
	{3497.0, 2.7441, 5753.3849},
	{3418.0, 2.8289, 3.5231},
	{3136.0, 3.6277, 77713.7715},
	{2676.0, 4.4181, 7860.4194},
	{2343.0, 6.1352, 3930.2097},
	{1324.0, 0.7425, 11506.7698},
	{1273.0, 2.0371, 529.691},
	{1199.0, 1.1096, 1577.3435},
	{990, 5.233, 5884.927},
	{902, 2.045, 26.298},
	{857, 3.508, 398.149},
	{780, 1.179, 5223.694},
	{753, 2.533, 5507.553},
	{505, 4.583, 18849.228},
	{492, 4.205, 775.523},
	{357, 2.92, 0.067},
	{317, 5.849, 11790.629},
	{284, 1.899, 796.298},
	{271, 0.315, 10977.079},
	{243, 0.345, 5486.778},
	{206, 4.806, 2544.314},
	{205, 1.869, 5573.143},
	{202, 2.458, 6069.777},
	{156, 0.833, 213.299},
	{132, 3.411, 2942.463},
	{126, 1.083, 20.775},
	{115, 0.645, 0.98},
	{103, 0.636, 4694.003},
	{102, 0.976, 15720.839},
	{102, 4.267, 7.114},
	{99, 6.21, 2146.17},
	{98, 0.68, 155.42},
	{86, 5.98, 161000.69},
	{85, 1.3, 6275.96},
	{85, 3.67, 71430.7},
	{80, 1.81, 17260.15},
	{79, 3.04, 12036.46},
	{75, 1.76, 5088.63},
	{74, 3.5, 3154.69},
	{74, 4.68, 801.82},
	{70, 0.83, 9437.76},
	{62, 3.98, 8827.39},
	{61, 1.82, 7084.9},
	{57, 2.78, 6286.6},
	{56, 4.39, 14143.5},
	{56, 3.47, 6279.55},
	{52, 0.19, 12139.55},
	{52, 1.33, 1748.02},
	{51, 0.28, 5856.48},
	{49, 0.49, 1194.45},
	{41, 5.37, 8429.24},
	{41, 2.4, 19651.05},
	{39, 6.17, 10447.39},
	{37, 6.04, 10213.29},
	{37, 2.57, 1059.38},
	{36, 1.71, 2352.87},
	{36, 1.78, 6812.77},
	{33, 0.59, 17789.85},
	{30, 0.44, 83996.85},
	{30, 2.74, 1349.87},
	{25, 3.16, 4690.48}},
	{{628331966747.0, 0, 0},
		{206059.0, 2.678235, 6283.07585},
		{4303.0, 2.6351, 12566.1517},
		{425.0, 1.59, 3.523},
		{119.0, 5.796, 26.298},
		{109.0, 2.966, 1577.344},
		{93, 2.59, 18849.23},
		{72, 1.14, 529.69},
		{68, 1.87, 398.15},
		{67, 4.41, 5507.55},
		{59, 2.89, 5223.69},
		{56, 2.17, 155.42},
		{45, 0.4, 796.3},
		{36, 0.47, 775.52},
		{29, 2.65, 7.11},
		{21, 5.34, 0.98},
		{19, 1.85, 5486.78},
		{19, 4.97, 213.3},
		{17, 2.99, 6275.96},
		{16, 0.03, 2544.31},
		{16, 1.43, 2146.17},
		{15, 1.21, 10977.08},
		{12, 2.83, 1748.02},
		{12, 3.26, 5088.63},
		{12, 5.27, 1194.45},
		{12, 2.08, 4694},
		{11, 0.77, 553.57},
		{10, 1.3, 6286.6},
		{10, 4.24, 1349.87},
		{9, 2.7, 242.73},
		{9, 5.64, 951.72},
		{8, 5.3, 2352.87},
		{6, 2.65, 9437.76},
		{6, 4.67, 4690.48}},
	{{52919.0, 0, 0},
		{8720.0, 1.0721, 6283.0758},
		{309.0, 0.867, 12566.152},
		{27, 0.05, 3.52},
		{16, 5.19, 26.3},
		{16, 3.68, 155.42},
		{10, 0.76, 18849.23},
		{9, 2.06, 77713.77},
		{7, 0.83, 775.52},
		{5, 4.66, 1577.34},
		{4, 1.03, 7.11},
		{4, 3.44, 5573.14},
		{3, 5.14, 796.3},
		{3, 6.05, 5507.55},
		{3, 1.19, 242.73},
		{3, 6.12, 529.69},
		{3, 0.31, 398.15},
		{3, 2.28, 553.57},
		{2, 4.38, 5223.69},
		{2, 3.75, 0.98}},
	{{289.0, 5.844, 6283.076},
		{35, 0, 0},
		{17, 5.49, 12566.15},
		{3, 5.2, 155.42},
		{1, 4.72, 3.52},
		{1, 5.3, 18849.23},
		{1, 5.97, 242.73}},
	{{114.0, 3.142, 0},
		{8, 4.13, 6283.08},
		{1, 3.84, 12566.15}},
	{{1, 3.14, 0}}}

// var B_TERMS = [B_COUNT][B_MAX_SUBCOUNT][TERM_COUNT]float64{{{280.0, 3.199, 84334.662},
var B_TERMS = [][][]float64{{{280.0, 3.199, 84334.662},
	{102.0, 5.422, 5507.553},
	{80, 3.88, 5223.69},
	{44, 3.7, 2352.87},
	{32, 4, 1577.34}},
	{{9, 3.9, 5507.55},
		{6, 1.73, 5223.69}}}

// var R_TERMS = [R_COUNT][R_MAX_SUBCOUNT][TERM_COUNT]float64{{{100013989.0, 0, 0},
var R_TERMS = [][][]float64{{{100013989.0, 0, 0},
	{1670700.0, 3.0984635, 6283.07585},
	{13956.0, 3.05525, 12566.1517},
	{3084.0, 5.1985, 77713.7715},
	{1628.0, 1.1739, 5753.3849},
	{1576.0, 2.8469, 7860.4194},
	{925.0, 5.453, 11506.77},
	{542.0, 4.564, 3930.21},
	{472.0, 3.661, 5884.927},
	{346.0, 0.964, 5507.553},
	{329.0, 5.9, 5223.694},
	{307.0, 0.299, 5573.143},
	{243.0, 4.273, 11790.629},
	{212.0, 5.847, 1577.344},
	{186.0, 5.022, 10977.079},
	{175.0, 3.012, 18849.228},
	{110.0, 5.055, 5486.778},
	{98, 0.89, 6069.78},
	{86, 5.69, 15720.84},
	{86, 1.27, 161000.69},
	{65, 0.27, 17260.15},
	{63, 0.92, 529.69},
	{57, 2.01, 83996.85},
	{56, 5.24, 71430.7},
	{49, 3.25, 2544.31},
	{47, 2.58, 775.52},
	{45, 5.54, 9437.76},
	{43, 6.01, 6275.96},
	{39, 5.36, 4694},
	{38, 2.39, 8827.39},
	{37, 0.83, 19651.05},
	{37, 4.9, 12139.55},
	{36, 1.67, 12036.46},
	{35, 1.84, 2942.46},
	{33, 0.24, 7084.9},
	{32, 0.18, 5088.63},
	{32, 1.78, 398.15},
	{28, 1.21, 6286.6},
	{28, 1.9, 6279.55},
	{26, 4.59, 10447.39}},
	{{103019.0, 1.10749, 6283.07585},
		{1721.0, 1.0644, 12566.1517},
		{702.0, 3.142, 0},
		{32, 1.02, 18849.23},
		{31, 2.84, 5507.55},
		{25, 1.32, 5223.69},
		{18, 1.42, 1577.34},
		{10, 5.91, 10977.08},
		{9, 1.42, 6275.96},
		{9, 0.27, 5486.78}},
	{{4359.0, 5.7846, 6283.0758},
		{124.0, 5.579, 12566.152},
		{12, 3.14, 0},
		{9, 3.63, 77713.77},
		{6, 1.87, 5573.14},
		{3, 5.47, 18849.23}},
	{{145.0, 4.273, 6283.076},
		{7, 3.92, 12566.15}},
	{{4, 2.56, 6283.08}}}

////////////////////////////////////////////////////////////////
///  Periodic Terms for the nutation in longitude and obliquity
////////////////////////////////////////////////////////////////

// var Y_TERMS = [Y_COUNT][TERM_Y_COUNT]int{{0, 0, 0, 0, 1},
var Y_TERMS = [][]int{{0, 0, 0, 0, 1},
	{-2, 0, 0, 2, 2},
	{0, 0, 0, 2, 2},
	{0, 0, 0, 0, 2},
	{0, 1, 0, 0, 0},
	{0, 0, 1, 0, 0},
	{-2, 1, 0, 2, 2},
	{0, 0, 0, 2, 1},
	{0, 0, 1, 2, 2},
	{-2, -1, 0, 2, 2},
	{-2, 0, 1, 0, 0},
	{-2, 0, 0, 2, 1},
	{0, 0, -1, 2, 2},
	{2, 0, 0, 0, 0},
	{0, 0, 1, 0, 1},
	{2, 0, -1, 2, 2},
	{0, 0, -1, 0, 1},
	{0, 0, 1, 2, 1},
	{-2, 0, 2, 0, 0},
	{0, 0, -2, 2, 1},
	{2, 0, 0, 2, 2},
	{0, 0, 2, 2, 2},
	{0, 0, 2, 0, 0},
	{-2, 0, 1, 2, 2},
	{0, 0, 0, 2, 0},
	{-2, 0, 0, 2, 0},
	{0, 0, -1, 2, 1},
	{0, 2, 0, 0, 0},
	{2, 0, -1, 0, 1},
	{-2, 2, 0, 2, 2},
	{0, 1, 0, 0, 1},
	{-2, 0, 1, 0, 1},
	{0, -1, 0, 0, 1},
	{0, 0, 2, -2, 0},
	{2, 0, -1, 2, 1},
	{2, 0, 1, 2, 2},
	{0, 1, 0, 2, 2},
	{-2, 1, 1, 0, 0},
	{0, -1, 0, 2, 2},
	{2, 0, 0, 2, 1},
	{2, 0, 1, 0, 0},
	{-2, 0, 2, 2, 2},
	{-2, 0, 1, 2, 1},
	{2, 0, -2, 0, 1},
	{2, 0, 0, 0, 1},
	{0, -1, 1, 0, 0},
	{-2, -1, 0, 2, 1},
	{-2, 0, 0, 0, 1},
	{0, 0, 2, 2, 1},
	{-2, 0, 2, 0, 1},
	{-2, 1, 0, 2, 1},
	{0, 0, 1, -2, 0},
	{-1, 0, 1, 0, 0},
	{-2, 1, 0, 0, 0},
	{1, 0, 0, 0, 0},
	{0, 0, 1, 2, 0},
	{0, 0, -2, 2, 2},
	{-1, -1, 1, 0, 0},
	{0, 1, 1, 0, 0},
	{0, -1, 1, 2, 2},
	{2, -1, -1, 2, 2},
	{0, 0, 3, 2, 2},
	{2, -1, 0, 2, 2}}

// var PE_TERMS = [Y_COUNT][TERM_PE_COUNT]float64{{-171996, -174.2, 92025, 8.9},
var PE_TERMS = [][]float64{{-171996, -174.2, 92025, 8.9},
	{-13187, -1.6, 5736, -3.1},
	{-2274, -0.2, 977, -0.5},
	{2062, 0.2, -895, 0.5},
	{1426, -3.4, 54, -0.1},
	{712, 0.1, -7, 0},
	{-517, 1.2, 224, -0.6},
	{-386, -0.4, 200, 0},
	{-301, 0, 129, -0.1},
	{217, -0.5, -95, 0.3},
	{-158, 0, 0, 0},
	{129, 0.1, -70, 0},
	{123, 0, -53, 0},
	{63, 0, 0, 0},
	{63, 0.1, -33, 0},
	{-59, 0, 26, 0},
	{-58, -0.1, 32, 0},
	{-51, 0, 27, 0},
	{48, 0, 0, 0},
	{46, 0, -24, 0},
	{-38, 0, 16, 0},
	{-31, 0, 13, 0},
	{29, 0, 0, 0},
	{29, 0, -12, 0},
	{26, 0, 0, 0},
	{-22, 0, 0, 0},
	{21, 0, -10, 0},
	{17, -0.1, 0, 0},
	{16, 0, -8, 0},
	{-16, 0.1, 7, 0},
	{-15, 0, 9, 0},
	{-13, 0, 7, 0},
	{-12, 0, 6, 0},
	{11, 0, 0, 0},
	{-10, 0, 5, 0},
	{-8, 0, 3, 0},
	{7, 0, -3, 0},
	{-7, 0, 0, 0},
	{-7, 0, 3, 0},
	{-7, 0, 3, 0},
	{6, 0, 0, 0},
	{6, 0, -3, 0},
	{6, 0, -3, 0},
	{-6, 0, 3, 0},
	{-6, 0, 3, 0},
	{5, 0, 0, 0},
	{-5, 0, 3, 0},
	{-5, 0, 3, 0},
	{-5, 0, 3, 0},
	{4, 0, 0, 0},
	{4, 0, 0, 0},
	{4, 0, 0, 0},
	{-4, 0, 0, 0},
	{-4, 0, 0, 0},
	{-4, 0, 0, 0},
	{3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0}}

///////////////////////////////////////////////

func rad2deg(radians float64) float64 {
	return (180.0 / PI) * radians
}

func deg2rad(degrees float64) float64 {
	return (PI / 180.0) * degrees
}

func integer(value float64) int {
	return int(value)
}

func limit_degrees(degrees float64) float64 {
	var limited float64

	degrees /= 360.0
	limited = 360.0 * (degrees - math.Floor(degrees))
	if limited < 0 {
		limited += 360.0
	}

	return limited
}

func limit_degrees180pm(degrees float64) float64 {
	var limited float64

	degrees /= 360.0
	limited = 360.0 * (degrees - math.Floor(degrees))
	if limited < -180.0 {
		limited += 360.0
	} else if limited > 180.0 {
		limited -= 360.0
	}

	return limited
}

func limit_degrees180(degrees float64) float64 {
	var limited float64

	degrees /= 180.0
	limited = 180.0 * (degrees - math.Floor(degrees))
	if limited < 0 {
		limited += 180.0
	}

	return limited
}

func limit_zero2one(value float64) float64 {
	var limited float64

	limited = value - math.Floor(value)
	if limited < 0 {
		limited += 1.0
	}

	return limited
}

func limit_minutes(minutes float64) float64 {
	limited := minutes

	if limited < -20.0 {
		limited += 1440.0
	} else if limited > 20.0 {
		limited -= 1440.0
	}

	return limited
}

func dayfrac_to_local_hr(dayfrac, timezone float64) float64 {
	return 24.0 * limit_zero2one(dayfrac+timezone/24.0)
}

func third_order_polynomial(a, b, c, d, x float64) float64 {
	return ((a*x+b)*x+c)*x + d
}

///////////////////////////////////////////////////////////////////////////////////////////////
func validate_inputs(spa *Spa_data) int {
	if (spa.Year < -2000) || (spa.Year > 6000) {
		return 1
	}
	if (spa.Month < 1) || (spa.Month > 12) {
		return 2
	}
	if (spa.Day < 1) || (spa.Day > 31) {
		return 3
	}
	if (spa.Hour < 0) || (spa.Hour > 24) {
		return 4
	}
	if (spa.Minute < 0) || (spa.Minute > 59) {
		return 5
	}
	if (spa.Second < 0) || (spa.Second >= 60) {
		return 6
	}
	if (spa.Pressure < 0) || (spa.Pressure > 5000) {
		return 12
	}
	if (spa.Temperature <= -273) || (spa.Temperature > 6000) {
		return 13
	}
	if (spa.Delta_ut1 <= -1) || (spa.Delta_ut1 >= 1) {
		return 17
	}
	if (spa.Hour == 24) && (spa.Minute > 0) {
		return 5
	}
	if (spa.Hour == 24) && (spa.Second > 0) {
		return 6
	}

	if math.Abs(spa.Delta_t) > 8000 {
		return 7
	}
	if math.Abs(spa.Timezone) > 18 {
		return 8
	}
	if math.Abs(spa.Longitude) > 180 {
		return 9
	}
	if math.Abs(spa.Latitude) > 90 {
		return 10
	}
	if math.Abs(spa.Atmos_refract) > 5 {
		return 16
	}
	if spa.Elevation < -6500000 {
		return 11
	}

	if (spa.Function == SPA_ZA_INC) || (spa.Function == SPA_ALL) {
		if math.Abs(spa.Slope) > 360 {
			return 14
		}
		if math.Abs(spa.Azm_rotation) > 360 {
			return 15
		}
	}

	return 0
}

///////////////////////////////////////////////////////////////////////////////////////////////
func julian_day(year, month, day, hour, minute int, second, dut1, tz float64) float64 {
	var day_decimal, julian_day, a float64

	day_decimal = float64(day) + (float64(hour)-tz+(float64(minute)+(second+dut1)/60.0)/60.0)/24.0

	if month < 3 {
		month += 12
		year--
	}

	// julian_day = int(365.25*(float64(year)+4716.0)) + int(30.6001*(float64(month)+1)) + day_decimal - 1524.5
	julian_day = float64(int(365.25*(float64(year)+4716.0))) + float64(int((30.6001 * (float64(month) + 1)))) + day_decimal - 1524.5

	if julian_day > 2299160.0 {
		// a = int(year / 100)
		a = float64(int(year / 100))
		// julian_day += (2 - a + int(a/4))
		julian_day += (2 - a + float64(int(a/4)))
	}

	return julian_day
}

func julian_century(jd float64) float64 {
	return (jd - 2451545.0) / 36525.0
}

func julian_ephemeris_day(jd, delta_t float64) float64 {
	return jd + delta_t/86400.0
}

func julian_ephemeris_century(jde float64) float64 {
	return (jde - 2451545.0) / 36525.0
}

func julian_ephemeris_millennium(jce float64) float64 {
	return (jce / 10.0)
}

func earth_periodic_term_summation(terms [][]float64, count int, jme float64) float64 { // terms[][TERM_COUNT]
	var i int
	sum := 0.0

	for i = 0; i < count; i++ {
		sum += terms[i][TERM_A] * math.Cos(terms[i][TERM_B]+terms[i][TERM_C]*jme)
	}

	return sum
}

func earth_values(term_sum []float64, count int, jme float64) float64 {
	var i int
	sum := 0.0

	for i = 0; i < count; i++ {
		sum += term_sum[i] * math.Pow(jme, float64(i))
	}

	sum /= 1.0 * math.Pow10(8) //maybe math.Exp

	return sum
}

func earth_heliocentric_longitude(jme float64) float64 {
	sum := make([]float64, L_COUNT)
	var i int

	for i = 0; i < L_COUNT; i++ {
		sum[i] = earth_periodic_term_summation(L_TERMS[i], l_subcount[i], jme)
	}

	return limit_degrees(rad2deg(earth_values(sum, L_COUNT, jme)))

}

func earth_heliocentric_latitude(jme float64) float64 {
	sum := make([]float64, B_COUNT)
	var i int

	for i = 0; i < B_COUNT; i++ {
		sum[i] = earth_periodic_term_summation(B_TERMS[i], b_subcount[i], jme)
	}

	return rad2deg(earth_values(sum, B_COUNT, jme))

}

func earth_radius_vector(jme float64) float64 {
	sum := make([]float64, R_COUNT)
	var i int

	for i = 0; i < R_COUNT; i++ {
		sum[i] = earth_periodic_term_summation(R_TERMS[i], r_subcount[i], jme)
	}

	return earth_values(sum, R_COUNT, jme)

}

func geocentric_longitude(l float64) float64 {
	theta := l + 180.0

	if theta >= 360.0 {
		theta -= 360.0
	}

	return theta
}

func geocentric_latitude(b float64) float64 {
	return -b
}

func mean_elongation_moon_sun(jce float64) float64 {
	return third_order_polynomial(1.0/189474.0, -0.0019142, 445267.11148, 297.85036, jce)
}

func mean_anomaly_sun(jce float64) float64 {
	return third_order_polynomial(-1.0/300000.0, -0.0001603, 35999.05034, 357.52772, jce)
}

func mean_anomaly_moon(jce float64) float64 {
	return third_order_polynomial(1.0/56250.0, 0.0086972, 477198.867398, 134.96298, jce)
}

func argument_latitude_moon(jce float64) float64 {
	return third_order_polynomial(1.0/327270.0, -0.0036825, 483202.017538, 93.27191, jce)
}

func ascending_longitude_moon(jce float64) float64 {
	return third_order_polynomial(1.0/450000.0, 0.0020708, -1934.136261, 125.04452, jce)
}

func xy_term_summation(i int, x []float64) float64 { // x[TERM_X_COUNT]
	var j int
	sum := 0.0

	for j = 0; j < TERM_Y_COUNT; j++ {
		sum += x[j] * float64(Y_TERMS[i][j])
	}

	return sum
}

func nutation_longitude_and_obliquity(jce float64, x []float64, del_psi, del_epsilon *float64) { // x[TERM_X_COUNT]
	var i int
	var xy_term_sum float64
	sum_psi := 0.0
	sum_epsilon := 0.0

	for i = 0; i < Y_COUNT; i++ {
		xy_term_sum = deg2rad(xy_term_summation(i, x))
		sum_psi += (float64(PE_TERMS[i][TERM_PSI_A]) + jce*float64(PE_TERMS[i][TERM_PSI_B])) * math.Sin(xy_term_sum)
		sum_epsilon += (float64(PE_TERMS[i][TERM_EPS_C]) + jce*float64(PE_TERMS[i][TERM_EPS_D])) * math.Cos(xy_term_sum)
	}

	*del_psi = sum_psi / 36000000.0
	*del_epsilon = sum_epsilon / 36000000.0
}

func ecliptic_mean_obliquity(jme float64) float64 {
	u := jme / 10.0

	return 84381.448 + u*(-4680.93+u*(-1.55+u*(1999.25+u*(-51.38+u*(-249.67+
		u*(-39.05+u*(7.12+u*(27.87+u*(5.79+u*2.45)))))))))
}

func ecliptic_true_obliquity(delta_epsilon, epsilon0 float64) float64 {
	return delta_epsilon + epsilon0/3600.0
}

func aberration_correction(r float64) float64 {
	return -20.4898 / (3600.0 * r)
}

func apparent_sun_longitude(theta, delta_psi, delta_tau float64) float64 {
	return theta + delta_psi + delta_tau
}

func greenwich_mean_sidereal_time(jd, jc float64) float64 {
	return limit_degrees(280.46061837 + 360.98564736629*(jd-2451545.0) +
		jc*jc*(0.000387933-jc/38710000.0))
}

func greenwich_sidereal_time(nu0, delta_psi, epsilon float64) float64 {
	return nu0 + delta_psi*math.Cos(deg2rad(epsilon))
}

func geocentric_right_ascension(lamda, epsilon, beta float64) float64 {
	lamda_rad := deg2rad(lamda)
	epsilon_rad := deg2rad(epsilon)

	return limit_degrees(rad2deg(math.Atan2(math.Sin(lamda_rad)*math.Cos(epsilon_rad)-
		math.Tan(deg2rad(beta))*math.Sin(epsilon_rad), math.Cos(lamda_rad))))
}

func geocentric_declination(beta, epsilon, lamda float64) float64 {
	beta_rad := deg2rad(beta)
	epsilon_rad := deg2rad(epsilon)

	return rad2deg(math.Asin(math.Sin(beta_rad)*math.Cos(epsilon_rad) +
		math.Cos(beta_rad)*math.Sin(epsilon_rad)*math.Sin(deg2rad(lamda))))
}

func observer_hour_angle(nu, longitude, alpha_deg float64) float64 {
	return limit_degrees(nu + longitude - alpha_deg)
}

func sun_equatorial_horizontal_parallax(r float64) float64 {
	return 8.794 / (3600.0 * r)
}

func right_ascension_parallax_and_topocentric_dec(latitude, elevation, xi, h, delta float64, delta_alpha, delta_prime *float64) {
	var delta_alpha_rad float64
	lat_rad := deg2rad(latitude)
	xi_rad := deg2rad(xi)
	h_rad := deg2rad(h)
	delta_rad := deg2rad(delta)
	u := math.Atan(0.99664719 * math.Tan(lat_rad))
	y := 0.99664719*math.Sin(u) + elevation*math.Sin(lat_rad)/6378140.0
	x := math.Cos(u) + elevation*math.Cos(lat_rad)/6378140.0

	delta_alpha_rad = math.Atan2(-x*math.Sin(xi_rad)*math.Sin(h_rad),
		math.Cos(delta_rad)-x*math.Sin(xi_rad)*math.Cos(h_rad))

	*delta_prime = rad2deg(math.Atan2((math.Sin(delta_rad)-y*math.Sin(xi_rad))*math.Cos(delta_alpha_rad),
		math.Cos(delta_rad)-x*math.Sin(xi_rad)*math.Cos(h_rad)))

	*delta_alpha = rad2deg(delta_alpha_rad)
}

func topocentric_right_ascension(alpha_deg, delta_alpha float64) float64 {
	return alpha_deg + delta_alpha
}

func topocentric_local_hour_angle(h, delta_alpha float64) float64 {

	return h - delta_alpha
}

func topocentric_elevation_angle(latitude, delta_prime, h_prime float64) float64 {
	lat_rad := deg2rad(latitude)
	delta_prime_rad := deg2rad(delta_prime)

	return rad2deg(math.Asin(math.Sin(lat_rad)*math.Sin(delta_prime_rad) +
		math.Cos(lat_rad)*math.Cos(delta_prime_rad)*math.Cos(deg2rad(h_prime))))
}

func atmospheric_refraction_correction(pressure, temperature, atmos_refract, e0 float64) float64 {
	del_e := 0.0

	PT := pressure/(temperature+273.15)

	if (e0 < 15.0 && e0 > -2.5) {
		del_e = PT * (0.1594+0.0196*e0+math.Pow(e0, 2)*2E-5) / (1+0.505*e0+0.0845*math.Pow(e0,2))
	}
	if (e0 >= 15.0 && e0 < 90.0) {
		del_e = 0.00452*PT/math.Tan(deg2rad(e0))
	}

	// if e0 >= -1*(SUN_RADIUS+atmos_refract) {
	// 	a := pressure * 283.0 * 1.02
	// 	b := 1010.0 * (temperature + 273.15) * 60.0 * math.Tan(deg2rad(e0+10.3/(e0+5.11)))
	// 	del_e = a / b
	// 	//del_e = (pressure / 1010.0) * (283.0 / (273.15 + temperature)) *
	// 	//	1.02 / (60.0 * math.Tan(deg2rad(e0+10.3/(e0+5.11))))
	// }
	//fmt.Println(del_e)
	return del_e
}

func topocentric_elevation_angle_corrected(e0, delta_e float64) float64 {
	return e0 + delta_e
}

func topocentric_zenith_angle(e float64) float64 {
	return 90.0 - e
}

func topocentric_azimuth_angle_astro(h_prime, latitude, delta_prime float64) float64 {
	h_prime_rad := deg2rad(h_prime)
	lat_rad := deg2rad(latitude)

	//
	// math.Cos(h_prime_rad)*math.Sin(lat_rad)-math.Tan(deg2rad(delta_prime))*math.Cos(lat_rad)))

	return limit_degrees(rad2deg(math.Atan2(math.Sin(h_prime_rad),
		math.Cos(h_prime_rad)*math.Sin(lat_rad)-math.Tan(deg2rad(delta_prime))*math.Cos(lat_rad))))
}

func topocentric_azimuth_angle(azimuth_astro float64) float64 {
	return limit_degrees(azimuth_astro + 180.0)
}

func surface_incidence_angle(zenith, azimuth_astro, azm_rotation, slope float64) float64 {
	zenith_rad := deg2rad(zenith)
	slope_rad := deg2rad(slope)

	return rad2deg(math.Acos(math.Cos(zenith_rad)*math.Cos(slope_rad) +
		math.Sin(slope_rad)*math.Sin(zenith_rad)*math.Cos(deg2rad(azimuth_astro-azm_rotation))))
}

func sun_mean_longitude(jme float64) float64 {
	return limit_degrees(280.4664567 + jme*(360007.6982779+jme*(0.03032028+
		jme*(1/49931.0+jme*(-1/15300.0+jme*(-1/2000000.0))))))
}

func eot(m, alpha, del_psi, epsilon float64) float64 {
	return limit_minutes(4.0 * (m - 0.0057183 - alpha + del_psi*math.Cos(deg2rad(epsilon))))
}

func approx_sun_transit_time(alpha_zero, longitude, nu float64) float64 {
	return (alpha_zero - longitude - nu) / 360.0
}

func sun_hour_angle_at_rise_set(latitude, delta_zero, h0_prime float64) float64 {
	h0 := -99999.0
	latitude_rad := deg2rad(latitude)
	delta_zero_rad := deg2rad(delta_zero)
	argument := (math.Sin(deg2rad(h0_prime)) - math.Sin(latitude_rad)*math.Sin(delta_zero_rad)) /
		(math.Cos(latitude_rad) * math.Cos(delta_zero_rad))

	if math.Abs(argument) <= 1 {
		h0 = limit_degrees180(rad2deg(math.Acos(argument)))
	}

	return h0
}

func approx_sun_rise_and_set(m_rts []float64, h0 float64) {
	h0_dfrac := h0 / 360.0

	m_rts[SUN_RISE] = limit_zero2one(m_rts[SUN_TRANSIT] - h0_dfrac)
	m_rts[SUN_SET] = limit_zero2one(m_rts[SUN_TRANSIT] + h0_dfrac)
	m_rts[SUN_TRANSIT] = limit_zero2one(m_rts[SUN_TRANSIT])
}

func rts_alpha_delta_prime(ad []float64, n float64) float64 {
	a := ad[JD_ZERO] - ad[JD_MINUS]
	b := ad[JD_PLUS] - ad[JD_ZERO]

	if math.Abs(a) >= 2.0 {
		a = limit_zero2one(a)
	}
	if math.Abs(b) >= 2.0 {
		b = limit_zero2one(b)
	}

	return ad[JD_ZERO] + n*(a+b+(b-a)*n)/2.0
}

func rts_sun_altitude(latitude, delta_prime, h_prime float64) float64 {
	latitude_rad := deg2rad(latitude)
	delta_prime_rad := deg2rad(delta_prime)

	return rad2deg(math.Asin(math.Sin(latitude_rad)*math.Sin(delta_prime_rad) +
		math.Cos(latitude_rad)*math.Cos(delta_prime_rad)*math.Cos(deg2rad(h_prime))))
}

func sun_rise_and_set(m_rts, h_rts, delta_prime []float64, latitude float64, h_prime []float64, h0_prime float64, sun int) float64 {
	return m_rts[sun] + (h_rts[sun]-h0_prime)/
		(360.0*math.Cos(deg2rad(delta_prime[sun]))*math.Cos(deg2rad(latitude))*math.Sin(deg2rad(h_prime[sun])))
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate required SPA parameters to get the right ascension (alpha) and declination (delta)
// Note: JD must be already calculated and in structure
////////////////////////////////////////////////////////////////////////////////////////////////
func calculate_geocentric_sun_right_ascension_and_declination(spa *Spa_data) {
	var x = make([]float64, TERM_X_COUNT)

	spa.jc = julian_century(spa.Jd)

	spa.jde = julian_ephemeris_day(spa.Jd, spa.Delta_t)
	spa.jce = julian_ephemeris_century(spa.jde)
	spa.jme = julian_ephemeris_millennium(spa.jce)

	spa.L = earth_heliocentric_longitude(spa.jme)
	spa.B = earth_heliocentric_latitude(spa.jme)
	spa.R = earth_radius_vector(spa.jme)

	spa.theta = geocentric_longitude(spa.L)
	spa.beta = geocentric_latitude(spa.B)

	x[TERM_X0], spa.x0 = mean_elongation_moon_sun(spa.jce), mean_elongation_moon_sun(spa.jce)
	x[TERM_X1], spa.x1 = mean_anomaly_sun(spa.jce), mean_anomaly_sun(spa.jce)
	x[TERM_X2], spa.x2 = mean_anomaly_moon(spa.jce), mean_anomaly_moon(spa.jce)
	x[TERM_X3], spa.x3 = argument_latitude_moon(spa.jce), argument_latitude_moon(spa.jce)
	x[TERM_X4], spa.x4 = ascending_longitude_moon(spa.jce), ascending_longitude_moon(spa.jce)

	nutation_longitude_and_obliquity(spa.jce, x, &(spa.Del_psi), &(spa.Del_epsilon))

	spa.epsilon0 = ecliptic_mean_obliquity(spa.jme)
	spa.Epsilon = ecliptic_true_obliquity(spa.Del_epsilon, spa.epsilon0)

	spa.del_tau = aberration_correction(spa.R)
	spa.lamda = apparent_sun_longitude(spa.theta, spa.Del_psi, spa.del_tau)
	spa.nu0 = greenwich_mean_sidereal_time(spa.Jd, spa.jc)
	spa.nu = greenwich_sidereal_time(spa.nu0, spa.Del_psi, spa.Epsilon)

	spa.alpha = geocentric_right_ascension(spa.lamda, spa.Epsilon, spa.beta)
	spa.delta = geocentric_declination(spa.beta, spa.Epsilon, spa.lamda)
}

////////////////////////////////////////////////////////////////////////
// Calculate Equation of Time (EOT) and Sun Rise, Transit, & Set (RTS)
////////////////////////////////////////////////////////////////////////

func calculate_eot_and_sun_rise_transit_set(spa *Spa_data) {
	var sun_rts Spa_data
	var nu, m, h0, n float64
	var alpha, delta = make([]float64, JD_COUNT), make([]float64, JD_COUNT)
	var m_rts, nu_rts, h_rts = make([]float64, SUN_COUNT), make([]float64, SUN_COUNT), make([]float64, SUN_COUNT)
	var alpha_prime, delta_prime, h_prime = make([]float64, SUN_COUNT), make([]float64, SUN_COUNT), make([]float64, SUN_COUNT)
	h0_prime := -1 * (SUN_RADIUS + spa.Atmos_refract)
	var i int

	sun_rts = *spa
	m = sun_mean_longitude(spa.jme)
	spa.eot = eot(m, spa.alpha, spa.Del_psi, spa.Epsilon)

	sun_rts.Hour, sun_rts.Minute, sun_rts.Second = 0, 0, 0
	sun_rts.Delta_ut1, sun_rts.Timezone = 0.0, 0.0

	sun_rts.Jd = julian_day(sun_rts.Year, sun_rts.Month, sun_rts.Day, sun_rts.Hour,
		sun_rts.Minute, sun_rts.Second, sun_rts.Delta_ut1, sun_rts.Timezone)

	calculate_geocentric_sun_right_ascension_and_declination(&sun_rts)
	nu = sun_rts.nu

	sun_rts.Delta_t = 0
	sun_rts.Jd--
	for i = 0; i < JD_COUNT; i++ {
		calculate_geocentric_sun_right_ascension_and_declination(&sun_rts)
		alpha[i] = sun_rts.alpha
		delta[i] = sun_rts.delta
		sun_rts.Jd++
	}

	m_rts[SUN_TRANSIT] = approx_sun_transit_time(alpha[JD_ZERO], spa.Longitude, nu)
	h0 = sun_hour_angle_at_rise_set(spa.Latitude, delta[JD_ZERO], h0_prime)

	if h0 >= 0 {

		approx_sun_rise_and_set(m_rts, h0)

		for i = 0; i < SUN_COUNT; i++ {

			nu_rts[i] = nu + 360.985647*m_rts[i]

			n = m_rts[i] + spa.Delta_t/86400.0
			alpha_prime[i] = rts_alpha_delta_prime(alpha, n)
			delta_prime[i] = rts_alpha_delta_prime(delta, n)

			h_prime[i] = limit_degrees180pm(nu_rts[i] + spa.Longitude - alpha_prime[i])

			h_rts[i] = rts_sun_altitude(spa.Latitude, delta_prime[i], h_prime[i])
		}

		spa.srha = h_prime[SUN_RISE]
		spa.ssha = h_prime[SUN_SET]
		spa.sta = h_rts[SUN_TRANSIT]

		spa.suntransit = dayfrac_to_local_hr(m_rts[SUN_TRANSIT]-h_prime[SUN_TRANSIT]/360.0,
			spa.Timezone)

		spa.Sunrise = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,
			spa.Latitude, h_prime, h0_prime, SUN_RISE), spa.Timezone)

		spa.Sunset = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,
			spa.Latitude, h_prime, h0_prime, SUN_SET), spa.Timezone)

	} else {
		spa.srha, spa.ssha, spa.sta, spa.suntransit, spa.Sunrise, spa.Sunset = -99999, -99999, -99999, -99999, -99999, -99999
	}

}

///////////////////////////////////////////////////////////////////////////////////////////
// Calculate all SPA parameters and put into structure
// Note: All inputs values (listed in header file) must already be in structure
///////////////////////////////////////////////////////////////////////////////////////////
func Spa_calculate(spa *Spa_data) int {
	var result int

	result = validate_inputs(spa)

	if result == 0 {
		spa.Jd = julian_day(spa.Year, spa.Month, spa.Day, spa.Hour,
			spa.Minute, spa.Second, spa.Delta_ut1, spa.Timezone)

		calculate_geocentric_sun_right_ascension_and_declination(spa)

		spa.H = observer_hour_angle(spa.nu, spa.Longitude, spa.alpha)

		spa.xi = sun_equatorial_horizontal_parallax(spa.R)

		right_ascension_parallax_and_topocentric_dec(spa.Latitude, spa.Elevation, spa.xi,
			spa.H, spa.delta, &(spa.del_alpha), &(spa.delta_prime))


		spa.alpha_prime = topocentric_right_ascension(spa.alpha, spa.del_alpha)
		spa.h_prime = topocentric_local_hour_angle(spa.H, spa.del_alpha)

		spa.e0 = topocentric_elevation_angle(spa.Latitude, spa.delta_prime, spa.h_prime)
		spa.del_e = atmospheric_refraction_correction(spa.Pressure, spa.Temperature,
			spa.Atmos_refract, spa.e0)
		spa.e = topocentric_elevation_angle_corrected(spa.e0, spa.del_e)

		spa.Zenith = topocentric_zenith_angle(spa.e)
		spa.azimuth_astro = topocentric_azimuth_angle_astro(spa.h_prime, spa.Latitude,
			spa.delta_prime)

		spa.Azimuth = topocentric_azimuth_angle(spa.azimuth_astro)

		if (spa.Function == SPA_ZA_INC) || (spa.Function == SPA_ALL) {
			spa.Incidence = surface_incidence_angle(spa.Zenith, spa.azimuth_astro,
				spa.Azm_rotation, spa.Slope)
		}

		if (spa.Function == SPA_ZA_RTS) || (spa.Function == SPA_ALL) {
			calculate_eot_and_sun_rise_transit_set(spa)
		}

	}

	return result
}

///////////////////////////////////////////////////////////////////////////////////////////
