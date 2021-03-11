* Project: 	GxInsurance
* Content:  HRS 2006 LB file
* Author: 	Pia Arce
* Date: 	November 27, 2020

/*Extract HRS LB file and save it as .dta 
NOTE: Set own wd
*/

*************************************************************************
******************************** PREAMBLE *******************************
*************************************************************************
global DIR "T:/econ/biroli/geighei/data/HRS/phenoraw/"
global SAVEHRSLB "T:/econ/biroli/geighei/code/GxInsurance/1_data/"


*2006
clear all
global HRSLB06 "$DIR/h06core"
infile using $HRSLB06/h06sta/H06LB_R.dct , using($HRSLB06/h06da/H06LB_R.da)
save "$SAVEHRSLB/H06LB_R.dta", replace

* 2008
clear all
global HRSLB08 "$DIR/h08core"
infile using $HRSLB08/h08sta/H08LB_R.dct  , using($HRSLB08/h08da/H08LB_R.da)
save "$SAVEHRSLB/H08LB_R.dta", replace

* 2010
clear all
global HRSLB10 "$DIR/h10core"
infile using $HRSLB10/h10sta/H10LB_R.dct  , using($HRSLB10/h10da/H10LB_R.da)
save "$SAVEHRSLB/H10LB_R.dta", replace

* 2012
clear all
global HRSLB12 "$DIR/h12core"
infile using $HRSLB12/h12sta/H12LB_R.dct  , using($HRSLB12/h12da/H12LB_R.da)
save "$SAVEHRSLB/H12LB_R.dta", replace

* 2014
clear all
global HRSLB14 "$DIR/h14core"
infile using $HRSLB14/h14sta/H14LB_R.dct  , using($HRSLB14/h14da/H14LB_R.da)
save "$SAVEHRSLB/H14LB_R.dta", replace

* 2016
clear all
global HRSLB16 "$DIR/h16core"
infile using $HRSLB16/h16sta/H16LB_R.dct  , using($HRSLB16/h16da/H16LB_R.da)
save "$SAVEHRSLB/H16LB_R.dta", replace
