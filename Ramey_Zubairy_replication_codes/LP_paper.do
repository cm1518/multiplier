*** Understanding the Size of the Government Spending Multiplier: It's in the Sign, (R. Barnichon, D. Debortoli and C. Matthes), Review of Economic Studies
*** This code estimates asymmetric effects using Ramey-Zubairy news shocks to defense spending (fig 4 of paper)
*** Code adapted from Ramey-Zubairy "Government Spending Multipliers in Good Times and in Bad: Evidence from U.S. Historical Data" JPE, April 2018
*** Requires:
*** rzdat.xlsx
********************************************************************************

*uncomment and adjust the line below to directly save the results in the correct folder
cd "C:\Users\Rage\Dropbox\replication_material_ReStud_final_for_regis\Ramey_Zubairy_replication_codes"

#delimit;

drop _all;
clear all;

/*
*install necessary packages
ssc install ivreg2;
ssc install ranktest;
*/

set more 1;
set matsize 800;

capture log close;
log using results.log, replace;


/*******************************************************************************
  SET PARAMETERS THAT GOVERN SPECIFICATION
*******************************************************************************/

local p = 8; /*number of lags of control variables*/

local hor = 30; /*horizon of irf */


*******************************************************************************;
** RAW DATA IMPORTATION AND DATA SETUP;
*******************************************************************************;

import excel rzdat.xls, sheet("rzdat") firstrow;

drop if quarter<1889;

gen qdate = q(1889q1) + _n-1;
tsset qdate, q;


gen nomit = 0;  /* indicator for no omit */


*** DEFINE QUARTIC TREND;

gen t = _n;
gen t2 = t^2;
gen t3 = t^3;
gen t4 = t^4;


*** NORMALIZATION;

/* choice of potential GDP for normalization: 
6th degree polynomial, fitted excluding Great Depression, WWII:  quarter>=1930 & quarter<1947 */
local ynorm rgdp_pott6; 


* BASIC VARIABLES;

gen newsy = news/(L.`ynorm'*L.pgdp);
gen rgov = ngov/pgdp;
gen rtax = nfedcurrreceipts_nipa/pgdp;
gen taxy = nfedcurrreceipts_nipa/ngdp;
gen debty = pubfeddebt_treas/L.ngdp;
gen lpgdp = ln(pgdp);
gen ly = ln(rgdp);
gen infl = 400*D.lpgdp;

* normalize variables and shorten names;
gen y = rgdp/`ynorm';
gen g = rgov/`ynorm';
 

local shock newsy;
gen neg=newsy<0;


local state neg; /* "state" is dummy equal to 1 if G shock is negative and 0 if G shock is positive */


*******************************************************************************;
** CUMULATIVE VARIABLES;
*******************************************************************************;

gen cumuly = 0;
gen cumulg = 0;
 
forvalues i = 0/`hor' {;

   gen f`i'cumuly = F`i'.y + cumuly;
   gen f`i'cumulg = F`i'.g + cumulg;
     
   replace cumuly = f`i'cumuly;
   replace cumulg = f`i'cumulg;
   
};


*******************************************************************************;
**  INTERACTION OF SHOCKS WITH STATE;
*******************************************************************************;

 foreach var in newsy { ;
 
   gen rec0`var' = `var'*`state'; /*shock is negative*/
   gen exp0`var' = `var'*(1-`state');  /*shock is positive*/
 
 };

*******************************************************************************;
** CREATE LISTS;
*******************************************************************************;

 gen h = t - 1;  /* h is the horizon for the irfs */
 global trendlist t t2 t3 t4;
		 
		 
forvalues i = 1/`p' {; 

  foreach var in newsy y g taxy debty infl{;

	gen `var'`i' = L`i'.`var';
 
  };
};


global newsyxlist L(1/`p').newsy L(1/`p').y L(1/`p').g L(1/`p').taxy $trendlist;
global newsylinshock newsy;

** INITIALIZE SUM OF EFFECTS TO 0 AND PARAMETERS SERIES TO MISSING;

gen sumliny = 0; gen sumling = 0;
gen sumexpy = 0; gen sumexpg = 0;
gen sumrecy = 0; gen sumrecg = 0;

foreach var in bylin byexp byrec bglin bgexp bgrec up90bylin up90byexp up90byrec up90bglin up90bgexp up90bgrec
  lo90bylin lo90byexp lo90byrec lo90bglin lo90bgexp lo90bgrec seylin seyexp seyrec seglin segexp segrec
  multlin multexp multrec {;
  
  quietly gen `var' = .;
  
}; 

 
*******************************************************************************;
** ESTIMATION OF IRFS
*******************************************************************************;


forvalues i = 0/`hor' {;

	* estimate linear effect;	
	ivreg2 F`i'.y $`shock'linshock $`shock'xlist, robust bw(auto);
	gen bylinh`i' = _b[$`shock'linshock];  
	gen seylinh`i' = _se[$`shock'linshock];
  
	ivreg2 F`i'.g $`shock'linshock $`shock'xlist, robust bw(auto);
	gen bglinh`i' = _b[$`shock'linshock];  
	gen seglinh`i' = _se[$`shock'linshock]; 

	* estimate asymmetric effects;		
	ivreg2 F`i'.y exp0`shock' rec0`shock' $`shock'xlist, robust bw(auto);	
	gen byexph`i' = _b[exp0`shock'];
	gen byrech`i' = _b[rec0`shock'];
    gen seyexph`i' = _se[exp0`shock'];
	gen seyrech`i' = _se[rec0`shock']; 

	ivreg2 F`i'.g exp0`shock' rec0`shock' $`shock'xlist, robust bw(auto);	
	gen bgexph`i' = _b[exp0`shock']; 
	gen bgrech`i' = _b[rec0`shock'];
	gen segexph`i' = _se[exp0`shock'];
	gen segrech`i' = _se[rec0`shock'];
  
  
	replace sumliny = bylinh`i' + sumliny;
	replace sumling = bglinh`i' + sumling;
  
	replace sumexpy = byexph`i' + sumexpy;
	replace sumexpg = bgexph`i' + sumexpg;
  
	replace sumrecy = byrech`i' + sumrecy;
	replace sumrecg = bgrech`i' + sumrecg;
  
	gen multlinh`i' = sumliny/sumling;
	gen multexph`i' = sumexpy/sumexpg;
	gen multrech`i' = sumrecy/sumrecg;
  
  foreach var in bylin byexp byrec bglin bgexp bgrec multlin multexp multrec seyexp seyrec segexp segrec {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	
  };
  
  foreach var in ylin glin yexp gexp yrec grec {;
  
    quietly replace up90b`var' = b`var'h`i' + 1.654*se`var'h`i' if h==`i';
	quietly replace lo90b`var' = b`var'h`i' - 1.654*se`var'h`i' if h==`i';
	
  };

};


display as text "MULTIPLIERS:  2 STEP";

rename multlin multlin2;
rename multexp multexp2;
rename multrec multrec2;

*save IRFs;
*outsheet h bglin up90bglin lo90bglin  bylin up90bylin lo90bylin bgexp up90bgexp lo90bgexp  byexp up90byexp lo90byexp bgrec up90bgrec lo90bgrec  byrec up90byrec lo90byrec  using RZ_irfs.csv if h<=`hor', comma replace;


label var bglin "Gov, linear model";
label var bylin "GDP, linear model";
label var bgexp "GOV, expansion";
label var byexp "GDP, expansion";
label var bgrec "Gov, recession";
label var byrec "GDP, recession";

drop mult???h*;

drop sey*;




foreach var in multlin1 multexp1 multrec1 Fkplin Fkpexp Fkprec seylin seyexp seyrec ptestdiff Fdifflin Fdiffexp Fdiffrec{;
  
  quietly gen `var' = .;
  
}; 

*******************************************************************************;
** ESTIMATION OF CUMULATIVE MULTIPLIER;
*******************************************************************************;

forvalues i = 0/`hor' {; 

  *estimate linear effect;
  ivreg2 f`i'cumuly (f`i'cumulg = $`shock'linshock) $`shock'xlist, robust bw(auto);
  gen Fkplinh`i'= e(widstat); /* Kleibergen-Paap rk Wald F statistic*/
  gen Fdifflinh`i'= Fkplinh`i'- 23.1085; 
  gen multlinh`i' = _b[f`i'cumulg];
  gen seylinh`i' = _se[f`i'cumulg]; /* HAC robust standard error*/

  *estimate asymmetric effects;
  ivreg2 f`i'cumuly (f`i'cumulg = exp0`shock') $`shock'xlist, robust bw(auto);
  gen Fkpexph`i'= e(widstat);
  gen Fdiffexph`i'= Fkpexph`i'- 23.1085; 
  gen multexph`i' = _b[f`i'cumulg];
  gen seyexph`i' = _se[f`i'cumulg];
  
  ivreg2 f`i'cumuly (f`i'cumulg = rec0`shock') $`shock'xlist, robust bw(auto);
  gen Fkprech`i'= e(widstat);  
  gen Fdiffrech`i'= Fkprech`i'- 23.1085; 
  gen multrech`i' = _b[f`i'cumulg];
  gen seyrech`i' = _se[f`i'cumulg];
  
  ivreg2 f`i'cumuly (f`i'cumulg f`i'cumulg = exp0`shock' rec0`shock') $`shock'xlist, robust bw(auto);
  test f`i'cumulg=f`i'cumulg;	

  gen ptestdiffh`i' = r(p);
  
  	
 foreach var in multlin multexp multrec {;
  
    quietly replace `var'1 = `var'h`i' if h==`i';
	
  };
  
  foreach var in seylin seyexp seyrec ptestdiff Fkplin Fkpexp Fkprec  {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	
  };
  
 foreach var in  Fdifflin Fdiffexp Fdiffrec {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	quietly replace `var' = 30 if `var'>30;
	
  };
};


display as text "First stage F-statistic (Kleibergen-Paap rk Wald F statistic): Linear, Expansion, Recession";

list h Fkplin Fkpexp Fkprec if h<=`hor';
outsheet h Fkplin Fkpexp Fkprec using junk.csv if h<=`hor', comma replace ;
outsheet h Fdifflin Fdiffexp Fdiffrec using junkfdiff.csv if h<=`hor', comma replace ;

display as text "Multipliers from 1 step and 3 step approaches: Linear, Expansion, Recession";

list h multlin1 multlin2 multexp1 multexp2 multrec1 multrec2 if h<=`hor';
outsheet h multlin1 multexp1 multrec1 using junk1step.csv if h<=`hor', comma replace ;

display as text "Multipliers and corresponding standard errors: Linear, Expansion, Recession";

list h multlin1 seylin multexp1 seyexp multrec1 seyrec ptestdiff if h<=`hor';
outsheet h multlin1 seylin multexp1 seyexp multrec1 seyrec ptestdiff using junkmultse.csv if h<=`hor', comma replace ;
outsheet h multlin1 seylin multexp1 seyexp multrec1 seyrec ptestdiff using RZ_M.csv if h<=`hor', comma replace ;



capture log close;

 
