    //			FERTILITY AND CONFLICT
    //			Working on DHS
    //			Duration Model
    //			First Version 	: 130710
    //          Last modify 	: 131020

    ///////////////////////////////////////////////////////////////////////////


clear all
set more off
set mem 1640m
set mat 800
set line 150
    // 	folder, initial set up
cap cd "/data/uctpjam/Stata/" // HPC
cap cd "E:\PhD Papers\Fertility decision during conflict\DataAnalysis\" // at uni
cap cd "/Volumes/MY WORK/PhD Papers/Fertility decision during conflict/DataAnalysis\"  //at mac
    ///////////////////////////////////////////////////////////////////////////



/* 
    Description:
		This file runs the duration analysis over different specifications.

		I got it from the file 130710DurationAnalysisShell.do, for previous
		description and changes see the description of that file.
   
	:::: update 131005 -> Insert the macros to be able to run the 
						conditional regressions 
 */

    
// Macros
////

loc ageOn   $AGE
loc fNum    $FILENUMBER
loc fName   $FILENAME
loc lA      $LOWAGE
loc uA      $UPPERAGE
loc poli    $POLI
loc set     $SET
loc dep     $DEPVAR
loc reps    $BSREPS
loc prefix 	$PREFIX 
loc esp		$ESP
loc cond 	$COND
loc intVar 	$INTER

loc inF "DataFiles/Working/"
loc inF1 "DataFiles/Working/MapsData/120922/"
loc doFld "DoFiles/JobMarketRuns/"

////////
cap log close
log using LogFiles/`prefix'Drt`dep'`lA'`uA'on`fNum'`ageOn'Int`intVar'InstSet`set'poli`poli'esp`esp'BS`reps'.log, replace

// 1. Declare the paramenters of this run

gl  fileOn  `fName'                 // Violence file
gl  ageCat  `ageOn'                 // Age Cathegory

gl  lAge    `lA'                    // Age lower bound
gl  uAge    `uA'                    // Age upper bound

gl  depVar  `dep'                   // outcome variable (has to be a diration variable)
gl  setOn   `set'                   // DD Network - instrument
gl  POLI    `poli'                  // Control Function polinomiom degree   
gl  BSR     `reps'                  // Bootstrap repetitions

gl control ///
    "i.Traf mCoca lPop uribe chavez mpioFamilias dptoGDPpk lfCens DHS_* yb_* urban etnia_* lPop wealth_* v136 *Dead y_* mpio_*"    
            


//  2. Prepare the data
    
    // Open durtation data
    u `inF'130628_DurationRegVars, clear
        keep if year>= 1990 & year <= 2009      // Availabilty of pUS pEU
        drop if DHS < 2000                      // I have not wealth before 2000 DHS
        drop if mig`lA' == 1                  // Dropping migrants after low age.


    // Municipality- year key   
    g code = codemun + "-" + string(year)
    la var code "Municipalty and year code"

    // Merge violence data
    so code
    mer m:1 code using `inF'$fileOn
        tab _merge
        keep if _merge == 3
        drop _merge


    // Data Set up
    qui do `doFld'130620DurationRenameHomicides 	// Rename the homicide variables
    qui do `doFld'130620DurationMergeIndVars 		// Women control vars
    qui do `doFld'130710DurationMergeMpioVars		// Municipality data
    qui do `doFld'130620DurationJenkins2005   		// Jenkins 2005 (Discrete time duration model)
    destring codemun, g(mun)					// I need this to run the Bootstrap    
        replace $depVar = . if $depVar == 2 | $depVar == -99  
        drop if $depVar == .    // Dropping the observations after first pregnancy or outside sample

    qui do `doFld'130710DurationCreateInstrument   // Merge and create the vars for the instrument


// 3. Run the estimations
	g uno = 1 		// if I use uno as conditional is the same of unconditional.
	foreach var in homRateTotal0 {
		gl VIOVAR `var'	
		do `doFld'DurationRegs_`esp'		// Running the model depending on the specification
											// The description in each regression do file.	
	}
	
   
    // drop instrument variables 
    cap drop firstVar - lastVar 
        
*/    
cap log close
exit, clear
//////
