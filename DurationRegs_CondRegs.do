
/* 
	This is a support file for the DurationAnalysis set.
	First Version	:	131006
	Last Modify		:	131019

		In this Do file I run the logistic discrete duration model condition
		the outcome variable to another time related event.

			e.g. 	I can run in this file the duration model of pregnancy conditional 
					on being married.

*/
 
// Import macros
	loc lAgeOn $lAge
	loc uAgeOn $uAge
	loc lowBound = `lAgeOn' - 1
	loc set $setOn
	loc yVar $depVar
	loc labelY : var label `yVar'
	loc vioVar $VIOVAR
	loc labelVio : var label `vioVar'
	loc cfpoli $POLI
	loc condVar $COND
	loc reps $BSR
	loc prefix $PREFIX

// Titles 
	loc tableTitle "`yVar' from `lAgeOn' to `uAgeOn' on `vioVar' - Violence and parents mortality"

// Period of analysis
    su year if `vioVar'  != .
    local yMin  = r(min)
    local yMax  = r(max)
    di "(`yMin' - `yMax')"

    // Min Age
    cap drop minAge
    egen minAge = min(age) if `vioVar' !=., by(newcaseid)
    table year if (age >= `lAgeOn')*(age <= `uAgeOn') == 1, ///
    c(min minAge max minAge)

    // Censoring
    so newcaseid
    foreach i in min max {
        cap drop `i'Age
        egen `i'Age = `i'(age) if `vioVar' != ., by(newcaseid)
    }
    cap drop lfCens
    cap drop rgCens
    g lfCens = (age == minAge)*(minAge > `lAgeOn')
    g rgCens = (age == maxAge)*(`yVar'==0)
   
    loc ageCond "(age >= `lAgeOn')*(age <= `uAgeOn') == 1" 
    loc yearCond "(year >= `yMin')*(year <= `yMax') == 1"

    loc dur1 lnj                      
                                
    loc esttabOpt "star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) se(%9.3f) nonum"
    loc fZones "Pacific Atlantic VenSouth VenNorth"                                                                  
    
    // Instrument lists
    des clst_*xln_pUS clst_*xln_pEU, varlist
    loc listIns1 = r(varlist)
	//clst_Pacificxln_pUS clst_Pacificxln_pEU clst_Atlanticxln_pUS clst_Atlanticxln_pEU clst_VenSouthxln_pUS clst_VenSouthxln_pEU clst_VenNorthxln_pUS clst_VenNorthxln_pEU
	//loc listIns1 "clst_Pacificxln_pUS  clst_Atlanticxln_pUS  clst_VenSouthxln_pEU clst_VenNorthxln_pUS clst_VenNorthxln_pEU"
	// loc listIns2 "CombiPacAtlln_pUS  CombiAtlVNln_pUS CombiVNVSln_pUS CombiVNVSln_pEU"
	
	des PI_*_US PI_*_EU, varlist
    loc listIns2 = r(varlist)
	//loc listIns2  "PI_Pacific_US PI_Atlantic_US PI_VenNorth_US PI_VenSouth_EU "
    loc listIns "`listIns1' `listIns2'"

// Prd * price * cluster

foreach b in `fZones' {
	foreach p in US EU {
		g mCoca_`b'_`p' = mCocaXln_p`p'*c_`b'
	}
}
g Traf = c_Pacific + c_Atlantic + c_VenSouth + c_VenNorth
g uribe = (year >= 2002)*(year <= 2010)
g chavez = (year >= 1999)*(year <= 2013)
/*
g PacAtl = c_Pacific*c_Atlantic
g PacVN = c_Pacific*c_VenNorth
g AtlVN = c_Atlantic*c_VenNorth
g AtlVS = c_Atlantic*c_VenSouth
g VNVS = c_VenNorth*c_VenSouth

foreach c in PacAtl PacVN AtlVN AtlVS VNVS {
	foreach p in US EU {
		g Combi`c'ln_p`p' = `c'*ln_p`p'
	}
}

foreach z in `fZones' {
	
}
*/
/// Drop observations outside my sample (see the notes in the description) 
preserve
keep if `ageCond' & `yearCond' & `condVar' >= 1

///////////// REGRESSIONS
cap est clear
// Standarized the vio variable
	cap drop zHom*
	cap drop z_`vioVar'
	sum `vioVar' if `ageCond' & `yearCond'
	loc varMean = r(mean)
	loc varSD = r(sd)
    zscore `vioVar' if `ageCond' & `yearCond'
    ren z_`vioVar' zHom
        
		
// For the bootstrap
	qui do DoFiles/130710ivDurLogit.do

	// 
    bootstrap diff = r(ppDiff), ///
	    reps(`reps') cluster(newcaseid)   : ///
            ivDurLogit $control if  `ageCond' & `yearCond' & `condVar' >= 1, ///
                minage(`lAgeOn') maxage(`uAgeOn') cfp(`cfpoli') z(`listIns') y(`yVar') w(zHom) 
                    loc bDiff = _b[diff]
                    loc seDiff = _se[diff]
 */


// First Stage - Vio = Drug Traffic
 	cap drop errors
    cap drop ControlF*
	eststo FS: ///
        reg zHom `listIns' `dur1' $control ///
    	[pw = awtvar] if `ageCond' & `yearCond' & `condVar' >= 1, cluster(codemun) 
    loc Mpios = e(N_clust)
	test `listIns'
    loc F_Inst = r(F)
	loc r2FS = e(r2)
	estadd scalar F_Inst = `F_Inst'

    predict errors, res
        forv p = 1/`cfpoli' {
			g ControlF`p' = errors^`p'
        }  
            
    des ControlF*, varlist
	loc ControlList = r(varlist)

// Second Stage
	eststo SS: ///
    	logit `yVar' zHom  ControlF* `dur1' $control ///
        [pw = awtvar] if  `ageCond' & `yearCond' & `condVar' >= 1, cluster(codemun) iter(100) 
    
	
	// N individuals
		cap drop id
    	bys newcaseid: g id = (_n == 1) if  `ageCond' & `yearCond' & `condVar' >= 1 & ControlF1 != .
        	replace id = sum(id) if  `ageCond' & `yearCond' & `condVar' >= 1 & ControlF1 != . 
        	sum id
		estadd scalar Ind = r(max)
		estadd scalar M = `Mpios'

	// First stage stats
		estadd scalar F_Inst = `F_Inst'

	// bootstrap results
		estadd scalar ppDiff = `bDiff'
        estadd scalar ppSE = `seDiff'
	
	// Test over Control Function  
        test `ControlList'
        loc ControlTest = r(chi2)
        estadd scalar CTest = `ControlTest' 
                    


// Publish
foreach out in csv tex {

    esttab FS using ///
	LatexFiles/`prefix'Drt`yVar'`lAgeOn'`uAgeOn'on`vioVar'DT`set'poli`cfpoli'Cond`condVar'.`out', r ///
    keep(`listIns') s(F_Inst r2, fmt(2 2))  `esttabOpt' ti("First Stage") nonotes
         
	esttab SS using ///
	LatexFiles/`prefix'Drt`yVar'`lAgeOn'`uAgeOn'on`vioVar'DT`set'poli`cfpoli'Cond`condVar'.`out' ///
	, a `esttabOpt'  keep(zHom ControlF* lnj _cons) ///
	s(ppDiff ppSE F_Inst CTest r2_p M Ind N, ///
	fmt(2 2 2 2 2 0 0 0)) ///
	addnotes( ///
		"Duration var : `labelY'" ///
		"Conditional on : `condVar'" ///
		"Age range : `lAgeOn' - `uAgeOn'" ///
		"Period : `yMin' - `yMax'" ///
		"Violence var : `labelVio'" ///
		"DT Set : `set'" ///
		"Control Funtion polynomial : `cfpoli'" ///
		"Bootstrap Repetitions : `reps'" ///
		"Control variables include : $control") ///
 		ti(`tableTitle') 
}

// Tab the event and the condition
tab `yVar' `condVar', missing
restore
/////
