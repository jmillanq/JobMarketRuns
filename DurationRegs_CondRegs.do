
/*
	This is a support file for the DurationAnalysis set.
	First Version	:	131006
	Last Modify		:	131006

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
    // Labels
    loc lab50 "All Ages" 
    loc lab51 "<=4 y.o."
    loc lab52 "5 to 14 y.o."
    loc lab53 "15 to 44 y.o."
    loc lab54 "45 to 64 y.o."
    loc lab55 ">=65 y.o."
    
    loc ageCond "(age >= `lAgeOn')*(age <= `uAgeOn') == 1" 
    loc yearCond "(year >= `yMin')*(year <= `yMax') == 1"

    loc dur1 lnj                      
                                
    loc esttabOpt "star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) se(%9.3f) nonum"
                                                                      
    
    // Instrument lists
    des clst_*xln_pUS clst_*xln_pEU, varlist
        loc listIns1 = r(varlist)
    des PI_*_US PI_*_EU, varlist
        loc listIns2 = r(varlist)
    loc listIns "`listIns1' `listIns2'"

/// Drop observations outside my sample (see the notes in the description) 
preserve


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
    loc F_Instr = r(F)
	loc r2FS = e(r2)
	estadd scalar F_Instr = `F_Instr'

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
		estadd scalar F_Inst = `F_Instr'

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
	LatexFiles/`prefix'Drt`dep'`lA'`uA'on`vioVar'DT`set'poli`cfpoli'Cond`condVar'.`out', r ///
    keep(`listIns') s(F_Instr r2, fmt(2 2))  `esttabOpt' ti("First Stage") nonotes
         
	esttab SS using ///
	LatexFiles/`prefix'Drt`dep'`lA'`uA'on`vioVar'DT`set'poli`cfpoli'Cond`condVar'.`out' ///
	, a `esttabOpt'  keep(zHom ControlF*) ///
	s(ppDiff ppSE F_Instr CTest r2_p M Ind N, ///
	fmt(2 2 2 2 2 0 0 0)) ///
	addnotes( ///
		"Duration var : `labelY'" ///
		"Age range : `lAgeOn' - `uAgeOn'" ///
		"Period : `yMin' - `yMax'" ///
		"Violence var : `labelVio'" ///
		"DT Set : `set'" ///
		"Control Funtion polynomial : `cfpoli'") ti(`tableTitle')
}

// Tab the even and the condition
tab `yVar' `condVar', missing

/////
