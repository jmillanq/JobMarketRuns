/*
	This is a support file for the DurationAnalysis set.
	First Version	:	131001
	Last Modify		:	131001

		In this Do file I run the logistic discrete duration model using 
		the interaction of the violence var with the dummies of 
		mother is dead (motherDead), father is dead (fatherDead) and both 
		parents are dead (bothDead).

		The exit file report table will include the F statistic on 
		vioVar + vioVar*'y'Dead


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
    
// Gen interaction parentsDead*zHom
	cap drop zHom_D*
    foreach x in mother father parents {
        g zHom_D`x' = zHom*`x'Dead        
	}
loc esp1 "zHom_D*"
loc esp2 "zHom_Dmother"
loc esp3 "zHom_Dfather"
loc esp4 "zHom_Dparents"


// First Stage - Vio = Drug Traffic
 	cap drop errors
    cap drop ControlF*
	eststo FS: ///
        reg zHom `listIns' `dur1' $control ///
    	[pw = awtvar] if `ageCond' & `yearCond', cluster(codemun) 
    loc Mpios = e(N_clust)
	test `listIns'
    loc F_Instr = r(F)
	loc r2FS = e(r2)
    predict errors, res
        forv p = 1/`cfpoli' {
			g ControlF`p' = errors^`p'
        }  
            
    des ControlF*, varlist
	loc ControlList = r(varlist)


forv mod = 1/4 {    
	eststo SS_`mod': ///
    	logit `yVar' zHom `esp`mod'' ControlF* `dur1' $control ///
        [pw = awtvar] if  `ageCond' & `yearCond' , cluster(codemun) iter(100) 
    
	// Test over the sum with the interaction
		foreach int in mother father parents {
			cap n test zHom + zHom_D`int' = 0
			cap n estadd scalar tot`int' = _b[zHom] + _b[zHom_D`int']
			cap n estadd scalar Chi2`int' = r(chi2)
		}
	
	// N individuals
		cap drop id
    	bys newcaseid: g id = (_n == 1) if  `ageCond' & `yearCond' & ControlF1 != .
        	replace id = sum(id) if  `ageCond' & `yearCond' & ControlF1 != . 
        	sum id
		estadd scalar Ind = r(max)
		estadd scalar M = `Mpios'
		
	// Test over Control Function  
        test `ControlList'
        loc ControlTest = r(chi2)
        estadd scalar CTest = `ControlTest' 
                    
}


// Publish
foreach out in csv tex {
	esttab SS_* using ///
	LatexFiles/`prefix'Drt`dep'`lA'`uA'on`vioVar'DT`set'poli`cfpoli'ParentsMort.`out' ///
	, r `esttabOpt'  keep(zHom zHom_D* *Dead ControlF*) ///
	s(totmother Chi2mother totfather Chi2father totparents Chi2parents CTest r2_p M Ind N, ///
	fmt(2 2 2 2 2 2 2 2 0 0 0)) ///
	addnotes( ///
		"Duration var : `labelY'" ///
		"Age range : `lAgeOn' - `uAgeOn'" ///
		"Period : `yMin' - `yMax'" ///
		"Violence var : `labelVio'" ///
		"DT Set : `set'" ///
		"Control Funtion polynomial : `cfpoli'") ti(`tableTitle')
}

/////
