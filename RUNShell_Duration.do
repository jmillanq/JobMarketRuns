 // Running shell, duration analysis

// 	folder, initial set up
cap cd "/data/uctpjam/Stata/" // HPC
cap cd "E:\PhD Papers\Fertility decision during conflict\DataAnalysis\" // at uni
cap cd "/Volumes/MY WORK/PhD Papers/Fertility decision during conflict/DataAnalysis\"  //at mac


//Macros

 gl AGE 		0
 gl FILENUMBER	5 
 gl FILENAME	homRatesEV_1985_2010
 gl LOWAGE		13
 gl UPPERAGE	19
 gl POLI		1
 gl SET			10
 gl DEPVAR		preg2
 gl BSREPS		2
 gl PREFIX		131020
 gl ESP			ContVarInter 		
 gl COND		uno 	// Only for ESP = CondRegs (uno = unconditional)
 gl INTER		v136  	// Only for ESP = ContVarInter	
// Run the do file§
// do DoFiles/131019nBirthsTest
do DoFiles/JobMarketRuns/DurationAnalysis_shell

exit, clear
