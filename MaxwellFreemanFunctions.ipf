#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <KBColorizeTraces>

//Version MPPF2025-09-22

//**************************************************************************************************
Function SetupAverageSpectra([ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts, NumParticles])
	// Created 01/22/20 by JAH and MPPF
	// Purpose: To perform an element-wise average of fluorescence spectra waves from CropAnalysis().
	
	Variable ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts, NumParticles
	Variable ParticleCounter = 0
	
	// Initialize all default input parameters.
	if (ParamIsDefault(ThresholdVar))
		ThresholdVar = 0
	endif
	if (ParamIsDefault(LowerBound))
		LowerBound = 200
	endif
	if (ParamIsDefault(UpperBound))
		UpperBound = 900
	endif
	if (ParamIsDefault(DependentVarStepSize))
		DependentVarStepSize = 0.5
	endif
	if (ParamIsDefault(NumFoundationPts))
		NumFoundationPts = 1401
	endif
	if (ParamIsDefault(NumParticles))
		NumParticles = 100
	endif
	

	// Make FoundationWv as the skeleton onto which to perform regridding of the dependent variable.
	linspace("FoundationWv", LowerBound, UpperBound, Spacing = DependentVarStepSize)
	Wave FoundationWv
	
	Make/O/N=(numpnts(FoundationWv)) Average280Wv, Average350Wv, Average405Wv, Average450Wv = 0
	Make/O/N=1 MeanPeakPosition280, MeanPeakPosition350, MeanPeakPosition405, MeanPeakPosition450 = 0
	Make/O/N=(10,4) MaxIntensities = 0
	
	// "While" loop regrids and filters null values for each particle spectrum supplied.
	Variable i = 0
	do
		i += 1
		String WavelengthsName = "p" + num2str(i) + "xRef"
		
		// Once all particles have been dealt with, break the "while" loop.
		if (!WaveExists($WavelengthsName))
			break
		endif
		ParticleCounter += 1
		
		
		String Intensity280Name = "p" + num2str(i) + "y280"
		String Intensity350Name = "p" + num2str(i) + "y350"
		String Intensity405Name = "p" + num2str(i) + "y405"
		String Intensity450Name = "p" + num2str(i) + "y450"

		//////////////////////////////////////////////////////////////////////////////////////////////
		// Core functionality in box is applied to all four excitations //////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////

		// For 280 nm
		
		// Following command forces dependent variable of interest onto a separate grid with defined
		// start, end, and number of points.
		Regrid_Param2($WavelengthsName,$Intensity280Name,FoundationWv,FoundationWv)
		
		
		// RemoveValue clears all y-values equal to last point in the regridded y-wave (named ThresholdVar).
		Wave ywave_regrid
		ThresholdVar = ywave_regrid[numpnts(ywave_regrid)-1]
		ReplaceValue(ywave_regrid, ThresholdVar, NaN)
		
		// Save regridded intensity wave to some uniquely named data structure (could try appending to a
		// matrix?).
		Average280Wv += ywave_regrid
		String DummyStr = "regridded_p" + num2str(i) + "y280"
		Make/O/N=(numpnts(ywave_regrid)) $DummyStr = ywave_regrid
		
		// Use CurveFit to get peak positions
		CurveFit/Q gauss $DummyStr((400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize) /D	// NOTE: the values in parentheses are indices to fit between, not wavelengths.
		String FitString = "fit_"  + DummyStr
		Wavestats/Q $FitString
		MeanPeakPosition280[numpnts(MeanPeakPosition280)-1] = V_maxloc*DependentVarStepSize + LowerBound
		InsertPoints numpnts(MeanPeakPosition280), 1, MeanPeakPosition280
		MaxIntensities[i-1][0] = V_max
		//////////////////////////////////////////////////////////////////////////////////////////////
		
		// For 350 nm
		Regrid_Param2($WavelengthsName,$Intensity350Name,FoundationWv,FoundationWv)
		Wave ywave_regrid
		ThresholdVar = ywave_regrid[numpnts(ywave_regrid)-1]
		ReplaceValue(ywave_regrid, ThresholdVar, NaN)
		Average350Wv += ywave_regrid
		DummyStr = "regridded_p" + num2str(i) + "y350"
		Make/O/N=(numpnts(ywave_regrid)) $DummyStr = ywave_regrid
		CurveFit/Q gauss $DummyStr((400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize) /D
		FitString = "fit_"  + DummyStr
		Wavestats/Q $FitString
		MeanPeakPosition350[numpnts(MeanPeakPosition350)-1] = V_maxloc*DependentVarStepSize + LowerBound
		InsertPoints numpnts(MeanPeakPosition350), 1, MeanPeakPosition350
		MaxIntensities[i-1][1] = V_max

		
		// For 405 nm
		Regrid_Param2($WavelengthsName,$Intensity405Name,FoundationWv,FoundationWv)
		ThresholdVar = ywave_regrid[numpnts(ywave_regrid)-1]
		ReplaceValue(ywave_regrid, ThresholdVar, NaN)
		Average405Wv += ywave_regrid
		DummyStr = "regridded_p" + num2str(i) + "y405"
		Make/O/N=(numpnts(ywave_regrid)) $DummyStr = ywave_regrid
		CurveFit/Q gauss $DummyStr((400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize) /D
		FitString = "fit_"  + DummyStr
		Wavestats/Q $FitString
		MeanPeakPosition405[numpnts(MeanPeakPosition405)-1] = V_maxloc*DependentVarStepSize + LowerBound
		InsertPoints numpnts(MeanPeakPosition405), 1, MeanPeakPosition405
		MaxIntensities[i-1][2] = V_max
		
		// For 450 nm
		Regrid_Param2($WavelengthsName,$Intensity450Name,FoundationWv,FoundationWv)
		ThresholdVar = ywave_regrid[numpnts(ywave_regrid)-1]
		ReplaceValue(ywave_regrid, ThresholdVar, NaN)
		Average450Wv += ywave_regrid
		DummyStr = "regridded_p" + num2str(i) + "y450"
		Make/O/N=(numpnts(ywave_regrid)) $DummyStr = ywave_regrid
		CurveFit/Q gauss $DummyStr((400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize) /D
		FitString = "fit_"  + DummyStr
		Wavestats/Q $FitString
		MeanPeakPosition450[numpnts(MeanPeakPosition450)-1] = V_maxloc*DependentVarStepSize + LowerBound
		InsertPoints numpnts(MeanPeakPosition450), 1, MeanPeakPosition450
		MaxIntensities[i-1][3] = V_max


	while(1)
	
	// Create waves of mean values for each excitation.
	Average280Wv /= ParticleCounter
	Average350Wv /= ParticleCounter
	Average405Wv /= ParticleCounter
	Average450Wv /= ParticleCounter
	
	
	Make/O/N=(numpnts(FoundationWv)) StdDevs280Wv, StdDevs350Wv, StdDevs405Wv, StdDevs450Wv = 0
	
	// "While" loop calculates quadratic sum (numerator)for standard deviation for each point 
	// in the average waves.
	Variable j = 0
	do
		j += 1
		
		WavelengthsName = "p" + num2str(j) + "xRef"
		
		// Once all particles have been dealt with, break the "while" loop.
		if (!WaveExists($WavelengthsName))
			break
		endif
		
		// Start std. dev. waves by calculating sum of squared differences between value and mean.
		DummyStr = "regridded_p" + num2str(j) + "y280"
		Wave DummyWv = $(DummyStr)
		StdDevs280Wv += (DummyWV - Average280WV)^2
		DummyStr = "regridded_p" + num2str(j) + "y350"
		Wave DummyWv = $(DummyStr)
		StdDevs350Wv += (DummyWv - Average350WV)^2
		DummyStr = "regridded_p" + num2str(j) + "y405"
		Wave DummyWv = $(DummyStr)
		StdDevs405Wv += (DummyWv - Average405WV)^2
		DummyStr = "regridded_p" + num2str(j) + "y450"
		Wave DummyWv = $(DummyStr)
		StdDevs450Wv += (DummyWv - Average450WV)^2
		
	while(1)
	
	// Compute full standard deviations waves.
	StdDevs280Wv = sqrt(StdDevs280Wv/(ParticleCounter-1))
	StdDevs350Wv = sqrt(StdDevs350Wv/(ParticleCounter-1))
	StdDevs405Wv = sqrt(StdDevs405Wv/(ParticleCounter-1))
	StdDevs450Wv = sqrt(StdDevs450Wv/(ParticleCounter-1))
	
	// Create average and std. dev. waves normalized to max value in the respective average wave.
	Wavestats/Q Average280Wv
	Make/O/N=(numpnts(Average280Wv)) NormAverage280Wv = Average280Wv/V_max
	Make/O/N=(numpnts(StdDevs280Wv)) NormStdDevs280Wv = StdDevs280Wv/V_max
	Wavestats/Q Average350Wv
	Make/O/N=(numpnts(Average350Wv)) NormAverage350Wv = Average350Wv/V_max
	Make/O/N=(numpnts(StdDevs350Wv)) NormStdDevs350Wv = StdDevs350Wv/V_max
	Wavestats/Q Average405Wv
	Make/O/N=(numpnts(Average405Wv)) NormAverage405Wv = Average405Wv/V_max
	Make/O/N=(numpnts(StdDevs405Wv)) NormStdDevs405Wv = StdDevs405Wv/V_max
	Wavestats/Q Average450Wv
	Make/O/N=(numpnts(Average450Wv)) NormAverage450Wv = Average450Wv/V_max
	Make/O/N=(numpnts(StdDevs450Wv)) NormStdDevs450Wv = StdDevs450Wv/V_max
	
	// Report mean peak positions and standard deviations of peak position.
	DeletePoints (numpnts(MeanPeakPosition280)-1), 1, MeanPeakPosition280	// Removes extra point at end.
	DeletePoints (numpnts(MeanPeakPosition350)-1), 1, MeanPeakPosition350	// Removes extra point at end.
	DeletePoints (numpnts(MeanPeakPosition405)-1), 1, MeanPeakPosition405	// Removes extra point at end.
	DeletePoints (numpnts(MeanPeakPosition450)-1), 1, MeanPeakPosition450	// Removes extra point at end.
	Wavestats/Q MeanPeakPosition280
	Variable/G SpeciesMeanPeakPosition280 = V_avg
	Variable/G SpeciesPeakPositionStdDev280 = V_sdev
	Wavestats/Q MeanPeakPosition350
	Variable/G SpeciesMeanPeakPosition350 = V_avg
	Variable/G SpeciesPeakPositionStdDev350 = V_sdev
	Wavestats/Q MeanPeakPosition405
	Variable/G SpeciesMeanPeakPosition405 = V_avg
	Variable/G SpeciesPeakPositionStdDev405 = V_sdev
	Wavestats/Q MeanPeakPosition450
	Variable/G SpeciesMeanPeakPosition450 = V_avg
	Variable/G SpeciesPeakPositionStdDev450 = V_sdev
	
	// Kill regridded waves for cleaner data browser.
	Variable k
	for (k = 1; k <= ParticleCounter; k += 1)
		DummyStr = "regridded_p" + num2str(j) + "y280"
		KillWaves $DummyStr
		DummyStr = "regridded_p" + num2str(j) + "y350"
		KillWaves $DummyStr
		DummyStr = "regridded_p" + num2str(j) + "y405"
		KillWaves $DummyStr
		DummyStr = "regridded_p" + num2str(j) + "y450"
		KillWaves $DummyStr
	endfor
	
	// Make an averages graph to see results.
	Display/K=1 NormAverage280Wv vs FoundationWv
	AppendtoGraph NormAverage350Wv vs FoundationWv
	AppendtoGraph NormAverage405Wv vs FoundationWv
	AppendtoGraph NormAverage450Wv vs FoundationWv
	
	ModifyGraph rgb(NormAverage350Wv)=(0,0,65535),rgb(NormAverage405Wv)=(2,39321,1),rgb(NormAverage450Wv)=(26411,1,52428)
	Label left "Normalized Average Intensity (a.u)";DelayUpdate
	Label bottom "Wavelengths (nm)"
	Legend/C/N=text0/J/F=0/A=LT "\\s(NormAverage280Wv) 280 nm\r\\s(NormAverage350Wv) 350 nm\r\\s(NormAverage405Wv) 405 nm\r\\s(NormAverage450Wv) 450 nm"

End Function

//**************************************************************************************************

Function ReplaceValue(InputWv, Value, NewValue)
// Created 01/22/20 by MPPF
//	Updated 2025-09-22 MPPF to address +-inf.
// Purpose: to convert all entries equal to some value in a wave to a new value. 1 line to replace 50 lines if you understand the language... :D
	Wave InputWv
	Variable Value, NewValue
	InputWv = ( NumType(Value) != 0 ) ? ( NumType(InputWv) == 2 ? (NumType(InputWv) == 2 ? NewValue : InputWv) : (NumType(InputWv) == 1 ? NewValue : InputWv) ) : ( InputWv == Value ? NewValue : InputWv )

//	InputWv = ( NumType(Value) == 2 ) ? ( NumType(InputWv) == 2 ? NewValue : InputWv ) : ( InputWv == Value ? NewValue : InputWv )

	
//	InputWv = InputWv == Value ? NewValue : InputWv
//	Variable Dims = WaveDims(InputWv)
//	Variable i,j,k,l
////Remove a value from a 1D wave.
//	If(Dims == 1)
//		For (i = 0 ; i < numpnts(InputWv) ; i += 1)
//			If (InputWv[i] == Value)
//				InputWv[i] = NaN
//			Endif
//		Endfor
//	EndIf
//	
////Remove a value from a 2D wave.
//	If(Dims == 2)
//		For (j = 0; j < DimSize(InputWv, 1); j += 1)
//			For (i = 0; i < DimSize(InputWv, 0); i += 1)
//				If (InputWv[i][j] == Value)
//					InputWv[i][j] = NaN
//				EndIf
//			EndFor
//		EndFor
//	EndIf
//	
////Remove a value from a 3D wave.
//	If(Dims == 3)
//		For (k = 0; k < DimSize(InputWv, 2); k += 1)
//			For (j = 0; j < DimSize(InputWv, 1); j += 1)
//				For (i = 0; i < DimSize(InputWv, 0); i += 1)
//					If (InputWv[i][j][k] == Value)
//						InputWv[i][j][k] = NaN
//					EndIf
//				EndFor
//			EndFor
//		EndFor
//	EndIf
//	
////Remove a value from a 3D wave.
//	If(Dims == 4)
//		For (l = 0; l< DimSize(InputWv, 3); l += 1)
//			For (k = 0; k < DimSize(InputWv, 2); k += 1)
//				For (j = 0; j < DimSize(InputWv, 1); j += 1)
//					For (i = 0; i < DimSize(InputWv, 0); i += 1)
//						If (InputWv[i][j][k][l] == Value)
//							InputWv[i][j][k][l] = NaN
//						EndIf
//					EndFor
//				EndFor
//			EndFor
//		EndFor
//	EndIf
End Function

//**************************************************************************************************

Function linspace(WvName, LowerBound, UpperBound, [Spacing, NumberOfPoints])
// Created 01/22/20 by MPPF
// Purpose: to create a linearly spaced wave between two numbers (inclusive) either by number of 
// points or step size. Default is to use step size of 1.
	String WvName
	Variable LowerBound, UpperBound, Spacing, NumberOfPoints
	Variable i
	
	Make/O/N=1 TempWv = 0
	if(!ParamIsDefault(Spacing) && !ParamIsDefault(NumberOfPoints))
		print("Error: only one optional parameter can be used. Specifying number of points and spacing simultaneously can lead to inconsistencies.")
	
	elseif(ParamIsDefault(Spacing) && !ParamIsDefault(NumberOfPoints))
		Make/O/N=2 BoundsWv = {LowerBound, UpperBound}
		Interpolate2/T=1/N=(NumberOfPoints)/Y=TempWv BoundsWv
	
	elseif(ParamIsDefault(NumberOfPoints) && !ParamIsDefault(Spacing))
		TempWv[0] = LowerBound
		for (i = LowerBound+Spacing; i <= UpperBound; i += Spacing)
			InsertPoints numpnts(TempWv), 1, TempWv
			TempWv[numpnts(TempWv)-1] = i
		endfor
		
	else
		Spacing = 1
		TempWv[0] = LowerBound
		for (i = LowerBound+Spacing; i <= UpperBound; i += Spacing)
			InsertPoints numpnts(TempWv), 1, TempWv
			TempWv[numpnts(TempWv)-1] = i
		endfor
	endif
	
	Make/O/N=(numpnts(TempWv)) $WvName = TempWv
	
	Killwaves BoundsWv, TempWv
	
End Function

//**************************************************************************************************
//
//Function  Regrid_Param2(xwave1, ywave1, xwave2, ywave2)
//// Purpose: Maps ywave1 (gridded to xwave1) to a new grid (xwave 2). Not sure what can be done with
//// ywave2, for now, just copy xwave2 name there. Output is called "ywave_regrid."
//	Wave xwave1, ywave1, xwave2, ywave2
//       duplicate/o ywave2 ywave_regrid
//       ywave_regrid = ywave1 [binarysearchinterp(xwave1, xwave2[p])]
//End Function

//**************************************************************************************************

Function BTFMStats([ExcitationWavelength, ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts])
// Created 02/05/2020 by MPPF
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: to provide some stats on ensemble of spectra retrieved from CropAnalysis(). Key outputs are
// NormAveragedSpectra [data, excitations, averages + std. devs.], 
// AveragedSpectra (same dimensions as NormAveragedSpectra), 
// PeakStats[Mean peak position + std. dev. in peak position + max. intensity of average wave + avg. peak intensity, excitations], 
// Peaks [peak position + peak intensity, excitations, particles].
//
// Notes on BTFMStats() in Maxwell Freeman Notebook #1, p. 51-54.
// Replaces SetupAverageSpectra().
	
	Wave/T ExcitationWavelength
	Variable ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts
	
	Variable ParticleCounter = 0 // Number of particles analyzed.
	
	// Initialize all default input parameters.
	if (ParamIsDefault(ExcitationWavelength))
		Make/O/T ExcitationWavelength = {"280", "350", "405", "450"}
	endif
	if (ParamIsDefault(ThresholdVar))
		ThresholdVar = 0
	endif
	if (ParamIsDefault(LowerBound))
		LowerBound = 200
	endif
	if (ParamIsDefault(UpperBound))
		UpperBound = 900
	endif
	if (ParamIsDefault(DependentVarStepSize))
		DependentVarStepSize = 0.5
	endif
	if (ParamIsDefault(NumFoundationPts))
		NumFoundationPts = 1401
	endif

	

	// Make FoundationWv as the skeleton onto which to perform regridding of the dependent variable.
	linspace("FoundationWv", LowerBound, UpperBound, Spacing = DependentVarStepSize)
	Wave FoundationWv
	
	// Make initial 3-tensor. Rows are points in spectra, columns are 4 excitation sources (ascending by wavelength), 
	// layers are different particles. Will be ix4xN where i is number of points in FoundationWv and N is number of particles.
	Make/O/N=(numpnts(FoundationWv),4,1) RegriddedSpectra
	
	// Make 3-tensor for averages and standard deviations where each row is a point in a spectrum,
	// each column is an excitation source, and the first and second layers are the averages and
	// and standard deviation waves respectively.
	Make/O/N=(numpnts(FoundationWv),4,2) AveragedSpectra
	
	// Make matrix called PeakStats that contains excitations as columns; first row is mean peak position, second
	// row is std. dev. of peak position, third row is max intensity of average spectrum, fourth row is average intensity
	// of CurveFit peaks. Compare rows 3 and 4 within a column to see if noise exceeds peak (if row 3 > row 4, then spectrum
	// is more intense outside of the peak--i.e. noise dominates). Temporarily add extra row for the peak location of the 
	// average spectra.
	Make/O/N=(5,4) PeakStats
	
	// Make 3-tensor named Peaks where first row is peak position, second row is peak intensity, columns are excitations,
	// and layers are for separate particles. Will be 2x4xN where N is number of particles.
	Make/O/N=(2,4,1) Peaks
	
	// Regrid spectra for each particle and add to RegriddedSpectra
	Variable ParticleNumber = 1
	String WavelengthsName = "p" + num2str(ParticleNumber) + "xRef"
	do

		// For each new particle, add an extra layer in the appropriate tensors to accomodate the data for the new particle.
		if (ParticleNumber > 1)
			InsertPoints/M=2 DimSize(RegriddedSpectra, 2), 1, RegriddedSpectra
			InsertPoints/M=2 DimSize(Peaks, 2), 1, Peaks
		endif

		
		// "For" loop to perform same thing for each excitation source.
		variable ExcitationID
		for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
		
			String IntensityName = "p" + num2str(ParticleNumber) + "y" + ExcitationWavelength[ExcitationID]
			
			// Following command forces dependent variable of interest onto a separate grid with defined
			// start, end, and number of points.
			Regrid_Param2($WavelengthsName,$IntensityName,FoundationWv,FoundationWv)
		
			// RemoveValue clears all y-values equal to the mode (named ThresholdVar) in the regridded y-wave.
			Wave ywave_regrid
			StatsQuantiles/Q/iNaN ywave_regrid
			ThresholdVar = V_mode
			ReplaceValue(ywave_regrid, ThresholdVar, NaN)
		
			// Store data in RegriddedSpectra 3-tensor
			RegriddedSpectra[0,numpnts(ywave_regrid)-1][ExcitationID][ParticleNumber-1] = ywave_regrid[p][0][0]
			
			// Use CurveFit to get peak positions
			// Note: Wavelength (in nm) = DependentVarStepSize*index + LowerBound is akin to Fahrenheit = 1.8*Celsius + 32.
			CurveFit/Q gauss ywave_regrid[(400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize] /D	// NOTE: the values in parentheses are indices to fit between, not wavelengths; intention is to fit between 400 nm and 580 nm.
			Wave fit_ywave_regrid
			Wavestats/Q fit_ywave_regrid
			Peaks[0][ExcitationID][ParticleNumber-1] = V_maxloc*DependentVarStepSize + LowerBound
			Peaks[1][ExcitationID][ParticleNumber-1] = V_max
			
		endfor

			
		ParticleNumber += 1
		WavelengthsName = "p" + num2str(ParticleNumber) + "xRef"
		
	while(WaveExists($WavelengthsName))
	
	ParticleCounter = ParticleNumber-1
	
	
	
	// "For" loop to iteratively generate standard deviations and averages at each point in spectra.
	variable j
	for (j=0; j < DimSize(RegriddedSpectra, 0); j+=1)
	
		for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
			
			Make/FREE/N=(DimSize(RegriddedSpectra, 2)) temp = RegriddedSpectra[j][ExcitationID][p]
			WaveStats/Q temp
			AveragedSpectra[j][ExcitationID][0] = V_avg
			AveragedSpectra[j][ExcitationID][1] = V_sdev
			
		endfor


	endfor
	
	
	// Create normalized spectra by dividing averages and std. devs. by maximum of average spectrum of 
	// given excitation source.
	// Also use Peaks to fill PeakStats with necessary data.
	Duplicate/O AveragedSpectra, NormAveragedSpectra
	
	for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
		
		Make/FREE/N=(DimSize(Peaks, 2)) temp = Peaks[0][ExcitationID][p]
		WaveStats/Q temp
		PeakStats[0][ExcitationID] = V_avg		// Mean peak position for each excitation source.
		PeakStats[1][ExcitationID] = V_sdev		// Standard deviation in peak position for each average above.
	
		Make/FREE/N=(DimSize(AveragedSpectra, 0)) temp = AveragedSpectra[p][ExcitationID][0]
		WaveStats/Q temp
		NormAveragedSpectra[][ExcitationID][] /= V_max
		PeakStats[2][ExcitationID] = V_max		// Maximum intensity of average spectrum.
		PeakStats[4][ExcitationID] = V_maxloc*DependentVarStepSize + LowerBound // Temporarily record the peak location in the average spectrum.
		
		Make/FREE/N=(DimSize(Peaks, 2)) temp = Peaks[1][ExcitationID][p]
		WaveStats/Q temp
		PeakStats[3][ExcitationID] = V_avg		// Average peak intensity of all curve fits to non-averaged spectra.
		
	endfor
	
	Killwaves fit_ywave_regrid, ywave_regrid
	
	// Make an averages graph to see results.
	Display/K=1 NormAveragedSpectra[][0][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][1][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][2][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][3][0] vs FoundationWv
	
	ModifyGraph rgb(NormAveragedSpectra#1)=(0,0,65535),rgb(NormAveragedSpectra#2)=(2,39321,1),rgb(NormAveragedSpectra#3)=(26411,1,52428) 
	Label left "Normalized Average Intensity (a.u)";DelayUpdate
	Label bottom "Wavelength (nm)"
	Legend/C/N=text0/J/F=0/A=LT "\\s(NormAveragedSpectra) 280 nm\r\\s(NormAveragedSpectra#1) 350 nm\r\\s(NormAveragedSpectra#2) 405 nm\r\\s(NormAveragedSpectra#3) 450 nm"

End Function

//**************************************************************************************************

Function BTFMStats2([SaveFolder, ExcitationWavelength, ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts])
// Created 02/05/2020 by MPPF
// Updated from version 1 on 07/14/2020: can handle nonlinear particle numbers to agree with CropAnalysis2().
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: to provide some stats on ensemble of spectra retrieved from ParticleCropAnalysis(). Key outputs are
// NormAveragedSpectra [data, excitations, averages + std. devs.], 
// AveragedSpectra (same dimensions as NormAveragedSpectra), 
// PeakStats[Mean peak position + std. dev. in peak position + max. intensity of average wave + mean peak intensity + position of average wave peak, excitations], 
// Peaks [peak position + peak intensity, excitations, particles].
// If a new folder name is supplied as a string in the SaveFolder parameter, then key data will be duplicated into a subfolder of the current data folder with the name supplied.
//
// Notes on BTFMStats() in Maxwell Freeman Notebook #1, p. 51-54.
// Replaces SetupAverageSpectra().

	String SaveFolder
	Wave/T ExcitationWavelength
	Variable ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts

	// Initialize all default input parameters.
	if (ParamIsDefault(SaveFolder))
		SaveFolder = ""
	endif
	if (ParamIsDefault(ExcitationWavelength))
		Make/O/T ExcitationWavelength = {"280", "350", "405", "450"}
	endif
	if (ParamIsDefault(ThresholdVar))
		ThresholdVar = 0
	endif
	if (ParamIsDefault(LowerBound))
		LowerBound = 200
	endif
	if (ParamIsDefault(UpperBound))
		UpperBound = 900
	endif
	if (ParamIsDefault(DependentVarStepSize))
		DependentVarStepSize = 0.5
	endif
	if (ParamIsDefault(NumFoundationPts))
		NumFoundationPts = 1401
	endif

	

	// Declare the wave that keeps track of number of particles and numerical particle prefixes (e.g. the 4 in 4y280).
	Wave/T ImageWavePrefixes

	// Make FoundationWv as the skeleton onto which to perform regridding of the dependent variable.
	linspace("FoundationWv", LowerBound, UpperBound, Spacing = DependentVarStepSize)
	Wave FoundationWv
	
	// Make initial 3-tensor. Rows are points in spectra, columns are 4 excitation sources (ascending by wavelength), 
	// layers are different particles. Will be ix4xN where i is number of points in FoundationWv and N is number of particles.
	Make/O/N=(numpnts(FoundationWv),4,NumPnts(ImageWavePrefixes)) RegriddedSpectra = NaN
	
	// Make 3-tensor for averages and standard deviations where each row is a point in a spectrum,
	// each column is an excitation source, and the first and second layers are the averages and
	// and standard deviation waves respectively.
	Make/O/N=(numpnts(FoundationWv),4,2) AveragedSpectra = NaN
	
	// Make matrix called PeakStats that contains excitations as columns; first row is mean peak position, second
	// row is std. dev. of peak position, third row is max intensity of average spectrum, fourth row is average intensity
	// of CurveFit peaks. Compare rows 3 and 4 within a column to see if noise exceeds peak (if row 3 > row 4, then spectrum
	// is more intense outside of the peak--i.e. noise dominates). Temporarily add extra row for the peak location of the 
	// average spectra.
	Make/O/N=(5,4) PeakStats = NaN
	
	// Make 3-tensor named Peaks where first row is peak position, second row is peak intensity, columns are excitations,
	// and layers are for separate particles. Will be 2x4xN where N is number of particles.
	Make/O/N=(2,4,NumPnts(ImageWavePrefixes)) Peaks = NaN
	
	
	
	// Regrid spectra for each particle and add to RegriddedSpectra
	Variable ParticleIterator

	for (ParticleIterator = 0; ParticleIterator < NumPnts(ImageWavePrefixes); ParticleIterator += 1)

		String WavelengthsName = "p" + ImageWavePrefixes[ParticleIterator] + "xRef"
		
		// "For" loop to perform same thing for each excitation source.
		variable ExcitationID
		for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
		
			String IntensityName = "p" + ImageWavePrefixes[ParticleIterator] + "y" + ExcitationWavelength[ExcitationID]
			
			// Following command forces dependent variable of interest onto a separate grid with defined
			// start, end, and number of points.
			Regrid_Param2($WavelengthsName, $IntensityName, FoundationWv, FoundationWv)
		
			// RemoveValue clears all y-values equal to the mode (named ThresholdVar) in the regridded y-wave.
			Wave ywave_regrid
			if (ParamIsDefault(ThresholdVar))
				StatsQuantiles/Q/iNaN ywave_regrid
				ThresholdVar = V_mode
				ReplaceValue(ywave_regrid, ThresholdVar, NaN)
			else
				RemovePoints(ywave_regrid, ThresholdVar, 0)
			endif
			// Store data in RegriddedSpectra 3-tensor
			RegriddedSpectra[0,numpnts(ywave_regrid)-1][ExcitationID][ParticleIterator] = ywave_regrid[p][0][0]
			
			// Use CurveFit to get peak positions
			// Note: Wavelength (in nm) = DependentVarStepSize*index + LowerBound is akin to Fahrenheit = 1.8*Celsius + 32.
			CurveFit/Q gauss ywave_regrid[(400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize] /D	// NOTE: the values in parentheses are indices to fit between, not wavelengths; intention is to fit between 400 nm and 580 nm.
			Wave fit_ywave_regrid
			Wavestats/Q fit_ywave_regrid
			Peaks[0][ExcitationID][ParticleIterator] = V_maxloc*DependentVarStepSize + LowerBound
			Peaks[1][ExcitationID][ParticleIterator] = V_max
			
		endfor
	endfor
	
	
	
	// "For" loop to iteratively generate standard deviations and averages at each point in spectra.
	variable i
	for (i = 0; i < DimSize(RegriddedSpectra, 0); i += 1)
	
		for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
			
			Make/FREE/N=(DimSize(RegriddedSpectra, 2)) temp = RegriddedSpectra[i][ExcitationID][p]
			WaveStats/Q temp
			AveragedSpectra[i][ExcitationID][0] = V_avg
			AveragedSpectra[i][ExcitationID][1] = V_sdev
			
		endfor


	endfor
	
	
	
	// Create normalized spectra by dividing averages and std. devs. by maximum of average spectrum of 
	// given excitation source.
	// Also use Peaks to fill PeakStats with necessary data.
	Duplicate/O AveragedSpectra, NormAveragedSpectra
	
	for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
		
		Make/FREE/N=(DimSize(Peaks, 2)) temp = Peaks[0][ExcitationID][p]
		WaveStats/Q temp
		PeakStats[0][ExcitationID] = V_avg		// Mean peak position for each excitation source.
		PeakStats[1][ExcitationID] = V_sdev		// Standard deviation in peak position for each average above.
	
		Make/FREE/N=(DimSize(AveragedSpectra, 0)) temp = AveragedSpectra[p][ExcitationID][0]
		WaveStats/Q temp
		NormAveragedSpectra[][ExcitationID][] /= V_max
		PeakStats[2][ExcitationID] = V_max		// Maximum intensity of average spectrum.
		PeakStats[4][ExcitationID] = V_maxloc*DependentVarStepSize + LowerBound // Temporarily record the peak location in the average spectrum.
		
		Make/FREE/N=(DimSize(Peaks, 2)) temp = Peaks[1][ExcitationID][p]
		WaveStats/Q temp
		PeakStats[3][ExcitationID] = V_avg		// Average peak intensity of all curve fits to non-averaged spectra.
		
	endfor
	
	
	
	Killwaves fit_ywave_regrid, ywave_regrid
	
	// Make an averages graph to see results.
	Display/K=1 NormAveragedSpectra[][0][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][1][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][2][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][3][0] vs FoundationWv
	
	ModifyGraph rgb(NormAveragedSpectra#1)=(0,0,65535),rgb(NormAveragedSpectra#2)=(2,39321,1),rgb(NormAveragedSpectra#3)=(26411,1,52428) 
	Label left "Normalized Average Intensity (a.u)";DelayUpdate
	Label bottom "Wavelength (nm)"
	Legend/C/N=text0/J/F=0/A=LT "\\s(NormAveragedSpectra) 280 nm\r\\s(NormAveragedSpectra#1) 350 nm\r\\s(NormAveragedSpectra#2) 405 nm\r\\s(NormAveragedSpectra#3) 450 nm"

	if (!StringMatch(SaveFolder,""))
		Wave AspectRatio, MajMom, MinMom
		NewDataFolder $SaveFolder
		Duplicate ImageWavePrefixes :$(SaveFolder):ImageWavePrefixes
		Duplicate FoundationWv :$(SaveFolder):FoundationWv
		Duplicate NormAveragedSpectra :$(SaveFolder):NormAveragedSpectra
		Duplicate AspectRatio :$(SaveFolder):AspectRatio
		Duplicate MinMom :$(SaveFolder):MajMom
	endif

End Function

//**************************************************************************************************

Function RPStats([SaveFolder, ExcitationWavelength, ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts])
// Created 02/05/2020 by MPPF
// Updated from version 1 on 07/14/2020: can handle nonlinear particle numbers to agree with CropAnalysis2().
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: to provide some stats on ensemble of spectra retrieved from RPCropAnalysis(). Key outputs are
// NormAveragedSpectra [data, excitations, averages + std. devs.], 
// AveragedSpectra (same dimensions as NormAveragedSpectra), 
// PeakStats[Mean peak position + std. dev. in peak position + max. intensity of average wave + mean peak intensity + position of average wave peak, excitations], 
// Peaks [peak position + peak intensity, excitations, particles].
// If a new folder name is supplied as a string in the SaveFolder parameter, then key data will be duplicated into a subfolder of the current data folder with the name supplied.
//
// Notes on BTFMStats(), the parent function to this, is in Maxwell Freeman Notebook #1, p. 51-54.
// Replaces SetupAverageSpectra().

	String SaveFolder
	Wave/T ExcitationWavelength
	Variable ThresholdVar, LowerBound, UpperBound, DependentVarStepSize, NumFoundationPts

	// Initialize all default input parameters.
	if (ParamIsDefault(SaveFolder))
		SaveFolder = ""
	endif
	if (ParamIsDefault(ExcitationWavelength))
		Make/O/T/FREE ExcitationWavelength = {"280", "405", "450", "532"}
	endif
	if (ParamIsDefault(ThresholdVar))
		ThresholdVar = 0
	endif
	if (ParamIsDefault(LowerBound))
		LowerBound = 200
	endif
	if (ParamIsDefault(UpperBound))
		UpperBound = 900
	endif
	if (ParamIsDefault(DependentVarStepSize))
		DependentVarStepSize = 0.5
	endif
	if (ParamIsDefault(NumFoundationPts))
		NumFoundationPts = 1401
	endif

	

	// Declare the wave that keeps track of number of particles and numerical particle prefixes (e.g. the 4 in 4y280).
	Wave/T ImageWavePrefixes

	// Make FoundationWv as the skeleton onto which to perform regridding of the dependent variable.
	If (!ParamIsDefault(NumFoundationPts))
		linspace("FoundationWv", LowerBound, UpperBound, NumberofPoints = NumFoundationPts)
	Else
		linspace("FoundationWv", LowerBound, UpperBound, Spacing = DependentVarStepSize)
	EndIf
	Wave FoundationWv
	
	// Make initial 3-tensor. Rows are points in spectra, columns are 4 excitation sources (ascending by wavelength), 
	// layers are different particles. Will be ix4xN where i is number of points in FoundationWv and N is number of particles.
	// RegridSpecRedim is used for clustering and will be identical to RegriddedSpectra except the columns and layers will be swapped.
	Make/O/N=(numpnts(FoundationWv),NumPnts(ExcitationWavelength),NumPnts(ImageWavePrefixes)) RegriddedSpectra = NaN
	Make/O/N=(numpnts(FoundationWv),NumPnts(ImageWavePrefixes),NumPnts(ExcitationWavelength)) RegridSpecRedim = NaN
	
	// Make 3-tensor for averages and standard deviations where each row is a point in a spectrum,
	// each column is an excitation source, and the first and second layers are the averages and
	// and standard deviation waves respectively.
	Make/O/N=(numpnts(FoundationWv),NumPnts(ExcitationWavelength),2) AveragedSpectra = NaN
	
	// Make matrix called PeakStats that contains excitations as columns; first row is mean peak position, second
	// row is std. dev. of peak position, third row is max intensity of average spectrum, fourth row is average intensity
	// of CurveFit peaks. Compare rows 3 and 4 within a column to see if noise exceeds peak (if row 3 > row 4, then spectrum
	// is more intense outside of the peak--i.e. noise dominates). Temporarily add extra row for the peak location of the 
	// average spectra.
	Make/O/N=(5,NumPnts(ExcitationWavelength)) PeakStats = NaN
	
	// Make 3-tensor named Peaks where first row is peak position, second row is peak intensity, columns are excitations,
	// and layers are for separate particles. Will be 2x4xN where N is number of particles.
	Make/O/N=(2,NumPnts(ExcitationWavelength),NumPnts(ImageWavePrefixes)) Peaks = NaN
	
	
	
	// Regrid spectra for each particle and add to RegriddedSpectra
	Variable ParticleIterator

	for (ParticleIterator = 0; ParticleIterator < NumPnts(ImageWavePrefixes); ParticleIterator += 1)

		String WavelengthsName = "p" + ImageWavePrefixes[ParticleIterator] + "xRef"
		
		// "For" loop to perform same regridding for each excitation source.
		variable ExcitationID
		for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
		
			String IntensityName = "p" + ImageWavePrefixes[ParticleIterator] + "_" + ExcitationWavelength[ExcitationID] + "nm"
			
			// Following command forces dependent variable of interest onto a separate grid with defined
			// start, end, and number of points.
			Regrid_Param2($WavelengthsName, $IntensityName, FoundationWv, FoundationWv)
		
			// RemoveValue clears all y-values equal to the mode (named ThresholdVar) in the regridded y-wave.
			Wave ywave_regrid
			if (ParamIsDefault(ThresholdVar))
				StatsQuantiles/Q/iNaN ywave_regrid
				ThresholdVar = V_mode
				ReplaceValue(ywave_regrid, ThresholdVar, NaN)
			else
				RemovePoints(ywave_regrid, ThresholdVar, 0)
			endif
			
			// Normalize to highest point and store data in RegriddedSpectra 3-tensor
//			Wavestats/Q/RMD=[0,*][ExcitationID][ParticleIterator,ParticleIterator] ywave_regrid
//			ywave_regrid /= V_max
			RegriddedSpectra[0,numpnts(ywave_regrid)-1][ExcitationID][ParticleIterator] = ywave_regrid[p][0][0]
			RegridSpecRedim[0,numpnts(ywave_regrid)-1][ParticleIterator][ExcitationID] = ywave_regrid[p][0][0]
			
			// Use CurveFit to get peak positions
			// Note: Wavelength (in nm) = DependentVarStepSize*index + LowerBound is akin to Fahrenheit = 1.8*Celsius + 32.
			CurveFit/Q gauss ywave_regrid[(400-LowerBound)/DependentVarStepSize,(580-LowerBound)/DependentVarStepSize] /D	// NOTE: the values in parentheses are indices to fit between, not wavelengths; intention is to fit between 400 nm and 580 nm.
			Wave fit_ywave_regrid
			Wavestats/Q fit_ywave_regrid
			Peaks[0][ExcitationID][ParticleIterator] = V_maxloc*DependentVarStepSize + LowerBound
			Peaks[1][ExcitationID][ParticleIterator] = V_max
			
		endfor
	endfor
	
	
	ReplaceValue(RegridSpecRedim, NaN, 0)
	
	
	// "For" loop to iteratively generate standard deviations and averages at each point in spectra.
	variable i
	for (i = 0; i < DimSize(RegriddedSpectra, 0); i += 1)
	
		for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
			
			Make/FREE/N=(DimSize(RegriddedSpectra, 2)) temp = RegriddedSpectra[i][ExcitationID][p]
			WaveStats/Q temp
			AveragedSpectra[i][ExcitationID][0] = V_avg
			AveragedSpectra[i][ExcitationID][1] = V_sdev
			
		endfor


	endfor
	
	
	
	// Create normalized spectra by dividing averages and std. devs. by maximum of average spectrum of 
	// given excitation source.
	// Also use Peaks to fill PeakStats with necessary data.
	Duplicate/O AveragedSpectra, NormAveragedSpectra
	
	for (ExcitationID = 0; ExcitationID < numpnts(ExcitationWavelength); ExcitationID += 1)
		
		Make/FREE/N=(DimSize(Peaks, 2)) temp = Peaks[0][ExcitationID][p]
		WaveStats/Q temp
		PeakStats[0][ExcitationID] = V_avg		// Mean peak position for each excitation source.
		PeakStats[1][ExcitationID] = V_sdev		// Standard deviation in peak position for each average above.
	
		Make/FREE/N=(DimSize(AveragedSpectra, 0)) temp = AveragedSpectra[p][ExcitationID][0]
		WaveStats/Q temp
		NormAveragedSpectra[][ExcitationID][] /= V_max
		PeakStats[2][ExcitationID] = V_max		// Maximum intensity of average spectrum.
		PeakStats[4][ExcitationID] = V_maxloc*DependentVarStepSize + LowerBound // Temporarily record the peak location in the average spectrum.
		
		Make/FREE/N=(DimSize(Peaks, 2)) temp = Peaks[1][ExcitationID][p]
		WaveStats/Q temp
		PeakStats[3][ExcitationID] = V_avg		// Average peak intensity of all curve fits to non-averaged spectra.
		
	endfor
	
	
	
	Killwaves fit_ywave_regrid, ywave_regrid
	
	// Make an averages graph to see results.
	Display/K=1 NormAveragedSpectra[][0][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][1][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][2][0] vs FoundationWv
	AppendtoGraph NormAveragedSpectra[][3][0] vs FoundationWv
	
	ModifyGraph lstyle(NormAveragedSpectra)=2,rgb(NormAveragedSpectra)=(26214,26214,26214),rgb(NormAveragedSpectra#1)=(26411,1,52428),rgb(NormAveragedSpectra#2)=(0,0,65535),rgb(NormAveragedSpectra#3)=(3,52428,1)
	Label left "Normalized Average Intensity (a.u)";DelayUpdate
	Label bottom "Wavelength (nm)"
	Legend/C/N=text0/J/F=0/A=LT "\\s(NormAveragedSpectra) 280 nm\r\\s(NormAveragedSpectra#1) 405 nm\r\\s(NormAveragedSpectra#2) 450 nm\r\\s(NormAveragedSpectra#3) 532 nm"

	if (!StringMatch(SaveFolder,""))
		Wave AspectRatio, MajMom, MinMom
		NewDataFolder $SaveFolder
		Duplicate ImageWavePrefixes :$(SaveFolder):ImageWavePrefixes
		Duplicate FoundationWv :$(SaveFolder):FoundationWv
		Duplicate NormAveragedSpectra :$(SaveFolder):NormAveragedSpectra
		Duplicate AspectRatio :$(SaveFolder):AspectRatio
		Duplicate MinMom :$(SaveFolder):MajMom
	endif

End Function

//**************************************************************************************************

Function CropAnalysis2([pathName, nameWaveAfterFileName, displayImages]) // ("", 1, 0)
 // IMPORT TIFFS FROM SUBFOLDER
 	String pathName
	Variable nameWaveAfterFileName, displayImages
	
	
	if (ParamIsDefault(pathName))											// Set default parameters to default values
		pathName = ""
	endif
	if (ParamIsDefault(nameWaveAfterFileName))
		nameWaveAfterFileName = 1
	endif
	if (ParamIsDefault(displayImages))
		displayImages = 0
	endif
	
	
	if (strlen(pathName) == 0 )														// see if there's a given path/folder.  Quotations ("") returns 
		NewPath/O/M="Choose a folder containing TIFFs" pathName
		If (V_flag != 0)
			Return -1
		EndIf
	EndIf
	
	Make/O/N=1000 Exposure280			// Make a wave for exposure time for 280
	Make/O/N=1000 Exposure350			// 350
	Make/O/N=1000 Exposure405 			// 405
	Make/O/N=1000 Exposure450			// 450
	
	Exposure280 = NaN					// Make them all NaN
	Exposure350 = NaN
	Exposure405 = NaN
	Exposure450 = NaN
	
	
	String fileName, fileLoc
	Make/N=0/T/O ImageWaveNames, ImageWavePrefixes
	
	fileLoc = IndexedFile(pathName, -1, ".tif")
	fileLoc = SortList(FileLoc, ";", 16)
	
	Variable n = 0, Count = 1
	
	do
		fileName = StringFromList(n, FileLoc, ";")		// find ".tif" files in the specified folder
		if (strlen(fileName) == 0)							// If a ".tif" image exists, keep going, if not, break the loop
			count -= 1
			break										// break the loop
		endif
		
		ImageLoad/Q/P=pathName/T=TIFF/N=image fileName	// Load the image into igor, specifically Tiffs
		
		if (V_Flag)
			string name = StringFromList(0, S_waveNames)	// 
			wave w = $name
			
			if (nameWaveAfterFileName)
				string desiredName = S_fileName								// desiredName is equal to the full image name + ".tif"
				desiredName = ParseFilepath(3, S_filename, ":", 0, 0)				// Parse the file's name and put it into S_Filename (option 3 removes file extension)
						
						String NoTiff = desiredName								
										
						String expr= "%s %e"									// Make expr a string EXPRESSION equal to "%s %e" // %s stores all text to the next white space, %e catches the number exp behind it
						String LaserExSt
						Variable ExpTime, ParticlePrefix
						
						sscanf NoTiff, Expr, LaserExSt, ExpTime				// Separate the string NoTiff by Expr (%s %e) and split it into LaserExSt string and ExpTime string
						sscanf NoTiff, "%d", ParticlePrefix
						
						If (!StringMatch(ImageWavePrefixes[inf],num2str(ParticlePrefix)))
							InsertPoints NumPnts(ImageWavePrefixes), 1, ImageWavePrefixes
							InsertPoints NumPnts(ImageWaveNames), 1, ImageWaveNames
							ImageWaveNames[inf] = LaserExSt	
							ImageWavePrefixes[inf] = num2str(ParticlePrefix)
						EndIf

						String ExcValue280 = num2str(ParticlePrefix) + "y280"		// Make strings for each wave type per particle
						String Excvalue350 = num2str(ParticlePrefix) + "y350"
						String Excvalue405 = num2str(ParticlePrefix) + "y405"
						String Excvalue450 = num2str(ParticlePrefix) + "y450"
						String ExcValueCalib = num2Str(ParticlePrefix) + "yCalib"						
						
						If (Stringmatch(LaserExSt, ExcValue280))				// Add the exposure time number to each exposure time wave.
							Exposure280[NumPnts(ImageWavePrefixes)-1] = ExpTime
						ElseIf (Stringmatch(LaserExSt, ExcValue350))
							Exposure350[NumPnts(ImageWavePrefixes)-1] = ExpTime
						ElseIf (Stringmatch(LaserExSt, ExcValue405))
							Exposure405[NumPnts(ImageWavePrefixes)-1] = ExpTime
						Elseif (Stringmatch(LaserExSt, ExcValue450))
							Exposure450[NumPnts(ImageWavePrefixes)-1] = ExpTime	
							Count +=1
						EndIf		
						
			
						if (Exists(LaserExSt))								// If the desired name exists, this wave gets an identical name with n1, n2, n... starting from 0
							LaserExSt = UniqueName(LaserExSt, 1, 0)
						EndIf
						
						Rename w, $LaserExSt
						
																
				EndIf
			
			if (displayImages)
				NewImage w
			endif
		endif
		n += 1 
		
	while (1)







 	DeleteNaNs(Exposure280) // Function defined in AlexTools_v16.5.
	DeleteNaNs(Exposure350)
	DeleteNaNs(Exposure405)
	DeleteNaNs(Exposure450)
	
 	// Calibrating Factors fix
 	Exposure280 *=1000 			// Change seconds to milliseconds
 	Exposure350 *=1000
 	Exposure405 *=1000
 	Exposure450 *=1000
 	
 	Make/O/N=4 PowerDensity
 	PowerDensity[0] = 0.004		// PD factor correction for 280
 	PowerDensity[1] = 0.5			// PD factor correction for 350
 	PowerDensity[2] = 4.5			// PD factor correction for 405
 	PowerDensity[3] = 0.3			// PD factor correction for 450
 	
 	Variable rn, t, r, y, j, b, o, h
	String wCalib, Calib, w450, w405, w350, w280, yCalib, y450, y405, y350, y280, pXRef, CurrentRef


	b = 0
	j = 1
	o = 1	
	n = 1
	

	For (n=0; n < NumPnts(ImageWavePrefixes); n+=1)

		yCalib = ImageWavePrefixes[n]+"yCalib"
		y450 = ImageWavePrefixes[n]+"y450"
		y405 = ImageWavePrefixes[n]+"y405"
		y350 = ImageWavePrefixes[n]+"y350"
		y280 = ImageWavePrefixes[n]+"y280"
	
		If (WaveExists($yCalib))
	
			Wave iCalib = $yCalib
			Wave i450 = $y450
			Wave i405 = $y405
			Wave i350 = $y350
			Wave i280 = $y280
		
			
			h = DimSize(iCalib,0)			
					
			Make/O /N=(h) XRef							// Make XRef wave equal to r (XRef is X dimension for spectra)
			t = 1
			

			linspace("xRef", 1, DimSize($yCalib, 0))
				
			CurrentRef = "p"+ImageWavePrefixes[n]+"xRef"
			Make/O /N=(h) $CurrentRef = xRef
		

			wCalib = "p"+ImageWavePrefixes[n]+"yCalib"
			w450 = "p"+ImageWavePrefixes[n]+"y450"
			w405 = "p"+ImageWavePrefixes[n]+"y405"
			w350 = "p"+ImageWavePrefixes[n]+"y350"
			w280 = "p"+ImageWavePrefixes[n]+"y280"
		
			MatrixOP/O AveragesWave = sumRows(iCalib)/numCols(iCalib)
			Make/O/N=(DimSize($CurrentRef,0)) $wCalib = AveragesWave
					
			MatrixOP/O AveragesWave = sumRows(i450)/numCols(i450)
			Make/O/N=(DimSize($CurrentRef,0)) $w450 = AveragesWave

			MatrixOP/O AveragesWave = sumRows(i405)/numCols(i405)
			Make/O/N=(DimSize($CurrentRef,0)) $w405 = AveragesWave
	
			MatrixOP/O AveragesWave = sumRows(i350)/numCols(i350)
			Make/O/N=(DimSize($CurrentRef,0)) $w350 = AveragesWave

			MatrixOP/O AveragesWave = sumRows(i280)/numCols(i280)
			Make/O/N=(DimSize($CurrentRef,0)) $w280 = AveragesWave
			
			b += 1			
			
		EndIf
	EndFor


	Wave Fit_CurveFitting
	Variable l
	n = 1
	j = 1
	
	Variable CountingUp, DelPts
     	Wave MajMom, MinMom, AngMom, AspectRatio, YSize, Degrees, XSize
     	     	
	CountingUp = 1
       Make/O/N=50 MajMom, MinMom, AngMom, AspectRatio, YSize, Degrees, XSize, BlueCalib, RedCalib
            
       Degrees = NaN
       MajMom = NaN
       MinMom = NaN
       AngMom = NaN
       AspectRatio = NaN
	YSize = NaN
	XSIze = NaN
	BlueCalib = NaN
	RedCalib = NaN
																								///////////////////////////////////////////////////////////////////
	Variable BrightnessThreshold = 25 // in grayscale intensity					//// Edit these to change thresholds of detection of particles,	//
	Variable ParticleSizeThreshold = 25 // in pixels									//// thereby changing calibrations.										//
       																						///////////////////////////////////////////////////////////////////
	For (n=0; n < NumPnts(ImageWavePrefixes); n+=1)

		String Measured, BinMeas
		Measured = ImageWavePrefixes[n] + "yCalib"

		Wave w0=$Measured
         												
		ImageThreshold/Q/T=(BrightnessThreshold)/I w0						
		Wave M_ImageThresh                                             
		Duplicate/O M_ImageThresh, RedPointWave
         
		ImageThreshold/Q/T=(BrightnessThreshold)/I w0
		Duplicate/O M_ImageThresh, BluePointWave
         
		Redimension/N=(-1,-1) RedPointWave
		Redimension/N=(-1,-1) BluePointWave
           
		Duplicate/O/R=(400,820)(0,104) RedPointWave, ParticleMeasure // (X)(Y) points
		Duplicate/O/R=(400,820)(0,104) RedPointWave, RedMeasure
		Duplicate/O/R=(0,350)(0,104) BluePointWave, BlueMeasure
        
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, BlueMeasure
		Wave M_Moments
		BlueCalib[n] = M_Moments[0]
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, RedMeasure
		RedCalib[n] = M_Moments[0]+400
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, ParticleMeasure
        
		MajMom[n] = M_Moments[2]*2
		MinMom[n] = M_Moments[3]*2
          
		If (M_Moments[4]>pi)
			M_Moments[4] = M_Moments[4]-pi
			AngMom[n] = M_Moments[4]
		Else
			AngMom[n] = M_Moments[4]
		EndIf
           
		AspectRatio[n] = MajMom[n]/MinMom[n]
     
		Wave W_YMax, W_YMin, W_XMax, W_XMin
		YSize[n] = W_YMax[0]-W_YMin[0]
		XSize[n] = W_XMax[0]-W_XMin[0]
        
		Degrees[n] = AngMom[n]*(180/pi)
	
		CountingUp +=1
		
   EndFor
   
   DeleteNaNs(MajMom)
   DeleteNaNs(BlueCalib)
   DeleteNaNs(RedCalib)
   DeleteNaNs(MinMom)
   DeleteNaNs(AngMom)
   DeleteNaNs(AspectRatio)
   DeleteNaNs(YSize)
   DeleteNaNs(XSize)
   DeleteNaNs(Degrees)
	
	
	String SizeMeas405
	Make/o/n=1000 SizeWave
	SizeWave = NaN
	
	For (n = 0; n < NumPnts(ImageWavePrefixes); n += 1)
		
		Make/O/N=2 DestWave
		Variable Mid
		
		Sizemeas405 = ImageWavePrefixes[n]+"y405"
		Wave MeasuredWave = $SizeMeas405
		String pSize, pSizeRef, SizeMeas2
		
		pSize = ImageWavePrefixes[n]+"pSize"		
			
		MatrixOP/O sizingwave = sumCols(MeasuredWave)/numRows(MeasuredWave)
		Make/O/N=(DimSize($Sizemeas405,1)) $pSize = SizingWave
		
		Wave CurveFitting = $pSize
		
		CurveFit/Q/NTHR=0 gauss CurveFitting /D
		String SizeMeas = "fit_" + pSize
		Wave Subtractor = $SizeMeas				
		Wavestats/Q Subtractor						// Measure the Guassians' characteristics	
		
		CurveFitting -= V_Min
		
		CurveFit/Q/NTHR=0 gauss CurveFitting /D
		SizeMeas2 = "fit_" + pSize
		Wave Subtractor = $SizeMeas2				
		Wavestats/Q Subtractor		
		
		Mid = V_max/2
		FindLevels/Q/D=destWave $SizeMeas, Mid
		
		SizeWave[n] = (DestWave[1]-DestWave[0])
		
	
	EndFor
	
	DeleteNaNs(SizeWave)
	
	
	String NameStr, p, Fit, Calibrate
		
	j = 1
	Make /O /N=1000 BluePts
	BluePts = NaN
	Make /O /N=1000 RedPts
	RedPts = NaN
	
	For (n = 0; n < NumPnts(ImageWavePrefixes); n += 1)					// Still need to try to recognize if the wave exists.
	
		NameStr = "p" + ImageWavePrefixes[n] + "yCalib"		// Create a string out of p[n]yCalib, n being the particle number
		Wave Calibrating = $NameStr				// Make the NameStr into a reference and make the Wave Calibrating that string
	
		
		if (WaveExists(Calibrating))
			CurveFit/Q/NTHR=0 gauss  Calibrating[0,350] /D // Measure the Gaussian around the blue calib point
			Calibrate = "fit_" + NameStr
			Wave CalibWave = $Calibrate
			Wavestats/Q CalibWave						// Measure the Guassians' characteristics	
			BluePts[n] = V_MaxLoc					// Add the Max Location Calib point to the BluePts wave 
		
			CurveFit/Q/NTHR=0 gauss  Calibrating[400,820]  /D // Measure the Gaussian around the red calib point
			Calib = "fit_" + NameStr
			Wave CalibWave = $Calib
			Wavestats/Q CalibWave
			RedPts[n] = V_MaxLoc
		EndIf
	EndFor
	

	DeleteNaNs(Bluepts)
	DeleteNaNs(Redpts)


	Variable c, e					// n, c, r
	Wave BluePts, RedPts
	String NameStrRef
	
	c = 632.8-403.457
	

	For (y = 0; y < DimSize(BluePts,0); y += 1)				
		
			NameStrRef = "p" + ImageWavePrefixes[y] + "xREF"
			
			Wave CalibRef = $NameStrRef

			e = c/(RedPts[y]-BluePts[y])
			
			CalibRef -= BlueCalib[y]
			
			CalibRef *= e
		
			CalibRef += 403.457
			
	EndFor

	
	String NameStrX, NameStrY, NewSizeWaveX, NewSizeF
	Wave BluePts, YSize, XSize, FSize

	//////////////////////////////////////////////////
	// Copy command for new sizing //
	//////////////////////////////////////////////////
	String NewSizeWave, NameStrXSize, NameStrFSize
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "y450"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "y450_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "y450_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "y450_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "y405"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "y405_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "y405_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "y405_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "y350"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "y350_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "y350_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "y350_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)

		NameStrY = "p" + ImageWavePrefixes[n] + "y280"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "y280_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "y280_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "y280_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	String PeakSW, PeakY, PeakX, PeakXRef, PeakF
	Wave Bluepts, AngMom
	
	Make/O/N=(DimSize(AngMom, 0)) PeaksSW, PeaksY, PeaksX, PeaksF
	
	PeaksSW = NaN
	PeaksY = NaN
	PeaksX = NaN
	PeaksF = NaN
	Variable up
	Up = 0
	
	
	Wave Degrees, YSize, XSize, FSize
		

		Up = 0
		
		Make/O FactorSplitY, FactorSplitX, FSize
	
		FactorSplitY = NaN
		FactorSplitX = NaN
		FSize = NaN
		
	For	(up = 0; up<DimSize(Degrees, 0); up +=1)
		If (Degrees[Up]<90)
			FactorSplitY[up] = Degrees[up]/90
			FactorSplitX[up]  = 1-FactorsplitY[up] 
		ElseIf (Degrees[Up]>90)
			FactorSplitX[up] = Degrees[up]/180
			FactorSplitY[up]  = 1-FactorsplitX[up]
		EndIf
		
		FSize[up] = abs((FactorSplitY[up]*YSize[up]))+((FactorSplitX[up]*XSize[up]))
	EndFor
	

	DeleteNaNs(FactorSplitX)
	DeleteNaNs(FactorSplitY)
	DeleteNaNs(FSize)
	
	////////////////////////////////////////////////
	// Display & Calibration Section //
	////////////////////////////////////////////////
	
	Display/K=1
	j = 1
	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
																	// Create a string out of p[n]yCalib, n being the particle number												
			NewSizeWave = "p" + ImageWavePrefixes[n] + "y450_New"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "y450_NewX"
			NewSizeF = "p" + ImageWavePrefixes[n] + "y450_NewF"
			NameStrY = "p" + ImageWavePrefixes[n] + "y450"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			
						
			if (!waveexists($NameStrY))
				Break
			EndIf
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[3]
			DisplayY *=100000
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure450[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[3]
			DisplayNew *=100000	
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure450[n]
			
			Wavestats/Q DisplayNewX
			DisplayNewx -=V_min
			DisplayNewx /=PowerDensity[3]
			DisplayNewx *=100000
			DisplayNewx /= XSize[n]
			DisplayNewx /= Exposure450[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[3]
			DisplayNewF *=100000
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure450[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)						// Still need to try to recognize if the wave exists.
																				// Create a string out of p[n]yCalib, n being the particle number
			NameStrY = "p" + ImageWavePrefixes[n] + "y405"
			NewSizeWave = "p" + ImageWavePrefixes[n] + "y405_New"
			Wave DisplayNew = $NewSizeWave
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "y405_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			
			NewSizeF = "p" + ImageWavePrefixes[n] + "y405_NewF"
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[2]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure405[n]

			Wavestats/Q DisplayNew
			
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[2]
			DisplayNew *=100000
				
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure405[n]	
			
			Wavestats/Q DisplayNewX
			
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[2]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure405[n]
			
			Wavestats/Q displaynewF
			
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[2]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure405[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	
	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
		
			NameStrY = "p" + ImageWavePrefixes[n] + "y350"		// Create a string out of p[n]yCalib, n being the particle number
			NewSizeWave = "p" + ImageWavePrefixes[n] + "y350_New"
			Wave DisplayNew = $NewSizeWave
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "y350_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			
			NewSizeF = "p" + ImageWavePrefixes[n] + "y350_NewF"
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[1]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure350[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[1]
			DisplayNew *=100000
			
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure350[n]	
			
			Wavestats/Q DisplayNewX
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[1]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure350[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[1]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure350[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	
	n = 1
	For(n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
																		// Create a string out of p[n]yCalib, n being the particle number
			NameStrY = "p" + ImageWavePrefixes[n] + "y280"		// Create a string out of p[n]yCalib, n being the particle number
			NewSizeWave = "p" + ImageWavePrefixes[n] + "y280_New"
			Wave DisplayNew = $NewSizeWave
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "y280_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			
			NewSizeF = "p" + ImageWavePrefixes[n] + "y280_NewF"
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[0]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure280[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[0]
			DisplayNew *=100000
			
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure280[n]	
			
			Wavestats/Q DisplayNewX
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[0]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure280[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[0]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure280[n]
			
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"

	Up = 0
	
	For (n = 0; n < DimSize(AngMom, 0); n+=1)

		PeakSW = "p" + ImageWavePrefixes[n] + "y450"
		PeakXRef = "p" + ImageWavePrefixes[n] + "xRef"
		PeakY = "p" + ImageWavePrefixes[n] + "y450_New"
		PeakX = "p" + ImageWavePrefixes[n] + "y450_NewX"
		PeakF = "p" + ImageWavePrefixes[n] + "y450_NewF"
		
		If (Waveexists($PeakSw)==0)
			Break
		EndIf
		
		Wavestats/Q $PeakSW
		PeaksSW[n] = V_Max
	
		Wavestats/Q $PeakY
		PeaksY[n] = V_Max
		
		Wavestats/Q $PeakX
		PeaksX[n] = V_Max
		
		Wavestats/Q $PeakF
		PeaksF[n] = V_max
		
		Wave Degrees
		
		Make/O FactorSplitY, FactorSplitX
	
		FactorSplitY = NaN
		FactorSplitX = NaN
			
		If (Degrees[Up]<90)
			FactorSplitY = Degrees/90
			FactorSplitX  = 1-FactorsplitY 
		ElseIf (Degrees[Up]>90)
			FactorSplitY = 90/Degrees
			FactorSplitX  = 1-FactorsplitY 
		EndIf
		up += 1
	EndFor
	
	DeleteNaNs(FactorSplitX)
	DeleteNaNs(FactorSplitY)
	DeleteNaNs(PeaksSW)
	DeleteNaNs(PeaksY)
	DeleteNaNs(PeaksX)
	
	Edit/K=1
	
	AppendToTable Majmom, Minmom, AspectRatio
End


//**************************************************************************************************

Function RPCropAnalysis([pathName, nameWaveAfterFileName, displayImages]) // ("", 1, 0)
 // IMPORT TIFFS FROM SUBFOLDER
 	String pathName
	Variable nameWaveAfterFileName, displayImages
	
	
	if (ParamIsDefault(pathName))											// Set default parameters to default values
		pathName = ""
	endif
	if (ParamIsDefault(nameWaveAfterFileName))
		nameWaveAfterFileName = 1
	endif
	if (ParamIsDefault(displayImages))
		displayImages = 0
	endif
	
//////////////// Initialize important variables. //////////////////////////

																									///////////////////////////////////////////////////////////////////
	Variable BrightnessThreshold = 50 // in (unsigned byte) grayscale intensity	//// Edit these to change thresholds of detection of particles,	//
	Variable ParticleSizeThreshold = 25 // in pixels										//// thereby changing calibrations.										//
																									///////////////////////////////////////////////////////////////////
	
	
/////////////// End variable declarations. ///////////////////////////////	
	
	if (strlen(pathName) == 0 )														// see if there's a given path/folder.  Quotations ("") returns 
		NewPath/O/M="Choose a folder containing TIFFs" pathName
		If (V_flag != 0)
			Return -1
		EndIf
	EndIf
	
	Make/O/N=1000 Exposure280			// Make a wave for exposure time for 280
	Make/O/N=1000 Exposure405 		// 405
	Make/O/N=1000 Exposure450			// 450
	Make/O/N=1000 Exposure532			// 532
	
	Exposure280 = NaN					// Make them all NaN
	Exposure405 = NaN
	Exposure450 = NaN
	Exposure532 = NaN
	
	
	String fileName, fileLoc
	Make/N=0/T/O ImageWavePrefixes//, ImageWaveNames
	
	fileLoc = IndexedFile(pathName, -1, ".tif")
	fileLoc = SortList(FileLoc, ";", 16)
	
	Variable n = 0
	
	do
		fileName = StringFromList(n, FileLoc, ";")		// find ".tif" files in the specified folder
		if (strlen(fileName) == 0)							// If a ".tif" image exists, keep going, if not, break the loop
			break										// break the loop
		endif
		
		ImageLoad/Q/O/P=pathName/T=TIFF/N=image fileName	// Load the image into igor, specifically Tiffs
		
		if (V_Flag)
			string name = StringFromList(0, S_waveNames)	// 
			wave w = $name
			
			if (nameWaveAfterFileName)
				string desiredName = S_fileName								// desiredName is equal to the full image name + ".tif"
				desiredName = ParseFilepath(3, S_filename, ":", 0, 0)				// Parse the file's name and put it into S_Filename (option 3 removes file extension)
						
						String NoTiff = desiredName								
										
//						String expr = "%s %e"									// Make expr a string EXPRESSION equal to "%s %e" // %s stores all text to the next white space, %e catches the number exp behind it
						String expr = "%s %*[#] %s %*s %*s %e %*[ISO] %e %*[sec]"		// This string EXPRESSION is for RaspPiSpec/FOBS
						String LaserExSt, SampleName
						Variable ExpTime, ParticlePrefix, Gain
						
//						sscanf NoTiff, Expr, LaserExSt, ExpTime				// Separate the string NoTiff by Expr (%s %e) and split it into LaserExSt string and ExpTime string
						sscanf NoTiff, Expr, SampleName, LaserExSt, Gain, ExpTime
//						sscanf NoTiff, "%d", ParticlePrefix
						
//						If (!StringMatch(ImageWavePrefixes[inf],num2str(ParticlePrefix)))
//							InsertPoints NumPnts(ImageWavePrefixes), 1, ImageWavePrefixes
//							InsertPoints NumPnts(ImageWaveNames), 1, ImageWaveNames
//							ImageWaveNames[inf] = LaserExSt	
//							ImageWavePrefixes[inf] = num2str(ParticlePrefix)
//						EndIf

						If (!StringMatch(ImageWavePrefixes[inf],SampleName))
							InsertPoints NumPnts(ImageWavePrefixes), 1, ImageWavePrefixes
//							InsertPoints NumPnts(ImageWaveNames), 1, ImageWaveNames
//							ImageWaveNames[inf] = LaserExSt	
							ImageWavePrefixes[inf] = SampleName
						EndIf

//						String ExcValue280 = num2str(ParticlePrefix) + "y280"		// Make strings for each wave type per particle
//						String Excvalue532 = num2str(ParticlePrefix) + "y532"
//						String Excvalue405 = num2str(ParticlePrefix) + "y405"
//						String Excvalue450 = num2str(ParticlePrefix) + "y450"
//						String ExcValueCalib = num2Str(ParticlePrefix) + "yCalib"						
						
						If (Stringmatch(LaserExSt, "280nm"))				// Add the exposure time number to each exposure time wave.
							Exposure280[NumPnts(ImageWavePrefixes)-1] = ExpTime
						ElseIf (Stringmatch(LaserExSt, "405nm"))
							Exposure405[NumPnts(ImageWavePrefixes)-1] = ExpTime
						Elseif (Stringmatch(LaserExSt, "450nm"))
							Exposure450[NumPnts(ImageWavePrefixes)-1] = ExpTime
						ElseIf (Stringmatch(LaserExSt, "532nm"))
							Exposure532[NumPnts(ImageWavePrefixes)-1] = ExpTime	
						EndIf		
						
			
//						if (Exists(LaserExSt))								// If the desired name exists, this wave gets an identical name with n1, n2, n... starting from 0
//							LaserExSt = UniqueName(LaserExSt, 1, 0)
//						EndIf

						String WaveNameString = SampleName + "_" + LaserExSt		// SampleName could be anything. LaserExSt will have values {"280nm","405nm","450nm","532nm","Calib"} May also contain "incandescent".
						ImageTransform/Q rgb2gray w
						Wave M_RGB2Gray
						Rename M_RGB2Gray, $WaveNameString
						
																
				EndIf
			
			if (displayImages)
				NewImage w
			endif
		endif
		n += 1 
		
	while (1)







 	DeleteNaNs(Exposure280) // Function defined in AlexTools_v16.5.
	DeleteNaNs(Exposure405)
	DeleteNaNs(Exposure450)
	DeleteNaNs(Exposure532)
	
 	// Calibrating Factors fix
 	Exposure280 *=1000 			// Change seconds to milliseconds
 	Exposure532 *=1000
 	Exposure405 *=1000
 	Exposure450 *=1000
 	
 	Make/O/N=4 PowerDensity			// These only apply to BTFM \\
 	PowerDensity[0] = 1//0.004		// PD factor correction for 280
 	PowerDensity[1] = 1//0.5			// PD factor correction for 350
 	PowerDensity[2] = 1//4.5			// PD factor correction for 405
 	PowerDensity[3] = 1//0.3			// PD factor correction for 450
 	
// 	Variable rn, t, r, y, j, b, o, h
	Variable h
	String wCalib, Calib, w450, w405, w532, w280, yCalib, y450, y405, y532, y280, pXRef, CurrentRef


//	b = 0
//	j = 1
//	o = 1
	

	For (n=0; n < NumPnts(ImageWavePrefixes); n+=1)

		yCalib = ImageWavePrefixes[n]+"_Calib"
		y532 = ImageWavePrefixes[n]+"_532nm"
		y450 = ImageWavePrefixes[n]+"_450nm"
		y405 = ImageWavePrefixes[n]+"_405nm"
		y280 = ImageWavePrefixes[n]+"_280nm"
	
		If (WaveExists($yCalib))
	
			Wave iCalib = $yCalib
			Wave i450 = $y450
			Wave i405 = $y405
			Wave i532 = $y532
			Wave i280 = $y280
		
			
			h = DimSize(iCalib,0)			
					
			Make/O/N=(h) XRef							// Make XRef wave equal to r (XRef is X dimension for spectra)
//			t = 1
			

			linspace("xRef", 1, DimSize($yCalib, 0))
				
			CurrentRef = "p"+ImageWavePrefixes[n]+"xRef"
			Make/O /N=(h) $CurrentRef = xRef
		

			wCalib = "p"+ImageWavePrefixes[n]+"_Calib"
			w450 = "p"+ImageWavePrefixes[n]+"_450nm"
			w405 = "p"+ImageWavePrefixes[n]+"_405nm"
			w532 = "p"+ImageWavePrefixes[n]+"_532nm"
			w280 = "p"+ImageWavePrefixes[n]+"_280nm"
		
			MatrixOP/O AveragesWave = sumRows(iCalib)/numCols(iCalib)
			Make/O/N=(DimSize($CurrentRef,0)) $wCalib = AveragesWave
					
			MatrixOP/O AveragesWave = sumRows(i450)/numCols(i450)
			Make/O/N=(DimSize($CurrentRef,0)) $w450 = AveragesWave

			MatrixOP/O AveragesWave = sumRows(i405)/numCols(i405)
			Make/O/N=(DimSize($CurrentRef,0)) $w405 = AveragesWave
	
			MatrixOP/O AveragesWave = sumRows(i532)/numCols(i532)
			Make/O/N=(DimSize($CurrentRef,0)) $w532 = AveragesWave

			MatrixOP/O AveragesWave = sumRows(i280)/numCols(i280)
			Make/O/N=(DimSize($CurrentRef,0)) $w280 = AveragesWave
			
//			b += 1			
			
		EndIf
	EndFor


	Wave Fit_CurveFitting
//	Variable l
//	n = 1
//	j = 1
	
//	Variable CountingUp, DelPts
// 	Wave MajMom, MinMom, AngMom, AspectRatio, YSize, Degrees, XSize
     	     	
//	CountingUp = 1
	Make/O/N=50 MajMom, MinMom, AngMom, AspectRatio, YSize, Degrees, XSize, BlueCalib, RedCalib
	    
	Degrees = NaN
	MajMom = NaN
	MinMom = NaN
	AngMom = NaN
	AspectRatio = NaN
	YSize = NaN
	XSIze = NaN
	BlueCalib = NaN
	RedCalib = NaN
																									///////////////////////////////////////////////////////////////////
//	Variable BrightnessThreshold = 50 // in (unsigned byte) grayscale intensity	//// Edit these to change thresholds of detection of particles,	//
//	Variable ParticleSizeThreshold = 25 // in pixels										//// thereby changing calibrations.										//
	Variable NumParticlesError = 1  		// Tracks successful particle counting;	///////////////////////////////////////////////////////////////////
       										// 1 is good, all other numbers are an error.													
	For (n=0; n < NumPnts(ImageWavePrefixes); n+=1)

		String Measured, BinMeas
		Measured = ImageWavePrefixes[n] + "_Calib"

		Wave w0=$Measured
         												
		ImageThreshold/Q/T=(BrightnessThreshold)/I w0						
		Wave M_ImageThresh                                             
		Duplicate/O M_ImageThresh, RedPointWave
         
		ImageThreshold/Q/T=(BrightnessThreshold)/I w0
		Duplicate/O M_ImageThresh, BluePointWave
         
		Redimension/N=(-1,-1) RedPointWave
		Redimension/N=(-1,-1) BluePointWave
           
//		Duplicate/O/R=(400,820)(0,104) RedPointWave, ParticleMeasure // (X)(Y) points
//		Duplicate/O/R=(400,820)(0,104) RedPointWave, RedMeasure
//		Duplicate/O/R=(0,350)(0,104) BluePointWave, BlueMeasure
        
 		Duplicate/O/R=(ceil(0.55*(DimSize(BluePointWave,0))), DimSize(BluePointWave,0)-1)(0, DimSize(BluePointWave,1)-1) BluePointWave, ParticleMeasure // Reformatting subsections of images to detect blue and red points for RaspPiSpec
		Duplicate/O/R=(0, floor(0.55*(DimSize(RedPointWave,0))))(0, DimSize(RedPointWave,1)-1) RedPointWave, RedMeasure											// 0-55% of image for red, 55%-100% for blue, use blue for particle measurement
		Duplicate/O/R=(ceil(0.55*(DimSize(BluePointWave,0))), DimSize(BluePointWave,0)-1)(0, DimSize(BluePointWave,1)-1) BluePointWave, BlueMeasure
		
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, BlueMeasure
		
		if(V_NumParticles > 1)
			print "Error: " + num2str(V_NumParticles) + " particles were found on blue side of " + Measured
			NumParticlesError = V_NumParticles
		elseif(V_NumParticles < 1)
			print "Error: no particles were found on blue side of " + Measured
			NumParticlesError = V_NumParticles
		endif
		
		Wave M_Moments
		BlueCalib[n] = ceil(0.55*(DimSize(BluePointWave,0))) + M_Moments[0]
		
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, RedMeasure
		
		if(V_NumParticles > 1)
			print "Error: " + num2str(V_NumParticles) + " particles were found on red side of " + Measured
			NumParticlesError = V_NumParticles
		elseif(V_NumParticles < 1)
			print "Error: no particles were found on red side of " + Measured
			NumParticlesError = V_NumParticles
		endif
		
		RedCalib[n] = M_Moments[0]
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, ParticleMeasure
        
		MajMom[n] = M_Moments[2]*2
		MinMom[n] = M_Moments[3]*2
          
		If (M_Moments[4]>pi)
			AngMom[n] = M_Moments[4]-pi
		Else
			AngMom[n] = M_Moments[4]
		EndIf
           

		AspectRatio[n] = MajMom[n]/MinMom[n]
     
		Wave W_YMax, W_YMin, W_XMax, W_XMin
		YSize[n] = W_YMax[0]-W_YMin[0]
		XSize[n] = W_XMax[0]-W_XMin[0]
        
		Degrees[n] = AngMom[n]*(180/pi)
	
//		CountingUp +=1
		
   EndFor
   
   DeleteNaNs(MajMom)
   DeleteNaNs(BlueCalib)
   DeleteNaNs(RedCalib)
   DeleteNaNs(MinMom)
   DeleteNaNs(AngMom)
   DeleteNaNs(AspectRatio)
   DeleteNaNs(YSize)
   DeleteNaNs(XSize)
   DeleteNaNs(Degrees)
	
	
	String SizeMeas450
	Make/o/n=1000 SizeWave
	SizeWave = NaN
	Variable SizeNormalizationError = 2
	// This for loop is used to get width of 450 spectral swath (i.e. particle size from 450 fluorescence); may not make sense for RaspPiSpec.
	For (n = 0; n < NumPnts(ImageWavePrefixes); n += 1)
		
		Make/O/N=2 DestWave
		Variable Mid
		
		Sizemeas450 = ImageWavePrefixes[n]+"_450nm"
		Wave MeasuredWave = $SizeMeas450
		String pSize, pSizeRef, SizeMeas2
		
		pSize = ImageWavePrefixes[n]+"pSize"		
			
		MatrixOP/O sizingwave = sumCols(MeasuredWave)/numRows(MeasuredWave) // maybe averages image in wrong dimension to get particle size--check later
		Make/O/N=(DimSize($Sizemeas450,1)) $pSize = SizingWave
		
		Wave CurveFitting = $pSize
		
		CurveFit/Q/NTHR=0 gauss CurveFitting /D
		String SizeMeas = "fit_" + pSize
		Wave Subtractor = $SizeMeas				
		Wavestats/Q Subtractor						// Measure the Guassians' characteristics	
		
		CurveFitting -= V_Min
		
		CurveFit/Q/NTHR=0 gauss CurveFitting /D
		SizeMeas2 = "fit_" + pSize
		Wave Subtractor = $SizeMeas2				
		Wavestats/Q Subtractor		
		
		Mid = V_max/2
		FindLevels/Q/D=destWave $SizeMeas, Mid
		if(V_LevelsFound > 2)
			SizeWave[n] = (DestWave[1]-DestWave[0])
			SizeNormalizationError = V_LevelsFound
		elseif (V_LevelsFound == 2)
			SizeWave[n] = (DestWave[1]-DestWave[0])
		else
			SizeWave[n] = 1
			SizeNormalizationError = V_LevelsFound
		endif
		
	
	EndFor
	
	DeleteNaNs(SizeWave)
	
	
	String NameStr, Calibrate
	
//	j = 1
	Make /O /N=1000 BluePts
	BluePts = NaN
	Make /O /N=1000 RedPts
	RedPts = NaN
	
	For (n = 0; n < NumPnts(ImageWavePrefixes); n += 1)					// Still need to try to recognize if the wave exists.
	
		NameStr = "p" + ImageWavePrefixes[n] + "_Calib"		// Create a string out of p[n]yCalib, n being the particle number
		Wave Calibrating = $NameStr				// Make the NameStr into a reference and make the Wave Calibrating that string
	
		
		if (WaveExists(Calibrating))
			CurveFit/Q/NTHR=0 gauss  Calibrating[ceil(0.55*(DimSize(Calibrating,0))), DimSize(Calibrating,0)-1] /D // Measure the Gaussian around the blue calib point
			Calibrate = "fit_" + NameStr
			Wave CalibWave = $Calibrate
			Wavestats/Q CalibWave						// Measure the Guassians' characteristics	
			BluePts[n] = V_MaxLoc					// Add the Max Location Calib point to the BluePts wave 
		
			CurveFit/Q/NTHR=0 gauss  Calibrating[0,floor(0.55*(DimSize(Calibrating,0)))]  /D // Measure the Gaussian around the red calib point
			Calib = "fit_" + NameStr
			Wave CalibWave = $Calib
			Wavestats/Q CalibWave
			RedPts[n] = V_MaxLoc
		EndIf
	EndFor
	

	DeleteNaNs(Bluepts)
	DeleteNaNs(Redpts)


	Variable c, e					
//	Wave BluePts, RedPts
	String NameStrRef
	
//	c = 632.8-403.457	// Using these wavelengths for BTFM
	c = 635-450			// Using these wavelengths for RaspPiSpec
	

	For (n = 0; n < DimSize(BluePts,0); n += 1)				
		
			NameStrRef = "p" + ImageWavePrefixes[n] + "xREF"
			
			Wave CalibRef = $NameStrRef

			e = c/(RedPts[n]-BluePts[n])
			
			CalibRef -= BlueCalib[n]		
			
			CalibRef *= e
		
			CalibRef += 450
			
	EndFor

	
	String NameStrY, NameStrX, NewSizeWaveX, NewSizeF
//	Wave BluePts, YSize, XSize, FSize

	//////////////////////////////////////////////////
	// Copy command for new sizing //
	//////////////////////////////////////////////////
	String NewSizeWave, NameStrXSize, NameStrFSize
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "_532nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_532nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_532nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_532nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "_450nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_450nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_450nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_450nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "_405nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_405nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_405nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_405nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)

		NameStrY = "p" + ImageWavePrefixes[n] + "_280nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_280nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_280nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_280nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	String PeakSW, PeakY, PeakX, PeakXRef, PeakF
//	Wave Bluepts, AngMom
	
	Make/O/N=(DimSize(AngMom, 0)) PeaksSW, PeaksY, PeaksX, PeaksF
	
	PeaksSW = NaN
	PeaksY = NaN
	PeaksX = NaN
	PeaksF = NaN
//	Variable up
//	Up = 0
	
	
//	Wave Degrees, YSize, XSize, FSize
		

//	Up = 0
	
	Make/O FactorSplitY, FactorSplitX, FSize

	FactorSplitY = NaN
	FactorSplitX = NaN
	FSize = NaN
		
	For	(n = 0; n < DimSize(Degrees, 0); n += 1)
		If (Degrees[n]<90)
			FactorSplitY[n] = Degrees[n]/90
			FactorSplitX[n]  = 1-FactorsplitY[n] 
		ElseIf (Degrees[n]>90)
			FactorSplitX[n] = Degrees[n]/180
			FactorSplitY[n]  = 1-FactorsplitX[n]
		EndIf
		
		FSize[n] = ((FactorSplitY[n]*YSize[n]))+((FactorSplitX[n]*XSize[n]))
	EndFor
	

	DeleteNaNs(FactorSplitX)
	DeleteNaNs(FactorSplitY)
	DeleteNaNs(FSize)
	
	////////////////////////////////////////////////
	// Display & Calibration Section //
	////////////////////////////////////////////////
	
	Display/K=1
//	j = 1
//	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
																	// Create a string out of p[n]yCalib, n being the particle number												
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_450nm_New"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_450nm_NewX"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_450nm_NewF"
			NameStrY = "p" + ImageWavePrefixes[n] + "_450nm"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			
						
			if (!waveexists($NameStrY))
				Break
			EndIf
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[3]
			DisplayY *=100000
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure450[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[3]
			DisplayNew *=100000	
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure450[n]
			
			Wavestats/Q DisplayNewX
			DisplayNewx -=V_min
			DisplayNewx /=PowerDensity[3]
			DisplayNewx *=100000
			DisplayNewx /= XSize[n]
			DisplayNewx /= Exposure450[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[3]
			DisplayNewF *=100000
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure450[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	
	Display/K=1
//	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)						// Still need to try to recognize if the wave exists.
																				// Create a string out of p[n]yCalib, n being the particle number
			NameStrY = "p" + ImageWavePrefixes[n] + "_405nm"
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_405nm_New"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_405nm_NewX"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_405nm_NewF"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[2]
			DisplayY *=100000
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure405[n]

			Wavestats/Q DisplayNew
			
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[2]
			DisplayNew *=100000	
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure405[n]	
			
			Wavestats/Q DisplayNewX
			
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[2]
			DisplayNewX *=100000
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure405[n]
			
			Wavestats/Q displaynewF
			
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[2]
			DisplayNewF *=100000
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure405[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	
//	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
		
			NameStrY = "p" + ImageWavePrefixes[n] + "_532nm"		// Create a string out of p[n]yCalib, n being the particle number
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_532nm_New"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_532nm_NewF"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_532nm_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[1]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure532[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[1]
			DisplayNew *=100000
			
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure532[n]	
			
			Wavestats/Q DisplayNewX
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[1]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure532[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[1]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure532[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	
//	n = 1
	For(n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
																		// Create a string out of p[n]yCalib, n being the particle number
			NameStrY = "p" + ImageWavePrefixes[n] + "_280nm"		// Create a string out of p[n]yCalib, n being the particle number
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_280nm_New"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_280nm_NewF"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_280nm_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNew = $NewSizeWave
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[0]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure280[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[0]
			DisplayNew *=100000
			
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure280[n]	
			
			Wavestats/Q DisplayNewX
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[0]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure280[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[0]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure280[n]
			
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"

//	Up = 0
	
	For (n = 0; n < DimSize(AngMom, 0); n+=1)

		PeakSW = "p" + ImageWavePrefixes[n] + "_450nm"
		PeakXRef = "p" + ImageWavePrefixes[n] + "xRef"
		PeakY = "p" + ImageWavePrefixes[n] + "_450nm_New"
		PeakX = "p" + ImageWavePrefixes[n] + "_450nm_NewX"
		PeakF = "p" + ImageWavePrefixes[n] + "_450nm_NewF"
		
		If (Waveexists($PeakSw)==0)
			Break
		EndIf
		
		Wavestats/Q $PeakSW
		PeaksSW[n] = V_Max
	
		Wavestats/Q $PeakY
		PeaksY[n] = V_Max
		
		Wavestats/Q $PeakX
		PeaksX[n] = V_Max
		
		Wavestats/Q $PeakF
		PeaksF[n] = V_max
		
//		Wave Degrees
		
		Make/O FactorSplitY, FactorSplitX
	
		FactorSplitY = NaN
		FactorSplitX = NaN
			
		If (Degrees[n]<90)
			FactorSplitY = Degrees/90
			FactorSplitX  = 1-FactorsplitY 
		ElseIf (Degrees[n]>90)
			FactorSplitY = 90/Degrees
			FactorSplitX  = 1-FactorsplitY 
		EndIf
//		up += 1
	EndFor
	
	DeleteNaNs(FactorSplitX)
	DeleteNaNs(FactorSplitY)
	DeleteNaNs(PeaksSW)
	DeleteNaNs(PeaksY)
	DeleteNaNs(PeaksX)
	
	Edit/K=1
	
	AppendToTable Majmom, Minmom, AspectRatio
	/////////////////////////// Final error messages. \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	If(NumParticlesError != 1)
		print "ERROR: ImageAnalyzeParticles did not detect the correct number of particles in one or more images. Results (especially particle sizes) may be skewed or incorrect."
	EndIf
	
	If(SizeNormalizationError != 2)
		print "ERROR: Particle fitting could not determine the size of one or more particles; normalization may be incorrect or skewed."
	EndIf
	
End

//**************************************************************************************************

Function MeasureIntensity()
// Created 12/07/2020 by MPPF
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: Measure the light intensity (brightness) of one or multiple images for comparison.
	


	/////////////////////////// Load image, work it up, and extract its name. /////////////////////////////////////////
	
	ImageLoad/O/Q/T=any/N=image	// Load the image into igor (any graphic file type)
	
	if (V_Flag)
	
		string name = StringFromList(0, S_waveNames)	// 
		wave ImageWave = $name
		
		String NoExtension = ParseFilepath(3, S_filename, ":", 0, 0)		// Parse the file's name and put it into S_Filename (option 3 removes file extension)									
		String TraceName
		sscanf NoExtension, "%s", TraceName
			
	endif
			
	ImageTransform/Q rgb2gray ImageWave						// Converts color to grayscale. Output is M_RGB2Gray.
	Wave M_RGB2Gray												// This is the wave output by ImageTransform for the rgb2gray operation.
	ImageStats M_RGB2Gray										// Computes stats on image. V_avg for average pixel value
	
	
	
	
	//////////////////////////////////////// Graph the intensity data. ////////////////////////////////////////////
	
	DoWindow IntensityGraph																													// Ask if IntensityGraph exists (i.e. if this function has been used recently).
	
	if(V_flag)																																	// If this is not the first time, go forth and append!
//		String IntensityGraphTraces = TraceNameList("IntensityGraph",";",1)															// Make a list of the trace names in IntensityGraph.
//		Variable LengthTracesList = ItemsInList(IntensityGraphTraces,";")															// Determine the length of the list of trace names.
//		String LastTraceName = stringfromlist(LengthTracesList-1, IntensityGraphTraces, ";")										// Pick out the most recent trace name.
//		
//		Variable TraceNumber
//		sscanf LastTraceName, "%*[HorizontalAverage]%d", TraceNumber																	// Extract the number suffix for the most recent trace name.
//		
//		TraceNumber += 1
//		String NewTraceName = "HorizontalAverage" + num2str(TraceNumber)																// Create the new trace name string with suffix iterated forward.
//		MatrixOP/O $NewTraceName = sumRows(M_RGB2Gray)/numCols(M_RGB2Gray)/255*100
//		AppendtoGraph $NewTraceName
		MatrixOP/O $TraceName = sumRows(M_RGB2Gray)/numCols(M_RGB2Gray)/255*100
		AppendtoGraph $TraceName
		Print "FILE NAME =", S_filename, "PLOTTED AS " + TraceName,  "AVERAGE = " + num2str(V_avg/255*100) + "%", "RMS = " + num2str(V_rms/255*100) +"%"
		
	else	
//		MatrixOP/O HorizontalAverage0 = sumRows(M_RGB2Gray)/numCols(M_RGB2Gray)/255*100																																		// If IntensityGraph doesn't exist, make it.
		MatrixOP/O $TraceName = sumRows(M_RGB2Gray)/numCols(M_RGB2Gray)/255*100
//		Display/N=IntensityGraph/K=0 HorizontalAverage0
		Display/N=IntensityGraph/K=0 $TraceName
		Label left "Intensity (% saturation)"
		Label bottom "Pixel Number"
		Legend
		Print "FILE NAME =", S_filename, "PLOTTED AS " + TraceName,  "AVERAGE = " + num2str(V_avg/255*100) + "%", "RMS = " + num2str(V_rms/255*100) +"%"
	endif
	
	QuickColorGraph()
	
End Function

//**************************************************************************************************

Function QuickColorGraph()
// Created 12/07/2020 by MPPF; idea found on WaveMetrics forum ("Automatically color traces on multi-trace graph") posted by user "chozo."
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: to quickly color the traces of the selected graph differently up to 21 unique traces.
// Colors are roughly (named by me): orange, blue, black, red, forest green, neon green, cyan, baby blue, purple, pink, magenta, yellow, gray, brown, 
// burgundy, lavender, hot pink, teal, khaki, cobalt, rusted yellow, tan.

	String Traces = TraceNameList("",";",1)			// These two lines determine the number of traces to be modified.
	Variable Items = ItemsinList(Traces,";")
	
	//Note the flag /N is unnecessary in the line below but specifies how many (rows,columns) I have.
	Make/FREE/N=(3,22) colors = { {65535,32764,16385}, {0,21845,65535}, {8738,8738,8738}, {65535,0,0}, {2,39321,1}, {0,65535,0}, {0,65535,65535}, {16385,49025,65535}, {26411,1,52428}, {65535,32768,45875}, {52428,1,41942}, {65535,65535,0}, {34952,34952,34952}, {26214,17479,0}, {39321,1,1}, {51664,44236,58982}, {65535,0,52428}, {16385,65535,41303}, {26214,26212,0}, {27242,42662,44204}, {52428,52425,1}, {54227,45232,36494} }
	Variable i,j
	for (i = 0; i < Items; i +=1)					// FOR loop iterates through all traces.
	
		j = i
		if(j >= DimSize(colors, 1))				// IF statement loops colors if there are more traces than available colors.
			j = mod(j,DimSize(colors,1))
		endif
		
		ModifyGraph rgb($StringFromList(i,Traces,";"))=(colors[0][j],colors[1][j],colors[2][j])			// This line does the actual modification.
		
	EndFor

End Function

//**************************************************************************************************

Function MakeModulairPMImagePlot([FileNamePath])
// Created 6/29/2021 by MPPF
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: To download raw Modulair-PM data and form an image plot from it.
	String FileNamePath
	
	if(ParamIsDefault(FileNamePath))
		FileNamePath = ""
	endif

// Load the downloaded Modulair-PM CSV into Igor as many 1D waves.
	LoadWave/J/D/A/W/O/Q/E=1/K=0/L={0,1,0,1,0}/R={English,2,2,2,2,"Year-Month-DayOfMonth",40} FileNamePath
	
//	Remove redundant lines of data or lines that are flagged.
	Wave timestamp, flag
	Variable NumWavesLoaded = V_flag
	Variable CheckIdx = 0, WaveToEdit = 0
	Do
		If((timestamp[CheckIdx] == timestamp[CheckIdx+1]) || (flag[CheckIdx] != 0))
			For(WaveToEdit = 0; WaveToEdit < NumWavesLoaded; WaveToEdit += 1)
				Wave TempWave = $StringFromList(WaveToEdit,S_WaveNames,";")
				Deletepoints CheckIdx, 1, TempWave
			EndFor
			CheckIdx -= 1
		ElseIf(timestamp[CheckIdx] < timestamp[CheckIdx+1])
			Print "Time is not monotonically increasing. Verify validity of downloaded data."
			Break
		EndIf
		
		CheckIdx += 1
	While(CheckIdx < numpnts(timestamp)-1)
	
// Convert 1D raw time-series size distribution waves into a single 2D wave. Units are particles per mL (averaged from twelve 5-second measurements over a minute).
	Wave bin0,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13,bin14,bin15,bin16,bin17,bin18,bin19,bin20,bin21,bin22,bin23
	Concatenate/O/DL {bin0,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13,bin14,bin15,bin16,bin17,bin18,bin19,bin20,bin21,bin22,bin23},ModPMTimeSizeDist
	
// Make a wave to bound the size bins recorded my the Modulair-PM isntrument. 
// Thus, 25 points in the wave bound each of 24 bins used by the instrument. Sizes in units of micrometers.
	Make/O/N=25/D SizeDistBounds = {0.35,0.46,0.66,1.0,1.3,1.7,2.3,3.0,4.0,5.2,6.5,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,25.0,28.0,31.0,34.0,37.0,40.0} 
	
// Normalize particle number concentrations by bin width.
	Make/N=(numpnts(SizeDistBounds))/O LogSizeDistBounds = log(SizeDistBounds)
	Variable i
	For(i = 0; i < DimSize(ModPMTimeSizeDist, 1); i += 1)
		ModPMTimeSizeDist[][i] /= (LogSizeDistBounds[i+1] - LogSizeDistBounds[i]) //The factor of 100 makes the result: number of particles per mL per 10 nm width of size bin. RECONSIDER THIS?
	EndFor
	
// Add an extra "dummy" point at the end of time waves to allow for image plot; this point is one minute earlier than data is recorded. 
// Units used are seconds since January 1st, 1904.
	Wave timestamp, timestamp_local
	Duplicate/O timestamp, timestamp_long
	Duplicate/O timestamp_local, timestamp_local_long
	InsertPoints numpnts(timestamp_long), 1, timestamp_long, timestamp_local_long
	timestamp_long[inf] = timestamp_long[numpnts(timestamp_long)-2] - 60
	timestamp_local_long[inf] = timestamp_local_long[numpnts(timestamp_long)-2] - 60
	

	
// Make the image plot itself using local time.
	DoWindow ModPMImagePlot
	If(V_flag)
		Print "ERROR: a Modulair-PM image plot window already exists."
	Else
		ReplaceValue(ModPMTimeSizeDist, 0, NaN)
		Display/K=1/N=ModPMImagePlot
		AppendImage/W=ModPMImagePlot ModPMTimeSizeDist vs {timestamp_local_long,SizeDistBounds}
		ModifyImage/W=ModPMImagePlot ModPMTimeSizeDist ctab= {*,*,Rainbow,1}
		Label/W=ModPMImagePlot left "Diameter (μm)"
		Label/W=ModPMImagePlot bottom "Time"
		ModifyGraph/W=ModPMImagePlot log(left)=1
		ColorScale/W=ModPMImagePlot/C/N=text0/A=MT/E vert=0,side=2,image=ModPMTimeSizeDist
		ColorScale/W=ModPMImagePlot/C/N=text0 "dN/dlogD\Bp\M (cm\\S-3\\M)"
	EndIf
	
	
	
End Function

//**************************************************************************************************

Function SizeFOBSParticles([Option,BrightnessThreshold,AreaThreshold])
// Created 05/02/2022 by MPPF
// Igor Version: Igor Pro 7.08 (Windows 10 64-bit)
//
// Purpose: To select either an image or a folder of images (Option = 1) of cropped particles from the FOBS
// in order to size the particle by fitting an ellipse to it. Each file/image should contain only one particle.
// The brightness and area thresholds for what defines a particle can be specified.
// Output is i x 3 wave called SizeInfo where i is number of images with single particles submitted (number in folder if Option = 1; otherwise i = 1)
// and the 3 columns are the major axes lengths, minor axes lengths, and ratio of major to minor axes, respectively.

	//* Declare necessary input variables *\\
	Variable Option, BrightnessThreshold, AreaThreshold
	
	//* Initialize defaults if variables not specified *\\
	If(ParamIsDefault(Option))
		Option = 0
	EndIf
	
	If(ParamIsDefault(BrightnessThreshold))
		BrightnessThreshold = 50
	EndIf
	
	If(ParamIsDefault(AreaThreshold))
		AreaThreshold = 10
	EndIf
	
	
	//* If Option = 1, evaluate for a folder of images *\\
	If(Option)
	
		NewPath/O/M="Choose a folder containing images." pathName				// Pop-up dialog to select folder containing images.
		String FileName, FileLoc
		FileLoc = IndexedFile(pathName, -1, "????")								// Index images in folder.
		FileLoc = SortList(FileLoc, ";", 16)										// Sort images in new indexed list.
		
		Variable n = 0
		do
			FileName = StringFromList(n, FileLoc, ";")		
			if (strlen(FileName) == 0)							
				break																		// 	Loop through all of the indexed images.
			endif
			
			ImageLoad/Q/O/T=any/N=ImportedParticle/P=pathName FileName			// Load one image.
			Wave ImportedParticle
			ImageTransform/Q rgb2gray ImportedParticle								// Convert color image to grayscale.
			Wave M_RGB2Gray
			ImageThreshold/Q/I/T=(BrightnessThreshold) M_RGB2Gray				// Create the threshold image based on BrightnessThreshold.
			Wave M_ImageThresh
			ImageAnalyzeParticles/E/Q/A=(AreaThreshold) stats M_ImageThresh	// Compute stats on particle based on AreaThreshold.
			Wave M_Moments
			
			//* If more than one particle in an image, report an error. If only one particle, save results in a temporary wave that will be concatenated into SizeInfo *\\
			If(V_NumParticles != 1)
				print("'Did not see exactly 1 particle for " + S_fileName + "! V_NumParticles = " + num2str(V_NumParticles) + "'")	// Report if more than one particle or no particles were found.
			Else
				Make/FREE/O/N=(1,3) TempSize = {{2*M_Moments[2]},{2*M_Moments[3]},{M_Moments[2]/M_Moments[3]}}
				Concatenate/NP=0 {TempSize}, SizeInfo
			EndIf
			
			n += 1
		while(1)
	
	//* If Option is not set to 1, default to evaluating for a single selected image *\\
	Else
	
		ImageLoad/Q/O/T=any/N=ImportedParticle ""			// Generate a dialog to select image and name resulting wave "ImportedParticle". The rest is same as above.
		Wave ImportedParticle
		ImageTransform/Q rgb2gray ImportedParticle
		Wave M_RGB2Gray
		ImageThreshold/Q/I/T=(BrightnessThreshold) M_RGB2Gray
		Wave M_ImageThresh
		ImageAnalyzeParticles/E/Q/A=(AreaThreshold) stats M_ImageThresh
		Wave M_Moments
		
		If(V_NumParticles != 1)
			print("'Did not see exactly 1 particle for " + S_fileName + "! V_NumParticles = " + num2str(V_NumParticles) + "'")
		Else
			Make/FREE/O/N=(1,3) TempSize = {{2*M_Moments[2]},{2*M_Moments[3]},{M_Moments[2]/M_Moments[3]}}
			Concatenate/NP=0 {TempSize}, SizeInfo
		EndIf
	
	EndIf
	
End Function

//**************************************************************************************************

Function RPCropAnalysis2([pathName, nameWaveAfterFileName, displayImages]) // ("", 1, 0)
 // IMPORT TIFFS FROM SUBFOLDER
 	String pathName
	Variable nameWaveAfterFileName, displayImages
	
	
	if (ParamIsDefault(pathName))											// Set default parameters to default values
		pathName = ""
	endif
	if (ParamIsDefault(nameWaveAfterFileName))
		nameWaveAfterFileName = 1
	endif
	if (ParamIsDefault(displayImages))
		displayImages = 0
	endif
	
//////////////// Initialize important variables. //////////////////////////

																									///////////////////////////////////////////////////////////////////
	Variable BrightnessThreshold = 20 // in (unsigned byte) grayscale intensity	//// Edit these to change thresholds of detection of particles,	//
	Variable ParticleSizeThreshold = 15 // in pixels										//// thereby changing calibrations.										//
																									///////////////////////////////////////////////////////////////////
	
	
/////////////// End variable declarations. ///////////////////////////////	
	
	if (strlen(pathName) == 0 )														// see if there's a given path/folder.  Quotations ("") returns 
		NewPath/O/M="Choose a folder containing TIFFs" pathName
		If (V_flag != 0)
			Return -1
		EndIf
	EndIf
	
	Make/O/N=1000 Exposure280			// Make a wave for exposure time for 280
	Make/O/N=1000 Exposure405 		// 405
	Make/O/N=1000 Exposure450			// 450
	Make/O/N=1000 Exposure532			// 532
	
	Exposure280 = NaN					// Make them all NaN
	Exposure405 = NaN
	Exposure450 = NaN
	Exposure532 = NaN
	
	
	String fileName, fileLoc
	Make/N=0/T/O ImageWavePrefixes//, ImageWaveNames
	Make/N=(0,4,3)/O ParticleSizes
	
	fileLoc = IndexedFile(pathName, -1, ".tif")
	fileLoc = SortList(FileLoc, ";", 16)
	
	Variable n = 0
	Variable LaserID
	
	do
		fileName = StringFromList(n, FileLoc, ";")		// find ".tif" files in the specified folder
		if (strlen(fileName) == 0)							// If a ".tif" image exists, keep going, if not, break the loop
			break										// break the loop
		endif
		
		ImageLoad/Q/O/P=pathName/T=TIFF/N=image fileName	// Load the image into igor, specifically Tiffs
		
		if (V_Flag)
			string name = StringFromList(0, S_waveNames)	// 
			wave w = $name
			
			if (nameWaveAfterFileName)
				string desiredName = S_fileName								// desiredName is equal to the full image name + ".tif"
				desiredName = ParseFilepath(3, S_filename, ":", 0, 0)				// Parse the file's name and put it into S_Filename (option 3 removes file extension)
						
						String NoTiff = desiredName								
										
//						String expr = "%s %e"									// Make expr a string EXPRESSION equal to "%s %e" // %s stores all text to the next white space, %e catches the number exp behind it
						String expr = "%s %*[#] %s %*s %*s %e %*[ISO] %e %*[sec]"		// This string EXPRESSION is for RaspPiSpec/FOBS
						String LaserExSt, SampleName
						Variable ExpTime, ParticlePrefix, Gain
						
//						sscanf NoTiff, Expr, LaserExSt, ExpTime				// Separate the string NoTiff by Expr (%s %e) and split it into LaserExSt string and ExpTime string
						sscanf NoTiff, Expr, SampleName, LaserExSt, Gain, ExpTime
//						sscanf NoTiff, "%d", ParticlePrefix
						
//						If (!StringMatch(ImageWavePrefixes[inf],num2str(ParticlePrefix)))
//							InsertPoints NumPnts(ImageWavePrefixes), 1, ImageWavePrefixes
//							InsertPoints NumPnts(ImageWaveNames), 1, ImageWaveNames
//							ImageWaveNames[inf] = LaserExSt	
//							ImageWavePrefixes[inf] = num2str(ParticlePrefix)
//						EndIf

						If (!StringMatch(ImageWavePrefixes[inf],SampleName))
							InsertPoints NumPnts(ImageWavePrefixes), 1, ImageWavePrefixes
//							InsertPoints NumPnts(ImageWaveNames), 1, ImageWaveNames
//							ImageWaveNames[inf] = LaserExSt	
							ImageWavePrefixes[inf] = SampleName
							InsertPoints inf, 1, ParticleSizes
						EndIf

//						String ExcValue280 = num2str(ParticlePrefix) + "y280"		// Make strings for each wave type per particle
//						String Excvalue532 = num2str(ParticlePrefix) + "y532"
//						String Excvalue405 = num2str(ParticlePrefix) + "y405"
//						String Excvalue450 = num2str(ParticlePrefix) + "y450"
//						String ExcValueCalib = num2Str(ParticlePrefix) + "yCalib"						
						
						If (Stringmatch(LaserExSt, "280nm"))				// Add the exposure time number to each exposure time wave.
							Exposure280[NumPnts(ImageWavePrefixes)-1] = ExpTime
							LaserID = 5
						ElseIf (Stringmatch(LaserExSt, "405nm"))
							Exposure405[NumPnts(ImageWavePrefixes)-1] = ExpTime
							LaserID = 0
						Elseif (Stringmatch(LaserExSt, "450nm"))
							Exposure450[NumPnts(ImageWavePrefixes)-1] = ExpTime
							LaserID = 1
						ElseIf (Stringmatch(LaserExSt, "532nm"))
							Exposure532[NumPnts(ImageWavePrefixes)-1] = ExpTime	
							LaserID = 2
						ElseIf (StringMatch(LaserExSt, "Calib"))
							LaserID = 3
						EndIf		
						
			
//						if (Exists(LaserExSt))								// If the desired name exists, this wave gets an identical name with n1, n2, n... starting from 0
//							LaserExSt = UniqueName(LaserExSt, 1, 0)
//						EndIf

						String WaveNameString = SampleName + "_" + LaserExSt		// SampleName could be anything. LaserExSt will have values {"280nm","405nm","450nm","532nm","Calib"} May also contain "incandescent".
						ImageTransform/Q rgb2gray w
						Wave M_RGB2Gray
						Duplicate/O M_RGB2Gray, $WaveNameString
						
//						Create intermediate mask wave for ImageAnalyzeParticles.
						ImageThreshold/Q/T=(BrightnessThreshold)/I M_RGB2Gray
						Wave M_ImageThresh
						Duplicate/O/RMD=[1500,*][] M_ImageThresh, M_ImageThreshCrop
						
						ImageAnalyzeParticles/E/Q/M=3/A=(ParticleSizeThreshold) stats, M_ImageThreshCrop
						Wave W_ImageObjArea, M_Moments
						If (LaserID != 5)
							If (DimSize(M_Moments, 0) != 0)
								ParticleSizes[inf][LaserID][0] = M_Moments[0][2]
								ParticleSizes[inf][LaserID][1] = M_Moments[0][3]
							Else
								ParticleSizes[inf][LaserID][0] = 0
								ParticleSizes[inf][LaserID][1] = 0
							EndIf
							ParticleSizes[inf][LaserID][2] = WaveMax(W_ImageObjArea)
						EndIf
															
				EndIf
			
			if (displayImages)
				NewImage w
			endif
		endif
		n += 1 
		
	while (1)
	
	ReplaceValue(ParticleSizes, NaN, 0)


 	DeleteNaNs(Exposure280) // Function defined in AlexTools_v16.5.
	DeleteNaNs(Exposure405)
	DeleteNaNs(Exposure450)
	DeleteNaNs(Exposure532)
	
 	// Calibrating Factors fix
// 	Exposure280 *=1000 			// Change seconds to milliseconds
// 	Exposure532 *=1000
// 	Exposure405 *=1000
// 	Exposure450 *=1000
 	
 	Make/O/N=4 PowerDensity			// These only apply to BTFM \\
 	PowerDensity[0] = 1//0.004		// PD factor correction for 280
 	PowerDensity[1] = 1//0.5			// PD factor correction for 350
 	PowerDensity[2] = 1//4.5			// PD factor correction for 405
 	PowerDensity[3] = 1//0.3			// PD factor correction for 450
 	
// 	Variable rn, t, r, y, j, b, o, h
	Variable h
	String wCalib, Calib, w450, w405, w532, w280, yCalib, y450, y405, y532, y280, pXRef, CurrentRef


//	b = 0
//	j = 1
//	o = 1
	

	For (n=0; n < NumPnts(ImageWavePrefixes); n+=1)
	
		yCalib = ImageWavePrefixes[n]+"_Calib"
		y532 = ImageWavePrefixes[n]+"_532nm"
		y450 = ImageWavePrefixes[n]+"_450nm"
		y405 = ImageWavePrefixes[n]+"_405nm"
		y280 = ImageWavePrefixes[n]+"_280nm"
	
		If (WaveExists($yCalib))
	
			Wave iCalib = $yCalib
			Wave i450 = $y450
			Wave i405 = $y405
			Wave i532 = $y532
			Wave i280 = $y280
		
			
			h = DimSize(iCalib,0)			
					
			Make/O/N=(h) XRef							// Make XRef wave equal to r (XRef is X dimension for spectra)
//			t = 1
			

			linspace("xRef", 1, DimSize($yCalib, 0))
				
			CurrentRef = "p"+ImageWavePrefixes[n]+"xRef"
			Make/O /N=(h) $CurrentRef = xRef
		

			wCalib = "p"+ImageWavePrefixes[n]+"_Calib"
			w450 = "p"+ImageWavePrefixes[n]+"_450nm"
			w405 = "p"+ImageWavePrefixes[n]+"_405nm"
			w532 = "p"+ImageWavePrefixes[n]+"_532nm"
			w280 = "p"+ImageWavePrefixes[n]+"_280nm"
		
			MatrixOP/O AveragesWave = sumRows(iCalib)/numCols(iCalib)
			Make/O/N=(DimSize($CurrentRef,0)) $wCalib = AveragesWave
					
			MatrixOP/O AveragesWave = sumRows(i450)/numCols(i450)
			Make/O/N=(DimSize($CurrentRef,0)) $w450 = AveragesWave

			MatrixOP/O AveragesWave = sumRows(i405)/numCols(i405)
			Make/O/N=(DimSize($CurrentRef,0)) $w405 = AveragesWave
	
			MatrixOP/O AveragesWave = sumRows(i532)/numCols(i532)
			Make/O/N=(DimSize($CurrentRef,0)) $w532 = AveragesWave

			MatrixOP/O AveragesWave = sumRows(i280)/numCols(i280)
			Make/O/N=(DimSize($CurrentRef,0)) $w280 = AveragesWave
			
//			b += 1			
			
		EndIf
	EndFor


	Wave Fit_CurveFitting
//	Variable l
//	n = 1
//	j = 1
	
//	Variable CountingUp, DelPts
// 	Wave MajMom, MinMom, AngMom, AspectRatio, YSize, Degrees, XSize
     	     	
//	CountingUp = 1
	Make/O/N=50 MajMom, MinMom, AngMom, AspectRatio, YSize, Degrees, XSize, BlueCalib, RedCalib
	    
	Degrees = NaN
	MajMom = NaN
	MinMom = NaN
	AngMom = NaN
	AspectRatio = NaN
	YSize = NaN
	XSIze = NaN
	BlueCalib = NaN
	RedCalib = NaN
																									///////////////////////////////////////////////////////////////////
//	Variable BrightnessThreshold = 50 // in (unsigned byte) grayscale intensity	//// Edit these to change thresholds of detection of particles,	//
//	Variable ParticleSizeThreshold = 25 // in pixels										//// thereby changing calibrations.										//
	Variable NumParticlesError = 1  		// Tracks successful particle counting;	///////////////////////////////////////////////////////////////////
       										// 1 is good, all other numbers are an error.													
	For (n=0; n < NumPnts(ImageWavePrefixes); n+=1)

		String Measured, BinMeas
		Measured = ImageWavePrefixes[n] + "_Calib"

		Wave w0=$Measured
         												
		ImageThreshold/Q/T=(BrightnessThreshold)/I w0						
		Wave M_ImageThresh                                             
		Duplicate/O M_ImageThresh, RedPointWave
         
		ImageThreshold/Q/T=(BrightnessThreshold)/I w0
		Duplicate/O M_ImageThresh, BluePointWave
         
		Redimension/N=(-1,-1) RedPointWave
		Redimension/N=(-1,-1) BluePointWave
           
//		Duplicate/O/R=(400,820)(0,104) RedPointWave, ParticleMeasure // (X)(Y) points
//		Duplicate/O/R=(400,820)(0,104) RedPointWave, RedMeasure
//		Duplicate/O/R=(0,350)(0,104) BluePointWave, BlueMeasure
        
 		Duplicate/O/R=(ceil(0.55*(DimSize(BluePointWave,0))), DimSize(BluePointWave,0)-1)(0, DimSize(BluePointWave,1)-1) BluePointWave, ParticleMeasure // Reformatting subsections of images to detect blue and red points for RaspPiSpec
		Duplicate/O/R=(0, floor(0.55*(DimSize(RedPointWave,0))))(0, DimSize(RedPointWave,1)-1) RedPointWave, RedMeasure											// 0-55% of image for red, 55%-100% for blue, use blue for particle measurement
		Duplicate/O/R=(ceil(0.55*(DimSize(BluePointWave,0))), DimSize(BluePointWave,0)-1)(0, DimSize(BluePointWave,1)-1) BluePointWave, BlueMeasure
		
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, BlueMeasure
		
		if(V_NumParticles > 1)
			print "Error: " + num2str(V_NumParticles) + " particles were found on blue side of " + Measured
			NumParticlesError = V_NumParticles
		elseif(V_NumParticles < 1)
			print "Error: no particles were found on blue side of " + Measured
			NumParticlesError = V_NumParticles
		endif
		
		Wave M_Moments
		BlueCalib[n] = ceil(0.55*(DimSize(BluePointWave,0))) + M_Moments[0]
		
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, RedMeasure
		
		if(V_NumParticles > 1)
			print "Error: " + num2str(V_NumParticles) + " particles were found on red side of " + Measured
			NumParticlesError = V_NumParticles
		elseif(V_NumParticles < 1)
			print "Error: no particles were found on red side of " + Measured
			NumParticlesError = V_NumParticles
		endif
		
		RedCalib[n] = M_Moments[0]
		ImageAnalyzeParticles /E/W/Q/M=3/A=(ParticleSizeThreshold) stats, ParticleMeasure
        
		MajMom[n] = M_Moments[2]*2
		MinMom[n] = M_Moments[3]*2
          
		If (M_Moments[4]>pi)
			AngMom[n] = M_Moments[4]-pi
		Else
			AngMom[n] = M_Moments[4]
		EndIf
           

		AspectRatio[n] = MajMom[n]/MinMom[n]
     
		Wave W_YMax, W_YMin, W_XMax, W_XMin
		YSize[n] = W_YMax[0]-W_YMin[0]
		XSize[n] = W_XMax[0]-W_XMin[0]
        
		Degrees[n] = AngMom[n]*(180/pi)
	
//		CountingUp +=1
		
   EndFor
   
   DeleteNaNs(MajMom)
   DeleteNaNs(BlueCalib)
   DeleteNaNs(RedCalib)
   DeleteNaNs(MinMom)
   DeleteNaNs(AngMom)
   DeleteNaNs(AspectRatio)
   DeleteNaNs(YSize)
   DeleteNaNs(XSize)
   DeleteNaNs(Degrees)
	
	
	String SizeMeas450
	Make/o/n=1000 SizeWave
	SizeWave = NaN
	Variable SizeNormalizationError = 2
	// This for loop is used to get width of 450 spectral swath (i.e. particle size from 450 fluorescence); may not make sense for RaspPiSpec.
	For (n = 0; n < NumPnts(ImageWavePrefixes); n += 1)
		
		Make/O/N=2 DestWave
		Variable Mid
		
		Sizemeas450 = ImageWavePrefixes[n]+"_450nm"
		Wave MeasuredWave = $SizeMeas450
		String pSize, pSizeRef, SizeMeas2
		
		pSize = ImageWavePrefixes[n]+"pSize"		
			
		MatrixOP/O sizingwave = sumCols(MeasuredWave)/numRows(MeasuredWave) // maybe averages image in wrong dimension to get particle size--check later
		Make/O/N=(DimSize($Sizemeas450,1)) $pSize = SizingWave
		
		Wave CurveFitting = $pSize
		
		CurveFit/Q/NTHR=0 gauss CurveFitting /D
		String SizeMeas = "fit_" + pSize
		Wave Subtractor = $SizeMeas				
		Wavestats/Q Subtractor						// Measure the Guassians' characteristics	
		
		CurveFitting -= V_Min
		
		CurveFit/Q/NTHR=0 gauss CurveFitting /D
		SizeMeas2 = "fit_" + pSize
		Wave Subtractor = $SizeMeas2				
		Wavestats/Q Subtractor		
		
		Mid = V_max/2
		FindLevels/Q/D=destWave $SizeMeas, Mid
		if(V_LevelsFound > 2)
			SizeWave[n] = (DestWave[1]-DestWave[0])
			SizeNormalizationError = V_LevelsFound
		elseif (V_LevelsFound == 2)
			SizeWave[n] = (DestWave[1]-DestWave[0])
		else
			SizeWave[n] = 1
			SizeNormalizationError = V_LevelsFound
		endif
		
	
	EndFor
	
	DeleteNaNs(SizeWave)
	
	
	String NameStr, Calibrate
	
//	j = 1
	Make /O /N=1000 BluePts
	BluePts = NaN
	Make /O /N=1000 RedPts
	RedPts = NaN
	
	For (n = 0; n < NumPnts(ImageWavePrefixes); n += 1)					// Still need to try to recognize if the wave exists.
	
		NameStr = "p" + ImageWavePrefixes[n] + "_Calib"		// Create a string out of p[n]yCalib, n being the particle number
		Wave Calibrating = $NameStr				// Make the NameStr into a reference and make the Wave Calibrating that string
	
		
		if (WaveExists(Calibrating))
			CurveFit/Q/NTHR=0 gauss  Calibrating[ceil(0.55*(DimSize(Calibrating,0))), DimSize(Calibrating,0)-1] /D // Measure the Gaussian around the blue calib point
			Calibrate = "fit_" + NameStr
			Wave CalibWave = $Calibrate
			Wavestats/Q CalibWave						// Measure the Guassians' characteristics	
			BluePts[n] = V_MaxLoc					// Add the Max Location Calib point to the BluePts wave 
		
			CurveFit/Q/NTHR=0 gauss  Calibrating[0,floor(0.55*(DimSize(Calibrating,0)))]  /D // Measure the Gaussian around the red calib point
			Calib = "fit_" + NameStr
			Wave CalibWave = $Calib
			Wavestats/Q CalibWave
			RedPts[n] = V_MaxLoc
		EndIf
	EndFor
	

	DeleteNaNs(Bluepts)
	DeleteNaNs(Redpts)


	Variable c, e					
//	Wave BluePts, RedPts
	String NameStrRef, CalibWaveName
	
//	c = 632.8-403.457	// Using these wavelengths for BTFM
//	c = 635-450			// Using these wavelengths for RaspPiSpec
	c = -450/768.7		// Using this factor when calibrating by distance from zeroth order image.
	

	For (n = 0; n < DimSize(BluePts,0); n += 1)				
		
			NameStrRef = "p" + ImageWavePrefixes[n] + "xREF"
			CalibWaveName = "p" + ImageWavePrefixes[n] + "_Calib"
			
			Wave CalibRef = $NameStrRef
			Wavestats/Q/R=[1300] $CalibWaveName
			CalibRef = (p-V_maxloc)*c
//			e = c/(RedPts[n]-BluePts[n])
//			
//			CalibRef -= BlueCalib[n]		
//			
//			CalibRef *= e
//		
//			CalibRef += 450
			
	EndFor

	
	String NameStrY, NameStrX, NewSizeWaveX, NewSizeF
//	Wave BluePts, YSize, XSize, FSize

	//////////////////////////////////////////////////
	// Copy command for new sizing //
	//////////////////////////////////////////////////
	String NewSizeWave, NameStrXSize, NameStrFSize
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "_532nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_532nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_532nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_532nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "_450nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_450nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_450nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_450nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)
	
		NameStrY = "p" + ImageWavePrefixes[n] + "_405nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_405nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_405nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_405nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	For (n = 0; n < DimSize(BluePts, 0); n += 1)

		NameStrY = "p" + ImageWavePrefixes[n] + "_280nm"
		NewSizeWave = "p" + ImageWavePrefixes[n] + "_280nm_New"
		NameStrXSize = "p" + ImageWavePrefixes[n] + "_280nm_NewX"
		NameStrFSize = "p" + ImageWavePrefixes[n] + "_280nm_NewF"
	
		Duplicate/O $NameStrY, $NewSizeWave
		Duplicate/O $NameStrY, $NameStrXSize
		Duplicate/O $NameStrY, $NameStrFSize

	EndFor
	
	String PeakSW, PeakY, PeakX, PeakXRef, PeakF
//	Wave Bluepts, AngMom
	
	Make/O/N=(DimSize(AngMom, 0)) PeaksSW, PeaksY, PeaksX, PeaksF
	
	PeaksSW = NaN
	PeaksY = NaN
	PeaksX = NaN
	PeaksF = NaN
//	Variable up
//	Up = 0
	
	
//	Wave Degrees, YSize, XSize, FSize
		

//	Up = 0
	
	Make/O FactorSplitY, FactorSplitX, FSize

	FactorSplitY = NaN
	FactorSplitX = NaN
	FSize = NaN
		
	For	(n = 0; n < DimSize(Degrees, 0); n += 1)
		If (Degrees[n]<90)
			FactorSplitY[n] = Degrees[n]/90
			FactorSplitX[n]  = 1-FactorsplitY[n] 
		ElseIf (Degrees[n]>90)
			FactorSplitX[n] = Degrees[n]/180
			FactorSplitY[n]  = 1-FactorsplitX[n]
		EndIf
		
		FSize[n] = ((FactorSplitY[n]*YSize[n]))+((FactorSplitX[n]*XSize[n]))
	EndFor
	

	DeleteNaNs(FactorSplitX)
	DeleteNaNs(FactorSplitY)
	DeleteNaNs(FSize)
	
	////////////////////////////////////////////////
	// Display & Calibration Section //
	////////////////////////////////////////////////
	
	Display/K=1
//	j = 1
//	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
																	// Create a string out of p[n]yCalib, n being the particle number												
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_450nm_New"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_450nm_NewX"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_450nm_NewF"
			NameStrY = "p" + ImageWavePrefixes[n] + "_450nm"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			
						
			if (!waveexists($NameStrY))
				Break
			EndIf
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[3]
			DisplayY *=100000
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure450[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[3]
			DisplayNew *=100000	
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure450[n]
			
			Wavestats/Q DisplayNewX
			DisplayNewx -=V_min
			DisplayNewx /=PowerDensity[3]
			DisplayNewx *=100000
			DisplayNewx /= XSize[n]
			DisplayNewx /= Exposure450[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[3]
			DisplayNewF *=100000
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure450[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	
	Display/K=1
//	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)						// Still need to try to recognize if the wave exists.
																				// Create a string out of p[n]yCalib, n being the particle number
			NameStrY = "p" + ImageWavePrefixes[n] + "_405nm"
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_405nm_New"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_405nm_NewX"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_405nm_NewF"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[2]
			DisplayY *=100000
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure405[n]

			Wavestats/Q DisplayNew
			
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[2]
			DisplayNew *=100000	
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure405[n]	
			
			Wavestats/Q DisplayNewX
			
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[2]
			DisplayNewX *=100000
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure405[n]
			
			Wavestats/Q displaynewF
			
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[2]
			DisplayNewF *=100000
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure405[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	
//	n = 1
	For (n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
		
			NameStrY = "p" + ImageWavePrefixes[n] + "_532nm"		// Create a string out of p[n]yCalib, n being the particle number
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_532nm_New"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_532nm_NewF"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_532nm_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayNew = $NewSizeWave
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[1]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure532[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[1]
			DisplayNew *=100000
			
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure532[n]	
			
			Wavestats/Q DisplayNewX
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[1]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure532[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[1]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure532[n]
		
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"
	
	Display/K=1
	
//	n = 1
	For(n = 0; n < DimSize(BluePts,0); n += 1)					// Still need to try to recognize if the wave exists.
																		// Create a string out of p[n]yCalib, n being the particle number
			NameStrY = "p" + ImageWavePrefixes[n] + "_280nm"		// Create a string out of p[n]yCalib, n being the particle number
			NewSizeWave = "p" + ImageWavePrefixes[n] + "_280nm_New"
			NewSizeF = "p" + ImageWavePrefixes[n] + "_280nm_NewF"
			NewSizeWaveX = "p" + ImageWavePrefixes[n] + "_280nm_NewX"
			NameStrX = "p" + ImageWavePrefixes[n] + "xRef"
			
			Wave DisplayY = $NameStrY	
			Wave DisplayX = $NameStrX
			Wave DIsplaynewx = $NewSizeWaveX
			Wave DisplayNew = $NewSizeWave
			Wave DisplayNewF = $NewSizeF
			
			//Quick baseline subtraction
			Wavestats/Q DisplayY
			DisplayY -=V_min
			DisplayY /=PowerDensity[0]
			DisplayY *=100000
			
			DisplayY /= SizeWave[n]
			DisplayY /= Exposure280[n]
			
			Wavestats/Q DisplayNew
			DisplayNew -=V_min
			DisplayNew /=PowerDensity[0]
			DisplayNew *=100000
			
			DisplayNew /= YSize[n]
			DisplayNew /= Exposure280[n]	
			
			Wavestats/Q DisplayNewX
			DisplayNewX -=V_min
			DisplayNewX /=PowerDensity[0]
			DisplayNewX *=100000
			
			DisplayNewX /= XSize[n]
			DisplayNewX /= Exposure280[n]
			
			Wavestats/Q DisplayNewF
			DisplayNewF -=V_min
			DisplayNewF /=PowerDensity[0]
			DisplayNewF *=100000
			
			DisplayNewF /= FSize[n]
			DisplayNewF /= Exposure280[n]
			
			AppendToGraph DisplayY vs DisplayX
	EndFor
	Label left "Intensity (a.u.)"
	Label bottom "Wavelength (nm)"

//	Up = 0
	
	For (n = 0; n < DimSize(AngMom, 0); n+=1)

		PeakSW = "p" + ImageWavePrefixes[n] + "_450nm"
		PeakXRef = "p" + ImageWavePrefixes[n] + "xRef"
		PeakY = "p" + ImageWavePrefixes[n] + "_450nm_New"
		PeakX = "p" + ImageWavePrefixes[n] + "_450nm_NewX"
		PeakF = "p" + ImageWavePrefixes[n] + "_450nm_NewF"
		
		If (Waveexists($PeakSw)==0)
			Break
		EndIf
		
		Wavestats/Q $PeakSW
		PeaksSW[n] = V_Max
	
		Wavestats/Q $PeakY
		PeaksY[n] = V_Max
		
		Wavestats/Q $PeakX
		PeaksX[n] = V_Max
		
		Wavestats/Q $PeakF
		PeaksF[n] = V_max
		
//		Wave Degrees
		
		Make/O FactorSplitY, FactorSplitX
	
		FactorSplitY = NaN
		FactorSplitX = NaN
			
		If (Degrees[n]<90)
			FactorSplitY = Degrees/90
			FactorSplitX  = 1-FactorsplitY 
		ElseIf (Degrees[n]>90)
			FactorSplitY = 90/Degrees
			FactorSplitX  = 1-FactorsplitY 
		EndIf
//		up += 1
	EndFor
	
	DeleteNaNs(FactorSplitX)
	DeleteNaNs(FactorSplitY)
	DeleteNaNs(PeaksSW)
	DeleteNaNs(PeaksY)
	DeleteNaNs(PeaksX)
	
	Edit/K=1
	
	AppendToTable Majmom, Minmom, AspectRatio
	/////////////////////////// Final error messages. \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	If(NumParticlesError != 1)
		print "ERROR: ImageAnalyzeParticles did not detect the correct number of particles in one or more images. Results (especially particle sizes) may be skewed or incorrect."
	EndIf
	
	If(SizeNormalizationError != 2)
		print "ERROR: Particle fitting could not determine the size of one or more particles; normalization may be incorrect or skewed."
	EndIf
	
End

//**************************************************************************************************

Function RPDCF() // ()
//Created 2025-05-13 by MPPF.
//Designed to interpret cropped images from VOBS focused on just RBG intensities that are not diffracted/dispersed.
 	String pathName
	Variable nameWaveAfterFileName, displayImages

	pathName = ""
	nameWaveAfterFileName = 1
	displayImages = 0
	
//	These will be the key output waves. 
//	ParticleNames will simply be used with NumParticles for debugging.
	Make/T/O/N=0 ParticleNames
	Make/O/N=0 NumParticles
//	ParticleSizes will be the lengths of major and minor axes from elliptical fit to particle in the layers dimension; the third layer will be pixel area of each particle.
//	Each particle is a new row. The 5 columns will be the 5 images (excitation source) as for ParticleColors with the calibration image as the last column (280->Calib).
	Make/O/N=(0,5,3) ParticleSizes
//	ParticleColors will be 3D, where first dimension is particles, second dimension is excitation source, third dimension is RGB.
	Make/O/N=(0,4,3) ParticleColors
//	ColorStats will store the average RGB intensities for each excitation; RGB will be rows, excitations as columns. Standard deviations will be second layer.
	Make/O/N=(3,4,2) ColorStats
//	SizeStats will store average and standard deviation (in layers) of: major axis, minor axis, total pixel areas from ParticleSizes (in rows). 
// 5 columns for each column of ParticleSizes.
	Make/O/N=(3,5,2) SizeStats
		
//////////////// Initialize important variables. //////////////////////////

																									///////////////////////////////////////////////////////////////////
	Variable BrightnessThreshold = 30 // in (unsigned byte) grayscale intensity	//// Edit these to change thresholds of detection of particles,	//
	Variable ParticleSizeThreshold = 15 // in pixels										//// thereby changing calibrations.										//																							///////////////////////////////////////////////////////////////////
																									///////////////////////////////////////////////////////////////////
	
/////////////// End variable declarations. ///////////////////////////////	
	
	if (strlen(pathName) == 0 )														// see if there's a given path/folder.  Quotations ("") returns 
		NewPath/O/M="Choose a folder containing TIFFs" pathName
		If (V_flag != 0)
			Return -1
		EndIf
	EndIf
	
	Make/O/N=0 Exposure280, Exposure405, Exposure450, Exposure532 = NaN
	Make/T/O/N=4 ExcitationWavelengths = {"280","405","450","532"}
	Variable ExposureVar
	
	
	String fileName, fileLoc
	Make/N=0/T/O ImageWavePrefixes//, ImageWaveNames
	
	fileLoc = IndexedFile(pathName, -1, ".tif")
	fileLoc = SortList(FileLoc, ";", 16)
	
	

	Variable n = 0	
	do
		fileName = StringFromList(n, FileLoc, ";")		// find ".tif" files in the specified folder
		if (strlen(fileName) == 0)							// If a ".tif" image exists, keep going, if not, break the loop
			break										// break the loop
		endif
		
		ImageLoad/Q/O/P=pathName/T=TIFF/N=image fileName	// Load the image into igor, specifically Tiffs
		
		if (V_Flag)
			string name = StringFromList(0, S_waveNames)	// 
			wave w = $name
			
			if (nameWaveAfterFileName)
				string desiredName = S_fileName								// desiredName is equal to the full image name + ".tif"
				desiredName = ParseFilepath(3, S_filename, ":", 0, 0)				// Parse the file's name and put it into S_Filename (option 3 removes file extension)
						
						String NoTiff = desiredName								
										
						String expr = "%s %*[#] %s %*s %*s %e %*[ISO] %e %*[sec]"		// This string EXPRESSION is for RaspPiSpec/FOBS. Make expr a string EXPRESSION equal to "%s %e" // %s stores all text to the next white space, %e catches the number exp behind it
						String LaserExSt, SampleName
						Variable ExpTime, ParticlePrefix, Gain
						
						sscanf NoTiff, Expr, SampleName, LaserExSt, Gain, ExpTime
						

						If (!StringMatch(ImageWavePrefixes[inf],SampleName))
							InsertPoints NumPnts(ImageWavePrefixes), 1, ImageWavePrefixes
							ImageWavePrefixes[inf] = SampleName
							InsertPoints inf, 1, ParticleColors
							InsertPoints inf, 1, ParticleSizes
						EndIf

						
//						Add the exposure time number to each exposure time wave and determine LaserID (which excitation source we are looking at).
						Variable LaserID
						If (Stringmatch(LaserExSt, "280nm"))
							LaserID = 0				
							InsertPoints inf, 1, Exposure280
							Exposure280[NumPnts(ImageWavePrefixes)-1] = ExpTime
						ElseIf (Stringmatch(LaserExSt, "405nm"))
							LaserID = 1
							InsertPoints inf, 1, Exposure405
							Exposure405[NumPnts(ImageWavePrefixes)-1] = ExpTime
						Elseif (Stringmatch(LaserExSt, "450nm"))
							LaserID = 2
							InsertPoints inf, 1, Exposure450
							Exposure450[NumPnts(ImageWavePrefixes)-1] = ExpTime
						ElseIf (Stringmatch(LaserExSt, "532nm"))
							LaserID = 3
							InsertPoints inf, 1, Exposure532
							Exposure532[NumPnts(ImageWavePrefixes)-1] = ExpTime	
						ElseIf (Stringmatch(LaserExSt, "Calib"))
							LaserID = 4
						EndIf		
						

//						SampleName could be anything. LaserExSt will have values {"280nm","405nm","450nm","532nm","Calib"} May also contain "incandescent".
						String WaveNameStringGray = SampleName + "_" + LaserExSt + "_Gray"
						String WaveNameString = SampleName + "_" + LaserExSt		
						InsertPoints inf, 1, ParticleNames
						ParticleNames[inf] = WaveNameString
						
//						Create a grayscale version  of image to facilitate masking and sizing.
						ImageTransform/Q rgb2gray w
						Wave M_RGB2Gray
						
//						Create intermediate mask wave for ImageAnalyzeParticles.
						ImageThreshold/Q/T=(BrightnessThreshold)/I M_RGB2Gray
						Wave M_ImageThresh
						
						Variable color

//						If (LaserID == 0)
//							InsertPoints inf, 1, ParticleColors
//							InsertPoints inf, 1, ParticleSizes
//						EndIf
						
						ImageAnalyzeParticles/E/Q/M=3/A=(ParticleSizeThreshold) stats, M_ImageThresh
						Wave W_ImageObjArea, M_Moments
						If (DimSize(M_Moments, 0) != 0)
							ParticleSizes[inf][LaserID][0] = M_Moments[0][2]
							ParticleSizes[inf][LaserID][1] = M_Moments[0][3]
						Else
							ParticleSizes[inf][LaserID][0] = 0
							ParticleSizes[inf][LaserID][1] = 0
						EndIf
						ParticleSizes[inf][LaserID][2] = WaveMax(W_ImageObjArea)


						If (LaserID != 4)
						
							String ExposureString = "Exposure" + ExcitationWavelengths[LaserID]
							Wave ExposureWave = $ExposureString
							ExposureVar = ExposureWave[inf]
							
							For (color = 0; color < 3; color += 1)
								Duplicate/O/R=[][][color] w, w_color
								ImageAnalyzeParticles/E/Q/M=3/A=(ParticleSizeThreshold)/D=w_color stats, M_ImageThresh
								Wave W_IntAvg
								ParticleColors[inf][LaserID][color] = WaveMax(W_IntAvg)/ExposureVar
							EndFor
							
						EndIf

						
						InsertPoints inf, 1, NumParticles
						NumParticles[inf] = V_numparticles
						
						
			
//						if(V_NumParticles > 1)
//							print "Error: " + num2str(V_NumParticles) + " particles were found on blue side of " + Measured
//							NumParticlesError = V_NumParticles
//						elseif(V_NumParticles < 1)
//							print "Error: no particles were found on blue side of " + Measured
//							NumParticlesError = V_NumParticles
//						endif

						
//						Rename M_RGB2Gray, $WaveNameStringGray
//						Rename w, $WaveNameString
						
																
				EndIf
			
			if (displayImages)
				NewImage w
			endif
		endif
		n += 1
	while (1)
	
	ReplaceValue(ParticleColors, Nan, 0)
	ReplaceValue(ParticleSizes, Nan, 0)
	
//	ParticleColorsTotals will contain total intensities instead of averages--this effectively pixel area of particle multiplied by averages.
	Duplicate/O ParticleColors, ParticleColorsTotals
	ParticleColorsTotals[][][] = ParticleColors[p][q][r]*ParticleSizes[p][q][3]
	
//	Save average and standard deviations of 12 pieces of information for each particle into ColorStats.
	Variable RGB, Excitations
	For (RGB = 0; RGB < 3; RGB += 1)
		For (Excitations = 0; Excitations < 4; Excitations += 1)
			WaveStats/Q/RMD=[][Excitations][RGB] ParticleColors
			ColorStats[RGB][Excitations][0] = V_avg
			ColorStats[RGB][Excitations][1] = V_sdev
		EndFor
	EndFor
	
//	Collecting SizeStats information and FinalSizes.
	Make/O/N=(Dimsize(ParticleSizes, 0), 3) FinalSizes
	Variable PartIdx
	For (PartIdx = 0; PartIdx < DimSize(ParticleSizes, 0); PartIdx += 1)
		WaveStats/Q/RMD=[PartIdx][1,3][2] ParticleSizes
		FinalSizes[PartIdx][] =  ParticleSizes[PartIdx][V_maxColLoc][q]
	EndFor
	
	Variable SizeMetric
	For (Excitations = 0; Excitations < 5; Excitations += 1)
		For (SizeMetric = 0; SizeMetric < 3; SizeMetric += 1)
			WaveStats/Q/RMD=[][Excitations][SizeMetric] ParticleSizes
			SizeStats[SizeMetric][Excitations][0] = V_avg
			SizeStats[SizeMetric][Excitations][1] = V_sdev
		EndFor
	EndFor
//	WaveStats/Q/RMD=[][1] ParticleSizes
//	SizeStats[1][0] = V_avg
//	SizeStats[1][1] = V_sdev
//	WaveStats/Q/RMD=[][2] ParticleSizes
//	SizeStats[2][0] = V_avg
//	SizeStats[2][1] = V_sdev
	
End Function
	
//**************************************************************************************************

Function PlotRPDCF([ParentFolder, GraphType])
//Created by MPPF; 2025-05-13.
//Function to plot data produced by RPDiscreteAnalysis() assuming that data for different particle types is stored
//in multiple subfolders stored in a parent folder and parent folder can't be "root."
//GraphType can be any of {0,1,2} corresponding to {averages of average intensities, sets of average intensities, sets of total intensities}.

	String ParentFolder
	Variable GraphType
	
//	Initialize input parameters. Full path to parent folder of folders containing data.
	If(ParamIsDefault(ParentFolder))
		ParentFolder = GetDataFolder(1)
		ParentFolder = ParentFolder[5,StrLen(ParentFolder)-2]
	EndIf
	
	If (ParamIsDefault(GraphType))
		GraphType = 1
	EndIf
	
	Make/T/N=3/O RGBLabels = {"Red", "Green", "Blue"}, SizeLabels = {"Major Axis", "Minor Axis", "Pixel Area"}
	Make/N=3/O LabelLocs = {0.5,1.5,2.5}
	
//	Create graph windows; overwrite if they already exist.
	Make/T/O GraphNames = {"Graph405", "Graph450", "Graph532", "GraphSize"}
	Variable GraphName
	For (GraphName = 0; GraphName < NumPnts(GraphNames); GraphName += 1)
		DoWindow/F $GraphNames[GraphName]
		If (V_flag)
			KillWindow $GraphNames[GraphName]
		EndIf
		Display/K=1/N=$GraphNames[GraphName]/I/W=(4*GraphName,GraphName,4*GraphName+7,GraphName+4)
	EndFor
	
	SetRandomSeed 0.999
	Variable NumFolders = CountObjects("root:" + ParentFolder + ":", 4)
	Variable i, j, r, g, b
	String Folder
	For (i = 0; i < NumFolders; i += 1)
		If (i == 0)
			Folder = StringFromList(1, StringFromList(i, DataFolderDir(1, $ParentFolder), ","),":")
		Else
			Folder = StringFromList(i, DataFolderDir(1, $ParentFolder), ",")
			If(i == (NumFolders-1))
				Folder = Folder[0,StrLen(Folder)-3]
			EndIf
		EndIf
		
		String TraceName405 = "Trace405" + Folder
		String TraceName450 = "Trace450" + Folder
		String TraceName532 = "Trace532" + Folder
		String TraceNameSize = "TraceSize" + Folder
		
		Make/O/T/N=4 TraceNames = {TraceName405, TraceName450, TraceName532, TraceNameSize}
		String ScatterXName = "ScatterX" + Folder
		Make/O/N=(3, DimSize(root:$(ParentFolder):$(Folder):ParticleColors, 0) ) $ScatterXName = NaN
		
		r=enoise(65535/2)+65535/2
		g=enoise(65535/2)+65535/2
		b=enoise(65535/2)+65535/2	
		
		
		For (GraphName = 0; GraphName < NumPnts(GraphNames); GraphName += 1)
			If(GraphType == 0)
				If (StringMatch(GraphNames[GraphName], "GraphSize"))
					AppendToGraph/W=$GraphNames[GraphName]/C=(r,g,b) root:$(ParentFolder):$(Folder):SizeStats[][4][0]/TN=$TraceNames[GraphName] vs SizeLabels
					ErrorBars/W=$GraphNames[GraphName] $TraceNames[GraphName] Y, wave=(root:$(ParentFolder):$(Folder):SizeStats[][1], root:$(ParentFolder):$(Folder):SizeStats[][4][1])
				Else
					AppendToGraph/W=$GraphNames[GraphName]/C=(r,g,b) root:$(ParentFolder):$(Folder):ColorStats[][GraphName+1][0]/TN=$TraceNames[GraphName] vs RGBLabels
					ErrorBars/W=$GraphNames[GraphName] $TraceNames[GraphName] Y, wave=(root:$(ParentFolder):$(Folder):ColorStats[*][GraphName+1][1], root:$(ParentFolder):$(Folder):ColorStats[*][GraphName+1][1])
				EndIf
				ModifyGraph/W=$GraphNames[GraphName] mode=3, marker=19
			ElseIf ( (GraphType == 1) || (GraphType == 2) )
				Wave ScatterX = $ScatterXName
				For (j = 0; j < 3; j += 1)
					If (GraphName == 3)
						If (j == 2)
							AppendToGraph/W=GraphSize/C=(r,g,b)/R root:$(ParentFolder):$(Folder):ParticleSizes[][4][j]/TN=$TraceNameSize vs ScatterX[j][]
//							ModifyGraph/W=$GraphNames[GraphName] mode=3, marker=16, msize = 1.5, userticks(bottom)={LabelLocs,SizeLabels}
//							Continue
						Else
							AppendToGraph/W=GraphSize/C=(r,g,b) root:$(ParentFolder):$(Folder):ParticleSizes[][4][j]/TN=$TraceNameSize vs ScatterX[j][]
						EndIf
						ModifyGraph/W=$GraphNames[GraphName] mode=3, marker=19, msize = 1.5, userticks(bottom)={LabelLocs,SizeLabels}
					Else
						ScatterX[j][] = j + ((i+1)/(NumFolders+1))
						If (GraphType == 1)
							AppendToGraph/W=$GraphNames[GraphName]/C=(r,g,b) root:$(ParentFolder):$(Folder):ParticleColors[][GraphName+1][j]/TN=$TraceNames[GraphName] vs ScatterX[j][]
						ElseIf (GraphType == 2)
							AppendToGraph/W=$GraphNames[GraphName]/C=(r,g,b) root:$(ParentFolder):$(Folder):ParticleColorsTotals[][GraphName+1][j]/TN=$TraceNames[GraphName] vs ScatterX[j][]
						EndIf
						ModifyGraph/W=$GraphNames[GraphName] mode=3, marker=19, msize = 1.5, userticks(bottom)={LabelLocs,RGBLabels}
					EndIf
				EndFor
			EndIf
		EndFor
	EndFor
	
	For (GraphName = 0; GraphName < NumPnts(GraphNames); GraphName += 1)
		If (GraphName == 3)
			SetAxis/W=$GraphNames[GraphName]/N=2 left 0,*
			SetAxis/W=$GraphNames[GraphName] bottom 0,3
			SetAxis/W=$GraphNames[GraphName]/N=2/Z right 0,*
		Else
			SetAxis/W=$GraphNames[GraphName] left 0,255
			SetAxis/W=$GraphNames[GraphName] bottom 0,3
			SetAxis/W=$GraphNames[GraphName]/Z right 0,*
			If (GraphType == 2)
				Label/W=$GraphNames[GraphName] left "Total Intensity (a.u)"
			Else
				Label/W=$GraphNames[GraphName] left "Average Intensity (a.u)"
			EndIf
		EndIf	
		Legend/W=$GraphNames[GraphName]/C/N=text0/F=0/A=RT/E=2
		SetDrawEnv/W=$GraphNames[GraphName] fillbgc=(0,0,0), fillfgc=(0,0,0), save
		DrawRRect/W=$GraphNames[GraphName] 0.332333333,0,0.334333333,1
		DrawRRect/W=$GraphNames[GraphName] 0.66566666,0,0.667666666,1

	EndFor
	
End Function

//**************************************************************************************************

Function pract([ParentFolder])
//Created by MPPF; 2025-05-13.
//Function to plot data produced by RPDiscreteAnalysis() assuming that data for different particle types is stored
//in multiple subfolders stored in a parent folder and parent folder can't be "root."

	String ParentFolder
	
//	Initialize input parameters. Full path to parent folder of folders containing data.
//	GraphType is an integer: 0 for averages + standard deviations, 1 for all data points.
	If(ParamIsDefault(ParentFolder))
		ParentFolder = GetDataFolder(1)
		ParentFolder = ParentFolder[5,StrLen(ParentFolder)-2]
	EndIf
	
	SetRandomSeed 0.999
	Variable NumFolders = CountObjects("root:" + ParentFolder + ":", 4)
	Variable i, j, r, g, b, GraphName
	String Folder
	Make/O/T/N=3 RGBLabels = {"Red","Green","Blue"}, Excitations = {"405","450","532"}
	For (i = 0; i < NumFolders; i += 1)

		If (i == 0)
			Folder = StringFromList(1, StringFromList(i, DataFolderDir(1, $ParentFolder), ","),":")
		Else
			Folder = StringFromList(i, DataFolderDir(1, $ParentFolder), ",")
			If(i == (NumFolders-1))
				Folder = Folder[0,StrLen(Folder)-3]
			EndIf
		EndIf
		
		DoWindow/F $Folder
		If (V_flag)
			KillWindow $Folder
		EndIf
		Display/K=1/N=$Folder/I/W=(4*i,i,4*i+7,i+4)
		
		r=enoise(65535/2)+65535/2
		g=enoise(65535/2)+65535/2
		b=enoise(65535/2)+65535/2	
		
		For (GraphName = 1; GraphName < 4; GraphName += 1)
			For (j=0;j<3;j+=1)
//				Display/K=1
				String TraceNameTotals = "TraceTotals" + Folder + RGBLabels[j] + Excitations[GraphName-1], LAxis = "L" + Folder + RGBLabels[j] + Excitations[GraphName-1], BAxis ="B" + Folder + RGBLabels[j] + Excitations[GraphName-1]
				
//				r=enoise(65535/2)+65535/2
//				g=enoise(65535/2)+65535/2
//				b=enoise(65535/2)+65535/2	
				
				AppendToGraph/L=$LAxis/B=$BAxis/C=(r,g,b) root:$(ParentFolder):$(Folder):ParticleColors[][GraphName][j]/TN=$TraceNameTotals vs root:$(ParentFolder):$(Folder):ParticleColorsTotals[][GraphName][j]
				ModifyGraph mode=3, marker=19, msize = 3, axisEnab($LAxis)={1-( (j+1)/3 )+0.07, 1-(j/3)-0.07}, axisEnab($BAxis)={0.07+(GraphName-1)/3, -0.07+GraphName/3}, freePos($LAxis) = {0,$BAxis}, freePos($BAxis) = {0,$LAxis}
			EndFor
		EndFor
		
		
	EndFor
	
//	For (GraphName = 0; GraphName < NumPnts(GraphNames); GraphName += 1)
//		SetAxis/W=$GraphNames[GraphName] left 0,*
//		SetAxis/W=$GraphNames[GraphName] bottom 0,3
//		SetAxis/W=$GraphNames[GraphName] right 0,*
//		Legend/W=$GraphNames[GraphName]/C/N=text0/F=0/A=RT/E=2
//		Label left "Intensity (a.u)"
//		SetDrawEnv/W=$GraphNames[GraphName] fillbgc=(0,0,0), fillfgc=(0,0,0), save
//		DrawRRect/W=$GraphNames[GraphName] 0.332333333,0,0.334333333,1
//		DrawRRect/W=$GraphNames[GraphName] 0.66566666,0,0.667666666,1
//	EndFor
	
End Function