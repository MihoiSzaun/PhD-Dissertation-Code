#pragma rtGlobals=1		// Use modern global access method.
//
// AlexTools_v14.ipf
// 	
// Created by J. Alex Huffman
//	University of Denver
//	Last updated: Nov. 24, 2013

// v13: Added AverageSizeDistMx 2/5/2013
// v14: Added AddExpinfo2Graph() 11/20/2013
//        AddedMakeDiurnalMx 11/24/2013
// v16.1: Changes to AverageSizeDistMx() to deal with NaNs 7/10/2015
// v16.5: Various changes prompted by Cyprus analysis, Feb. 2019
//		e.g. CollectBatchAveragesManyWaves ..

//***************************************************************************
//***************************************************************************

Function  Regrid()
       Wave xwave1 = bg_time
       wave ywave1 = org_bg
       wave xwave2 = root:index:t_Series
       wave ywave2 = org
       duplicate/o ywave2 ywave_regrid
       ywave_regrid = ywave1 [binarysearchinterp(xwave1, xwave2[p])]
       //appendtotable ywave_regrid
End Function
//***************************************************************************
Function  Regrid_Param2(xwave1, ywave1, xwave2, ywave2)
	Wave xwave1, ywave1, xwave2, ywave2
//       Wave xwave1 = bg_time
//       wave ywave1 = org_bg
//       wave xwave2 = root:index:t_Series
//       wave ywave2 = org
       duplicate/o ywave2 ywave_regrid
       ywave_regrid = ywave1 [binarysearchinterp(xwave1, xwave2[p])]
       //appendtotable ywave_regrid
End Function
//***************************************************************************
Function  Regrid_Matrix(xwave1, MxWv1, xwave2)
	// Created 2/5/2019 by JAH
	Wave xwave1, MxWv1, xwave2
	Variable iCol, iRow
	Make/o/n=(numpnts(xwave2), DimSize(MxWv1,1)) NewMx_regrid
	For (iCol = 0 ; iCol <= DimSize(MxWv1, 1) ; iCol += 1)
	       NewMx_regrid [][iCol] = MxWv1 [binarysearchinterp(xwave1, xwave2[p])] [iCol]
	Endfor
End Function
//***************************************************************************
Function RemoveZeros(InputWv)
	Wave InputWv
	Variable n
	For (n = 0 ; n <= numpnts(InputWv) ; n+= 1)
		If (InputWv[n] == 0)
			InputWv[n] = nan
		Elseif (InputWv[n] < -0.1)
			InputWv[n] = 0
		Endif
	Endfor
End Function
//***************************************************************************
Function DeleteNaNs(InputWv)
	Wave InputWv
	Variable n
	For (n = 0 ; n <= numpnts(InputWv) ; n+= 1)
		If (numtype(InputWv[n]) == 2)
			Deletepoints n,1,InputWv
			n-=1
		Endif
	Endfor
End Function
//***************************************************************************
Function RemoveInfs(InputWv,MatrixVar)
	Wave InputWv
	Variable MatrixVar
	Variable n, m
	
	If (MatrixVar == 0)
		For (n = 0 ; n <= numpnts(InputWv) ; n+= 1)
			If (numtype(InputWv[n]) == 1)
				InputWv[n] = nan
			Endif
		Endfor
	ElseIf (MatrixVar == 1)
		For (m = 0 ; m <= DimSize(InputWv,0) ; m += 1)	
			For (n = 0 ; n <= DimSize(InputWv,1) ; n+= 1)
				If (numtype(InputWv[m][n]) == 1)
					InputWv[m][n] = nan
				Endif
			Endfor
		Endfor
	Endif
End Function
//***************************************************************************

Function AddZeros(InputWv)
// Created 07/20/21 by MPPF as modified version of previous fxn (by JAH) to accommodate multi-dimensional waves.
// Purpose: to convert all entries equal to NaN in a wave to zeros.
	Wave InputWv
	InputWv = NumType(InputWv) == 2 ? 0 : InputWv
//	Variable Dims = WaveDims(InputWv)
//	Variable i,j,k,l
////Replace NaNs with zeros for a 1D wave.
//	If(Dims == 1)
//		For (i = 0 ; i < numpnts(InputWv) ; i += 1)
//			If (numtype(InputWv[i]) == 2)
//				InputWv[i] = 0
//			Endif
//		Endfor
//	EndIf
//	
////2D wave.
//	If(Dims == 2)
//		For (j = 0; j < DimSize(InputWv, 1); j += 1)
//			For (i = 0; i < DimSize(InputWv, 0); i += 1)
//				If (numtype(InputWv[i][j]) == 2)
//					InputWv[i][j] = 0
//				EndIf
//			EndFor
//		EndFor
//	EndIf
//	
////3D wave.
//	If(Dims == 3)
//		For (k = 0; k < DimSize(InputWv, 2); k += 1)
//			For (j = 0; j < DimSize(InputWv, 1); j += 1)
//				For (i = 0; i < DimSize(InputWv, 0); i += 1)
//					If (numtype(InputWv[i][j][k]) == 2)
//						InputWv[i][j][k] = 0
//					EndIf
//				EndFor
//			EndFor
//		EndFor
//	EndIf
//	
////4D wave.
//	If(Dims == 4)
//		For (l = 0; l< DimSize(InputWv, 3); l += 1)
//			For (k = 0; k < DimSize(InputWv, 2); k += 1)
//				For (j = 0; j < DimSize(InputWv, 1); j += 1)
//					For (i = 0; i < DimSize(InputWv, 0); i += 1)
//						If (numtype(InputWv[i][j][k][l]) == 2)
//							InputWv[i][j][k][l] = 0
//						EndIf
//					EndFor
//				EndFor
//			EndFor
//		EndFor
//	EndIf
End Function

//***************************************************************************
Function RemoveTimePoints(InputWv, TimeWave)
	Wave InputWv, TimeWave
	Variable n
	For (n = numpnts(InputWv) ; n >= 0 ; n -= 1)
		If (numtype(InputWv[n]) == 2)
			Deletepoints n,1, InputWv, TimeWave
		Endif
	Endfor
End Function
//***************************************************************************
Function RemovePoints(ImportWv,Value,State)// State = 0 means '<', State = 1 means '>'
	// Created 08.13.2010
	Wave ImportWv
	Variable Value
	Variable State
	Variable n
	For (n = 0 ; n < numpnts(ImportWv) ; n += 1)
		If (State == 0)
			If (ImportWv[n] < Value)
				ImportWv[n] = nan
			Endif
		Elseif (State == 1)
			If (ImportWv[n] > Value)
				ImportWv[n] = nan
			Endif
		Endif
			endfor
End Function
//***************************************************************************
Function PrintAvgFromCursors()
	// Created 4.29.09
	print (vcsr(A) + vcsr(B)) / 2
End Function
//***************************************************************************
Function FillPoints(InputWv)
	Wave InputWv
	Variable n, xStore
	
	For (n = 0 ; n <= numpnts(InputWv) ; n += 1)
		If (numtype(InputWv[n]) == 0)
			xStore = InputWv[n]
		Endif
		InputWv[n] = xStore
	Endfor
End Function
//***************************************************************************
Function CheckMonotonic(InputWv)
	Wave InputWv
	Variable n
	
	Make/o/n=(numpnts(InputWv)) CheckMonotonicWv
	For (n = 0 ; n <= numpnts(InputWv) ; n += 1)
		CheckMonotonicWv[n] = InputWv[n+1] - InputWv[n]
	Endfor
	Display CheckMonotonicWv
	SetAxis left -1,1
End Function
//***************************************************************************
Function PointOutNaNs(InputWv, DisplayVar)
	Wave InputWv
	Variable DisplayVar
	Variable n
	
	Make/o/n=(numpnts(InputWv)) PointOutNaNWv
	PointOutNaNWv = nan
	For (n = 0 ; n <= numpnts(InputWv) ; n += 1)
		If (NumType(InputWv[n]) == 2)
			PointOutNaNWv[n] = 1
		Elseif (NumType(InputWv[n]) == 1)
			PointOutNaNWv[n] = 1
		Endif
	Endfor
	
	If (DisplayVar == 1)
		Display /W=(372,62.75,821.25,323.75) PointOutNaNWv
		ModifyGraph mode=3
		ModifyGraph msize=8
		ModifyGraph mirror=2
		ModifyGraph fSize=14
		Label bottom "Point Number"
		SetAxis left 1,1
		ShowInfo
		TextBox/C/N=text0/F=0/A=MC/X=-15.78/Y=-44.86 "\\Z14NaN values located at markers"
	Endif
End Function
//***************************************************************************
Function RoundWave(InputWv, InputStr, RndVar)
	Wave InputWv
	String InputStr
	Variable RndVar	// How much to round the wave by (in minutes)
	
	String NewStr = InputStr + "_rnd"
	Duplicate/o InputWv $NewStr
	Wave LocalWv = $NewStr
	LocalWv /= (RndVar * 60)
	LocalWv = round(LocalWv)
	LocalWv *= (RndVar * 60)
	
End Function

//***************************************************************************
Function TimeAverage(InputDateWv, InputYWv, TimeSpan, SampleInterval, OffsetVar, DeleteNaNVar)
	Wave InputDateWv, InputYWv
	Variable TimeSpan	// Time period to average (in minutes)
	Variable SampleInterval // Frequency of instrument sampling (in minutes)
	Variable OffsetVar // # positive offset points (in minutes)
	Variable DeleteNaNVar // Removes all rows where no data if = 1
	Variable n, PointerVar, TempSumVar, NumAvgPntsVar, NewPositionVar, StatsPointerVar
	Duplicate/d/o InputDateWv LocalXWv
	Duplicate/d/o InputYWv LocalYWv
	LocalXWv += (OffsetVar*60)	
	LocalXWv /= (TimeSpan * 60)	// Added 2014.10.07 to offset if linking with regrid in next step
	LocalXWv = round(LocalXWv)
	LocalXWv *= (TimeSpan * 60)
	LocalXWv -= (OffsetVar*60)
	
	Make/d/o/n=(TimeSpan/SampleInterval + 4) StatsCache
	StatsCache = nan
	StatsPointerVar = 0 + OffSetVar
	Make/d/o/n=0 AveragedYWv, AveragedXWv
	n = 0 
	PointerVar = LocalXWv[n]
	TempSumVar = 0
	NumAvgPntsVar = 0
	NewPositionVar = 0
	
	For (n = 0 ; n < numpnts(LocalXWv) ; n += 1)
		If (LocalXWv[n] == PointerVar)
			StatsCache[StatsPointerVar] = InputYWv[n]
			StatsPointerVar += 1
		Else
			InsertPoints numpnts(AveragedYWv),1,AveragedYWv
			InsertPoints numpnts(AveragedXWv),1,AveragedXWv
			Wavestats/q/z StatsCache
			AveragedYWv[NewPositionVar] = V_avg
			AveragedXWv[NewPositionVar] = LocalXWv[n-1]		
			PointerVar = LocalXWv[n]
			StatsPointerVar = 0
			wave wave0 = nan
			StatsCache = nan
			TempSumVar = 0
			NumAvgPntsVar = 0
			NewPositionVar += 1
			n -= 1
		Endif
	Endfor
	
	// Forces new wave to start with data (instead of long section of nans)
	If (DeleteNaNVar == 1)
		For (n = 0 ; n < numpnts(AveragedXWv) ; n += 1)
			If (NumType(AveragedYWv[n]) == 2)
				Deletepoints n,1,AveragedXWv, AveragedYWv
				n -= 1
			Endif
		Endfor
	Endif
	
//	// Cleans up waves by removing extra nans at end
//	For (n = numpnts(AveragedXWv) ; n >= 0 ; n -= 1)
//		If (AveragedXWv[n] == 0)
//			Deletepoints n,1,AveragedXWv, AveragedYWv
//		Endif
//		//n-=1
//	Endfor

	SetScale d 0, 0, "dat", AveragedXWv
	Killwaves/z LocalXWv, LocalYWv, StatsCache
	
End Function

//***************************************************************************
Function TimeAverageMx(InputDateWv, InputYMx, TimeSpan, SampleInterval, OffsetVar, DeleteNaNVar)
	// Built 2/12/2019 by JAH from TimeAverage(..) as of v16.5
	Wave InputDateWv, InputYMx
	Variable TimeSpan	// Time period to average (in minutes)
	Variable SampleInterval // Frequency of instrument sampling (in minutes)
	Variable OffsetVar // # positive offset points (in minutes)
	Variable DeleteNaNVar // Removes all rows where no data if = 1
	Variable n, PointerVar, TempSumVar, NumAvgPntsVar, NewPositionVar, StatsPointerVar, iCol
	Duplicate/d/o InputDateWv LocalXWv
	Duplicate/d/o InputYMx LocalYMx
	LocalXWv += (OffsetVar*60)	
	LocalXWv /= (TimeSpan * 60)	// Added 2014.10.07 to offset if linking with regrid in next step
	LocalXWv = round(LocalXWv)
	LocalXWv *= (TimeSpan * 60)
	LocalXWv -= (OffsetVar*60)
	
	Make/d/o/n=((TimeSpan/SampleInterval + 4),DimSize(InputYMx,1)) StatsCache	// now a matrix
	Make/d/o/n=(TimeSpan/SampleInterval + 4) TempCacheExtractWv
	StatsCache = nan
	StatsPointerVar = 0 + OffSetVar
		//Temp:  ..Mx Function Creation, 2/12/19
		//Make/d/o/n=0 AveragedYMx, AveragedXWv
	Make/d/o/n=0 AveragedXWv
	Make/d/o/n=(0,DimSize(InputYMx,1)) AveragedYMx
	
	n = 0 
	PointerVar = LocalXWv[n]
	TempSumVar = 0
	NumAvgPntsVar = 0
	NewPositionVar = 0
	
	For (n = 0 ; n < numpnts(LocalXWv) ; n += 1)
		If (LocalXWv[n] == PointerVar)
			For (iCol = 0 ; iCol <= (DimSize(InputYMx,1)) ; iCol += 1)
				//StatsCache[StatsPointerVar] = InputYWv[n]
				StatsCache[StatsPointerVar][iCol] = InputYMx[n][iCol]		// STATS CACHE
			Endfor
			StatsPointerVar += 1
		Else
			InsertPoints numpnts(AveragedYMx),1,AveragedYMx
			InsertPoints numpnts(AveragedXWv),1,AveragedXWv
			AveragedXWv[NewPositionVar] = LocalXWv[n-1]		
			PointerVar = LocalXWv[n]
			For (iCol = 0 ; iCol <= (DimSize(InputYMx,1)) ; iCol += 1)
				TempCacheExtractWv = StatsCache[p][iCol]		
				Wavestats/q/z TempCacheExtractWv						// STATS CACHE
				AveragedYMx[NewPositionVar][iCol] = V_avg
			Endfor
			
			StatsPointerVar = 0
			wave wave0 = nan
			StatsCache = nan
			TempSumVar = 0
			NumAvgPntsVar = 0
			NewPositionVar += 1
			n -= 1
		Endif
	Endfor
	
	SetScale d 0, 0, "dat", AveragedXWv
	Duplicate/o AveragedXWv AveragedXWv_img
	Wave AveragedXWv_img
	InsertPoints numpnts(AveragedXWv_img),1,AveragedXWv_img
	Variable iPoint = numpnts(AveragedXWv_img)
	AveragedXWv_img[iPoint-1] = AveragedXWv_img[iPoint-2] + TimeSpan*60
	Killwaves/z LocalXWv, LocalYMx, StatsCache, TempCacheExtractWv
	
End Function
//***************************************************************************

Function CalcDiurnalAvg(NameStr, InputYWv, InputXWv, PointA, PointB, PointC, PointD, PlotVar, SDevVar)
	String NameStr	// Name of Y-wave input by user (using  " ")
	Wave InputYWv, InputXWv	// Wave names of y-, x- (with path if necessary?)
	Variable PointA, PointB // Gives the range over which to look at the data
	Variable PointC, PointD // Gives a second range over which to look at the data
	Variable PlotVar, SDevVar // PlotVar = 1 will make a plot, SDevVar = 1 will add error bars
	
	Variable n, i, n0, n1
	String ToDStr, InputYWv_HrStr, MeanStr, MedianStr, SDevStr
	

	Duplicate/o InputXWv InputXWv_ToD
	For (n = 0 ; n <= numpnts(InputXWv) ; n += 1)
		InputXWv_ToD[n] = mod(InputXWv[n],(60*60*24))
	Endfor
	
	Make/o/n=24 DiurnalOutput_Mean, DiurnalOutput_Median, DiurnalOutput_SDev
	Make/o/n=24 DiurnalOutput_Mean, DiurnalOutput_5th, DiurnalOutput_25th, DiurnalOutput_75th, DiurnalOutput_95th
	
	Duplicate/o InputYWv RangeMask 
	DiurnalOutput_Mean = nan
	DiurnalOutput_Median = nan
	RangeMask = nan
	RangeMask[PointA, PointB] = 1
	RangeMask[PointC, PointD] = 1
	
	For (n = 0 ; n <= 23 ; n += 1)
		ToDStr = "Msk_Hr" + num2str(n)
		Make/o/n=(numpnts(InputXWv)) $ToDStr
		Wave ToDMask = $ToDStr
		ToDMask = nan
		
		For (i = 0 ; i <= numpnts(InputXWv) ; i += 1)
			n0 = n*3600
			n1 = (n+1)*3600
			If (InputXWv_ToD[i] >= n0 && InputXWv_ToD[i] < n1)
				ToDMask[i] = 1
			Endif
		Endfor
		
		InputYWv_HrStr =  "InputYWv_" + num2str(n)
		Duplicate/o InputYWv $InputYWv_HrStr
		Wave InputYWv_Hr = $InputYWv_HrStr
		InputYWv_Hr = ToDMask * InputYWv * RangeMask
		Duplicate/o InputYWv_Hr TempWv
		Sort TempWv, TempWv
		Wavestats/q TempWv
		Deletepoints V_npnts, 100000, TempWv
				
		DiurnalOutput_Mean[n] = V_Avg
		DiurnalOutput_Median[n] = TempWv[(V_npnts/2) -1]
		DiurnalOutput_5th[n] = TempWv[(V_npnts*1/20) -1]
		DiurnalOutput_25th[n] = TempWv[(V_npnts*1/4) -1]
		DiurnalOutput_75th[n] = TempWv[(V_npnts*3/4) -1]
		DiurnalOutput_95th[n] = TempWv[(V_npnts*19/20) -1]
		DiurnalOutput_SDev[n] = V_SDev
		
		Killwaves ToDMask, InputYWv_Hr
	Endfor
	
	Make/o/n=24 TimeOfDay
	TimeOfDay = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
	 
	MeanStr = NameStr + "_DiurnalMean"
	MedianStr = NameStr + "_DiurnalMedian"
	SDevStr = NameStr + "_DiurnalSDev"
	Duplicate/o DiurnalOutput_Mean $MeanStr
	Duplicate/o DiurnalOutput_Median $MedianStr
	Duplicate/o DiurnalOutput_SDev $SDevStr
	Wave MeanWv = $MeanStr
	Wave MedianWv = $MedianStr
	Wave SDevWv = $SDevStr
	 
	 If (PlotVar == 1)
	 	Display /W=(35.25,42.5,523.5,360.5) $MeanStr,MedianWv vs TimeOfDay
		ModifyGraph mode=4
		ModifyGraph marker($MeanStr)=19,marker($MedianStr)=7
		ModifyGraph lSize=1.5
		ModifyGraph rgb($MeanStr)=(0,0,0),rgb($MedianStr)=(0,15872,65280)
		ModifyGraph msize=5
		ModifyGraph mrkThick($MedianStr)=2
		ModifyGraph offset($MedianStr)={0.25,0}
		ModifyGraph mirror=2
		ModifyGraph fSize=14
		Label left NameStr + " (Number Concentration (L\S-1\M))"
		Label bottom "Time of Day"
		TextBox/C/N=text1/X=1.00/Y=1.00/F=0 NameStr
		If (SDevVar == 1)
			ErrorBars/T=1.5/L=1.5 $MeanStr Y,wave=($SDevStr,$SDevStr)
			ErrorBars/T=1.5/L=1.5 $MedianStr Y,wave=($SDevStr,$SDevStr)	
			Legend/C/N=text0/J/F=0/B=1/X=75.05/Y=82.69 "\\Z14\\s(" + MeanStr + ") Mean\r\\s(" + MedianStr + ") Median\r\\Z12(Error Bars = sdev)"
		Else
			Legend/C/N=text0/J/F=0/B=1/X=79/Y=87 "\\Z14\\s(" + MeanStr + ") Mean\r\\s(" + MedianStr + ") Median"
		Endif
	Endif
	
	Killwaves TempWv, InputXWv_ToD, DiurnalOutput_Mean, DiurnalOutput_Median, DiurnalOutput_SDev	
		
End Function

//***************************************************************************

Function ExampleLoadFile()
	svar PathString, FileString, ImportName
	nvar TableVar
	Variable n

	String PathandFileString = "D:RESEARCH-MPI:Projects:" + PathString + ":" + FileString + ".txt"
	
	String FolderName = "root:" + ImportName
	Wave CheckWv = root:$ImportName
	If (DataFolderExists(FolderName) == 1)
		abort "Aborting Calculation.\r \rThis folder already exists.  Please check file names and try again.\rAborted from function: Load()."
	Endif
	If (TableVar == 0)
		LoadWave/J/D/Q/A=$ImportName/K=0/L={0,7,0,0,68} PathandFileString
	Elseif (TableVar ==1)
		LoadWave/J/D/Q/A=$ImportName/E=1/K=0/L={0,7,0,0,68} PathandFileString
	Endif

   // This part deals with the fact that the diameter wave is imported as a text wave
	String FolderName2 = ImportName + "3"
	Make/o/n=(numpnts($FolderName2)) TempWv3
	Wave/t SourceWv = $FolderName2
	For (n = 0 ; n <= numpnts(TempWv3) ; n += 1)
		TempWv3[n] = str2num(SourceWv[n])
	Endfor
	TempWv3[0,1] = nan
	Killwaves SourceWv
	Rename TempWv3 $FolderName2
	
End Function
//***************************************************************************
Function MakeMatrixOutOfWaves(ColNumber, RowNumber)
	Variable ColNumber, RowNumber
	Variable n
	String NameStr
	
	Make/o/n=(RowNumber,ColNumber) NewNameMx	
	NewNameMx = nan
	For (n = 0 ; n <= ColNumber ; n += 1)
		NameStr = "wave" + num2str(n)
		Wave NameWv = $NameStr
		NewNameMx[][n] = NameWv[p]
	Endfor
	
End Function
//***************************************************************************
Function WeightedAvg(ImportWvToAvg, ImportWvForWeight, FromPt, ToPt)	// See Notebook MPIC-2, p. 29
	Wave ImportWvToAvg
	Wave ImportWvForWeight
	Variable FromPt
	Variable ToPt
	Variable n, WeightVar
	
	If (ToPt == 0)
		ToPt = numpnts(ImportWvToAvg)
	Endif

	WeightVar = sum(ImportWvForWeight,FromPt,ToPt)

	Duplicate/d/o ImportWvForWeight WeightTempWv
	WeightTempWv /= sum(ImportWvForWeight,FromPt,ToPt)
	
	Duplicate/d/o ImportWvToAvg WeightedAvgTempWv
	WeightedAvgTempWv *= WeightTempWv
	
	print "Weighted Average is: " + num2str(sum(WeightedAvgTempWv,FromPt,ToPt))
		
	Killwaves/z WeightTempWv, WeightedAvgTempWv
	
End Function
//***************************************************************************
Function IntegrateSizeDistMx(DataMatrix, DiameterWv, TimeWv, TimePtA, TimePtB, SizePtA, SizePtB)
	Wave DataMatrix, DiameterWv, TimeWv
	Variable TimePtA, TimePtB,SizePtA, SizePtB 
	Variable IntA, IntB
	Variable n, m

	Make/o/n=(numpnts(TimeWv)) IntegratedWv, TempMxSlice
	Duplicate/o DiameterWv LogDiameterWv
	LogDiameterWv = log(DiameterWv)
	IntegratedWv = nan
	IntA = LogDiameterWv[SizePtA]
	IntB = LogDiameterWv[SizePtB]
	
	For (n = TimePtA ; n < TimePtB ; n += 1)
		For (m = 0 ; m < numpnts(DiameterWv) ; m += 1)
			TempMxSlice[m] = DataMatrix[n][m]
		Endfor
		IntegratedWv[n] = areaXY(LogDiameterWv, TempMxSlice, IntA, IntB)
	Endfor
	
	killwaves/z LogDiameterWv, TempMxSclice

End Function
//***************************************************************************
Function IntegrateSizeDistMxOption(DataMatrix, DiameterWv, TimeWv, TimePtA, TimePtB, SizePtA, SizePtB, OptionVar)
	// Duplicated version of IntegrateSizeDistMx
	// Purpose: To integrate matrices whether provided as dN/dlogDp, dN/Dp, or otherwise ...
	// Written by JAH (& MIG), 8/17/2015
	Wave DataMatrix, DiameterWv, TimeWv
	Variable TimePtA, TimePtB,SizePtA, SizePtB
	Variable OptionVar // OptionVar: 1 = dN/dlogDp, 2 = dN/Dp, 3 = dN 
	Variable IntA, IntB
	Variable n, m

	Make/o/n=(numpnts(TimeWv)) IntegratedWv, TempMxSlice
	Duplicate/o DiameterWv LogDiameterWv
	LogDiameterWv = log(DiameterWv)
	IntegratedWv = nan

	If (OptionVar == 1)
		IntA = LogDiameterWv[SizePtA]
		IntB = LogDiameterWv[SizePtB]
	ElseIf (OptionVar == 2)
		IntA = DiameterWv[SizePtA]
		IntB = DiameterWv[SizePtB]
	ElseIf (OptionVar == 3)
		IntA = DiameterWv[SizePtA]
		IntB = DiameterWv[SizePtB]
	Endif
	
	For (n = TimePtA ; n < TimePtB ; n += 1)
		For (m = 0 ; m < numpnts(DiameterWv) ; m += 1)
			TempMxSlice[m] = DataMatrix[n][m]
		Endfor
		If (OptionVar == 1)
			IntegratedWv[n] = areaXY(LogDiameterWv, TempMxSlice, IntA, IntB)
		Elseif (OptionVar == 2)
			IntegratedWv[n] = areaXY(DiameterWv, TempMxSlice, IntA, IntB)
		Elseif (OptionVar == 3)
			IntegratedWv[n] = sum(TempMxSlice, IntA,IntB)
		Endif
	Endfor

End Function
//***************************************************************************

Function AverageSizeDistMx(DataMatrix, DiameterWv, TimePtA, TimePtB, CursorVar, PlotVar)
	// Created 2.5.2013 by JAH
	
	Wave DataMatrix		// Matrix: Time=Rows, Diameter=Columns
	Wave DiameterWv 
	Variable TimePtA, TimePtB
	Variable CursorVar	// CursorVar = 1 means use cursors instead of manual point input
							// Just make sure to highlight plot where cursors are highlighted
	Variable PlotVar		// PlotVar = 1 means to auto-plot
	Variable iRow, iCol, iCount, n
	
	Duplicate/o DataMatrix MatrixTemp
	Duplicate/o DiameterWv DiamTemp

	Make/o/n=(numpnts(DiamTemp)) AvgSizeDist
	AvgSizeDist = 0
	iCount = 0

	If (CursorVar == 1)
		For (iRow = pcsr(A) ; iRow <= pcsr(B) ; iRow += 1)	
			AvgSizeDist += MatrixTemp[iRow][p]
			iCount += 1
		Endfor	
	Else
		For (iRow = TimePtA ; iRow <= TimePtB ; iRow += 1)
			variable probe1 = MatrixTemp[irow][0]
			If (numtype(MatrixTemp[iRow][0]) == 0)		// Ignores if first column of row is nan
				AvgSizeDist += MatrixTemp[iRow][p]
				iCount += 1
			Endif
		Endfor
	Endif

	AvgSizeDist /= iCount


	If (PlotVar == 1)
		Display/W=(789.75,66.5,1239.75,356.75) AvgSizeDist vs DiamTemp
		ModifyGraph mode=4
		ModifyGraph marker=19
		ModifyGraph lSize=1.5
		ModifyGraph rgb=(36864,14592,58880)
		ModifyGraph msize=3
		ModifyGraph log(bottom)=1
		ModifyGraph mirror=2
		ModifyGraph fSize=15
		Label bottom "Diameter"
		Legend/C/N=text0/J/F=0/B=1/X=3.10/Y=4.75 "\\Z15\\s(AvgSizeDist) AvgSizeDist"
	EndIf
	
	Killwaves/z MatrixTemp, DiamTemp
	
End Function

//***************************************************************************
//Function IntegrateSizeDist(DataWv, DiameterWv, IntA, IntB)
//	Wave DataWv, DiameterWv
//	Variable IntA, IntB
//	Variable n, m
//
////	Make/o/n=(numpnts(TimeWv)) IntegratedWv, TempMxSlice
//	Duplicate/o DiameterWv LogDiameterWv
//	LogDiameterWv = log(DiameterWv)
//	IntegratedWv = nan
//	
////	For (n = TimePtA ; n <= TimePtB ; n += 1)
//		For (m = 0 ; m < numpnts(DiameterWv) ; m += 1)
//			TempMxSlice[m] = DataWv[n][m]
//		Endfor
//		IntegratedWv[n] = areaXY(LogDiameterWv, TempMxSlice, IntA, IntB)
////	Endfor
//
//End Function
//***************************************************************************
Function MakeWhiskers(SubFolderName)	
	// IMPORTANT NOTE: Sub-folder must be in root:Whiskers:$SubFolderName and must include the following waves, named as:
	//	MeanWv, MedianWv, p5Wv, p25Wv, p75Wv, p95Wv (pXX = XXth percentile)
	String SubFolderName
	Variable xVar, n, RangeVarLow, RangeVarHigh
	
	String fldrSav0= GetDataFolder(1)	
	String TempStr = "root:Whiskers:" + SubFolderName
	If (DataFolderExists("root:Whiskers") == 0)
		abort "root:Whiskers and appropriate Sub-Folders do not exist."
	ElseIf (DataFolderExists(TempStr) == 0)
		abort "root:Whiskers:Sub-Folder does not exist."
	Endif
	SetDataFolder root:Whiskers:$SubFolderName
	xVar = numpnts(MeanWv) // Number of boxes
	Wave MeanWv, MedianWv, p5Wv, p25Wv, p75Wv, p95Wv

	Duplicate/o MeanWv p25BarWv, p75BarWv
	p25BarWv = MedianWv - p25Wv
	p75BarWv = p75Wv - MedianWv
	
	Make/o/n=(xVar) WhiskerXWave, WhiskerXPosition, WhiskerYPosition
	For (n = 0 ; n < xVar ; n += 1)
		WhiskerXWave[n] = n+1
	Endfor

	RangeVarLow = WhiskerXWave[0] - 0.4
	RangeVarHigh = WhiskerXWave[numpnts(WhiskerXWave)] + 0.4

	If (numpnts(WhiskerXWave) == 1)
		WhiskerXPosition = {0}
		WhiskerYPosition = {-58.8}
	Elseif (numpnts(WhiskerXWave) == 2)
		WhiskerXPosition = {-27.4, 27.4}
		WhiskerYPosition = {-58.8, -58.8}	
	Elseif (numpnts(WhiskerXWave) == 3)
		WhiskerXPosition = {-35, 0,35}
		WhiskerYPosition = {-55,-61,-55}	
	Elseif (numpnts(WhiskerXWave) == 4)
		WhiskerXPosition = {-38,-12.8,13.6,39.5}
		WhiskerYPosition = {-55,-61,-55,-61}	
	Elseif (numpnts(WhiskerXWave) == 5)
		WhiskerXPosition = {-42.5,-21,0,21.7,42.2}
		WhiskerYPosition = {-55,-61,-55,-61,-55}
	Elseif (numpnts(WhiskerXWave) == 6)
		WhiskerXPosition = {-42,-25.5,-8,9,27,44}
		WhiskerYPosition = {-55,-61,-55,-61,-55,-61}
	Else
		WhiskerXPosition = 0
		WhiskerYPosition = 58.8	
	Endif

	Display /W=(29.25,46.25,377.25,329.75) MedianWv,MedianWv vs WhiskerXWave
	AppendToGraph p95Wv, p75Wv, p5Wv, p25Wv, MeanWv vs WhiskerXWave
	ModifyGraph margin(bottom)=36
	ModifyGraph mode(MedianWv)=2,mode(MedianWv#1)=2,mode(p95Wv)=1
	ModifyGraph mode(p75Wv)=2,mode(p5Wv)=1,mode(p25Wv)=2,mode(MeanWv)=3
	ModifyGraph marker(MeanWv)=8, opaque(MeanWv)=1, rgb(MedianWv#1)=(47872,47872,47827)
	ModifyGraph rgb(MedianWv)=(0,0,0),rgb(p95Wv)=(0,0,0),rgb(p75Wv)=(0,0,0)
	ModifyGraph rgb(p5Wv)=(0,0,0),rgb(p25Wv)=(0,0,0),rgb(MeanWv)=(0,0,0)
	ModifyGraph msize(MeanWv)=4,mrkThick(MeanWv)=1.5
	ModifyGraph toMode(p95Wv)=1,toMode(p5Wv)=1
	ModifyGraph tick(bottom)=3
	ModifyGraph mirror=2
	ModifyGraph noLabel(bottom)=1
	ModifyGraph fSize=14
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=1.5
	ModifyGraph ZisZ(left)=1
	ModifyGraph manTick=0
	ModifyGraph margin(left)=65
	SetAxis left 0,*
	SetAxis bottom RangeVarLow, RangeVarHigh
	ErrorBars/T=2/L=0 MedianWv BOX,const=0.2,wave=(p75BarWv,p25BarWv)
	ErrorBars/T=0/L=2 MedianWv#1 X,const=0.2
	ErrorBars/T=0/L=1.5 p95Wv X,const=0.1
	ErrorBars/T=0/L=1.5 p5Wv X,const=0.1
	ShowInfo
	SetDataFolder fldrSav0

End Function
//***************************************************************************
// Created 2011.09.14
Function MakeBarDateFn(InputStartWv, InputStopWv)
	Wave InputStartWv, InputStopWv
	variable n, m
	
	Duplicate/o InputStartWv TempWv1
	Make/o/n=(numpnts(InputStartWv)) TempWv2
	TempWv2 = 1
	m=0
	
	For (n = 1 ; n <= numpnts(TempWv1) ; n += 2)
		InsertPoints n,1, TempWv1, TempWv2
		TempWv1[n] = InputStopWv[m]
		m+= 1
	Endfor
	
End Function
//***************************************************************************
// Created 2011.09.23

Function MultiplyMaskThroughMatrix(InputMx, InputMskWv)
	Wave InputMx, InputMskWv
	Variable iRow, iCol
	
	
	Duplicate/o InputMx MaskedMx
	For (iRow = 0 ; iRow < numpnts(InputMskWv) ; iRow += 1)
		For (iCol = 0 ; iCol < DimSize(InputMx,1) ; iCol += 1)
			MaskedMx[iRow][iCol] = InputMx[iRow][iCol] * InputMskWv[iRow]
		Endfor
	Endfor
End Function

//***************************************************************************
// Created 2012.06.26 by JAH

Function BinCorrelationFn(InputXWv, InputYWv, LowerXBound, UpperXBound, NumXBins, PrintOptionVar)
	// Takes a correlation (scatter) plot and bins with respect to the x-values so that the trends can be more easily seen
	
	// Simply needs input x and y waves in the folder currently selected
	// Critical output waves: YAvgWv, XWv		(re-name before running again)
	//	Other output waves: YSDYWv, NumSumWv, NumSumTxtWv
	
	Wave InputXWv, InputYWv
	Variable LowerXBound, UpperXBound, NumXBins, PrintOptionVar
	
	Variable n, xCursor
	Variable binwidth, localcursor, shiftnumvar, shiftvar, localLowerXBound, returnvar
	Variable iCursor
	
	Duplicate/o InputXWv LocalXWv
	Duplicate/o InputYWv LocalYWv
	Sort LocalXWv, LocalXWv, LocalYWv
		
	binwidth= (UpperXBound - LowerXBound) / NumXBins
	localcursor = 0
	Make/o/n=1 LocalSumWv, LocalXSumWv, YSumWv, YAvgWv, YSDWv, NumSumWv, XWv, YMedianWv
	Make/o/t/n=1 NumSumTxtWv

	n = 0
	ReturnVar = 0
	
	Do															// This deals with x-waves whose first point is outside of the first requested bin
		If ( (LowerXBound + binwidth*n) < LocalXWv[0])
			If (LocalXWv[0] < LowerXBound + binwidth*(n+1) )
				ReturnVar = 1
				LocalLowerXBound = LowerXBound + (n) * binwidth
			Endif
		Endif
		n += 1
	While (ReturnVar == 0) 

	xCursor =LocalLowerXBound
	
	For (n = 0 ; n < numpnts(LocalXWv) ; n += 1)					// Routine goes through all data points
		If (LocalXWv[n] >= xCursor)									//If loop #1: adds values to LocalSumWv if in one particular bin
			LocalSumWv[localcursor] = LocalYWv[n]
			LocalXSumWv[localcursor] = LocalXWv[n]					
			localcursor += 1
			InsertPoints numpnts(LocalSumWv), 1, LocalSumWv, LocalXSumWv
		Endif
		
		
		If (LocalXWv[n+1] > UpperXBound)								// If loop #2: once the bin is finished, the routine enters here to do stats and keep values
		ElseIf (LocalXWv[n+1] >= XCursor + binwidth)						// the counter (XCursor) is increased and the LocalSumWv is flushed
			Deletepoints numpnts(LocalSumWv)-1, 1, LocalSumWv, LocalXSumWv
			Wavestats/q LocalSumWv
			BinCorrelation_MedianSubFn()
			nvar LocalMedianVar
			YSumWv[numpnts(YSumWv)] = v_avg * v_npnts
			YAvgWv[numpnts(YSumWv)] = v_avg
			YMedianWv[numpnts(YSumWv)] = LocalMedianVar
			YSDWv[numpnts(YSumWv)] = v_sdev
			NumSumWv[numpnts(NumSumWv)] = v_npnts
			XWv[numpnts(NumSumWv)] = XCursor 
			NumSumTxtWv[numpnts(NumSumWv)] = num2str(v_npnts)
			InsertPoints numpnts(YSumWv), 1, YSumWv, YAvgWv, YSDWv, NumSumWv, XWv, NumSumTxtWv, YMedianWv
			Make/o/n=1 LocalSumWv, LocalXSumWv
			XCursor += binwidth
		ElseIf (n == (numpnts(LocalXWv)-1) )							// This part of the loop provides identical commands to the immediatley previous part
			Deletepoints numpnts(LocalSumWv)-1, 1, LocalSumWv, LocalXSumWv			// This section allows the last point to be calculated
			Wavestats/q LocalSumWv
			YSumWv[numpnts(YSumWv)] = v_avg * v_npnts
			YAvgWv[numpnts(YSumWv)] = v_avg
			YMedianWv[numpnts(YSumWv)] = LocalMedianVar
			YSDWv[numpnts(YSumWv)] = v_sdev
			NumSumWv[numpnts(NumSumWv)] = v_npnts
			XWv[numpnts(NumSumWv)] = XCursor 
			NumSumTxtWv[numpnts(NumSumWv)] = num2str(v_npnts)
			InsertPoints numpnts(YSumWv), 1, YSumWv, YAvgWv, YSDWv, NumSumWv, XWv, NumSumTxtWv, YMedianWv
			Make/o/n=1 LocalSumWv, LocalXSumWv
			XCursor += binwidth
		Endif
	Endfor
	
	XWv[numpnts(XWv)] = XWv[numpnts(XWv)-2] + binwidth				// Without this, the bar is not created correctly (missing last bar)
	YSumWv[numpnts(XWv)] = nan
	YAvgWv[numpnts(XWv)] = nan
	YSDWv[numpnts(XWv)] = nan
	YMedianWv[numpnts(XWv)] = nan
	NumSumTxtWv[numpnts(XWv)] = ""
	
	Duplicate/o YAvgWv HalfYAvgWv
	HalfYAvgWv /= 4
	
	Variable PlotOffsetVar = binwidth / 2

	If (PrintOptionVar == 1)
		BinCorrelationPrintFn(NumSumWv, YSumWv, YAvgWv, YSDYWv, XWv, NumSumTxtWv, PlotOffsetVar)		// YMedianWv not yet updated to print here
	Endif
	Killwaves/z LocalXWv, LocalYWv, LocalSumWv, LocalXSumWv, YSumWv

End Function
//***************************************************************************
Function BinCorrelation_MedianSubFn()
	Wave LocalSumWv, LocalXSumWv
	Variable n
	Variable/g LocalMedianVar = nan

	Duplicate/o LocalSumWv TempYWv
	Duplicate/o LocalXSumWv TempXWv
	
	Sort TempXWv, TempXWv, TempYWv
	For (n = 0 ; n <= numpnts(TempXWv) ; n += 1)
		If (numtype(TempXWv[n]) == 2)
			Deletepoints n, 1e9, TempXWv, TempYWv
		Endif
	Endfor
	
	Sort TempYWv, TempYWv, TempXWv
	For (n = 0 ; n <= numpnts(TempYWv) ; n += 1)
		If (numtype(TempYWv[n]) == 2)
			Deletepoints n, 1e9, TempYWv, TempXWv
		Endif
	Endfor
	
	LocalMedianVar = TempYWv[numpnts(TempYWv)/2]
	
	Killwaves/z LocalTempSumWv, TempLocalXWv, TempLocalYWv
End Function

//***************************************************************************
Function BinCorrelationPrintFn(NumSumWv, YSumWv, YAvgWv, YSDYWv, XWv, NumSumTxtWv, PlotOffsetVar)
	Wave NumSumWv, YSumWv, YAvgWv, YSDYWv, XWv, NumSumTxtWv
	Variable PlotOffsetVar
	
	Display /W=(218.25,64.25,789,497.75) YAvgWv vs XWv
	AppendToGraph HalfYAvgWv[0,106] vs XWv
	AppendToGraph YAvgWv vs XWv
	ModifyGraph mode(YAvgWv)=5,mode(HalfYAvgWv)=3,mode(YAvgWv#1)=2
	ModifyGraph rgb(YAvgWv)=(0,39168,39168),rgb(HalfYAvgWv)=(0,0,0),rgb(YAvgWv#1)=(0,0,0)
	ModifyGraph msize(HalfYAvgWv)=2
	ModifyGraph hbFill(YAvgWv)=2
	ModifyGraph useNegPat(YAvgWv)=1
	ModifyGraph hBarNegFill(YAvgWv)=4
	ModifyGraph useBarStrokeRGB(YAvgWv)=1
	ModifyGraph offset(HalfYAvgWv)={PlotOffsetVar,0},offset(YAvgWv#1)={PlotOffsetVar,0}
	ModifyGraph textMarker(HalfYAvgWv)={NumSumTxtWv,"default",0,0,5,0.00,0.00}
	ModifyGraph mirror=2
	ModifyGraph fSize=14
	Label left "Y-Axis Label"
	Label bottom "X-Axis Label"
	ErrorBars YAvgWv#1 Y,wave=(YSDWv,YSDWv)
	Legend/C/N=text0/J/F=0/B=1/A=LT/X=-1.15/Y=0.00 "\\Z14\\s(YAvgWv) YAvgWv"
End Function



//***************************************************************************

Function MakeThreshMask(InputWv, threshVar)
	// 2012.07.26
	// JAH
	
	Wave InputWv
	Variable threshVar
	Variable n
	
	Duplicate/o InputWv LocalWv, ThreshMsk
	
	threshMsk = 1
	
	//For (n = 0 ; n < numpnts(LocalWv) ; n += 1)
	For (n = 65988 ; n < 70219 ; n += 1)
		If (LocalWv[n] < threshVar)
			threshMsk[n] = 0
		Endif
	Endfor		
	
	Killwaves/z LocalWv
	
End Function

//***************************************************************************

Function DewPointFn(InputRHWv, InputTempWv)
	// 2012.12.27
	// JAH
	
	Wave InputRHWv, InputTempWv
	Variable a, b, c, d
	Variable n
	
	a = 6.1121
	b = 18.564
	c = 255.57
	d = 254.4
	
	Duplicate/o InputRHWv TempGammaWv, DewPoint, AboveDewPoint
	
	TempGammaWv = nan
	DewPoint = nan
	
	TempGammaWv = ln( (InputRHWv/100) * exp ( (b-(InputTempWv/d)) * (InputTempWv/(c+InputTempWv)) )   )
	DewPoint = c*TempGammaWv/(b-TempGammaWv)
	AboveDewPoint = InputTempWv-DewPoint
		
	//Killwaves/z TempGammaWv
	
End Function

//***************************************************************************

Function AddExpinfo2Graph()
	// Written by Harald Stark, Univ. Colorado - Boulder (Group J. Jimenez), Oct. 2013
  pathinfo home
  String message = "\\Z06Location:"+S_path
  message +="\rname: \\{igorinfo(1)}"
  message +="\rdate: \\{date()}, \\{time()}"
  message +="\rgraphname: \\{winlist(\"*\",\"\",\"WIN:\")}"
  TextBox/C/N=expinfo/F=0/A=RT/E=2 message
End Function 

//***************************************************************************

Function CalcDiurnalMx(InputMx, InputXWv, InputYWv, ExtraMskWv, PointA, PointB, KeepDiurnalHrSDWavesVar, PlotVar)

	Wave InputMx, InputXWv, InputYWv, ExtraMskWv
	Variable PointA, PointB, KeepDiurnalHrSDWavesVar, PlotVar
	String ToDStr, DiurnalMeanStr, DiurnalMedianStr, MxHrStr, MeanStr, MedianStr
	Variable n, m, iDex, nCol, i, n0, n1
	
//	Duplicate/o InputXWv DiurnalPlotXWv
//	InsertPoints (numpnts(InputXWv)),1, DiurnalPlotXWv 
//	DiurnalPlotXWv[numpnts(DiurnalPlotXWv)] = 2*InputXWv[numpnts(InputXWv)-1] - InputXWv[numpnts(InputXWv)-2]
	
	Duplicate/o InputYWv DiurnalPlotYWv
	InsertPoints (numpnts(InputYWv)),1, DiurnalPlotYWv 
	DiurnalPlotYWv[numpnts(DiurnalPlotYWv)] = 2*InputYWv[numpnts(InputYWv)-1] - InputYWv[numpnts(InputYWv)-2]
	
	Make/o/n=24 TimeOfDay
	TimeOfDay = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
	
	Make/o/n=25 DiurnalPlotTODWv
	 DiurnalPlotTODWv = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
	

	Make/o/n=(numpnts(InputXWv)) RangeMask
	RangeMask = nan
	RangeMask[PointA, PointB] = 1
	
	If (WaveExists(ExtraMskWv) == 1)
		RangeMask *= ExtraMskWv
	Endif

	Duplicate/o InputXWv InputXWv_ToD
	For (n = 0 ; n <= numpnts(InputXWv) ; n += 1)
		InputXWv_ToD[n] = mod(InputXWv[n],(60*60*24))
	Endfor
	
	For (n = 0 ; n <= 23 ; n += 1)
		ToDStr = "Msk_Hr" + num2str(n)
		Make/o/n=(numpnts(InputXWv)) $ToDStr
		Wave ToDMask = $ToDStr
		ToDMask = nan
		
		For (i = 0 ; i <= numpnts(InputXWv) ; i += 1)
			n0 = n*3600
			n1 = (n+1)*3600
			If (InputXWv_ToD[i] >= n0 && InputXWv_ToD[i] < n1)
				ToDMask[i] = 1 * RangeMask[i]
			Endif
		Endfor
	Endfor
	
		DiurnalMeanStr = "DiurnalMean"
		DiurnalMedianStr = "DiurnalMedian"
		nCol = DimSize(InputMx,1)
		Make/o/n=(24,nCol) DiurnalMeanMx, DiurnalMedianMx
			Wave DiurnalMeanMx 
			Wave DiurnalMedianMx 
			
		DiurnalMeanMx = nan
		DiurnalMedianMx = nan
		
	For (n = 0 ; n <= 23 ; n += 1)	// goes through daily hours
			MxHrStr = "Mx_" + num2str(n) + "Hr"
			Duplicate/o InputMx $MxHrStr
				Wave MatrixWvName = $MxHrStr

			ToDStr = "Msk_Hr" + num2str(n)
				Wave ToDMask = $ToDStr
				variable probe1 = dimsize(inputmx,0)
			For (m = 0 ; m < DimSize(InputMx,0) ; m += 1)			// For 0 - num(samples)
				MatrixWvName[m][] = InputMx[m][q] * ToDMask[m]	
			Endfor
			
			// Now average across time
			Variable t = DimSize(InputMx,0) - 1		//@change to '_mp'
			For (m = t ; m >= 0 ; m -= 1)
				If (numtype(MatrixWvName[m][1]) == 2)
					Deletepoints/M=0 m,1,MatrixWvName
				Endif
			Endfor

			MeanStr = "ToD_Hr"+ num2str(n) + "_Mean"
			MedianStr = "ToD_Hr"+ num2str(n) + "_Median"
				Make/o/n=(DimSize(MatrixWvName,1)) $MeanStr, $MedianStr
				Wave MeanWv = $MeanStr
				Wave MedianWv = $MedianStr
			
			For (m = 0 ; m <= DimSize(MatrixWvName,1) ; m += 1)	// Number of diameter points
				Make/o/n=(DimSize(MatrixWvName,0)) Temp_ExtractedRow
					Wave Temp_ExtractedRow = Temp_ExtractedRow
				Temp_ExtractedRow = MatrixWvName[p][m]	
				Wavestats/z/q Temp_ExtractedRow
				Sort Temp_ExtractedRow, Temp_ExtractedRow
				MeanWv[m] = V_avg
				MedianWv[m] = Temp_ExtractedRow[round((V_npnts/2) -1)]
			Endfor
			
			// Here it adds the SDs back into one matrix each
//			DiurnalMeanStr = "DiurnalMean"
//			DiurnalMedianStr = "DiurnalMedian"
//			Wave DiurnalMeanMx = $DiurnalMeanStr
//			Wave DiurnalMedianMx = $DiurnalMedianStr
			
			For (m = 0 ; m <= DimSize(MatrixWvName,1) ; m += 1)	// For every diameter point
				DiurnalMeanMx[n][m] = MeanWv[m]	// DiurnalMeanMx(Hours,Diameter)
				DiurnalMedianMx[n][m] = MedianWv[m]
			Endfor
			
			Killwaves/z Temp_ExtractedRow
			//Killwaves/z MatrixWvName
				
	Endfor
	
	Killwaves/z InputXWv_ToD
	
	For (n = 0 ; n <= 23 ; n += 1)
		ToDStr = "Msk_Hr" + num2str(n)
		Killwaves/z $ToDStr
	Endfor
	
	If (KeepDiurnalHrSDWavesVar != 1)
		For (n = 0 ; n <= 23 ; n += 1)
			String KillStrA, KillStrB
			KillStrA = "ToD_Hr" + num2str(n) + "_Mean"
			KillStrB = "ToD_Hr" + num2str(n) + "_Median"
			Killwaves/z $KillStrA, $KillStrB
		Endfor
	Endif
	Killwaves/z RangeMask
	
	If (PlotVar == 1)
		Display /W=(35.25,42.5,520.5,346.25)
		AppendImage DiurnalMeanMx vs {DiurnalPlotTODWv,DiurnalPlotYWv}
		ModifyImage DiurnalMeanMx ctab= {*,*,Rainbow,1}
		ModifyGraph log(left)=1
		ModifyGraph mirror=2
		ModifyGraph fSize=16
		Label left "Diameter (\\F'Symbol'm\\F'Arial'm)"
		Label bottom "Time of Day (Hr)"
		TextBox/C/N=text0/F=0/A=MC/X=-43.55/Y=43.63 "\\Z16Mean"
		AddExpinfo2Graph()

		Display /W=(551.25,37.25,1036.5,341)
		AppendImage DiurnalMedianMx vs {DiurnalPlotTODWv,DiurnalPlotYWv}
		ModifyImage DiurnalMedianMx ctab= {*,*,Rainbow,1}
		ModifyGraph log(left)=1
		ModifyGraph mirror=2
		ModifyGraph fSize=16
		Label left "Diameter (\\F'Symbol'm\\F'Arial'm)"
		Label bottom "Time of Day (Hr)"
		TextBox/C/N=text0/F=0/A=MC/X=-42.20/Y=44.27 "\\Z16Median"
		AddExpinfo2Graph()
	Endif
	
End Function

//***************************************************************************
Function CalcWhiskerValues(InputWv, InputStr)
	// Created by JAH 2014.04.28
	Wave InputWv
	String InputStr
	String WvNameStr = InputStr + "Wv"
	Make/o/n=7 NewWv
	
	Duplicate/o InputWv SortedWv
	Sort SortedWv, SortedWv
	Wavestats/q SortedWv
	Deletepoints V_npnts, 1e7, SortedWv
				
	NewWv[0] = V_avg
	NewWv[1] = SortedWv[round((V_npnts/2) -1)]
	NewWv[2] = SortedWv[round((v_npnts*1/20) -1)]
	NewWv[3] = SortedWv[round((v_npnts*19/20) -1)]
	NewWv[4] = SortedWv[round((v_npnts*1/4) -1)]
	NewWv[5] = SortedWv[round((v_npnts*3/4) -1)]
	NewWv[6] = V_sdev
	
	Make/o/t/n=7 WhiskerKeyWv
	WhiskerKeyWv[0] = "Mean"
	WhiskerKeyWv[1] = "Median"
	WhiskerKeyWv[2] = "5th Percentile"
	WhiskerKeyWv[3] = "95th Percentile"
	WhiskerKeyWv[4] = "25th Percentile"
	WhiskerKeyWv[5] = "75th Percentile"
	WhiskerKeyWv[6] =  "Standard Dev."
	
	Duplicate/o NewWv $WvNameStr
	Killwaves/z NewWv, SortedWv

End Function

//***************************************************************************
Function CalcWhiskerValuesRange(InputWv, InputStr, FromPtA, ToPtB)
	// Created by JAH 2014.04.28
	Wave InputWv
	String InputStr
	Variable FromPtA, ToPtB
	String WvNameStr = InputStr + "Wv"
	Make/o/n=7 NewWv
	
	Duplicate/o InputWv SortedWv
	Deletepoints ToPtB+1, 1e7, SortedWv
	Deletepoints 0,FromPtA, SortedWv
	Sort SortedWv, SortedWv
	Wavestats/q/z SortedWv
	//Deletepoints V_npnts, 1e7, SortedWv
				
	NewWv[0] = V_avg
	NewWv[1] = SortedWv[round((V_npnts/2) -1)]
	NewWv[2] = SortedWv[round((v_npnts*1/20) -1)]
	NewWv[3] = SortedWv[round((v_npnts*19/20) -1)]
	NewWv[4] = SortedWv[round((v_npnts*1/4) -1)]
	NewWv[5] = SortedWv[round((v_npnts*3/4) -1)]
	NewWv[6] = V_sdev
	
	Make/o/t/n=7 WhiskerKeyWv
	WhiskerKeyWv[0] = "Mean"
	WhiskerKeyWv[1] = "Median"
	WhiskerKeyWv[2] = "5th Percentile"
	WhiskerKeyWv[3] = "95th Percentile"
	WhiskerKeyWv[4] = "25th Percentile"
	WhiskerKeyWv[5] = "75th Percentile"
	WhiskerKeyWv[6] =  "Standard Dev."
	
	Duplicate/o NewWv $WvNameStr
	Killwaves/z NewWv, SortedWv

End Function
//***************************************************************************
Function ExtractNumbersfromWIBS(SizeRangeParam, FromPt, ToPt) 
	// Created 12.16.2014 by JAH
		// Function - To 'automate' calculation of integrated number concentration values from WIBS for comparison
			// Accesses toDoList from WASP panel
			// ** IMPORTANT: Each size distribution must already have been plotted at least once or will return error for that ToDo
			// ** ALSO needs existing wave called root:LogSizeWV

	String SizeRangeParam // Input string to name resultant wave with size range
	Variable FromPt // = Cell/point within LogSizeWv to serve as beginning point of integration
	Variable ToPt // = Cell/point within LogSizeWv to serve as end point of integration
	
	svar ListStr = root:panel:toDoList
	String LocalToDoStr
	
	Variable LoopLengthVar = ItemsInList(ListStr)
	Variable LoopCounterVar
	String AllConcWv = "All_" + SizeRangeParam + "_Conc" // These are the name of the key waves
	String FlConcWv = "Fl_" + SizeRangeParam + "_Conc"
	Make/o/n=(LoopLengthVar) $AllConcWv, $FlConcWV
	Make/o/t/n=(LoopLengthVar) ToDoStrWv
	Wave LocalAllConcWv = $AllConcWv
		LocalAllConcWv = nan
	Wave LocalFlConcWv = $FlConcWv
		LocalFlConcWv = nan
	Make/o/n=2 PointerWv	// This serves as the vehicle to bring info back from the sub-function
		PointerWv = nan
		
	For (LoopCounterVar = 0 ; LoopCounterVar < LoopLengthVar ; LoopCounterVar += 1)
		LocalToDoStr = stringfromlist(LoopCounterVar,ListStr)	 // Grabs value from ToDo list
		ToDoStrWv[LoopCounterVar] = LocalToDoStr
		ExtractNumbersfromWIBS_sub(LocalToDoStr, SizeRangeParam, FromPt, ToPt) // *Calls sub-function 
		LocalAllConcWv[LoopCounterVar] = PointerWv[0] // Stores values from sub-function in final wave
		LocalFlConcWv[LoopCounterVar] = PointerWv[1]
	Endfor
	
	Edit/W=(42,47.75,379.5,624.5) ToDoStrWv, LocalAllConcWv, LocalFlConcWv
	Killwaves/z TempAllAvgColumn, TempFlAvgColumn
	Killwaves/z PointerWv
End Function
//*******************
Function ExtractNumbersfromWIBS_sub(ToDoWvParam, SizeRangeParam, FromPt, ToPt) 
		// Called by function ExtractNumbersfromWIBS(..)
		// Main calculations done in this function
			// Hard-wired to run twice -> Once for "all" once for "Fl"
			// To hard-wire to do more times, add additional string/wave declarations and also increase length of PointerVar
	String ToDoWvParam	// Input string that is just the name of a given ToDo wave that already exists
	String SizeRangeParam // Input string to name resultant wave with size range
	Variable FromPt // = Cell/point within LogSizeWv to serve as beginning point of integration
	Variable ToPt // = Cell/point within LogSizeWv to serve as end point of integration
	
	Wave LogSizeWv // Needs LogSizeWv = Log of beginning point of each wave 
	String LocalAllWvStr = "All_SD_" + ToDoWvParam + "avg"
	String LocalFlWvStr = "Fl_SD_" + ToDoWvParam + "avg"
	Wave LocalAllWv = root:DataProducts:$LocalAllWvStr
	Wave LocalFlWv = root:DataProducts:$LocalFlWvStr
	
	Make/o/d/n=(DimSize(LocalAllWv, 0)) TempAllAvgColumn, TempFlAvgColumn
	TempAllAvgColumn = LocalAllWv[p][0] // Takes first column of 2-D wave created by WASP
	TempFlAvgColumn = LocalFlWv[p][0]
		
	//Variable PointerVarAll, PointerVarFl
	Wave PointerWv
	PointerWv[0] = areaxy(LogSizeWv, TempAllAvgColumn, LogSizeWv[FromPt], LogSizeWv[ToPt])	// This is the integration step
	PointerWv[1] = areaxy(LogSizeWv, TempFlAvgColumn, LogSizeWv[FromPt], LogSizeWv[ToPt])
	
	//Print ToDoWvParam + ", " + SizeRangeParam + ": All = " + num2str(PointerVarAll) + ", Fl = " + num2str(PointerVarFl) + " #/liter"
		
End Function
//***************************************************************************
Function ConvertNumberMxToMass(InputNumberMx, InputDiamWv, InputDensityVar)
	// Written 2015.07.31 by JAH
	Wave InputNumberMx // Must be: time = columns, size = rows [see #2 below for opposite]
	Wave InputDiamWv
	Variable InputDensityVar
	Variable rowVar, colVar
	

	Duplicate/o InputNumberMx NewMassMx
	
	For (colVar = 0 ; colVar < dimsize(InputNumberMx, 1) ; colVar += 1)
		For (rowVar = 0 ; rowVar < numpnts(InputDiamWv) ; rowVar += 1)
			variable probe1 = inputdensityvar
			variable probe2 = inputdiamwv[rowvar]
			variable probe3 = inputnumbermx[rowvar][colvar] 
			NewMassMx[rowVar][colVar] = (4*pi*InputDensityVar/3) * (InputDiamWv[rowVar]/2)^3 * (InputNumberMx[rowVar][colVar])
				// Previous line is the crux: 4/3 pi r^3 N
			variable probe4 = newmassmx[rowvar][colvar]
			variable probe5 = (4*pi*probe1/3) * (probe2/2)^3 * (probe3)
		Endfor
	Endfor 
	
	
End Function

//***************************************************************************
Function ConvertNumberMxToMass2(InputNumberMx, InputDiamWv, InputDensityVar)
	// Written 2015.07.31 by JAH --> modified 2019.02.12 by JAH to opposite orientation
	Wave InputNumberMx // Must be: time = rows, diam/size = columns
	Wave InputDiamWv
	Variable InputDensityVar
	Variable rowVar, colVar
	

	Duplicate/o InputNumberMx NewMassMx
	
	For (rowVar = 0 ; rowVar < dimsize(InputNumberMx, 0) ; rowVar += 1)
		For (colVar = 0 ; colVar < numpnts(InputDiamWv) ; colVar += 1)
			//variable probe1 = inputdensityvar
			//variable probe2 = inputdiamwv[rowvar]
			//variable probe3 = inputnumbermx[rowvar][colvar] 
			NewMassMx[rowVar][colVar] = (4*pi*InputDensityVar/3) * (InputDiamWv[colVar]/2)^3 * (InputNumberMx[rowVar][colVar])
				// Previous line is the crux: 4/3 pi r^3 N
			//variable probe4 = newmassmx[rowvar][colvar]
			//variable probe5 = (4*pi*probe1/3) * (probe2/2)^3 * (probe3)
		Endfor
	Endfor 
	
	
End Function

//***************************************************************************
Function ConvertNumberMxToSurfArea2(InputNumberMx, InputDiamWv)
	// Written 2019.02.12 by JAH --> nearly equiv. to and modified from ConvertNumberMxtoMass2
	Wave InputNumberMx // Must be: time = rows, diam/size = columns
	Wave InputDiamWv
	Variable rowVar, colVar
	

	Duplicate/o InputNumberMx NewSurfAreaMx
	
	For (rowVar = 0 ; rowVar < dimsize(InputNumberMx, 0) ; rowVar += 1)
		For (colVar = 0 ; colVar < numpnts(InputDiamWv) ; colVar += 1)
			//variable probe1 = inputdensityvar
			//variable probe2 = inputdiamwv[rowvar]
			//variable probe3 = inputnumbermx[rowvar][colvar] 
			NewSurfAreaMx[rowVar][colVar] = (4*pi) * (InputDiamWv[colVar]/2)^2 * (InputNumberMx[rowVar][colVar])
				// Previous line is the crux: 4/3 pi r^3 N
			//variable probe4 = newmassmx[rowvar][colvar]
			//variable probe5 = (4*pi*probe1/3) * (probe2/2)^3 * (probe3)
		Endfor
	Endfor 
	
	
End Function

//***************************************************************************

Function BatchAverage(WaveToAverage, FromPtWv, ToPtWv, SampleNameWv)
	// Written 2015.07.31 by JAH
	// Purpose: To average a given wave (e.g. time series) into discreet values, but do some many times automatically
	Wave WaveToAverage, FromPtWv, ToPtWv
	Wave/t SampleNameWv
	Variable n

	For (n = 0 ; n <= numpnts(FromPtWv) ; n += 1)
		CalcWhiskerValuesRange(WaveToAverage, SampleNameWv[n], FromPtWv[n], ToPtWv[n])
	Endfor
	
End Function
//***************************************************************************

Function BatchAverage2(WaveToAverage, FromPtWv, ToPtWv, SampleNameWv)
	// Written 2015.07.31 by JAH
	//	modified 2019.02.11 by JAH
	// Purpose: To average a given wave (e.g. time series) into discreet values, but do some many times automatically
	Wave WaveToAverage, FromPtWv, ToPtWv
	Wave/t SampleNameWv
	Variable n

	For (n = 0 ; n <= numpnts(FromPtWv) ; n += 1)
		If (numtype(FromPtWv[n]) == 0)
				CalcWhiskerValuesRange(WaveToAverage, SampleNameWv[n], FromPtWv[n], ToPtWv[n])
		Endif
	Endfor
	
End Function


//***************************************************************************

Function CollectBatchAverages(WaveToAverage, FromPtWv, ToPtWv, SampleNameWv, SubFolderStr)
	// Written 2015.07.31 by JAH
	//	Formerly called "MarieFunction"
	// Purpose: To collect data run from BatchAverage(..) into one wave
	Wave WaveToAverage, FromPtWv, ToPtWv
	Wave/t SampleNameWv
	String SubFolderStr
	Variable n

	BatchAverage2(WaveToAverage, FromPtWv, ToPtWv, SampleNameWv)
	
	Duplicate/o FromPtWv MeanWv, MedianWv, StDevWv
	For (n = 0 ; n < numpnts(FromPtWv) ; n += 1)
		If (numtype(FromPtWv[n]) == 0)
			String GoGetStr = SampleNameWv[n] + "Wv"
			Wave GotGetWv = $GoGetStr
			MeanWv[n] = GotGetWv[0]
			MedianWv[n] = GotGetWv[1]
			StDevWv[n] = GotGetWv[6]
			Killwaves/z GotGetWv
		Endif
	Endfor	
	
 	Duplicate/o MeanWv root:$(SubFolderStr):MeanWv
 	Duplicate/o MedianWv root:$(SubFolderStr):MedianWv
 	Duplicate/o StDevWv root:$(SubFolderStr):StDevWv
 	
 	//This part can be un-commented if run for only one set of waves (i.e. not "CollectBatchAveragesManyWaves")
 	//Killwaves/z MeanWv, MedianWv, StDevWv, WhiskerKeyWv

End Function

//***************************************************************************

Function CollectBatchAveragesManyWaves(FromPtWv, ToPtWv, SampleNameWv, SubFolderStr, ListWv)
	// Written 2019.02.11 by JAH
		// Purpose: To run BatchAverage(..) [many time-averages per waves] on many different waves
	Wave FromPtWv, ToPtWv
	Wave/t SampleNameWv
	String SubFolderStr
	Wave/t ListWv // This is a text wave with a list of all the time series to be averaged
	String ListWvStr
	//Wave PointerWv
	Variable n
	
	For (n = 0 ; n <= numpnts(ListWv) ; n += 1)
		ListWvStr = ListWv[n] // only for diagnostics
		//PointerWv = $ListWvStr // THIS IS THE PROBLEM
		CollectBatchAverages($ListWv[n], FromPtWv, ToPtWv, SampleNameWv, SubFolderStr)
		String MedianWvStr, MeanWvStr, StDevWvStr
		MeanWvStr = ListWv[n] + "_MeanWv"
		MedianWvStr = ListWv[n] + "_MedianWv"
		StDevWvStr = ListWv[n] + "_StDevWv"
		Duplicate/o MeanWv root:$(SubFolderStr):$MeanWvStr
		Duplicate/o MedianWv root:$(SubFolderStr):$MedianWvStr
		Duplicate/o StDevWv root:$(SubFolderStr):$StDevWvStr
	Endfor
	
	Killwaves/z MeanWv, MedianWv, StDevWv, WhiskerKeyWv
		
End Function

//***************************************************************************
Function AverageAroundSpecificPts(XWvInput, YWvInput, RefPtsWv, AvgWindowVar)
	// Created 2018.01.05 by JAH
	// Function created for Ben to use to average certain points in spectra for use as cluster inputs
	// The idea is for this to do y-averages around a few discrete x-points
	//	User can decide which x-points to use and which y-waves to input, by communicating these as parameters
	//	Note that the names of the wave parameters need not be the names of the waves themselves
	//	New waves (see below) will be overwritten each new run, so be careful to rename as necessary
	
	Wave XWvInput, YWvInput // These are input x and y waves
	Wave RefPtsWv // These are the values that you want as x-centers of an average
					// Values in RefPtsWv should be wave points, not x-values
	Variable AvgWindowVar // This dictates the number of points on either side to average
	Variable n, startpt, endpt
	
	Make/o/n=(numpnts(RefPtsWv)) AvgPtsWv, StDevPtsWv	// These will be the names of the output waves, re-name as necessary
	AvgPtsWv = nan
	StDevPtsWv = nan
	
	For (n = 0 ; n < numpnts(RefPtsWv) ; n += 1)
		startpt = RefPtsWv[n] - AvgWindowVar
		endpt = RefPtsWv[n] + AvgWindowVar
		WaveStats/q/r=(startpt, endpt) YWvInput
		AvgPtsWv[n] = V_avg
		StDevPtsWv[n] = V_sdev
	Endfor
	
End Function	


//***************************************************************************
