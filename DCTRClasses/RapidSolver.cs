using FunctionLibrary;
using Serilog;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    /// <summary>
    /// Solver which implements the DCTR algorithm
    /// </summary>
    public class RapidSolver : IRHTSolver
    {
        #region Public Methods

        public void Echo()
        {
            Log.Information($"Num width levels={Settings.NWidthLevels}");
            Log.Information($"Num bins split={Settings.NShifts}");
        }

        public void FinalizeCalcs(bool requiresOverlap)
        {
            // requiresOverlap is true when the full spectral interval is too large to handle in one piece
            // In this case, the code will split the full interval into overlapping sub-intervals
            // Finally, these overlapping sub-intervals are combined into one output array
            if (requiresOverlap) generateFinalOutputArray();

            // This step writes the final output array to disk
            Write(Parent.OutputFolder);

            _outArray = null;
        }

        public void Initialize(string baseFolder)
        {
            _currentRunID = Guid.NewGuid();

            _outArray = null;

            _dopplerEngine.SetMolecularWeights(baseFolder);

            _widthRanges = new double[Settings.NWidthLevels];
            List<int> shiftList = [];
            for (int i = -Settings.NShifts / 2 + 1; i <= Settings.NShifts / 2; i++)
            {
                shiftList.Add(i);
            }

            // Geometric width ratio based on minimum and maximum Lorentz half-widths
            // in the inhomogeneous path, and the number of discretization levels
            _widthRatio = Math.Pow(MaxWidth / MinWidth, 1.0 / Settings.NWidthLevels);

            for (int i = 0; i < Settings.NWidthLevels; i++)
            {
                // Actual representative Lorentz bin half-widths
                _widthRanges[i] = i == 0 ? MinWidth : _widthRanges[i - 1] * WidthRatio;

                // Weights for splitting the line intensity among neighboring bins
                // The dictionary defines one set of weights for each Lorentz half-width bin
                _splitWeightsDict.Add(i, new LorentzSplitWeights());
                _splitWeightsDict[i].Lambda = _widthRanges[i] / Parent.WaveNumberResolution;

                // Bin indices - 0 refers to binned line center (see DCTR paper)
                _splitWeightsDict[i].ShiftList = shiftList;
            }

            TotalLineCount = 0;
            TotalLineCalculationCountEstimate = 0;
        }

        // Initializes from the run specification file, which contains
        // all the run parameters, such as inhomogeneous path specification, database
        // specification, etc.
        // This is usually "DCTRInput.dat"
        public void Initialize(List<string> inpLineList)
        {
            foreach (string curLine in inpLineList)
            {
                string line = curLine;

                while (line.Contains("  ")) line = line.Replace("  ", " ");
                line = line.Replace(" ", "\t");
                string[] lineList = line.Split('\t');

                if (lineList[0] == "NWidths")
                {
                    Settings.NWidthLevels = int.Parse(lineList[1]);

                    if (Settings.NWidthLevels < 2) Settings.NWidthLevels = 2;
                    if (Settings.NWidthLevels > 1024) Settings.NWidthLevels = 1024;
                }

                if (lineList[0] == "NBinSplit")
                {
                    Settings.NShifts = int.Parse(lineList[1]);

                    if (Settings.NShifts % 2 == 1) Settings.NShifts++;
                    if (Settings.NShifts < 2) Settings.NShifts = 2;
                    if (Settings.NShifts > 8) Settings.NShifts = 8;
                }
            }
        }

        // Workhorse function, performs the actual run
        public void Run(double startWNum, double endWNum, bool requiresOverlap, StateMachine stateMachine,
            Dictionary<string, DataBuffer> bufferDict)
        {
            TimingFunctions.InitializeTime("RHTC: Rapid calcs");

            startWNum += Parent.StartWNumOffset;
            endWNum += Parent.StartWNumOffset;

            _outArray = new double[Settings.NWidthLevels][];

            Log.Information($"DCTR: Setting initial parameters");

            // Initialization
            setInitialParameters(startWNum, endWNum, requiresOverlap);

            int[] indexList = [-1, -1];

            foreach (string species in bufferDict.Keys)
            {
                if (!stateMachine.StateDefinitions.ContainsKey(species)) continue;

                List<RunEngineState> stateList = stateMachine.StateDefinitions[species];

                DataBuffer dBuffer = bufferDict[species];

                indexList[0] = -1;
                indexList[1] = -1;

                dBuffer.GetNextIndexList(indexList, startWNum, endWNum);
                while (indexList[1] != -1)
                {
                    foreach (RunEngineState state in stateList)
                    {
                        // "state" represents a homogeneous slab within the inhomogeneous path
                        state.Initialize(species, dBuffer.SpeciesDatabaseDef);

                        // Update line parameters (such as shifted line center, scaled intensity, etc.)
                        // based on the current state (homogeneous slab)
                        setLineCurrentStates(dBuffer.LineEnsemble, indexList[0] + 1, indexList[1] + 1, state, species);

                        // For debugging only
                        if (dBuffer.DebugMode)
                            dBuffer.LineEnsemble.AppendToFile(dBuffer.DebugFileName, indexList[0], indexList[1],
                            Parent.IntensityCutoff);

                        // Gets list of indices of lines which fall within each discretized Lorentz half-width bin
                        Dictionary<int, List<int>> lineDict = updateLineDictionary(
                            dBuffer.LineEnsemble, indexList[0] + 1,
                            indexList[1] + 1, out long lCount);

                        performUpdate(dBuffer.LineEnsemble, lineDict, state, species);

                        int NWNums = (int)(Parent.WaveNumberSpread / Parent.WaveNumberResolution);
                        if (NWNums < 0) NWNums = Parent.WaveNumberRange.Length <= 2 ? _NBaseArray : _NArray;

                        TotalLineCount += lCount;
                        TotalLineCalculationCountEstimate += lCount * NWNums;
                    }

                    dBuffer.GetNextIndexList(indexList, startWNum, endWNum);
                }
            }

            TimingFunctions.AddTime("RHTC: Rapid calcs");
        }

        // This could be the final output array (for smaller spectral interval runs)
        // Or it could be a sub-array (for larger spectral interval runs, which are
        // broken down into sub-intervals)
        public void SaveCurrentArray(bool requiresOverlap)
        {
            performFreqDomainConvolution();

            if (!requiresOverlap || _usedWidthIndices.Count == 0 || _outArray[0] == null)
            {
                return;
            }

            string fName = $"{Parent.BaseFolder}\\Temp\\{_currentRunID}";
            fName += $"_{Parent.CurWNumRangeIndex}.dat";

            using BinaryWriter sw = new(
                new FileStream(fName, FileMode.Create));
            sw.Write(_arrayStartWNum);
            sw.Write(Parent.WaveNumberResolution);
            sw.Write(_NArray);

            for (int i = 0; i < _NArray; i++)
            {
                sw.Write(_outArray[0][i]);
            }
        }

        public void SetWaveNumberRange()
        {
            double NN = (Parent.EndWaveNumber - Parent.StartWaveNumber) / Parent.WaveNumberResolution;

            int NMax = 1 << 20;

            double NLevels = NN / NMax;

            int NLevelsInt = NLevels < 1 ? 1 : (int)(2 * NLevels + 1);

            Parent.WaveNumberRange = new double[NLevelsInt + 1];

            double dWNums = (int)((Parent.EndWaveNumber - Parent.StartWaveNumber) / NLevelsInt);

            for (int i = 0; i <= NLevelsInt; i++)
            {
                Parent.WaveNumberRange[i] = Parent.StartWaveNumber + dWNums * i;

                if (Parent.WaveNumberRange[i] > Parent.EndWaveNumber) Parent.WaveNumberRange[i] = Parent.EndWaveNumber;
            }
        }

        // Updates the output arrays with either the line intensity value (for Lorentz profile runs)
        // or the Doppler profile scaled by the line intensity value (for Voigt profile runs)
        public void Update(LineEnsemble lEnsemble, int ix, RunEngineState state, string species)
        {
            int ixSingle = ix * LineEnsemble.NDataColumns;
            int ixIndex = ix * LineEnsemble.NIndexColumns;

            int index = lEnsemble.IndexArray[ixIndex];
            int widthIndex = lEnsemble.IndexArray[ixIndex + 1];

            double wCenter = lEnsemble.DataArray[ixSingle + 8];
            double curIntensity = lEnsemble.DataArray[ixSingle + 9];
            double curWidth = lEnsemble.DataArray[ixSingle + 10];
            int isotopeID = lEnsemble.IndexArray[ixIndex + 2];

            double[]? dArray = null;

            double repWNum = _arrayStartWNum + index * Parent.WaveNumberResolution;

            double kk = (wCenter - repWNum) / Parent.WaveNumberResolution;

            double[] weights = _splitWeightsDict[widthIndex].GetWeights(kk);
            List<int> shiftList = _splitWeightsDict[widthIndex].ShiftList;

            if (double.IsNaN(weights.Sum()))
            {
                weights = [1 - kk, kk];
                shiftList = [0, 1];
            }

            // For Doppler and Voigt profiles, the Doppler array is
            // calculated and inserted into the appropriate locations in
            // the profile array
            // dEngine being null implies that the Doppler array is a stick function
            if (Parent.Profile != Profile.Lorentz && _dopplerEngine != null)
            {
                Result<double[]?> dArrayRes = _dopplerEngine.GetDistributionArray(state.Temperature,
                wCenter, species, isotopeID, Parent.WaveNumberResolution);

                if (dArrayRes.TryGetValue(out double[]? dVal))
                    dArray = dVal;
            }

            if (index >= _NArray) return;
            if (index < 0) return;

            if (_outArray[widthIndex] == null) _outArray[widthIndex] = new double[_NArray];

            // dArray can be null if the profile is Lorentz, or
            // if the Doppler width is too small
            if (dArray != null)
            {
                double mult = curIntensity;

                //Result<double[]?> dArrayRes = _dopplerEngine.GetDistributionArray(state.Temperature,
                //wCenter, species, isotopeID, Parent.WaveNumberResolution);

                if (Parent.Profile == Profile.Voigt)
                {
                    mult *= Parent.WaveNumberResolution;
                }

                double v;
                int ixx;

                int wNumBinShiftIx = 0;
                foreach (int wNumBinShift in shiftList)
                {
                    if (index + wNumBinShift < _NArray)
                    {
                        _outArray[widthIndex][index + wNumBinShift] +=
                            dArray[0] * mult * weights[wNumBinShiftIx];
                    }

                    wNumBinShiftIx++;
                }

                // Splits the Doppler intensity array into several arrays centered
                // over a number of bins (see DCTR paper)
                for (int i = 1; i < dArray.Length; i++)
                {
                    v = dArray[i] * mult;

                    wNumBinShiftIx = 0;
                    foreach (int wNumBinShift in shiftList)
                    {
                        ixx = index - i + wNumBinShift;

                        if (ixx >= 0 && ixx < _NArray) _outArray[widthIndex][ixx] += v * weights[wNumBinShiftIx];

                        ixx = index + i + wNumBinShift;
                        if (ixx >= 0 && ixx < _NArray) _outArray[widthIndex][ixx] += v * weights[wNumBinShiftIx];

                        wNumBinShiftIx++;
                    }
                }
            }
            else
            {
                // In this case, there is no Doppler array, just a single intensity value
                // This intensity value is split among several bins neighboring the bin which
                // is closest to the line center
                int wNumBinShiftIx = 0;
                foreach (int wNumBinShift in shiftList)
                {
                    int ixx = index + wNumBinShift;
                    if (ixx >= 0 && ixx < _NArray) _outArray[widthIndex][ixx] += curIntensity * weights[wNumBinShiftIx];

                    wNumBinShiftIx++;
                }
            }

            // Currently not used
            _intensitySumList[widthIndex] += curIntensity;
            _intensityWidthSumList[widthIndex] += curIntensity / curWidth;
        }

        /// <summary>
        /// Time reporting, etc.
        /// </summary>
        /// <param name="reportDict"></param>
        public void UpdateReportDictionary(Dictionary<string, string> reportDict)
        {
            Dictionary<string, double> cpuTimes = TimingFunctions.GetCPUTimes();

            Log.Information(string.Join(Environment.NewLine, cpuTimes.Keys));

            if (!cpuTimes.TryGetValue("RHTC: DataBuffer: GetData", out double t1))
            {
                t1 = 0;
            }

            double t2 = cpuTimes["RHTC: Rapid calcs"] - t1;
            cpuTimes["RHTC: Rapid calcs"] = t2;
            TimingFunctions.UpdateTime("RHTC: Rapid calcs", t2);
            double t3 = cpuTimes["RHTC: Rapid solver finalization"];

            List<RunEngineState> stateList = Parent.StateMachine.StateDefinitions.FirstOrDefault().Value;
            int NSlabs = stateList.Count;

            reportDict.Add("Total calculation time", (t1 + t2 + t3) + " s");
            reportDict.Add("Wavenumber array count", _NBaseArray.ToString("N0"));
            reportDict.Add("Calculation array count", _NArray.ToString("N0"));
            reportDict.Add($"Total line count (for {NSlabs} homogeneous slabs)",
                TotalLineCount.ToString("N0"));
            reportDict.Add("Average line count per homogeneous slab",
                (TotalLineCount / (double)NSlabs).ToString("N0"));
            reportDict.Add("Total line calculation count (estimate)",
                TotalLineCalculationCountEstimate.ToString("N0"));

            double rate1 = TotalLineCalculationCountEstimate / (t1 + t2 + t3);
            double rate2 = TotalLineCalculationCountEstimate / (t2 + t3);
            double rate3 = TotalLineCalculationCountEstimate / t3;

            reportDict.Add("End time", DateTime.Now.ToString());

            reportDict.Add("Processing rate (with file read and array set-up)", rate1.ToString("E"));
            reportDict.Add("Processing rate (with array set-up)", rate2.ToString("E"));
            reportDict.Add("Processing rate (FDCT)", rate3.ToString("E"));
        }

        public void Write(string folderName)
        {
            string fName = folderName + "\\absorptances.dat";

            using StreamWriter sw = new(fName);
            int nOffset = (_NArray - _NBaseArray) / 2;

            // When spectrum is broken down into regions and recombined,
            // there is no offset
            if (_outArray[0] == null || _outArray[0].Length == _NBaseArray) nOffset = 0;

            _arrayStartWNum = Parent.StartWaveNumber - nOffset * Parent.WaveNumberResolution;

            for (int i = 0; i < _NBaseArray; i++)
            {
                //if (_arrayStartWNum < 0) _arrayStartWNum = 0 + Parent.StartWNumOffset;

                double wNum = _arrayStartWNum + (i + nOffset) * Parent.WaveNumberResolution;
                double aVal = _outArray[0] != null ? _outArray[0][i + nOffset] : 0;

                sw.WriteLine($"{wNum}\t{aVal}");
            }
        }

        public void WriteDetails(string folder)
        {
            using StreamWriter sw = new($"{folder}\\Widths.dat");
            sw.WriteLine($"Nominal number of widths\t{Settings.NWidthLevels}");
            sw.WriteLine($"Widths\tMin={MinWidth}; Max={MaxWidth}");
            sw.WriteLine($"Actual number of widths\t{_usedWidthIndices.Count}");

            List<string> wList = [.. _actualWidths.Select((item, ix) =>
                ix.ToString() + " " + item.ToString()).Where((item, ix) =>
                       _usedWidthIndices.Contains(ix))];

            sw.WriteLine($"Actual widths\t{string.Join(Environment.NewLine + "\t", wList)}");
        }

        public void WriteSettings(string folder)
        {
            Settings.WriteToFolder(folder);
        }

        #endregion Public Methods

        #region Public Properties

        // Maximum value of Lorentz half-width of all lines within the inhomogeneous gas path
        // This can be pre-determined and tabulated (see DCTR paper)
        public double MaxWidth
        {
            get { return _maxWidth; }
            set { _maxWidth = value; }
        }

        // Minimum value of Lorentz half-width of all lines within the inhomogeneous gas path
        // This can be pre-determined and tabulated (see DCTR paper)
        public double MinWidth
        {
            get { return _minWidth; }
            set { _minWidth = value; }
        }

        // The solver instance resides within a parent RunEngine instance
        // The Parent property below provides a link to the run engine, which
        // contains the solver
        public RunEngine? Parent
        {
            get { return _parent; }
            set { _parent = value; }
        }

        public RapidSolverSettings Settings
        {
            get { return _settings; }
            set { _settings = value; }
        }

        // Keeps track of total number of line-locations handled during the run
        // (See DCTR paper)
        public long TotalLineCalculationCountEstimate
        {
            get { return _totalLineCalculationCountEstimate; }
            set { _totalLineCalculationCountEstimate = value; }
        }

        // Keeps track of total number of lines handled during the run
        public long TotalLineCount
        {
            get { return _totalLineCount; }
            set { _totalLineCount = value; }
        }

        // Geometric binning width ratio (see DCTR paper)
        public double WidthRatio
        {
            get { return _widthRatio; }
        }

        #endregion Public Properties

        #region Private Methods

        /// <summary>
        /// For large spectral intervals, the code splits the overall interval into sub-intervals
        /// Each sub-interval generates a separate output file
        /// In this method, the individual output files are combined into an overall array for the full
        /// spectral interval, and this final output is saved
        /// </summary>
        void generateFinalOutputArray()
        {
            _NBaseArray = RHTFunctions.GetArrayLength(Parent.StartWaveNumber,
                Parent.EndWaveNumber, Parent.WaveNumberResolution);

            // Final output array
            _outArray = new double[1][];
            _outArray[0] = new double[_NBaseArray];

            // List of temporary files with sub-interval output arrays
            string[] fileList = Directory.GetFiles(Parent.BaseFolder + "\\Temp", _currentRunID.ToString() + "*.dat");

            double minWNum = double.PositiveInfinity, maxWNum = double.NegativeInfinity;

            // Combine to final output array
            foreach (string fName in fileList)
            {
                using (BinaryReader br = new(
                    new FileStream(fName, FileMode.Open)))
                {
                    double wNumStart = br.ReadDouble();
                    double wRes = br.ReadDouble();

                    int NVals = br.ReadInt32();

                    if (wNumStart < minWNum) minWNum = wNumStart;
                    if (wNumStart + NVals * wRes > maxWNum) maxWNum = wNumStart + NVals * wRes;

                    // Starting index in final output array
                    int ix = (int)Math.Round((wNumStart - Parent.StartWaveNumber) / Parent.WaveNumberResolution);

                    for (int i = 0; i < NVals; i++)
                    {
                        if (ix < 0)
                        {
                            double v = br.ReadDouble();

                            ix++;

                            continue;
                        }

                        if (ix > _NBaseArray - 1) break;

                        _outArray[0][ix++] += br.ReadDouble();
                    }
                }

                File.Delete(fName);
            }

            _NArray = RHTFunctions.GetArrayLength(minWNum,
                maxWNum, Parent.WaveNumberResolution);
        }

        void performFreqDomainConvolution()
        {
            TimingFunctions.InitializeTime("RHTC: Rapid solver finalization");

            // Contains the actual binned Lorentz half-widths used in the computation
            _actualWidths = new double[Settings.NWidthLevels];

            // _usedWidthIndices is a hash set of all the Lorentz half-width indices which
            // were actually encountered during the computation of the inhomogeneous gas path
            List<int> usedWidthIndices = [.. _usedWidthIndices];
            int NIndices = usedWidthIndices.Count;

            // _NArray * 2
            int NArray2 = _NArray << 1;

            // _NArray / 2
            int NArray_2 = _NArray >> 1;

            // Placeholder for later improvement
            string widthCalc = "Geometric";

            foreach (int ix in usedWidthIndices)
            {
                // Widths are calculated either geometrically or harmonically
                if (widthCalc == "Geometric")
                    _actualWidths[ix] = MinWidth * Math.Pow(WidthRatio, ix + 0.5);
                else
                    _actualWidths[ix] = _intensitySumList[ix] / _intensityWidthSumList[ix];
            }

            // Placeholder for later improvement (the flag which decides whether or not
            // to parallelize the calculations, can be set based on run parameters)
            bool parallelize = true;

            if (parallelize)
            {
                // At most 4 parallel threads will be spawned
                int nGroups = 4;

                Parallel.For(0, nGroups, k =>
                {
                    FastDCT dCont1 = new(NArray2);

                    for (int j = k; j < NIndices; j += nGroups)
                    {
                        int ix = usedWidthIndices[j];

                        if (Parent.Profile != Profile.Doppler)
                        {
                            // Representative Lorentz profile for the specified binned half-width
                            double[] lProf = RHTFunctions.GetLorentzProfile(_actualWidths[ix],
                                Parent.WaveNumberResolution, NArray2, Parent.WaveNumberSpread);

                            // Zero-padding to the right of the Lorentz profile
                            for (int i = _NArray; i < NArray2; i++) lProf[i] = 0;

                            double[] newArray = new double[NArray2];
                            // For the output array, the zero-padding is half to the left,
                            // and half to the right of the array
                            for (int i = NArray_2; i < _NArray + NArray_2; i++)
                                newArray[i] = _outArray[ix][i - NArray_2];

                            // Inverse DCT
                            _outArray[ix] = dCont1.IDCT(newArray);

                            // Inverse DCT
                            lProf = dCont1.IDCT(lProf);

                            // Element-wise product of IDCTs
                            for (int i = 0; i < NArray2; i++) _outArray[ix][i] *= lProf[i] * NArray2;
                        }
                    }
                });
            }
            else
            {
                // Same calculations as above, just not parallelized
                for (int j = 0; j < NIndices; j++)
                {
                    int ix = usedWidthIndices[j];

                    if (Parent.Profile != Profile.Doppler)
                    {
                        double[] lProf = RHTFunctions.GetLorentzProfile(_actualWidths[ix],
                            Parent.WaveNumberResolution, NArray2, Parent.WaveNumberSpread);

                        for (int i = _NArray; i < NArray2; i++) lProf[i] = 0;

                        double[] newArray = new double[NArray2];
                        for (int i = NArray_2; i < _NArray + NArray_2; i++)
                        {
                            newArray[i] = _outArray[ix][i - NArray_2];
                        }

                        _outArray[ix] = dCont.IDCT(newArray);

                        lProf = dCont.IDCT(lProf);

                        for (int i = 0; i < NArray2; i++) _outArray[ix][i] *= lProf[i] * NArray2;
                    }
                }
            }

            // _outArray is a 2-d array
            // The following code sums up all the rows, and puts the result in the first row
            for (int j = 1; j < NIndices; j++)
            {
                int ix = usedWidthIndices[j];
                for (int i = 0; i < NArray2; i++)
                    _outArray[usedWidthIndices[0]][i]
                        += _outArray[ix][i];
            }

            if (usedWidthIndices.Count != 0)
            {
                // For Doppler profile runs, there is no deconvolution
                if (Parent.Profile != Profile.Doppler)
                {
                    // Deconvolution consists of taking the DCT of the row sum of the
                    // output array
                    double[] newArray = dCont.DCT(_outArray[usedWidthIndices[0]]);

                    // Take the mid-part of the deconvolved array (ignore the first 1/4th and last 1/4th,
                    // take only the middle 1/2 of the array)
                    // This compensates for the earlier zero-padding
                    _outArray[0] = new double[_NArray];
                    for (int i = NArray_2; i < NArray_2 + _NArray; i++) _outArray[0][i - NArray_2] = newArray[i];
                }
                else
                    _outArray[0] = _outArray[usedWidthIndices[0]];
            }

            TimingFunctions.AddTime("RHTC: Rapid solver finalization");
        }

        void performUpdate(LineEnsemble lEnsemble, Dictionary<int, List<int>> lineDict,
            RunEngineState state, string species)
        {
            List<int> keyList = [.. lineDict.Keys];

            TimingFunctions.InitializeTime("RHTC: Update lines");

            bool parallelize = true;

            if (parallelize)
            {
                int nGroups = 16;

                //Parallel.For(0, nGroups, i =>
                for (int i = 0; i < nGroups; i++)
                {
                    for (int j = i; j < keyList.Count; j += nGroups)
                    {
                        //DopplerEngine dEngine = new();

                        // The full range of Doppler half-widths over the inhomogeneous path is
                        // discretized into 100,000 bins (see DCTR paper)
                        _dopplerEngine.Initialize(Parent.StartWaveNumber, Parent.EndWaveNumber, 100000);

                        int lIx = keyList[j], ixIndex, isotopeID;

                        foreach (int line in lineDict[lIx])
                        {
                            ixIndex = line * LineEnsemble.NIndexColumns;

                            isotopeID = lEnsemble.IndexArray[ixIndex + 2];

                            Update(lEnsemble, line, state, species);
                        }
                    }
                }//);
            }
            else
            {
                for (int i = 0; i < keyList.Count; i++)
                {
                    //DopplerEngine dEngine = new();

                    // The full range of Doppler half-widths over the inhomogeneous path is
                    // discretized into 100,000 bins (see DCTR paper)
                    _dopplerEngine.Initialize(Parent.StartWaveNumber, Parent.EndWaveNumber, 100000);

                    int lIx = keyList[i], ixIndex, isotopeID;

                    foreach (int line in lineDict[lIx])
                    {
                        ixIndex = line * LineEnsemble.NIndexColumns;

                        isotopeID = lEnsemble.IndexArray[ixIndex + 2];

                        Update(lEnsemble, line, state, species);
                    }
                }
            }

            TimingFunctions.AddTime("RHTC: Update lines");
        }

        void setInitialParameters(double startWNum, double endWNum, bool requiresOverlap)
        {
            // Nearest higher power of 2 to the actual number of wavenumber bins
            _NArray = RHTFunctions.Get2PowerLength(startWNum,
                endWNum, Parent.WaveNumberResolution);

            // Actual number of wavenumber bins
            _NBaseArray = RHTFunctions.GetArrayLength(startWNum,
                endWNum, Parent.WaveNumberResolution);

            if (requiresOverlap) _NArray <<= 1;

            // Allocating for twice the actual size, to allow for zero-padding
            dCont = new FastDCT(_NArray << 1);

            int nOffset = (_NArray - _NBaseArray) / 2;

            _arrayStartWNum = startWNum - nOffset * Parent.WaveNumberResolution;
            if (_arrayStartWNum < 0) _arrayStartWNum = 0 + Parent.StartWNumOffset;

            _intensitySumList = new double[Settings.NWidthLevels];
            _intensityWidthSumList = new double[Settings.NWidthLevels];

            _usedWidthIndices = [];
        }

        void setLineCurrentStates(LineEnsemble lEnsemble, int startIndex,
            int endIndex, RunEngineState state, string species)
        {
            TimingFunctions.InitializeTime("RHTC: Set Current State");

            bool parallelize = true;

            if (parallelize)
            {
                int nGroups = 16;

                Parallel.For(0, nGroups, i =>
                {
                    for (int j = i + startIndex; j < endIndex; j += nGroups)
                    {
                        // Updates line parameters for line "j" based on the current
                        // homogeneous slab within the inhomogeneous path
                        lEnsemble.SetCurrentState(j, state, species, _arrayStartWNum,
                            _widthRanges, WidthRatio, Parent.WaveNumberResolution,
                            Parent.IntensityCutoff);
                    }
                });
            }
            else
            {
                for (int i = startIndex; i < endIndex; i++)
                {
                    // Updates line parameters for line "i" based on the current
                    // homogeneous slab within the inhomogeneous path
                    lEnsemble.SetCurrentState(i, state, species, _arrayStartWNum,
                        _widthRanges, WidthRatio, Parent.WaveNumberResolution,
                        Parent.IntensityCutoff);
                }
            }

            TimingFunctions.AddTime("RHTC: Set Current State");
        }

        Dictionary<int, List<int>> updateLineDictionary(
    LineEnsemble lEnsemble, int startIndex,
    int endIndex, out long lCount)
        {
            TimingFunctions.InitializeTime("RHTC: Update dictionary");

            lCount = 0;

            // Key: index of discretized Lorentz half-width bin
            // Value: list of line indices whose half-widths fall within
            // the discretized half-width bin
            Dictionary<int, List<int>> lineDict = [];

            int wNumIndex, widthIndex;

            // For each line in the line ensemble whose index lies between
            // startIndex and endIndex, the below code calculates its half-width
            // The code then determines the index of the discretized Lorentz half-width
            // and adds the line index to the list of line indices falling within the
            // discretized half-width index
            for (int i = startIndex; i < endIndex; i++)
            {
                if (lEnsemble.DataArray[i * LineEnsemble.NDataColumns + 11] < Parent.IntensityCutoff)
                    continue;

                wNumIndex = lEnsemble.IndexArray[i * LineEnsemble.NIndexColumns];
                widthIndex = lEnsemble.IndexArray[i * LineEnsemble.NIndexColumns + 1];

                if (!lineDict.TryGetValue(widthIndex, out List<int>? value))
                {
                    value = [];

                    lineDict[widthIndex] = value;
                }

                value.Add(i);

                if (wNumIndex > 0 &&
                    wNumIndex < _NArray)
                    _usedWidthIndices.Add(widthIndex);

                lCount++;
            }

            TimingFunctions.AddTime("RHTC: Update dictionary");

            return lineDict;
        }

        #endregion Private Methods

        #region Private Properties

        // Actual line widths encountered during the run
        double[] _actualWidths = null;

        double _arrayStartWNum = 0;

        Guid _currentRunID = Guid.Empty;

        // Container class for the Fast (I)DCT algorithm
        FastDCT? dCont = null;

        readonly DopplerEngine _dopplerEngine = new();

        double[]? _intensitySumList = null;
        double[]? _intensityWidthSumList = null;

        // _NBaseArray is the actual array length which goes from the start wavenumber to
        // the end wavenumber with the desired resolution
        // _NArray is the nearest higher power of 2 to _NBaseArray
        int _NArray = 0, _NBaseArray = 0;

        double[][]? _outArray = null;

        readonly Dictionary<int, LorentzSplitWeights> _splitWeightsDict = [];

        // Indices of binned Lorentz half-widths which are actually encountered during the run
        // For example, the run might define 64 width bins geometrically, between 0.01 cm-1 and
        // 1 cm-1. However, lines might not fall into every one of these 64 width bins, and some
        // width bins would have no corresponding lines.
        // The hashset ensures that there are no repeats of width bin indices.
        HashSet<int> _usedWidthIndices = [];

        double[]? _widthRanges = null;

        // Maximum value of Lorentz half-width of all lines within the inhomogeneous gas path
        // This can be pre-determined and tabulated (see DCTR paper)
        double _maxWidth = 0.8;

        // Minimum value of Lorentz half-width of all lines within the inhomogeneous gas path
        // This can be pre-determined and tabulated (see DCTR paper)
        double _minWidth = 0.005;

        // The solver instance resides within a parent RunEngine instance
        // The Parent property below provides a link to the run engine, which
        // contains the solver
        RunEngine? _parent = null;

        RapidSolverSettings _settings = new();

        // Keeps track of total number of line-locations handled during the run
        // (See DCTR paper)
        long _totalLineCalculationCountEstimate = 0;

        // Keeps track of total number of lines handled during the run
        long _totalLineCount = 0;

        // Geometric binning width ratio (see DCTR paper)
        double _widthRatio = 1;

        #endregion Private Properties
    }
}
