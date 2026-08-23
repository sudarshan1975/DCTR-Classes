using Serilog;
using FunctionLibrary;
using System.Diagnostics.CodeAnalysis;

namespace DCTRClasses
{
    /// <summary>
    /// LBL: Line-by-line
    /// DCTR: Discrete cosine transform rapid (algorithm)
    /// The current code does not contain an LBL implementation
    /// So LBL is redundant
    /// </summary>
    public enum RunType
    {
        LBL,
        DCTR
    }

    public enum Profile
    {
        Lorentz,
        Doppler,
        Voigt
    }

    /// <summary>
    /// An engine class which implements DCTR runs
    /// </summary>
    public class RunEngine
    {
        #region Public Methods

        public void Echo()
        {
            update();

            if (!IsValid)
            {
                Log.Warning($"Invalid run engine: Echo");

                return;
            }

            Log.Information($"WNumStart={StartWaveNumber} cm-1");
            Log.Information($"WNumEnd={EndWaveNumber} cm-1");
            Log.Information($"WNum resolution={WaveNumberResolution} cm-1");
            Log.Information($"WNum spread={WaveNumberSpread} cm-1");
            Log.Information($"Profile={Profile}");
            Log.Information($"Intensity cutoff={IntensityCutoff} cm-1");
            Log.Information($"Species count={_speciesDefinitions.Count}");
            Log.Information($"Species={string.Join(",",
                _speciesDefinitions.Select(x => $"{x.Value} {x.Key}"))}");

            Solver.Echo();

            Log.Information($"Path specification:");

            _stateMachine.Echo();
        }

        Dictionary<string, string> generateReportDictionary()
        {
            Dictionary<string, string> reportDict = [];

            reportDict.Add("Run description", Description);
            reportDict.Add("Starting wavenumber", StartWaveNumber + " cm-1");
            reportDict.Add("Ending wavenumber", EndWaveNumber + " cm-1");
            reportDict.Add("Wavenumber resolution", WaveNumberResolution + " cm-1");

            return reportDict;
        }

        public void Initialize(string baseFolder)
        {
            BaseFolder = baseFolder;

            //_dopplerEngine.SetMolecularWeights(BaseFolder);

            ReadFromFile(BaseFolder + @"\DCTRInput.dat");

            if (StateMachine is null)
            {
                Log.Warning($"Could not initialize run engine: null state machine");

                return;
            }

            Solver.Parent = this;

            StateMachine.Parent = this;

            //_folderDefinitions.Add("HITRAN CO2", BaseFolder + @"\Data\HITRAN");
            //_folderDefinitions.Add("HITRAN H2O", BaseFolder + @"\Data\HITRAN");

            _folderDefinitions.Add("HITRAN2020 CO2", BaseFolder + @"\Data\HITRAN2020\CO2");
            _folderDefinitions.Add("HITRAN2020 H2O", BaseFolder + @"\Data\HITRAN2020\H2O");

            _folderDefinitions.Add("HITEMP2010 CO2", BaseFolder + @"\Data\HITEMP2010\CO2");
            _folderDefinitions.Add("HITEMP2010 H2O", BaseFolder + @"\Data\HITEMP2010\H2O");

            _folderDefinitions.Add("CDSD CO2", BaseFolder + @"\Data\CDSD");

            //_folderDefinitions.Add("TEST CO2", BaseFolder + @"\Data\TestDatabase");
            //_folderDefinitions.Add("TEST H2O", BaseFolder + @"\Data\TestDatabase");
        }

        public void Initialize(IList<string> inpList)
        {
            SpeciesDefinitions = [];

            StateMachine = new();
            List<string> stateLines = [];

            List<string> leftoverLines = [];

            int lineIx = 0;

            while (lineIx < inpList.Count)
            {
                string line = inpList[lineIx];

                while (line.Contains("  ")) line = line.Replace("  ", " ");
                line = line.Replace(" ", "\t");
                //line = line.Replace("\t", "=");
                string[] lineList = line.Split('=');

                bool doneFlag = false;

                if (lineList.Length < 2)
                {
                    lineIx++;

                    continue;
                }

                lineList[0] = lineList[0].ToLower();

                if (lineList[0] == "description")
                {
                    Description = lineList[1];

                    doneFlag = true;
                }

                if (lineList[0] == "wavenumberlevels" ||
                    lineList[0] == "wavenumberranges")
                {
                    List<double> WaveNumberRanges = [..
                        lineList.Where((item, ix) => ix > 0).Select(double.Parse)];

                    StartWaveNumber = WaveNumberRanges[0];
                    EndWaveNumber = WaveNumberRanges[^1];

                    doneFlag = true;
                }

                if (lineList[0] == "wnumstart")
                {
                    StartWaveNumber = double.Parse(lineList[1]);

                    doneFlag = true;
                }

                if (lineList[0] == "wnumend")
                {
                    EndWaveNumber = double.Parse(lineList[1]);

                    doneFlag = true;
                }

                if (lineList[0] == "wavenumberresolution")
                {
                    WaveNumberResolution = double.Parse(lineList[1]);

                    if (WaveNumberResolution < 0.0001) WaveNumberResolution = 0.0001;
                    if (WaveNumberResolution > 1) WaveNumberResolution = 1;

                    doneFlag = true;
                }

                if (lineList[0] == "wavenumberspread")
                {
                    WaveNumberSpread = double.Parse(lineList[1]);
                    if (WaveNumberSpread <= 0) WaveNumberSpread = -1;

                    doneFlag = true;
                }

                if (lineList[0] == "wavenumberoffset")
                {
                    StartWNumOffset = double.Parse(lineList[1]);

                    doneFlag = true;
                }

                if (lineList[0] == "intensitycutoff")
                {
                    IntensityCutoff = double.Parse(lineList[1]);

                    doneFlag = true;
                }

                if (lineList[0] == "runtype")
                {
                    RunType = RunType.DCTR;

                    bool b = Enum.TryParse(lineList[1], out RunType result);
                    if (b) RunType = result;

                    doneFlag = true;
                }

                if (lineList[0] == "profile")
                {
                    if (lineList[1].Equals("lorentz", StringComparison.CurrentCultureIgnoreCase))
                        Profile = Profile.Lorentz;
                    else Profile = Profile.Voigt;

                    doneFlag = true;
                }

                if (lineList[0] == "nwidths")
                {
                    Solver.Settings.NWidthLevels = int.Parse(lineList[1]);
                }

                if (lineList[0] == "nbinsplit")
                {
                    Solver.Settings.NShifts = int.Parse(lineList[1]);
                }

                if (lineList[0] == "species")
                {
                    Log.Information($"Species: \"{lineList[1]}\"");
                    lineList[1] = lineList[1].Replace(",", "\t");
                    lineList[1] = lineList[1].Replace(" ", "\t");
                    while (lineList[1].Contains("\t\t")) lineList[1] = lineList[1].Replace("\t\t", "\t");

                    lineList = lineList[1].ToUpper().Split('\t');

                    for (int i = 0; i < lineList.Length; i += 2)
                    {
                        SpeciesDefinitions[lineList[i + 1].Trim()] =
                            lineList[i].Trim();
                    }

                    doneFlag = true;
                }

                if (lineList[0] == "path file" || lineList[0] == "pathfile")
                {
                    if (lineList.Length > 1)
                    {
                        string fName = lineList[1].Trim();

                        fName = $"{BaseFolder}\\PathSpecifications\\{fName}";

                        Log.Information($"Reading path from file: {fName}");

                        if (fName.Length > 0)
                        {
                            StateMachine.ReadFromFile(fName);// string.Join(" ", lineList.Where((item, ix) => ix > 0)));
                        }
                    }

                    doneFlag = true;
                }

                if (lineList[0].StartsWith("slab"))// && !lineList[0].StartsWith("statefile"))
                {
                    stateLines.Add(line);

                    while (line != "}")
                    {
                        lineIx++;

                        line = inpList[lineIx];

                        stateLines.Add(line);
                    }

                    doneFlag = true;
                }

                if (!doneFlag && line != "") leftoverLines.Add(line);

                lineIx++;
            }

            if (StartWaveNumber < 0) StartWaveNumber = 0;
            if (EndWaveNumber < StartWaveNumber + 1) EndWaveNumber = StartWaveNumber + 1;

            Solver.Initialize(leftoverLines);

            if (stateLines.Count > 0) StateMachine.Initialize(stateLines);
        }

        public static void InitializeDataBuffer(DataBuffer dBuffer, string species, double startWNum, double endWNum)
        {
            dBuffer.Initialize(_folderDefinitions[species], _filePatternDefinitions[species], startWNum, endWNum);
        }

        public void ReadFromFile(string fileName)
        {
            List<string> initList = [];

            using (StreamReader sr = new(fileName))
            {
                initList = [.. sr.ReadToEnd().Split([Environment.NewLine],
                    StringSplitOptions.None)];
            }

            Initialize(initList);
        }

        static void report(string fileName, Dictionary<string, string> reportDict)
        {
            List<string> descList =
            [
                "Line count is evaluated over the full optical path (all homogeneous slabs).",
                "Thus, each line could be evaluated multiple times, and the line count increases accordingly.",
                "Lines are filtered by the intensity cutoff, after adjusting the intensity for",
                "temperature for each specific homogeneous slab.",
                "Therefore, the average line count per slab need not be an integer value.",
            ];

            string descString = string.Join(Environment.NewLine, descList);

            TimingFunctions.ReportTime(fileName);

            using StreamWriter sw = new(fileName, true);
            sw.WriteLine();

            sw.WriteLine(descString);

            sw.WriteLine();

            foreach (string key in reportDict.Keys)
            {
                sw.WriteLine($"{key}\t{reportDict[key]}");
            }
        }

        public void Run()
        {
            if (!IsValid)
            {
                Log.Warning($"Invalid run engine: could not run");

                return;
            }

            if (!Directory.Exists(OutputFolder)) Directory.CreateDirectory(OutputFolder);

            // Gets the minimum possible and maximum possible Lorentz half-width that
            // would be encountered in the full inhomogeneous gas path
            // This is determined based on pre-tabulated values, which are tabulated
            // based on temperature and species mole fraction
            Result<double[]?> minMaxWidthsRes = StateMachine.GetMinMaxWidths(
                SpeciesDefinitions, StartWaveNumber, EndWaveNumber);

            if (!minMaxWidthsRes.TryGetValue(out double[]? minMaxWidths))
            {
                Log.Warning($"RunEngine: Could not run: null min/max width output");

                return;
            }

            Solver.MinWidth = minMaxWidths[0];
            Solver.MaxWidth = minMaxWidths[1];

            Solver.Initialize(BaseFolder);

            List<string> inputFileList = [];

            TimingFunctions.Clear();

            Dictionary<string, DataBuffer> bufferDict = [];

            if (StateMachine.StateDefinitions is null)
            {
                Log.Warning($"Run engine: state definitions is null: could not run");

                return;
            }

            foreach (string speciesKey in SpeciesDefinitions.Keys)
            {
                if (!StateMachine.StateDefinitions.ContainsKey(speciesKey)) continue;

                string speciesDef = SpeciesDefinitions[speciesKey] + " " + speciesKey;

                string InputFolder = _folderDefinitions[speciesDef];
                string InputFilePattern = _filePatternDefinitions[speciesDef];

                Log.Information($"Initializing data buffer: {speciesDef}");

                DataBuffer dataBufferNew = new()
                {
                    DebugMode = false,
                    DebugFileName = $@"{BaseFolder}\DebugFolder\{speciesKey}.dbg",
                    SpeciesDatabaseDef = speciesDef
                };

                dataBufferNew.Initialize(InputFolder, InputFilePattern, StartWaveNumber, EndWaveNumber);

                // Set up one data buffer per species
                bufferDict.Add(speciesKey, dataBufferNew);
            }

            Solver.SetWaveNumberRange();

            if (WaveNumberRange is null)
            {
                Log.Warning($"Run engine: wavenumber range is null: could not run");

                return;
            }

            // Large spectral intervals are split into sub-intervals
            // (depending on the desired resolution) to facilitate
            // calculations
            bool requiresOverlap = WaveNumberRange.Length > 2;

            if (!requiresOverlap)
            {
                Log.Information($"Running wavenumber range: {StartWaveNumber} - {EndWaveNumber} cm-1");

                Solver.Run(StartWaveNumber, EndWaveNumber, requiresOverlap, StateMachine, bufferDict);

                Solver.SaveCurrentArray(requiresOverlap);
            }
            else
            {
                for (int i = 0; i < WaveNumberRange.Length - 1; i++)
                {
                    Log.Information($"Running wavenumber range:" +
                        $" {WaveNumberRange[i]} - {WaveNumberRange[i + 1]} cm-1");

                    CurWNumRangeIndex = i;

                    Solver.Run(WaveNumberRange[i], WaveNumberRange[i + 1], true, StateMachine, bufferDict);

                    Solver.SaveCurrentArray(true);
                }
            }

            Solver.FinalizeCalcs(requiresOverlap);

            Dictionary<string, string> reportDict = generateReportDictionary();

            Solver.UpdateReportDictionary(reportDict);

            writeFileList(OutputFolder + "\\InputFiles.dat", inputFileList);

            report(OutputFolder + "\\TimeReport.dat", reportDict);

            Log.Information($"Writing output data");

            Solver.WriteDetails(OutputFolder);

            Solver.WriteSettings(OutputFolder);

            WriteToFile(OutputFolder + "\\RunDefinition.dat");

            if (StateMachine.FileName is null)
            {
                Log.Warning($"Run engine: state machine has no file name" +
                    $" specification: could not complete run");

                return;
            }

            File.Copy(StateMachine.FileName, OutputFolder + "\\PathFile.dat");

            StateMachine.WriteStateArray(OutputFolder + "\\PathArray.dat");

            StateMachine.WriteStateSignature(OutputFolder + "\\PathSignature.dat");

            Log.Information($"Run complete: Success");
        }

        static void writeFileList(string fileName, List<string> inputFileList)
        {
            using StreamWriter sw = new(fileName);
            sw.WriteLine($"File list{Environment.NewLine}{string.Join(Environment.NewLine + "\t",
                inputFileList)}");
        }

        public void WriteToFile(string fName)
        {
            update();

            if (!IsValid)
            {
                Log.Warning($"Invalid run engine: could not write to file");

                return;
            }

            using StreamWriter sw = new(fName);
            sw.WriteLine($"Run type:\t{RunType}");
            sw.WriteLine($"Profile:\t{Profile}");
            foreach (string key in _speciesDefinitions.Keys)
                sw.WriteLine($"{_speciesDefinitions[key]} {key}");
            sw.Write("Path definition:\t");
            _stateMachine.WriteToStream(sw);
            sw.WriteLine($"StartWNum (cm-1):\t{StartWaveNumber}");
            sw.WriteLine($"EndWNum (cm-1):\t{EndWaveNumber}");
            sw.WriteLine($"WNum Resolution (cm-1):\t{WaveNumberResolution}");
            sw.WriteLine($"Wavenumber Spread  (cm-1):\t{WaveNumberSpread}");
            sw.WriteLine($"Wavenumber Offset (cm-1):\t{StartWNumOffset}");
            sw.WriteLine($"Intensity Cutoff (cm-1):\t{IntensityCutoff}");
        }

        #endregion Public Methods

        #region Public Properties

        public string BaseFolder
        {
            get { return _baseFolder; }
            set { _baseFolder = value; }
        }

        public int CurWNumRangeIndex
        {
            get { return _curWNumRangeIndex; }
            set { _curWNumRangeIndex = value; }
        }

        public string Description
        {
            get { return _description; }
            set { _description = value; }
        }

        // Cm-1
        public double EndWaveNumber
        {
            get { return _endWaveNumber; }
            set { _endWaveNumber = value; }
        }

        public double IntensityCutoff
        {
            get { return _intensityCutoff; }
            set { _intensityCutoff = value; }
        }

        [MemberNotNullWhen(true, [nameof(_speciesDefinitions),
            nameof(_stateMachine), nameof(_outputFolder),
            nameof(SpeciesDefinitions), nameof(StateMachine), nameof(OutputFolder)])]
        public bool IsValid
        {
            get
            {
                update();

                return _isValid;
            }
        }

        public string? OutputFolder
        {
            get { return _outputFolder; }
            set { _outputFolder = value; }
        }

        public Profile Profile
        {
            get { return _profile; }
            set { _profile = value; }
        }

        public RunType RunType
        {
            get { return _runType; }
            set { _runType = value; }
        }

        public RapidSolver Solver
        {
            get { return _solver; }
            set { _solver = value; }
        }

        // Key - species (CO2, H2O, etc.)
        // Value - source (HITEMP, HITRAN, CDSD)
        public Dictionary<string, string>? SpeciesDefinitions
        {
            get { return _speciesDefinitions; }
            set { _speciesDefinitions = value; }
        }

        // Cm-1
        public double StartWaveNumber
        {
            get { return _startWaveNumber; }
            set { _startWaveNumber = value; }
        }

        public double StartWNumOffset
        {
            get { return _startWNumOffset; }
            set { _startWNumOffset = value; }
        }

        public StateMachine? StateMachine
        {
            get { return _stateMachine; }
            set { _stateMachine = value; }
        }

        public double[]? WaveNumberRange
        {
            get { return _waveNumberRange; }
            set { _waveNumberRange = value; }
        }

        // Cm-1
        public double WaveNumberResolution
        {
            get { return _waveNumberResolution; }
            set { _waveNumberResolution = value; }
        }

        // Cm-1 - specifies Lorentz line spread
        public double WaveNumberSpread
        {
            get { return _waveNumberSpread; }
            set { _waveNumberSpread = value; }
        }

        #endregion Public Properties

        void update()
        {
            if (_isUpdated) return;

            _isUpdated = true;

            _isValid = _outputFolder is not null;
            _isValid = _isValid && _stateMachine is not null;
            _isValid = _isValid && _speciesDefinitions is not null;
        }

        #region Private Methods

        static RunEngine()
        {
            _filePatternDefinitions.Add("HITRAN CO2", "02_*");
            _filePatternDefinitions.Add("HITRAN H2O", "01_*");

            _filePatternDefinitions.Add("HITRAN2020 CO2", "02_*");
            _filePatternDefinitions.Add("HITRAN2020 H2O", "01_*");

            _filePatternDefinitions.Add("HITEMP2010 CO2", "02_*");
            _filePatternDefinitions.Add("HITEMP2010 H2O", "01_*");

            _filePatternDefinitions.Add("CDSD CO2", "cdsd_*");

            _filePatternDefinitions.Add("TEST CO2", "02_*");
            _filePatternDefinitions.Add("TEST H2O", "01_*");
        }

        #endregion Private Methods

        #region Private Properties

        string _baseFolder = "";

        int _curWNumRangeIndex = -1;

        string _description = "";

        //readonly DopplerEngine _dopplerEngine = new();

        // Cm-1
        double _endWaveNumber = 0;

        static readonly Dictionary<string, string> _filePatternDefinitions = [];

        static readonly Dictionary<string, string> _folderDefinitions = [];

        double _intensityCutoff = 0;

        bool _isUpdated = false;

        bool _isValid = false;

        string? _outputFolder = null;

        Profile _profile = Profile.Lorentz;

        RunType _runType = RunType.DCTR;

        RapidSolver _solver = new();

        // Key - species (CO2, H2O, etc.)
        // Value - source (HITEMP, HITRAN, CDSD)
        Dictionary<string, string>? _speciesDefinitions = null;

        double _startWNumOffset = 0;

        // Cm-1
        double _startWaveNumber = 0;

        StateMachine? _stateMachine = null;

        double[]? _waveNumberRange = null;

        // Cm-1
        double _waveNumberResolution = 0.01;

        // Cm-1 - specifies Lorentz line spread
        double _waveNumberSpread = -1;

        #endregion Private Properties
    }
}
