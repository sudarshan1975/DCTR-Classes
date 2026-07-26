using Serilog;

namespace DCTRClasses
{
    /// <summary>
    /// This class allows loading of line data in a buffered fashion
    /// Line data are sometimes spread across many files (for ex. - the CDSD
    /// database contains more than 2,500 files, with data for around 700 million lines)
    /// The DataBuffer class keeps track of files within each spectral sub-interval and
    /// loads files upon demand
    /// </summary>
    public class DataBuffer
    {
        #region Public Methods

        // Convenience method to get the names of files represented by the buffer
        public List<string> GetFileList()
        {
            return [.. _fileList.Select(item => item.FileName)];
        }

        /// <summary>
        /// Gives the starting and ending indices of files corresponding
        /// to the specified start and end wavenumbers
        /// </summary>
        /// <param name="indexList"></param>
        /// <param name="startWNum"></param>
        /// <param name="endWNum"></param>
        public void GetNextIndexList(int[] indexList, double startWNum, double endWNum)
        {
            if (LineEnsemble == null || indexList[1] >= LineEnsemble.LineCount - 1)
            {
                UpdateData();
            }
            else
            {
                if (LineEnsemble.DataArray[(indexList[1] + 1) * LineEnsemble.NDataColumns] >= endWNum)
                {
                    indexList[0] = -1;
                    indexList[1] = -1;

                    return;
                }
            }

            indexList[0] = -1;
            indexList[1] = -1;

            if (LineEnsemble == null) return;

            indexList[0] = binarySearch(startWNum);
            indexList[1] = binarySearch(endWNum);

            //if (!DebugMode) return;

            //if (indexList[0] < 0 || indexList[1] < 0) return;
        }

        // Convenience method, currently deprecated
        public static Tuple<double[,], double[,], double[,], double[,]> GetWidthTables(string species)
        {
            List<double> TempList = [];
            for (double T = 300; T <= 6000; T += 300)
            {
                TempList.Add(T);
            }

            TempList.Insert(0, 100);
            TempList.Insert(1, 200);

            List<double> wNumList = [];
            for (double wNum = 0; wNum < 30001; wNum += 100)
            {
                wNumList.Add(wNum);
            }

            double[,] minSelfWidthArray = new double[TempList.Count, wNumList.Count - 1];
            double[,] maxSelfWidthArray = new double[TempList.Count, wNumList.Count - 1];
            double[,] minAirWidthArray = new double[TempList.Count, wNumList.Count - 1];
            double[,] maxAirWidthArray = new double[TempList.Count, wNumList.Count - 1];

            for (int i = 0; i < TempList.Count; i++)
            {
                for (int j = 0; j < wNumList.Count - 1; j++)
                {
                    minSelfWidthArray[i, j] = double.PositiveInfinity;

                    maxSelfWidthArray[i, j] = double.NegativeInfinity;

                    minAirWidthArray[i, j] = double.PositiveInfinity;

                    maxAirWidthArray[i, j] = double.NegativeInfinity;
                }
            }

            DataBuffer dBuffer = new();

            RunEngine.InitializeDataBuffer(dBuffer, species, 0, 30000);

            dBuffer.UpdateData();

            int wNumIx = 0;

            while (dBuffer.LineEnsemble != null)
            {
                for (int i = 0; i < dBuffer.LineEnsemble.LineCount; i++)
                {
                    int lineCenterIx = i * LineEnsemble.NDataColumns;

                    int TIx = 0;
                    foreach (double T in TempList)
                    {
                        double lineCenter = dBuffer.LineEnsemble.DataArray[lineCenterIx];

                        while (lineCenter > wNumList[wNumIx + 1])
                        {
                            wNumIx++;
                        }

                        double w = dBuffer.LineEnsemble.GetHalfWidth(i, 0, T);

                        if (minAirWidthArray[TIx, wNumIx] > w) minAirWidthArray[TIx, wNumIx] = w;
                        if (maxAirWidthArray[TIx, wNumIx] < w) maxAirWidthArray[TIx, wNumIx] = w;

                        w = dBuffer.LineEnsemble.GetHalfWidth(i, 1, T);

                        if (minSelfWidthArray[TIx, wNumIx] > w) minSelfWidthArray[TIx, wNumIx] = w;
                        if (maxSelfWidthArray[TIx, wNumIx] < w) maxSelfWidthArray[TIx, wNumIx] = w;

                        TIx++;
                    }
                }

                dBuffer.UpdateData();
            }

            return new Tuple<double[,], double[,], double[,], double[,]>(minSelfWidthArray, maxSelfWidthArray,
                minAirWidthArray, maxAirWidthArray);
        }

        public void Initialize(string folder, string pattern, double startWNum, double endWNum)
        {
            FolderName = folder;
            FilePattern = pattern;
            StartWavenumber = startWNum;
            EndWavenumber = endWNum;

            Initialize();
        }

        public void Initialize()
        {
            if (DebugMode)
            {
                if (File.Exists(DebugFileName)) File.Delete(DebugFileName);
            }

            _fileList = [];

            // Get list of files which match the specified pattern, from the specified folder
            string[] fileList = Directory.GetFiles(FolderName, FilePattern);

            foreach (string f in fileList)
            {
                FileData newFData = new();
                newFData.SetFileName(f);

                _fileList.Add(newFData);
            }

            // Truncate file list to only contain files relevant to the specified wavenumber range
            // The wavenumber range is augmented by the neighbor line cutoff
            _fileList = [.._fileList.Where(item => item.StartWavenumber < EndWavenumber + NeighborLinesWNumCutoff
            && item.EndWavenumber > StartWavenumber - NeighborLinesWNumCutoff)];

            // List of starting wavenumbers for each data file
            // For example: for CDSD, each file contains lines for a wavenumber range of
            // several hundred cm-1; binary files are named according to the start wavenumber
            // and end wavenumber of lines contained in the file
            List<double> sWNumList = [.. _fileList.Select(item => item.StartWavenumber)];

            // Reorder the file list by start wavenumbers in ascending order
            // This step is usually redundant, since files are anyway presented
            // in the order of start wavenumbers by the operating system
            _fileList = [.. sWNumList.Zip(_fileList, (val1, val2) => new { val1, val2 }).
                OrderBy(pair => pair.val1).Select(pair => pair.val2)];

            _fileIndex = 0;
        }

        public override string ToString()
        {
            return FolderName + " (" + StartWavenumber + " to " + EndWavenumber + " cm-1)";
        }

        /// <summary>
        /// This method is periodically called to pull data from the next file
        /// Data are stored in the buffer within the "LineEnsemble" property
        /// </summary>
        public void UpdateData()
        {
            if (_fileList is null || _fileIndex >= _fileList.Count)
            {
                LineEnsemble = null;

                return;
            }

            TimingFunctions.InitializeTime("RHTC: DataBuffer: GetData");

            Log.Information($"Reading file: {Path.GetFileName(_fileList[_fileIndex].FileName)}");

            LineEnsemble = new();

            LineEnsemble.ReadFromBinaryFile(_fileList[_fileIndex].FileName);

            _fileIndex++;

            TimingFunctions.AddTime("RHTC: DataBuffer: GetData");
        }

        #endregion Public Methods

        #region Public Properties

        // For debugging only
        public bool DebugMode
        {
            get { return _debugMode; }
            set { _debugMode = value; }
        }

        public string? DebugFileName
        {
            get { return _debugFileName; }
            set { _debugFileName = value; }
        }

        // Ending wavenumber of spectral interval (cm-1)
        public double EndWavenumber
        {
            get { return _endWavenumber; }
            set { _endWavenumber = value; }
        }

        // File name pattern for loading binary line data files
        public string FilePattern
        {
            get { return _filePattern; }
            set { _filePattern = value; }
        }

        // Base folder containing line data binary files
        public string? FolderName
        {
            get { return _folderName; }
            set { _folderName = value; }
        }

        // Holds line data in memory as flattened arrays
        public LineEnsemble? LineEnsemble
        {
            get { return _lineEnsemble; }
            set { _lineEnsemble = value; }
        }

        // For each spectral sub-interval, some lines outside the sub-interval
        // are also considered, since these lines have overlap with the sub-interval
        public double NeighborLinesWNumCutoff
        {
            get { return _neighborLinesWNumCutoff; }
            set { _neighborLinesWNumCutoff = value; }
        }

        // Ex.: CDSD CO2 or HITRAN H2O
        public string SpeciesDatabaseDef
        {
            get { return _speciesDatabaseDef; }
            set { _speciesDatabaseDef = value; }
        }

        // Starting wavenumber of spectral interval (cm-1)
        public double StartWavenumber
        {
            get { return _startWavenumber; }
            set { _startWavenumber = value; }
        }

        #endregion Public Properties

        #region Private Methods

        /// <summary>
        /// Implements a binary search to get the index of the file
        /// corresponding to the given wavenumber location
        /// </summary>
        /// <param name="wNum">Wavenumber location in cm-1</param>
        /// <returns></returns>
        int binarySearch(double wNum)
        {
            if (LineEnsemble is null || LineEnsemble.DataArray is null) return -1;

            double[] dataArray = LineEnsemble.DataArray;

            int N = dataArray.Length / LineEnsemble.NDataColumns;
            if (N == 0) return -1;

            if (dataArray[0] >= wNum)
            {
                return -1;
            }

            if (dataArray[(N - 1) * LineEnsemble.NDataColumns] <= wNum) return N - 1;

            int ix1 = 0, ix2 = N, ix = (ix1 + ix2) >> 1;

            while (true)
            {
                int ixSingle = ix * LineEnsemble.NDataColumns;

                if (dataArray[ixSingle] <= wNum && dataArray[(ix + 1) * LineEnsemble.NDataColumns] >= wNum) return ix;

                if (dataArray[ixSingle] > wNum)
                {
                    ix2 = ix;

                    ix = (ix2 + ix1) >> 1;
                }
                else if (dataArray[ixSingle] < wNum)
                {
                    ix1 = ix;

                    ix = (ix1 + ix2) >> 1;
                }
            }
        }

        #endregion Private Methods

        #region Private Properties

        bool _debugMode = false;

        string? _debugFileName = null;

        // Ending wavenumber of spectral interval (cm-1)
        double _endWavenumber = 10000;

        // Internal book-keeping
        int _fileIndex = 0;

        // Description of files to be considered for the run
        List<FileData>? _fileList = null;

        // File name pattern for loading binary line data files
        string _filePattern = "";

        // Base folder containing line data binary files
        string? _folderName = null;

        // Holds line data in memory as flattened arrays
        LineEnsemble? _lineEnsemble = null;

        // For each spectral sub-interval, some lines outside the sub-interval
        // are also considered, since these lines have overlap with the sub-interval
        double _neighborLinesWNumCutoff = 0;

        // Ex.: CDSD CO2 or HITRAN H2O
        string _speciesDatabaseDef = "";

        // Starting wavenumber of spectral interval (cm-1)
        double _startWavenumber = 200;

        #endregion Private Properties
    }
}
