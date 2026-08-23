using Serilog;
using System.Diagnostics.CodeAnalysis;

namespace DCTRClasses
{
    public class LineEnsemble
    {
        #region Public Methods

        public void AppendToFile(string fileName, int startIx, int endIx, double intensityCutoff)
        {
            if (!IsValid)
            {
                Log.Warning($"Could not append to file in line ensemble");

                return;
            }

            int curIx = startIx * NDataColumns;

            using StreamWriter sw = new(fileName, true);
            for (int i = startIx; i < endIx; i++)
            {
                if (curIx >= 0 && DataArray[curIx + 1] >= intensityCutoff)
                {
                    for (int j = 0; j < NDataColumns; j++)
                    {
                        sw.Write(DataArray[curIx + j]);
                        if (j < NDataColumns - 1) sw.Write("\t");
                    }

                    sw.WriteLine();
                }

                curIx += NDataColumns;
            }
        }

        /// <summary>
        /// Gets the composite half-width (self and air-broadened, weighted by mole fraction)
        /// </summary>
        /// <param name="ix">Row index (specifies the line which is being considered)</param>
        /// <param name="moleFraction"></param>
        /// <param name="T">In Deg. K</param>
        /// <returns></returns>
        public double GetHalfWidth(int ix, double moleFraction, double T)
        {
            if (!IsValid)
            {
                Log.Warning($"Invalid line ensemble: could not get half width");

                return double.NaN;
            }

            int ixSingle = ix * NDataColumns;

            double value = 0, v, f;
            T = TRef / T;

            v = (1 - moleFraction) * DataArray[ixSingle + 4];
            f = Math.Pow(T, DataArray[ixSingle + 6]);
            value += v * f;

            v = moleFraction * DataArray[ixSingle + 5];
            f = Math.Pow(T, DataArray[ixSingle + 7]);
            value += v * f;

            return value;
        }

        public static int GetIndex(int lineIndex, string parameter)
        {
            return lineIndex * NDataColumns + DataColumnDict[parameter];
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="N">Total number of lines</param>
        public void Initialize(int N)
        {
            // 0: Center
            // 1: Intensity
            // 2: Lower energy
            // 3: Pressure shift
            // 4: wAir
            // 5: wSelf
            // 6: nAir
            // 7: nSelf
            // 8: Current center
            // 9: Current intensity
            // 10: Current half-width
            // 11: Raw scaled intensity
            DataArray = new double[N * NDataColumns];

            // 0: Wavenumber index (calculated)
            // 1: Width index (calculated)
            // 2: Isotope ID (from source)
            IndexArray = new int[N * NIndexColumns];

            LineCount = N;
        }

        public void ReadFromBinaryFile(string fileName)
        {
            if (fileName.Contains("CDSD")) TRef = 300;
            else TRef = 296;

            LineCount = (int)(new FileInfo(fileName).Length / 68);

            using BinaryReader br = new(File.Open(fileName, FileMode.Open));
            Initialize(LineCount);

            if (!IsValid)
            {
                Log.Warning($"Invalid line ensemble: could not read from binary file");

                return;
            }

            int lCount = LineCount * NDataColumns;

            int ix = 0;
            for (int i = 0; i < lCount; i += NDataColumns)
            {
                DataArray[i] = br.ReadDouble();
                DataArray[i + 1] = br.ReadDouble();
                DataArray[i + 2] = br.ReadDouble();
                DataArray[i + 3] = br.ReadDouble();
                DataArray[i + 4] = br.ReadDouble();
                DataArray[i + 5] = br.ReadDouble();
                DataArray[i + 6] = br.ReadDouble();
                DataArray[i + 7] = br.ReadDouble();

                IndexArray[ix + 2] = br.ReadInt32();

                ix += NIndexColumns;
            }
        }

        public void ReadFromFile(string fileName)
        {
            if (!IsValid)
            {
                Log.Warning($"Invalid line ensemble: could not read from file");

                return;
            }

            // HITEMP
            List<int> countList = [2, 1, 12, 10, 10, 5];
            countList.AddRange([5, 10, 4, 8, 4]);

            List<double[]> arrayList = [];
            List<int> indexArrayList = [];

            string? l;
            using (StreamReader sr = new(fileName))
            {
                l = sr.ReadLine();

                while (l != null)
                {
                    List<string> splitList = StringFunctions.SplitString(l, countList);

                    double lineCenter = Convert.ToDouble(splitList[2]);
                    double lineIntensity = Convert.ToDouble(splitList[3]);

                    double wAir = Convert.ToDouble(splitList[5]);
                    double wSelf = Convert.ToDouble(splitList[6]);
                    double lowerEnergy = Convert.ToDouble(splitList[7]);
                    double nAir = Convert.ToDouble(splitList[8]);
                    double pressureShift = Convert.ToDouble(splitList[9]);

                    double nSelf;

                    int isotopeID;

                    if (splitList[1] == "0") isotopeID = 10;
                    else if (splitList[1] == "A") isotopeID = 11;
                    else if (splitList[1] == "B") isotopeID = 12;
                    else isotopeID = Convert.ToInt32(splitList[1]);

                    if (splitList.Count > 10)
                    {
                        bool b = double.TryParse(splitList[10], out double v);
                        nSelf = b ? v : nAir;
                    }
                    else
                    {
                        nSelf = nAir;
                    }

                    arrayList.Add([lineCenter, lineIntensity, lowerEnergy, pressureShift, wAir, wSelf, nAir, nSelf]);
                    indexArrayList.Add(isotopeID);

                    l = sr.ReadLine();
                }
            }

            int N = arrayList.Count;

            Initialize(N);

            Log.Information($"Read {N} lines");

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    DataArray[i * NDataColumns + j] = arrayList[i][j];
                }

                IndexArray[i * NIndexColumns + 2] = indexArrayList[i];
            }
        }

        /// <summary>
        /// Sets parameters such as scaled intensity, composite half-width, etc.
        /// </summary>
        /// <param name="ix">Row index (represents individual lines)</param>
        /// <param name="state">Contains specifications of each homogeneous slab within the inhomogeneous path</param>
        /// <param name="species">Gas species</param>
        /// <param name="startWNum"></param>
        /// <param name="widthRanges"></param>
        /// <param name="widthRatio"></param>
        /// <param name="wNumRes"></param>
        /// <param name="intensityCutoff">Line intensity threshold</param>
        public void SetCurrentState(int ix, RunEngineState state, string species, double startWNum
    , double[] widthRanges, double widthRatio,
    double wNumRes, double intensityCutoff)
        {
            if (!IsValid)
            {
                Log.Warning($"Invalid line ensemble: could not set current state");

                return;
            }

            int ixSingle = ix * NDataColumns, ixIndex = ix * NIndexColumns;

            int isotopeID = IndexArray[ixIndex + 2];

            double LEC2 = DataArray[ixSingle + 2] * c2;
            double LCC2 = -DataArray[ixSingle] * c2;
            double OnesComplementCNuTRefInv = DataArray[ixSingle + 1] / (1 - Math.Exp(LCC2 / TRef));

            double v = (1 - Math.Exp(LCC2 / state.Temperature)) * OnesComplementCNuTRefInv;
            v *= Math.Exp(LEC2 * (1 / TRef - 1 / state.Temperature));

            int pfIsotopeID = isotopeID;
            if (!state.PartitionFunctionValueDictionary[species].ContainsKey(pfIsotopeID)) { pfIsotopeID = 1; }

            v /= state.PartitionFunctionValueDictionary[species][pfIsotopeID];

            // Raw scaled intensity
            DataArray[ixSingle + 11] = v;

            if (v < intensityCutoff) return;

            // Scaled intensity
            DataArray[ixSingle + 9] = v * state.IntensityScalingFactor;

            // Pressure shift of line center
            double Delta0 = DataArray[ixSingle + 3] * state.TotalPressure * (1 - state.MoleFraction);

            // Pressure-shifted line center
            DataArray[ixSingle + 8] = DataArray[ixSingle] + Delta0;

            // Line half-width
            DataArray[ixSingle + 10] = state.TotalPressure * GetHalfWidth(ix, state.MoleFraction, state.Temperature);

            // This calculates the index of the wavenumber location in the output array
            IndexArray[ixIndex] = (int)((DataArray[ixSingle + 8] -
                startWNum) / wNumRes);

            // This calculates the representative Lorentz half-width bin index
            if (widthRanges == null)
            {
                IndexArray[ixIndex + 1] = -1;
            }
            else
            {
                int w = (int)(Math.Log(DataArray[ixSingle + 10] / widthRanges[0])
                    / Math.Log(widthRatio));

                if (w >= widthRanges.Length)
                {
                    w = widthRanges.Length - 1;
                }

                IndexArray[ixIndex + 1] = w;

                if (IndexArray[ixIndex + 1] < 0) IndexArray[ixIndex + 1] = 0;
            }
        }

        public void WriteBinaryFile(string fileName)
        {
            if (!IsValid)
            {
                Log.Warning($"Invalid line ensemble: could not write binary file");

                return;
            }

            using BinaryWriter bw = new(File.Open(fileName, FileMode.Create));
            for (int i = 0; i < LineCount; i++)
            {
                int ix = i * NDataColumns;

                for (int j = 0; j < 8; j++) bw.Write(DataArray[ix + j]);
                bw.Write(IndexArray[i * NIndexColumns + 2]);
            }
        }

        #endregion Public Methods

        #region Public Properties

        // Dictionary which specifies the column index of each line parameter
        // Key represents the line parameter; value represents the column index
        readonly static Dictionary<string, int> DataColumnDict = [];

        [MemberNotNullWhen(true, [nameof(DataArray), nameof(IndexArray)])]
        public bool IsValid
        {
            get
            {
                update();

                return _isValid;
            }
        }

        public int LineCount = 0;

        // Static array which contains line data; rows represent individual lines,
        // columns represent line parameters; the 2-d array is flattened into a 1-d array
        // for fast access
        public double[]? DataArray = null;

        // Static array which contains line indices (see description of "NIndexColumns" above)
        // Rows represent individual lines, columns represent line indices; the 2-d array is flattened
        // into a 1-d array for fast access
        public int[]? IndexArray = null;

        // NDataColumns is the total number of columns which represents each line
        // Each column contains one particular line parameter, such as intensity, line center, etc.
        // NIndexColumns contains useful calculated indices, such as wavenumber location,
        // or source indices, such as isotope ID
        public static readonly int NDataColumns = 12, NIndexColumns = 3;

        // Reference temperature (K)
        public double TRef = 296;

        #endregion Public Properties

        #region Private Methods

        static LineEnsemble()
        {
            DataColumnDict.Add("LineCenter", 0);
            DataColumnDict.Add("LineIntensity", 1);
            DataColumnDict.Add("LowerEnergy", 2);
            DataColumnDict.Add("PressureShift", 3);
            DataColumnDict.Add("AirWidth", 4);
            DataColumnDict.Add("SelfWidth", 5);
            DataColumnDict.Add("AirIndex", 6);
            DataColumnDict.Add("SelfIndex", 7);
            DataColumnDict.Add("ShiftedCenter", 8);
            DataColumnDict.Add("ScaledIntensity", 9);
            DataColumnDict.Add("ScaledWidth", 10);
            DataColumnDict.Add("RawScaledIntensity", 11);
        }

        void update()
        {
            if (_isUpdated) return;

            _isUpdated = true;

            _isValid = DataArray is not null && IndexArray is not null;
        }

        #endregion Private Methods

        #region Private Properties

        // Radiation constant
        const double c2 = 1.4388028496642257;

        // Avogadro number (molecules/mol)
        //const double NA = 6.02214076e23;

        // Universal gas constant (J/mol-K)
        //static readonly double R = 8.314;

        bool _isUpdated = false;

        bool _isValid = false;

        #endregion Private Properties
    }
}
