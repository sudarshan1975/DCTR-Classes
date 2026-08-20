using FunctionLibrary;

namespace DCTRClasses
{
    public class PartitionFunction
    {
        // The partition function is represented as a series of cubic splines
        // The string key specifies the species (eg.: "HITRAN CO2" or "HITEMP H2O")
        // Within the inner dictionary, the integer key specifies the isotope ID
        // while the spline list represents a list of splines for each isotope
        // As an example, the spline list might contain 10 spline representations,
        // each of which is valid within a limited temperature range
        static Dictionary<string, Dictionary<int, List<CubicSpline>>> _splineDict = [];

        // See description for _splineDict
        // The string key specifies the species (eg.: "HITRAN CO2" or "HITEMP H2O")
        // Within the inner dictionary, the integer key specifies the isotope ID
        // while the spline list represents a list of temperature ranges for each isotope
        // For example, the full temperature range of (say) 0 to 6000 K may be divided
        // into smaller ranges of 50 K each, with a representative spline for each range
        static Dictionary<string, Dictionary<int, List<double>>> _temperatureDict = [];

        // Specifies the base temperature for each species (eg.: "HITRAN CO2" or "HITEMP H2O")
        // For example - HITRAN and HITEMP have a base temperature of 296 K,
        // while CDSD has 300 K
        static Dictionary<string, double> _baseTemperatureDict = [];

        // The string key specifies the species (eg.: "HITRAN CO2" or "HITEMP H2O")
        // Within the inner dictionary, the integer key specifies the isotope ID
        // The base value is a normalization value (the value of the partition function
        // for that species/ isotope at its base temperature)
        static Dictionary<string, Dictionary<int, double>> _baseValueDict = [];

        static PartitionFunction()
        {

        }

        public static void Initialize(string inpFile)
        {
            ReadFromBinaryFile(inpFile);
        }

        public static void WriteToBinaryFile(string fileName)
        {
            using BinaryWriter bw = new(File.Open(fileName, FileMode.Create));
            bw.Write(_splineDict.Count);

            foreach (string key in _splineDict.Keys)
            {
                bw.Write(key);

                bw.Write(_splineDict[key].Count);

                foreach (int isotopeID in _splineDict[key].Keys)
                {
                    List<CubicSpline> splineList = _splineDict[key][isotopeID];

                    bw.Write(isotopeID);

                    bw.Write(splineList.Count);

                    foreach (CubicSpline spline in splineList)
                    {
                        spline.WriteToBinaryStream(bw);
                    }
                }
            }

            bw.Write(_temperatureDict.Count);

            foreach (string key in _temperatureDict.Keys)
            {
                bw.Write(key);

                bw.Write(_temperatureDict[key].Count);

                foreach (int isotopeID in _temperatureDict[key].Keys)
                {
                    List<double> TList = _temperatureDict[key][isotopeID];

                    bw.Write(isotopeID);

                    bw.Write(TList.Count);

                    foreach (double T in TList)
                    {
                        bw.Write(T);
                    }
                }
            }

            bw.Write(_baseTemperatureDict.Count);

            foreach (string key in _baseTemperatureDict.Keys)
            {
                bw.Write(key);

                bw.Write(_baseTemperatureDict[key]);
            }

            bw.Write(_baseValueDict.Count);

            foreach (string key in _baseValueDict.Keys)
            {
                bw.Write(key);

                bw.Write(_baseValueDict[key].Count);

                foreach (int isotopeID in _baseValueDict[key].Keys)
                {
                    bw.Write(isotopeID);

                    bw.Write(_baseValueDict[key][isotopeID]);
                }
            }
        }

        public static void ReadFromBinaryFile(string fileName)
        {
            _splineDict = [];
            _temperatureDict = [];
            _baseTemperatureDict = [];
            _baseValueDict = [];

            using BinaryReader br = new(File.Open(fileName, FileMode.Open));
            int keyCount = br.ReadInt32();

            for (int i = 0; i < keyCount; i++)
            {
                string key = br.ReadString();

                _splineDict.Add(key, []);

                int isotopeCount = br.ReadInt32();

                for (int j = 0; j < isotopeCount; j++)
                {
                    int isotopeID = br.ReadInt32();

                    _splineDict[key].Add(isotopeID, []);

                    List<CubicSpline> splineList = _splineDict[key][isotopeID];

                    int splineCount = br.ReadInt32();

                    for (int splineIndex = 0; splineIndex < splineCount; splineIndex++)
                    {
                        CubicSpline spline = new();

                        splineList.Add(spline);

                        spline.ReadFromBinaryStream(br);
                    }
                }
            }

            keyCount = br.ReadInt32();

            for (int i = 0; i < keyCount; i++)
            {
                string key = br.ReadString();

                _temperatureDict.Add(key, []);

                int isotopeCount = br.ReadInt32();

                for (int j = 0; j < isotopeCount; j++)
                {
                    int isotopeID = br.ReadInt32();

                    _temperatureDict[key].Add(isotopeID, []);

                    int TCount = br.ReadInt32();

                    for (int TIx = 0; TIx < TCount; TIx++)
                    {
                        _temperatureDict[key][isotopeID].Add(br.ReadDouble());
                    }
                }
            }

            keyCount = br.ReadInt32();

            for (int i = 0; i < keyCount; i++)
            {
                string key = br.ReadString();

                double T = br.ReadDouble();

                _baseTemperatureDict.Add(key, T);
            }

            keyCount = br.ReadInt32();

            for (int i = 0; i < keyCount; i++)
            {
                string key = br.ReadString();

                _baseValueDict.Add(key, []);

                int TCount = br.ReadInt32();

                for (int j = 0; j < TCount; j++)
                {
                    int isotopeID = br.ReadInt32();

                    double T = br.ReadDouble();

                    _baseValueDict[key].Add(isotopeID, T);
                }
            }
        }

        public static string GetDescriptionString(string desc)
        {
            if (_temperatureDict.ContainsKey(desc)) return desc;

            if (desc.EndsWith("H2O"))
            {
                return "HITRAN H2O";
            }

            return "HITRAN CO2";
        }

        static List<double> getTemperatureList(string desc, int isotopeID)
        {
            desc = GetDescriptionString(desc);

            if (!_temperatureDict[desc].TryGetValue(isotopeID, out List<double>? value))
            {
                return _temperatureDict[desc][1];
            }

            return value;
        }

        static List<CubicSpline> getSplineList(string desc, int isotopeID)
        {
            desc = GetDescriptionString(desc);

            if (!_splineDict[desc].TryGetValue(isotopeID, out List<CubicSpline>? value))
            {
                return _splineDict[desc][1];
            }

            return value;
        }

        static double getBaseValue(string desc, int isotopeID)
        {
            desc = GetDescriptionString(desc);

            if (!_baseValueDict[desc].TryGetValue(isotopeID, out double value))
            {
                return _baseValueDict[desc][1];
            }

            return value;
        }

        public static Dictionary<int, double> GetValues(double T, string desc)
        {
            Dictionary<int, double> outDict = [];

            foreach (int isotopeID in _temperatureDict[desc].Keys)
            {
                outDict.Add(isotopeID, GetValue(desc, isotopeID, T));
            }

            return outDict;
        }

        public static double GetRawValue(string desc, int isotopeID, double T)
        {
            List<double> TList = getTemperatureList(desc, isotopeID);
            List<CubicSpline> splineList = getSplineList(desc, isotopeID);

            int ix = TList.BinarySearch(T);

            if (ix < 0) ix = ~ix;
            ix--;

            if (ix < 0) return double.NaN;
            if (ix >= TList.Count) return double.NaN;

            Result<double?> result = splineList[ix].Query(T);

            if (result.TryGetValue(out double? value))
            {
                return value.Value;
            }

            return double.NaN;
        }

        /// <summary>
        /// This is the workhorse function, which yields the value of the
        /// partition function for the specified species, isotope, and temperature
        /// </summary>
        /// <param name="desc">Species description (for example: "HITRAN CO2")</param>
        /// <param name="isotopeID">Isotope ID</param>
        /// <param name="T">Temperature in Kelvin</param>
        /// <returns></returns>
        public static double GetValue(string desc, int isotopeID, double T)
        {
            List<double> TList = getTemperatureList(desc, isotopeID);
            List<CubicSpline> splineList = getSplineList(desc, isotopeID);

            int ix = TList.BinarySearch(T);

            if (ix < 0) ix = ~ix;
            ix--;

            if (ix < 0) return double.NaN;
            if (ix >= TList.Count) return double.NaN;

            Result<double?> result = splineList[ix].Query(T);

            Logging.LogResult(result);

            if (result.TryGetValue(out double? value))
            {
                return value.Value / getBaseValue(desc, isotopeID);
            }

            return double.NaN;
        }
    }
}
