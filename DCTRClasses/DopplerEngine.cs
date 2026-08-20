using Newtonsoft.Json.Linq;
using System.Diagnostics.CodeAnalysis;
using FLIB = FunctionLibrary;

namespace DCTRClasses
{
    public class DopplerEngine
    {
        #region Public Methods

        /// <summary>
        /// One-sided Doppler distribution
        /// </summary>
        /// <param name="TInKelvin"></param>
        /// <param name="wNum">Wavenumber in cm-1</param>
        /// <param name="species"></param>
        /// <param name="wNumRes">Wavenumber resolution of output array in cm-1</param>
        /// <returns></returns>
        public FLIB.Result<double[]?> GetDistributionArray(double TInKelvin,
            double wNum, string species, int isotopeID, double wNumRes)
        {
            // This relies on the fact that line databases present lines in increasing
            // order of line center. Thus, an initial Doppler profile for the first wavenumber
            // bin (see the comments on the "Initialize" method) may be calculated, and cached
            // This cached profile is used, until the center of a line is best represented by
            // the next wavenumber bin, at which point a new Doppler profile (corresponding to
            // the next bin) is calculated, and cached
            if (_curWNum > 0 && wNum < _curWNum * _WNumMult)
            {
                if (_cachedArray is null) return new(null, false, $"Invalid cached array", FLIB.Severity.WARNING);

                return new(_cachedArray);
            }

            // Update the cached representative bin wavenumber
            if (_curWNum <= 0) _curWNum = wNum;
            else while (_curWNum * _WNumMult < wNum) _curWNum *= _WNumMult;

            double stdDev = GetStdDev(TInKelvin, wNum, species, isotopeID);

            // Returning null tells the calling code that the Doppler profile is
            // essentially an impulse function
            if (double.IsNaN(stdDev)) return new(null);

            // Check whether 4 Doppler standard deviations lies
            // within one wavenumber location in the output array
            int N = (int)(4 * stdDev / wNumRes + 0.5);

            // If the above condition is true, the Doppler profile
            // may be approximated as an impulse function
            if (N < 2) return new(null);

            double[] outArray = new double[N];
            double s2 = 2 * stdDev * stdDev;
            double sdTerm = stdDev * dTerm;

            for (int i = 0; i < N; i++)
            {
                double wNumCur = i * wNumRes;

                outArray[i] = Math.Exp(-wNumCur * wNumCur / s2) / sdTerm;
            }

            // Ensure a normalized profile, even when the profile width is
            // not much larger than the desired output resolution
            if (N < 5)
            {
                sdTerm = 0;
                for (int i = 1; i < N; i++) sdTerm += outArray[i];

                outArray[0] = 2 * (0.5 / wNumRes - sdTerm);
            }

            _cachedArray = outArray;

            return new(outArray);
        }

        /// <summary>
        /// Gets the half-width at half-height (HWHH) for a Doppler profile
        /// at the given temperature and wavenumber location, for the specified
        /// isotope of the specified species/ molecule
        /// </summary>
        /// <param name="TInKelvin">Temperature</param>
        /// <param name="wNumInInvCm">Wavenumber location</param>
        /// <param name="species"></param>
        /// <param name="isotopeID"></param>
        /// <returns></returns>
        public double GetHWHH(double TInKelvin, double wNumInInvCm, string species, int isotopeID)
        {
            update();

            if (!IsValid) return double.NaN;

            if (!_molWtDict.TryGetValue(species, out Dictionary<int, double>? molWtVal)) return double.NaN;
            if (!molWtVal.TryGetValue(isotopeID, out double mWt)) return double.NaN;

            double d = CDoppler * Math.Sqrt(TInKelvin / mWt) * wNumInInvCm;

            return d;
        }

        /// <summary>
        /// Gets the standard deviation of the representative normal distribution
        /// for the given Doppler profile, specified by the temperature, wavenumber location,
        /// species, and isotope
        /// </summary>
        /// <param name="TInKelvin"></param>
        /// <param name="wNumInInvCm"></param>
        /// <param name="species"></param>
        /// <param name="isotopeID"></param>
        /// <returns></returns>
        public double GetStdDev(double TInKelvin, double wNumInInvCm, string species, int isotopeID)
        {
            // Doppler profiles are essentially Gaussian/ normal distribution profiles
            // The Doppler half-width at half-height is directly proportional to the
            // standard deviation of the normal distribution (the proportionality factor
            // is C2)
            return GetHWHH(TInKelvin, wNumInInvCm, species, isotopeID) / C2;
        }

        /// <summary>
        /// Similar to GetValueFromStdDev(...), except that it accepts the HWHH as input
        /// </summary>
        /// <param name="wNumDiff"></param>
        /// <param name="HWHH"></param>
        /// <returns></returns>
        public static double GetValueFromHWHH(double wNumDiff, double HWHH)
        {
            double stdDev = HWHH / C2;// Math.Sqrt(2 * Math.Log(2));

            return GetValueFromStdDev(wNumDiff, stdDev);
        }

        /// <summary>
        /// Gets the height of the Doppler curve at the given wavenumber location, for
        /// the given standard deviation (which is related to the HWHH by the factor C2)
        /// </summary>
        /// <param name="wNumDiff">Wavenumber difference from line center, in cm-1</param>
        /// <param name="stdDev">Doppler standard deviation</param>
        /// <returns></returns>
        public static double GetValueFromStdDev(double wNumDiff, double stdDev)
        {
            double w2 = wNumDiff * wNumDiff, s2 = 2 * stdDev * stdDev;

            double v = Math.Exp(-w2 / s2);
            return v / stdDev / dTerm;
        }

        /// <summary>
        /// Similar to GetValuesFromStdDev(...), except that it accepts the HWHH as input
        /// </summary>
        /// <param name="wNumDiffs"></param>
        /// <param name="HWHH"></param>
        /// <returns></returns>
        public static double[] GetValuesFromHWHH(double[] wNumDiffs, double HWHH)
        {
            double stdDev = HWHH / Math.Sqrt(2 * Math.Log(2));

            return GetValuesFromStdDev(wNumDiffs, stdDev);
        }

        /// <summary>
        /// Same as GetValueFromStdDev(...), except for an array of wavenumber locations
        /// </summary>
        /// <param name="wNumDiffs">Wavenumber differences from line center, in cm-1</param>
        /// <param name="stdDev"></param>
        /// <returns></returns>
        public static double[] GetValuesFromStdDev(double[] wNumDiffs, double stdDev)
        {
            double s2 = 2 * stdDev * stdDev;
            double sdTerm = stdDev * dTerm;

            double[] v = new double[wNumDiffs.Length];

            for (int i = 0; i < wNumDiffs.Length; i++)
            {
                v[i] = Math.Exp(-wNumDiffs[i] * wNumDiffs[i] / s2) / sdTerm;
            }

            return v;
        }

        /// <summary>
        /// Sets parameters to implement low-resolution Doppler half-width binning
        /// The range of Doppler half-widths between the specified wavenumber range
        /// is discretized into "NLevels" number of bins
        /// For example: During a run from 1 to 10,000 cm-1, the range is geometrically
        /// divided into (say) 100,000 bins (NLevels would be specified as 100,000)
        /// Then, for a line whose center is at a certain wavenumber, the Doppler profile
        /// corresponding to the nearest bin is used, rather than the Doppler profile at
        /// the actual line center
        /// This greatly reduces the number of Doppler profiles which are calculated during the run
        /// </summary>
        /// <param name="startWNum">Start wavenumber in cm-1</param>
        /// <param name="endWNum">End wavenumber in cm-1</param>
        /// <param name="NLevels"></param>
        public void Initialize(double startWNum, double endWNum, int NLevels)
        {
            if (startWNum <= 0) startWNum = 1;

            _WNumMult = Math.Exp((Math.Log(endWNum / startWNum)) / NLevels);

            _cachedArray = null;

            _curWNum = -1;
        }

        public void SetMolecularWeights(string baseFolder)
        {
            _isUpdated = false;

            _molWtDict = null;

            string fName = $@"{baseFolder}\Data\molecularWeights.dat";

            if (!File.Exists(fName)) return;

            _molWtDict = [];

            string? curKey = null;

            using StreamReader sr = new(fName);
            string? line = sr.ReadLine();

            while (line != null)
            {
                string[] lineList = line.Split('\t');

                if (lineList.Length < 2)
                {
                    _molWtDict.Add(line, []);

                    curKey = line;
                }
                else
                {
                    if (curKey is null) continue;

                    int isotopeID = int.Parse(lineList[0]);
                    double mWt = double.Parse(lineList[1]) / 1000;

                    _molWtDict[curKey].Add(isotopeID, mWt);
                }

                line = sr.ReadLine();
            }
        }

        #endregion Public Methods

        #region Public Properties

        [MemberNotNullWhen(true, nameof(_molWtDict))]
        public bool IsValid
        {
            get { return _isValid; }
        }

        #endregion Public Properties

        #region Private Methods

        void update()
        {
            if (_isUpdated) return;

            _isValid = _molWtDict is not null;

            _isUpdated = true;
        }

        #endregion Private Methods

        #region Private Properties

        double[]? _cachedArray = [];

        // Constant property for thread safety
        const double C2 = 1.1774100225155; // Math.Sqrt(2 * Math.Log(2));

        const double CDoppler = 1.13246322E-08;

        double _curWNum = -1;

        // Const property for thread safety
        const double dTerm = 2.5066282746310; // Math.Sqrt(2 * Math.PI);

        bool _isUpdated = false;

        bool _isValid = false;

        // First key: molecule name; second key: isotope ID; value: molecular weight
        Dictionary<string, Dictionary<int, double>>? _molWtDict = null;

        double _WNumMult = 1;

        #endregion Private Properties
    }
}
