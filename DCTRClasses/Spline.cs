using MNum = MathNet.Numerics;
using FunctionLibrary;

namespace DCTRClasses
{
    /// <summary>
    /// Template class
    /// </summary>
    public class Spline
    {
        #region Public Methods

        /// <summary>
        /// Returns y value for given x
        /// </summary>
        /// <param name="x">Goes from minVal to maxVal (i.e. - x is not normalized)</param>
        /// <returns></returns>
        public Result<double?> Query(double x)
        {
            Result<double[]?> coeffsResult = getCoefficients(x);

            if (!coeffsResult.TryGetValue(out double[]? coeffs, out string? message))
            {
                return new(null, false, message, Severity.WARNING);
            }

            if (coeffs.Length != 3)
            {
                return new(null, false, $"Invalid spline coefficient length:" +
                    $" {coeffs.Length} instead of 3", Severity.WARNING);
            }

            double x2 = x * x;
            double x3 = x2 * x;

            double y = coeffs[0] * x3 + coeffs[1] * x2 + coeffs[2] * x + coeffs[3];

            return new(y);
        }

        public Result<List<double>?> Query(List<double> xList)
        {
            List<double> yList = [];

            foreach (double x in xList)
            {
                Result<double?> yRes = Query(x);

                if (!yRes.TryGetValue(out double? y, out string? message))
                {
                    return new(null, false, message, Severity.WARNING);
                }

                yList.Add(y.Value);
            }

            return new(yList);
        }

        public Result<double?> QuerySlope(double x)
        {
            Result<double[]?> coeffsRes = getCoefficients(x);

            if (!coeffsRes.TryGetValue(out double[]? coeffs, out string? message))
            {
                return new(null, false, message, Severity.WARNING);
            }

            double x2 = x * x;

            double y = 3 * coeffs[0] * x2 + 2 * coeffs[1] * x + coeffs[2];

            return new(y);
        }

        public Result<List<double>?> QuerySlope(List<double> xList)
        {
            List<double> yList = [];

            foreach (double x in xList)
            {
                Result<double?> yRes = QuerySlope(x);

                if (!yRes.TryGetValue(out double? y, out string? message))
                {
                    return new(null, false, message, Severity.WARNING);
                }

                yList.Add(y.Value);
            }

            return new(yList);
        }

        public Result<double?> QuerySlope2(double x)
        {
            Result<double[]?> coeffsRes = getCoefficients(x);

            if (!coeffsRes.TryGetValue(out double[]? coeffs, out string? message))
            {
                return new(null, false, message, Severity.WARNING);
            }

            double y = 6 * coeffs[0] * x + 2 * coeffs[1];

            return new(y);
        }

        public Result<List<double>?> QuerySlope2(List<double> xList)
        {
            List<double> yList = [];

            foreach (double x in xList)
            {
                Result<double?> yRes = QuerySlope2(x);

                if (!yRes.TryGetValue(out double? y, out string? message))
                {
                    return new(null, false, message, Severity.WARNING);
                }

                yList.Add(y.Value);
            }

            return new(yList);
        }

        public Result<double?> QuerySlope3(double x)
        {
            Result<double[]?> coeffsRes = getCoefficients(x);

            if (!coeffsRes.TryGetValue(out double[]? coeffs, out string? message))
            {
                return new(null, false, message, Severity.WARNING);
            }

            double y = 6 * coeffs[0];

            return new(y);
        }

        public Result<List<double>?> QuerySlope3(List<double> xList)
        {
            List<double> yList = [];

            foreach (double x in xList)
            {
                Result<double?> yRes = QuerySlope3(x);

                if (!yRes.TryGetValue(out double? y, out string? message))
                {
                    return new(null, false, message, Severity.WARNING);
                }

                yList.Add(y.Value);
            }

            return new(yList);
        }

        #endregion Public Methods

        #region Public Properties

        public List<double> XList
        {
            get { return _xList; }
            set { _xList = value; }
        }

        public List<double> YList
        {
            get { return _yList; }
            set { _yList = value; }
        }

        public int Count
        {
            get { return _count; }
            set { _count = value; }
        }

        #endregion Public Properties

        #region Private Methods

        bool checkList(List<double> inpList)
        {
            if (inpList == null) return false;

            if (inpList.Count < 2) return false;

            int NInc = 0, NDec = 0;

            for (int i = 0; i < inpList.Count - 1; i++)
            {
                double diffVal = inpList[i + 1] - inpList[i];

                if (diffVal < -1e-12) NDec++;

                if (diffVal > 1e-12) NInc++;
            }

            if (NInc + NDec != inpList.Count - 1) return false;

            if (NInc != 0 && NDec != 0) return false;

            return true;
        }

        Result<double[]?> getCoefficients(double x)
        {
            update();

            if (!_isValid) return new(null, false, $"Invalid spline", Severity.WARNING);

            if (double.IsNaN(x) || double.IsInfinity(x))
            {
                return new(null, false, $"Invalid X value: [{x}]", Severity.WARNING);
            }

            if (x < XList[0] || x > XList[Count - 1])
            {
                return new(null, false, $"X value out of range: [{x}]", Severity.WARNING);
            }

            int index = _xList.BinarySearch(x);

            if (index < 0) index = ~index;
            index -= 1;
            if (index > Count - 2) index = Count - 2;

            double a = _coeffList[0, index], b = _coeffList[1, index], c = _coeffList[2, index], d = _coeffList[3, index];

            return new([a, b, c, d]);
        }

        protected void set(List<double> xList, List<double> yList)
        {
            if (xList == null || yList == null) return;

            if (xList.Count != yList.Count) throw new Exception("Mismatched number of entries");

            if (xList.Count < 2) throw new Exception("Lists have too few entries");

            for (int i = 0; i < XList.Count; i++)
            {
                if (double.IsNaN(xList[i]) || double.IsInfinity(xList[i])) throw new Exception("Invalid x value at index " + (i + 1));
                if (double.IsNaN(yList[i]) || double.IsInfinity(yList[i])) throw new Exception("Invalid y value at index " + (i + 1));
            }

            bool b = checkList(xList);

            if (!b) throw new Exception("X values should be monotonically increasing or decreasing");

            XList = new List<double>(xList);
            YList = new List<double>(yList);

            Count = XList.Count;
        }

        void update()
        {
            if (_isUpdated) return;

            _isUpdated = true;

            _isValid = _coeffList is not null && _xList is not null && _yList is not null;
        }

        protected static void writeMatrixToFile(MNum.LinearAlgebra.Matrix<double> matrix, string fName)
        {
            int M = matrix.RowCount, N = matrix.ColumnCount;

            using (StreamWriter sw = new StreamWriter(fName))
            {
                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (j != 0) sw.Write("\t");

                        sw.Write(matrix[i, j]);
                    }

                    sw.WriteLine();
                }
            }
        }

        #endregion Private Methods

        #region Private Properties

        List<double>? _xList = null;

        List<double>? _yList = null;

        int _count = 0;

        protected double[,]? _coeffList = new double[4, 0];

        protected bool _isUpdated = false;

        protected bool _isValid = false;

        #endregion Private Properties
    }
}
