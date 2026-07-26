using MNUM = MathNet.Numerics;

namespace DCTRClasses
{
    public class LorentzSplitWeights
    {
        #region Public Methods

        public LorentzSplitWeights()
        {
            recalculate();
        }

        public double GetMaxErrorEstimate(double k)
        {
            double[] weights = GetWeights(k);

            int NShifts = ShiftList.Count;

            List<double> xArray = [];

            List<List<double>> lArrayList = [];
            List<double> lArrayNominal = [];

            int xix = 0;
            for (int x = -100; x < 101; x++)
            {
                xArray.Add(x);
                lArrayList.Add([]);

                for (int i = 0; i < NShifts; i++)
                {
                    lArrayList[xix].Add(weights[i] / (_lambda * _lambda + (x + ShiftList[i]) * (x + ShiftList[i])));
                }

                lArrayNominal.Add(1 / (_lambda * _lambda + (x + k) * (x + k)));

                xix++;
            }

            double[] lArraySum = [.. lArrayList.Select(item => item.Sum())];

            double maxError = lArrayNominal.Select((item, ix) => Math.Abs(lArraySum[ix] / item - 1)).Max();

            return maxError;
        }

        public double GetSumSquareErrorEstimate(double k)
        {
            double[] weights = GetWeights(k);

            int NShifts = ShiftList.Count;

            List<double> xArray = [];

            List<List<double>> lArrayList = [];
            List<double> lArrayNominal = [];

            double dx = 0.1;

            int xix = 0;
            for (double x = -100; x < 100 + dx / 2; x += dx)
            {
                xArray.Add(x);
                lArrayList.Add([]);

                for (int i = 0; i < NShifts; i++)
                {
                    lArrayList[xix].Add(weights[i] / (_lambda * _lambda + (x + ShiftList[i]) * (x + ShiftList[i])));
                }

                lArrayNominal.Add(1 / (_lambda * _lambda + (x + k) * (x + k)));

                xix++;
            }

            double[] lArraySum = [.. lArrayList.Select(item => item.Sum())];

            double sumSquareError = lArrayNominal.Select((item, ix) => Math.Pow(lArraySum[ix] - item, 2)).Sum();

            return sumSquareError;
        }

        public double[] GetWeights(double k)
        {
            int NShifts = ShiftList.Count;
            int N = NShifts - 1;

            double[] r = getRHSArray(k);
            double[] outWeights = new double[NShifts];

            if (NShifts == 2)
            {
                outWeights[0] = 1 - k;
                outWeights[1] = k;

                return outWeights;
            }

            outWeights[N] = 1;

            for (int i = 0; i < N; i++)
            {
                int iN = i * N;

                for (int j = 0; j < N; j++)
                {
                    // _flattenedInverseMatrix represents a flattened symmetric matrix
                    // So the index can be either i*N+j, or i+j*N
                    outWeights[i] += _flattenedInverseMatrix[iN + j] * r[j];
                }

                outWeights[N] -= outWeights[i];
            }

            return outWeights;
        }

        public double S(double m, double n)
        {
            return 1 / (S00_Inv + (m - n) * (m - n));
        }

        #endregion Public Methods

        #region Public Properties

        // Line half-width to spectral resolution ratio
        public double Lambda
        {
            get { return _lambda; }
            set
            {
                _lambda = value;

                recalculate();
            }
        }

        // Inverse of S00: 1/S00=4*Lambda^2
        public double S00_Inv { get; set; } = 4;

        public List<int> ShiftList
        {
            get { return _shiftList; }
            set
            {
                _shiftList = value;

                recalculate();
            }
        }

        #endregion Public Properties

        #region Private Methods

        double[,] getErrorCoefficientsMatrix()
        {
            int NShifts = ShiftList.Count;

            double[,] polyCoeffsMatrix = new double[NShifts - 1, NShifts - 1];

            for (int p = 0; p < NShifts - 1; p++)
            {
                for (int i = 0; i < NShifts - 1; i++)
                {
                    polyCoeffsMatrix[p, i] = getMatrixCoefficient(p, i);
                }
            }

            return polyCoeffsMatrix;
        }

        /// <summary>
        /// For a polynomial represented as: y=a_0+a_1*x+a_2*x^2+....+a_Order*x^Order
        /// </summary>
        /// <param name="p">Equation number</param>
        /// <param name="i">Polynomial coefficient index (0 to Order)</param>
        /// <returns></returns>
        double getMatrixCoefficient(int p, int i)
        {
            int NShifts = ShiftList.Count;

            double Xi = ShiftList[i], Xp = ShiftList[p], XNs = ShiftList[NShifts - 1];

            double outVal = S(Xi, Xp) - S(XNs, Xp) - S(Xi, XNs) + S(0, 0);

            return outVal;
        }

        double[] getRHSArray(double k)
        {
            int NShifts = ShiftList.Count;

            double[] outArray = new double[NShifts - 1];

            double S00 = S(0, 0), XNs = ShiftList[NShifts - 1];

            for (int p = 0; p < NShifts - 1; p++)
            {
                double Xp = ShiftList[p];

                outArray[p] = S00 + S(k, Xp) - S(XNs, Xp) - S(k, XNs);
            }

            return outArray;
        }

        void recalculate()
        {
            S00_Inv = 4 * _lambda * _lambda;

            int NShifts = ShiftList.Count;
            int N = NShifts - 1;

            _flattenedInverseMatrix = new double[N * N];

            double[,] p = getErrorCoefficientsMatrix();

            MNUM.LinearAlgebra.Matrix<double> m = MNUM.LinearAlgebra.Matrix<double>.Build.Dense(N, N);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    m[i, j] = p[i, j];
                }
            }

            MNUM.LinearAlgebra.Matrix<double> mInv = m.Inverse();
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    // m and mInv are both symmetric, so it doesn't matter if the
                    // flattened inverse matrix is ordered by rows first, or columns first
                    _flattenedInverseMatrix[i + j * N] = mInv[i, j];
                }
            }
        }

        #endregion Private Methods

        #region Private Properties

        double[]? _flattenedInverseMatrix = null;

        double _lambda = 1;

        List<int> _shiftList = [0, 1];

        #endregion Private Properties
    }
}
