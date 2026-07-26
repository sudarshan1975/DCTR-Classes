using MathNet.Numerics.LinearAlgebra;
using MNum = MathNet.Numerics;

namespace DCTRClasses
{
    public class CubicSpline : Spline
    {
        #region Public Methods

        public void ReadFromBinaryStream(BinaryReader br)
        {
            _isUpdated = false;

            Count = br.ReadInt32();

            XList = [];
            YList = [];
            _coeffList = new double[4, Count - 1];

            for (int i = 0; i < Count; i++)
            {
                XList.Add(br.ReadDouble());
                YList.Add(br.ReadDouble());
            }

            for (int i = 0; i < Count - 1; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    _coeffList[j, i] = br.ReadDouble();
                }
            }
        }

        public void SetData(List<double> xList, List<double> yList, double startSlope = double.NaN, double endSlope = double.NaN)
        {
            set(xList, yList);

            updateCoefficients(startSlope, endSlope);
        }

        public void WriteToBinaryStream(BinaryWriter bw)
        {
            bw.Write(Count);

            for (int i = 0; i < Count; i++)
            {
                bw.Write(XList[i]);
                bw.Write(YList[i]);
            }

            for (int i = 0; i < Count - 1; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    bw.Write(_coeffList[j, i]);
                }
            }
        }

        #endregion Public Methods

        #region Public Properties
        #endregion Public Properties

        #region Private Methods

        void updateCoefficients(double startSlope = double.NaN, double endSlope = double.NaN)
        {
            if (Count < 2) return;

            int NCoeffs = 4 * (Count - 1);

            bool clampStart = !double.IsNaN(startSlope) && !double.IsInfinity(startSlope);
            bool clampEnd = !double.IsNaN(endSlope) && !double.IsInfinity(endSlope);

            Matrix<double> CoeffsMatrix = Matrix<double>.Build.Dense(NCoeffs, NCoeffs);

            Matrix<double> RHSMatrix = Matrix<double>.Build.Dense(NCoeffs, 1);

            for (int i = 0; i < Count - 1; i++)
            {
                int i4 = 4 * i;

                RHSMatrix[i4, 0] = YList[i];
                RHSMatrix[i4 + 1, 0] = YList[i + 1];
                RHSMatrix[i4 + 2, 0] = 0;
                RHSMatrix[i4 + 3, 0] = 0;

                double x = XList[i];
                double x2 = x * x;
                double x3 = x2 * x;

                CoeffsMatrix[i4, i4] = x3;
                CoeffsMatrix[i4, i4 + 1] = x2;
                CoeffsMatrix[i4, i4 + 2] = x;
                CoeffsMatrix[i4, i4 + 3] = 1;

                double x1 = XList[i + 1];
                double x12 = x1 * x1;
                double x13 = x12 * x1;

                CoeffsMatrix[i4 + 1, i4] = x13;
                CoeffsMatrix[i4 + 1, i4 + 1] = x12;
                CoeffsMatrix[i4 + 1, i4 + 2] = x1;
                CoeffsMatrix[i4 + 1, i4 + 3] = 1;

                if (i != Count - 2)
                {
                    CoeffsMatrix[i4 + 2, i4] = 3 * x12;
                    CoeffsMatrix[i4 + 2, i4 + 1] = 2 * x1;
                    CoeffsMatrix[i4 + 2, i4 + 2] = 1;
                    CoeffsMatrix[i4 + 2, i4 + 3] = 0;
                    CoeffsMatrix[i4 + 2, i4 + 4] = -3 * x12;
                    CoeffsMatrix[i4 + 2, i4 + 5] = -2 * x1;
                    CoeffsMatrix[i4 + 2, i4 + 6] = -1;
                    CoeffsMatrix[i4 + 2, i4 + 7] = 0;

                    CoeffsMatrix[i4 + 3, i4] = 6 * x12;
                    CoeffsMatrix[i4 + 3, i4 + 1] = 2;
                    CoeffsMatrix[i4 + 3, i4 + 2] = 0;
                    CoeffsMatrix[i4 + 3, i4 + 3] = 0;
                    CoeffsMatrix[i4 + 3, i4 + 4] = -6 * x12;
                    CoeffsMatrix[i4 + 3, i4 + 5] = -2;
                    CoeffsMatrix[i4 + 3, i4 + 6] = -0;
                    CoeffsMatrix[i4 + 3, i4 + 7] = 0;
                }
                else
                {
                    if (clampStart)
                    {
                        CoeffsMatrix[i4 + 2, 0] = 3 * XList[0] * XList[0];
                        CoeffsMatrix[i4 + 2, 1] = 2 * XList[0];
                        CoeffsMatrix[i4 + 2, 2] = 1;

                        RHSMatrix[i4 + 2, 0] = startSlope;
                    }
                    else
                    {
                        CoeffsMatrix[i4 + 2, i4] = 6 * XList[0];
                        CoeffsMatrix[i4 + 2, i4 + 1] = 2;
                    }

                    if (clampEnd)
                    {
                        CoeffsMatrix[i4 + 3, i4] = 3 * XList[Count - 1] * XList[Count - 1];
                        CoeffsMatrix[i4 + 3, i4 + 1] = 2 * XList[Count - 1];
                        CoeffsMatrix[i4 + 3, i4 + 2] = 1;

                        RHSMatrix[i4 + 3, 0] = endSlope;
                    }
                    else
                    {
                        CoeffsMatrix[i4 + 3, i4] = 6 * XList[Count - 1];
                        CoeffsMatrix[i4 + 3, i4 + 1] = 2;
                    }
                }
            }

            Matrix<double> invCoeffsMatrix = CoeffsMatrix.Inverse();

            Matrix<double> outMatrix = invCoeffsMatrix * RHSMatrix;

            _coeffList = new double[4, Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    _coeffList[j, i] = outMatrix[4 * i + j, 0];
                }
            }
        }

        #endregion Private Methods

        #region Private Properties
        #endregion Private Properties
    }
}
