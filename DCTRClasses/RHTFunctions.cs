using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    public class RHTFunctions
    {
        #region Public Methods

        public static int Get2PowerLength(double sWNum, double eWNum, double wNumRes)
        {
            int N = (int)((eWNum - sWNum) / wNumRes + 0.5);
            int NPower2 = 1;
            while (NPower2 < N) NPower2 <<= 1;

            return NPower2;
        }

        public static int GetArrayLength(double sWNum, double eWNum, double wNumRes)
        {
            int N = (int)((eWNum - sWNum) / wNumRes + 0.5);

            return N;
        }

        /// <summary>
        /// Array of transmissivities for each wavenumber interval
        /// </summary>
        /// <param name="absCoeffs">Of size NX2, where N is the number of wave
        /// number bins. The first row contains wave number values, while the
        /// second row contains absorption coefficients (cm-1).</param>
        /// <param name="wNumAvg">Averaging wavenumber interval (cm-1)</param>
        /// <param name="opticalPathLength">Product of mole fraction and path length (cm)</param>
        /// <returns></returns>
        public static double[,] GetAveragedTransmissivities(double[,] absCoeffs, double wNumAvg,
            double opticalPathLength)
        {
            int M = absCoeffs.GetLength(0);

            double wRes = absCoeffs[1, 0] - absCoeffs[0, 0];

            int NAvg = (int)(wNumAvg / wRes + 0.5);

            int MAvg = M / NAvg;

            double[,] outArray = new double[MAvg, 2];

            int start, end;

            for (int i = 0; i < MAvg; i++)
            {
                start = i * NAvg;
                end = start + NAvg;

                for (int j = start; j < end; j++)
                {
                    outArray[i, 0] += absCoeffs[j, 0];
                    outArray[i, 1] += Math.Exp(-absCoeffs[j, 1] * opticalPathLength);
                }

                outArray[i, 0] /= NAvg;
                outArray[i, 1] /= NAvg;
            }

            return outArray;
        }

        public static double[] GetAveragedTransmissivities(double[] transmissivities,
            double startWNum, double wNumRes,
            double avgIntvl, out double[] outWNums, int startIndex = 0)
        {
            int N = transmissivities.Length;
            int NAvg = (int)(avgIntvl / wNumRes);
            int NShift = NAvg / 2;

            List<double> outList = [];
            List<double> outWNumList = [];

            int i = startIndex;
            double wNumCur = (startIndex + NShift) * wNumRes + startWNum;
            while (i < N - NAvg)
            {
                double tSum = 0;
                for (int j = i; j < i + NAvg; j++)
                {
                    tSum += transmissivities[j];
                }

                outWNumList.Add(wNumCur);
                wNumCur += NAvg * wNumRes;

                outList.Add(tSum / NAvg);

                i += NAvg;
            }

            outWNums = [.. outWNumList];

            return [.. outList];
        }

        public static double GetLorentzIntensityAtWavenumber(double intensity,
            double wNumDiff, double halfWidth)
        {
            double f = halfWidth / (halfWidth * halfWidth + wNumDiff * wNumDiff) / Math.PI;

            return intensity * f;
        }

        // Profile starts from 0 wavenumbers and is one-sided
        public static double[] GetLorentzProfile(double halfWidth,
            double wnRes, int N, double LorentzWavenumberRange)
        {
            double[] outArray = new double[N];
            double w2 = halfWidth * halfWidth;

            double prevVal = 0, curVal;

            for (int i = 0; i < N; i++)
            {
                double wnDiff = wnRes * i;

                if (LorentzWavenumberRange > 0 && wnDiff > LorentzWavenumberRange) break;

                if (wnRes > 3 * halfWidth)
                {
                    if (i == 0)
                    {
                        prevVal = Math.Atan(0.5 * wnRes / halfWidth) /
                            (wnRes * Math.PI);

                        outArray[i] = 2 * prevVal;
                    }
                    else
                    {
                        curVal = Math.Atan((i + 0.5) * wnRes / halfWidth) /
                            (wnRes * Math.PI);

                        outArray[i] = curVal - prevVal;

                        prevVal = curVal;
                    }

                    continue;
                }

                outArray[i] = halfWidth / (w2 + wnDiff * wnDiff) / Math.PI;
            }

            double oSum = outArray[0] + 2 * outArray.Where((item, ix) => ix > 0).Sum();

            return outArray;
        }

        /// <summary>
        /// Array of transmissivities at each wavenumber
        /// </summary>
        /// <param name="absCoeffs">Array of absorption coefficients
        /// First row - wavenumbers (cm-1)
        /// Second row - absorption coefficients (cm-1)</param>
        /// <param name="opticalPathLen">Product of mole fraction and path length (cm)</param>
        /// <returns></returns>
        public static double[,] GetTransmissivities(double[,] absCoeffs, double opticalPathLen)
        {
            int M = absCoeffs.GetLength(0), N = 2;

            double[,] outArray = new double[M, N];

            for (int i = 0; i < M; i++)
            {
                outArray[i, 0] = absCoeffs[i, 0];
                outArray[i, 1] = Math.Exp(-absCoeffs[i, 1] * opticalPathLen);
            }

            return outArray;
        }

        public static double[] GetWidthRange(double minTemp,
            double maxTemp, string species)
        {
            double[] minMaxWidths = new double[2];

            double min = Interpolation.Interpolate(_TempDict[species],
                _minWidthDict[species], minTemp);
            minMaxWidths[0] = Interpolation.Interpolate(_TempDict[species],
                _minWidthDict[species], maxTemp);
            minMaxWidths[0] = min < minMaxWidths[0] ? min : minMaxWidths[0];

            double max = Interpolation.Interpolate(_TempDict[species],
                _maxWidthDict[species], minTemp);
            minMaxWidths[1] = Interpolation.Interpolate(_TempDict[species],
                _maxWidthDict[species], maxTemp);
            minMaxWidths[1] = max > minMaxWidths[1] ? max : minMaxWidths[1];

            return minMaxWidths;
        }

        public static double[] GetWidthRange(double Temp, string species)
        {
            double[] minMaxWidths =
            [
                Interpolation.Interpolate(_TempDict[species],
                    _minWidthDict[species], Temp),
                Interpolation.Interpolate(_TempDict[species],
                    _maxWidthDict[species], Temp),
            ];
            return minMaxWidths;
        }

        #endregion Public Methods

        #region Public Properties
        #endregion Public Properties

        #region Private Methods

        static RHTFunctions()
        {
            using StreamReader sr = new(
                new MemoryStream(Properties.Resources.Min_Max_Line_Widths));
            string? line = sr.ReadLine();
            while (line != null)
            {
                string species = line;

                if (!_TempDict.ContainsKey(species))
                {
                    _TempDict.Add(species, []);

                    _minWidthDict.Add(species, []);

                    _maxWidthDict.Add(species, []);
                }

                while (line != "}")
                {
                    line = sr.ReadLine();

                    if (line is not null) line = line.Trim();

                    if (line == "{" || line == "}")
                    {
                        continue;
                    }

                    string[] lineList = line is not null ? line.Split('\t') : [];

                    double T = double.Parse(lineList[0]);
                    double minV = double.Parse(lineList[1]);
                    double maxV = double.Parse(lineList[2]);

                    _TempDict[species].Add(T);
                    _minWidthDict[species].Add(minV);
                    _maxWidthDict[species].Add(maxV);
                }

                line = sr.ReadLine();
            }
        }

        #endregion Private Methods

        #region Private Properties

        static readonly Dictionary<string, List<double>> _TempDict = [];

        static readonly Dictionary<string, List<double>> _minWidthDict = [];

        static readonly Dictionary<string, List<double>> _maxWidthDict = [];

        #endregion Private Properties
    }
}
