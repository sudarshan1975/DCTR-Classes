using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    public class TimingFunctions
    {
        #region Public Methods

        public static void AddTime(string desc)
        {
            if (!_stopwatchDict.TryGetValue(desc, out Stopwatch? sw)) return;

            _timeDict[desc] += sw.ElapsedTicks / _swFreq;

            _countDict[desc]++;

            sw.Stop();
        }

        public static void Clear()
        {
            _stopwatchDict = [];
            _timeDict = [];
            _countDict = [];
        }

        public static int GetCount(string desc)
        {
            if (!_countDict.TryGetValue(desc, out int outVal)) return -1;

            return outVal;
        }

        public static Dictionary<string, int> GetCounts()
        {
            Dictionary<string, int> outDict = [];

            foreach (string key in _countDict.Keys) outDict[key] = _countDict[key];

            return outDict;
        }

        public static double GetCPUTime(string desc)
        {
            if (!_timeDict.TryGetValue(desc, out double outVal)) return double.NaN;

            return outVal;
        }

        public static Dictionary<string, double> GetCPUTimes()
        {
            Dictionary<string, double> outDict = [];

            foreach (string key in _timeDict.Keys) outDict[key] = _timeDict[key];

            return outDict;
        }

        public static string GetTimeString()
        {
            List<string> reportStr = [];

            foreach (string dsc in _timeDict.Keys)
            {
                double t = _timeDict[dsc];
                int count = _countDict[dsc];

                reportStr.Add(dsc + ": Time: " + t + " s; Count: " + count);
            }

            return string.Join(Environment.NewLine, reportStr);
        }

        public static void InitializeTime(string desc, bool reset = false)
        {
            if (!_stopwatchDict.ContainsKey(desc))
            {
                _stopwatchDict[desc] = new();

                _timeDict[desc] = 0;

                _countDict[desc] = 0;
            }
            else if (reset)
            {
                _timeDict[desc] = 0;
                _countDict[desc] = 0;
            }

            _stopwatchDict[desc].Restart();
        }

        public static void ReportTime(string fileName = null)
        {
            string reportStr = GetTimeString();

            if (fileName != null)
            {
                using StreamWriter sw = new(fileName);
                sw.WriteLine(reportStr);

                return;
            }
        }

        public static void UpdateTime(string timeStr, double t)
        {
            _timeDict[timeStr] = t;
        }

        #endregion Public Methods

        #region Public Properties
        #endregion Public Properties

        #region Private Methods
        #endregion Private Methods

        #region Private Properties

        static Dictionary<string, double> _timeDict = [];
        static Dictionary<string, Stopwatch> _stopwatchDict = [];
        static Dictionary<string, int> _countDict = [];

        static readonly double _swFreq = Stopwatch.Frequency;

        #endregion Private Properties
    }
}
