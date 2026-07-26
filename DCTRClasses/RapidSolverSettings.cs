using Serilog;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    public class RapidSolverSettings
    {
        #region Public Methods

        public void WriteToFolder(string folder)
        {
            using StreamWriter sw = new(folder + "\\Settings.dat");
            sw.WriteLine("NWidthLevels:\t" + NWidthLevels);

            sw.WriteLine("NShifts:\t" + NShifts);
        }

        #endregion Public Methods

        #region Public Properties

        // Number of wavenumber bins, over which to split the line intensity
        // (See the DCTR paper)
        public int NShifts
        {
            get { return _nShifts; }
            set
            {
                if (value < 2 || value > 8)
                {
                    Log.Warning($"Attempt to set invalid number of wavenumber bin splits ({value})");

                    return;
                }

                _nShifts = value;
            }
        }

        // Number of width levels for Lorentz half-width discretization
        // (See the DCTR paper)
        public int NWidthLevels
        {
            get { return _nWidthLevels; }
            set
            {
                if (_nWidthLevels < 2)
                {
                    Log.Warning($"Attempt to set invalid number of width levels ({value})");

                    return;
                }

                _nWidthLevels = value;
            }
        }

        #endregion Public Properties

        #region Private Properties

        // Number of wavenumber bins, over which to split the line intensity
        // (See the DCTR paper)
        int _nShifts = 2;

        // Number of width levels for Lorentz half-width discretization
        // (See the DCTR paper)
        int _nWidthLevels = 64;

        #endregion Private Properties
    }
}
