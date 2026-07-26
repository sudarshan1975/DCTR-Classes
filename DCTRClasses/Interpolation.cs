using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    public class Interpolation
    {
        public static double Interpolate(List<double> x, List<double> y, double xValue)
        {
            int ix = x.BinarySearch(xValue);// SortFunctions.BinarySearchListIndex(x, xValue);
            if (ix < 0) ix = ~ix;
            ix--;

            if (ix < 0) return y[0];
            if (ix >= x.Count - 1) return y[x.Count - 1];

            double f = (xValue - x[ix]) / (x[ix + 1] - x[ix]);

            return y[ix] * (1 - f) + y[ix + 1] * f;
        }
    }
}
