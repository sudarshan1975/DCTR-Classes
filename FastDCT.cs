using System.Numerics;

namespace DCTRClasses
{
    /// <summary>
    /// A Fast Computational Algorithm for the Discrete Cosine Transform
    /// WEN-HSIUNG CHEN, C. HARRISON SMITH, AND S.C.FRALICK
    /// IEEE Transactions on Communications 1 September 1977
    /// </summary>
    public class FastDCT
    {
        #region Public Methods

        // Assumes N is an integer power of 2 greater than 2
        public FastDCT(int N)
        {
            if (N < 4) throw new Exception($"Input size needs to be at least 4");

            bool isPowerOfTwo = BitOperations.IsPow2(N);

            if (!isPowerOfTwo) throw new Exception($"Input size needs to be an integer power of 2");

            Size = N;

            cosArray = new double[N];
            sinArray = new double[N];

            cosArrayBitReversed = new double[N];
            sinArrayBitReversed = new double[N];

            double c = Math.PI / (2 * N);

            for (int i = 0; i < N; i++)
            {
                cosArray[i] = Math.Cos(c * i);
                sinArray[i] = Math.Sin(c * i);
            }

            indexArray = new int[N];

            int n = 1, m = 2;

            while (n < N)
            {
                for (int i = n; i < m; i++)
                {
                    indexArray[i - n] = indexArray[i - n] << 1;
                    indexArray[i] = indexArray[i - n] + 1;
                }

                n = m;

                m <<= 1;
            }

            for (int i = 0; i < indexArray.Length; i++)
            {
                cosArrayBitReversed[i] = cosArray[indexArray[i]];
                sinArrayBitReversed[i] = sinArray[indexArray[i]];
            }
        }

        public double[] DCT(double[] inpArray, bool bitReverse = true)
        {
            double[] outArray = new double[Size];

            for (int i = 0; i < Size; i++) outArray[i] = inpArray[i];

            int f = Size;
            while (f > 2)
            {
                InPlaceButterfly(outArray, 0, f - 1);

                f >>= 1;
            }

            // Original paper multiplies both by cosArray[Size>>1]
            FlipInPlace(outArray, 0, 1, 1, 1, cosArray[Size >> 1], -cosArray[Size >> 1]);

            /*************************/

            int n1 = 4, n2 = 2, n3 = 1, n4, n5, n6, n7, n8;
            double mult = cosArray[Size >> 1];
            bool nFlag;

            int ix1, ix2;

            while (n1 < Size)
            {
                n4 = n1 + n2 + n3 - 1;
                n5 = n1 + n3;

                for (int i = 0; i < n3; i++)
                {
                    FlipInPlace(outArray, n4 - i, n5 + i, mult, mult, mult, -mult);
                }

                n4 = n1;
                n5 = n2;

                ix1 = 2;
                ix2 = 0;

                double c1, c2, s1, s2, wUpper1, wUpper2, wLower1, wLower2;

                while (n5 > 1)
                {
                    nFlag = false;

                    while (n4 < n1 << 1)
                    {
                        InPlaceButterfly(outArray, n4, n4 + n5 - 1, nFlag);

                        n4 += n5;
                        nFlag = !nFlag;
                    }

                    if (n5 > 2)
                    {
                        n4 = n5 >> 2;

                        n6 = n1 + n4;

                        while (n6 < n1 + n2)
                        {
                            c1 = cosArrayBitReversed[ix1 + ix2];
                            s1 = sinArrayBitReversed[ix1 + ix2];
                            c2 = cosArrayBitReversed[(ix1 << 1) - ix2 - 1];
                            s2 = sinArrayBitReversed[(ix1 << 1) - ix2 - 1];

                            wUpper1 = -c1;
                            wUpper2 = c2;
                            wLower1 = s1;
                            wLower2 = s2;

                            for (int i = 0; i < n4 << 1; i++)
                            {
                                n7 = n6 + i;
                                n8 = (n1 << 1) + n1 - n7 - 1;

                                if (i == n4)
                                {
                                    wUpper1 = -s1;
                                    wUpper2 = -s2;
                                    wLower1 = -c1;
                                    wLower2 = c2;
                                }

                                FlipInPlace(outArray, n7, n8, wUpper1, wLower1, wUpper2, wLower2);
                            }

                            n6 += (n4 << 2);
                            ix2++;
                        }
                    }

                    n4 = n1;
                    n5 >>= 1;

                    ix1 <<= 1;
                    ix2 = 0;
                }

                n1 <<= 1;
                n2 <<= 1;
                n3 <<= 1;
            }

            /*************************/

            int n;

            int k = 1, m = 4;
            n = 2;

            while (n < Size)
            {
                for (int i = 0; i < k; i++)
                {
                    ix1 = n + i;
                    ix2 = m - 1 - i;

                    FlipInPlace(outArray, ix1, ix2, sinArrayBitReversed[ix1],
                        cosArrayBitReversed[ix1], -sinArrayBitReversed[ix2], cosArrayBitReversed[ix2]);
                }

                k <<= 1;
                n <<= 1;
                m <<= 1;
            }

            if (bitReverse) ReorderByBitReverse(outArray);

            return outArray;
        }

        public static void FlipInPlace(double[] outArray, int ixUpper, int ixLower,
            double wUpperToUpper, double wLowerToUpper, double wUpperToLower, double wLowerToLower)
        {
            double d = outArray[ixUpper] * wUpperToUpper + outArray[ixLower] * wLowerToUpper;
            outArray[ixLower] = outArray[ixUpper] * wUpperToLower + outArray[ixLower] * wLowerToLower;
            outArray[ixUpper] = d;
        }

        public double[] IDCT(double[] inpArray, bool bitReverse = true)
        {
            double[] outArray = new double[Size];

            for (int i = 0; i < Size; i++) outArray[i] = inpArray[i] * 2 / Size;
            outArray[0] /= Math.Sqrt(2); // This is different from the original paper

            if (bitReverse) ReorderByBitReverse(outArray);

            /***********************************/

            int n, ix1, ix2;

            int k = 1, m = 4;
            n = 2;

            while (n < Size)
            {
                for (int i = 0; i < k; i++)
                {
                    ix1 = n + i;
                    ix2 = m - 1 - i;

                    FlipInPlace(outArray, ix1, ix2, sinArrayBitReversed[ix1], -sinArrayBitReversed[ix2],
                        cosArrayBitReversed[ix1], cosArrayBitReversed[ix2]);
                }

                k <<= 1;
                n <<= 1;
                m <<= 1;
            }

            /***********************************/

            double c = cosArray[Size >> 1];

            FlipInPlace(outArray, 0, 1, c, c, c, -c);

            /***********************************/

            int n1 = 4, n2 = 2, n3 = 1, n4, n5, n6, n7, n8;
            double mult;
            bool nFlag;

            double c1, c2, s1, s2, wUpper1, wUpper2, wLower1, wLower2;

            while (n1 < Size)
            {
                n5 = 2;
                n6 = n3;
                n7 = n6;
                n8 = (n6 << 1) - 1;

                c1 = cosArrayBitReversed[n7];
                s1 = sinArrayBitReversed[n7];
                c2 = cosArrayBitReversed[n8];
                s2 = sinArrayBitReversed[n8];

                wUpper1 = -c1;
                wUpper2 = s1;
                wLower1 = c2;
                wLower2 = s2;

                while (n5 <= n2)
                {
                    n4 = n1;
                    nFlag = false;

                    while (n4 < n1 << 1)
                    {
                        InPlaceButterfly(outArray, n4, n4 + n5 - 1, nFlag);

                        n4 += n5;
                        nFlag = !nFlag;
                    }

                    if (n5 < n2)
                    {
                        ix1 = n1 + (n5 >> 1);

                        while (ix1 < n1 + n2)
                        {
                            for (int i = 0; i < n5; i++)
                            {
                                ix2 = ix1 + i;

                                int j = n1 * 3 - ix2 - 1;

                                if (i == n5 >> 1)
                                {
                                    wUpper1 = -s1;
                                    wUpper2 = -c1;
                                    wLower1 = -s2;
                                    wLower2 = c2;
                                }

                                FlipInPlace(outArray, ix2, j, wUpper1, wLower1, wUpper2, wLower2);
                            }

                            n7++;
                            n8--;

                            c1 = cosArrayBitReversed[n7];
                            s1 = sinArrayBitReversed[n7];
                            c2 = cosArrayBitReversed[n8];
                            s2 = sinArrayBitReversed[n8];

                            wUpper1 = -c1;
                            wUpper2 = s1;
                            wLower1 = c2;
                            wLower2 = s2;

                            ix1 += (n5 << 1);
                        }

                        n6 >>= 1;
                        n7 = n6;
                        n8 = (n6 << 1) - 1;

                        c1 = cosArrayBitReversed[n7];
                        s1 = sinArrayBitReversed[n7];
                        c2 = cosArrayBitReversed[n8];
                        s2 = sinArrayBitReversed[n8];

                        wUpper1 = -c1;
                        wUpper2 = s1;
                        wLower1 = c2;
                        wLower2 = s2;
                    }

                    n5 <<= 1;
                }

                n1 <<= 1;
                n2 <<= 1;
                n3 <<= 1;
            }

            /*******************************/

            n1 = 4;
            n2 = 2;
            n3 = 1;
            mult = cosArray[Size >> 1];

            while (n1 < Size)
            {
                n4 = n1 + n2 + n3 - 1;
                n5 = n1 + n3;

                for (int i = 0; i < n3; i++)
                {
                    FlipInPlace(outArray, n4 - i, n5 + i, mult, mult, mult, -mult);
                }

                n1 <<= 1;
                n2 <<= 1;
                n3 <<= 1;
            }

            /*******************************/

            int f = 4;
            while (f <= Size)
            {
                InPlaceButterfly(outArray, 0, f - 1);

                f <<= 1;
            }

            return outArray;
        }

        public static void InPlaceButterfly(double[] inpArray, int start, int finish, bool negateFlag = false)
        {
            int s = start, f = finish, inc = 1;
            if (negateFlag)
            {
                s = finish;
                f = start;
                inc = -1;
            }

            double d;

            while (s - inc != f)
            {
                d = inpArray[s] + inpArray[f];
                inpArray[f] = inpArray[s] - inpArray[f];
                inpArray[s] = d;

                s += inc;
                f -= inc;
            }
        }

        #endregion Public Methods

        #region Public Properties

        public int Size
        {
            get { return _size; }
            private set { _size = value; }
        }

        #endregion Public Properties

        #region Private Methods

        void ReorderByBitReverse(double[] outArray)
        {
            double d;

            for (int i = 0; i < Size; i++)
            {
                if (indexArray[i] >= i) continue;

                d = outArray[indexArray[i]];
                outArray[indexArray[i]] = outArray[i];
                outArray[i] = d;
            }
        }

        #endregion Private Methods

        #region Private Properties

        // First array - cos(i*pi/(2*N)); second array - sin(i*pi/(2*N))
        readonly double[] cosArray = [];
        readonly double[] sinArray = [];

        readonly double[] cosArrayBitReversed = [];
        readonly double[] sinArrayBitReversed = [];

        // Bit-reversed indices
        readonly int[] indexArray = [];

        int _size = 4;

        #endregion Private Properties
    }
}
