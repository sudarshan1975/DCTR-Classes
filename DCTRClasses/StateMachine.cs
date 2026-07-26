using Serilog;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    /// <summary>
    /// "State" refers to a homogeneous slab within the inhomogeneous path
    /// </summary>
    public class StateMachine
    {
        #region Public Methods

        public void Echo()
        {
            List<RunEngineState> stateList = StateDefinitions.FirstOrDefault().Value;

            if (stateList.Count == 1) Log.Information("Homogeneous path");
            else Log.Information("Inhomogeneous path with {stateList.Count} slabs");
        }

        public double[] GetMinMaxWidths(Dictionary<string, string> speciesDefinitions,
            double wNumStart, double wNumEnd)
        {
            double minWidth = double.PositiveInfinity,
                maxWidth = double.NegativeInfinity;

            double[] TempArray = [100, 200, 300, 600, 900, 1200, 1500, 1800, 2100,
                2400, 2700, 3000, 3300, 3600, 3900, 4200, 4500, 4800, 5100, 5400, 5700, 6000];

            int NTemp = 22, NWNumBins = 300;

            foreach (string species in speciesDefinitions.Keys)
            {
                if (!StateDefinitions.ContainsKey(species))
                {
                    Log.Information($"Ignoring species: {species}");

                    continue;
                }

                string dBase = speciesDefinitions[species];

                string fileName = Parent.BaseFolder + @"\Data\" + dBase + " " + species + " min self widths.dat";
                var minSelfArray = DataFunctions.ReadArrayFromFile(fileName);

                fileName = Parent.BaseFolder + @"\Data\" + dBase + " " + species + " max self widths.dat";
                var maxSelfArray = DataFunctions.ReadArrayFromFile(fileName);

                fileName = Parent.BaseFolder + @"\Data\" + dBase + " " + species + " min air widths.dat";
                var minAirArray = DataFunctions.ReadArrayFromFile(fileName);

                fileName = Parent.BaseFolder + @"\Data\" + dBase + " " + species + " max air widths.dat";
                var maxAirArray = DataFunctions.ReadArrayFromFile(fileName);

                for (int i = 0; i < NTemp - 1; i++)
                {
                    double T = TempArray[i];

                    foreach (RunEngineState state in StateDefinitions[species])
                    {
                        if (state.Temperature < T || state.Temperature > TempArray[i + 1]) continue;

                        for (int j = 0; j < NWNumBins; j++)
                        {
                            double wNum = j * 100;

                            if (wNumStart > wNum + 100 || wNumEnd < wNum) continue;

                            double minCurWidth = minSelfArray[j, i] * state.MoleFraction + minAirArray[j, i] * (1 - state.MoleFraction);
                            double maxCurWidth = maxSelfArray[j, i] * state.MoleFraction + maxAirArray[j, i] * (1 - state.MoleFraction);

                            minCurWidth *= state.TotalPressure;
                            maxCurWidth *= state.TotalPressure;

                            if (minWidth > minCurWidth)
                            {
                                minWidth = minCurWidth;
                            }

                            if (maxWidth < maxCurWidth)
                            {
                                maxWidth = maxCurWidth;
                            }
                        }
                    }
                }
            }

            double[] outArray = [minWidth, maxWidth];

            return outArray;
        }

        public double[,] GetStateArrays()
        {
            List<string> stateDefList = [];
            List<double[]> outData = [];
            List<string> speciesList = ["CO2", "H2O"];

            int speciesIx = 0;
            foreach (string species in speciesList)
            {
                if (StateDefinitions.ContainsKey(species))
                {
                    foreach (RunEngineState state in StateDefinitions[species])
                    {
                        if (!stateDefList.Contains(state.StateDefinition))
                        {
                            stateDefList.Add(state.StateDefinition);
                        }

                        int ix = stateDefList.IndexOf(state.StateDefinition);

                        // Path distance, temperature, pressure,
                        // CO2 mole fraction, H2O mole fraction
                        while (outData.Count <= ix) outData.Add(new double[5]);

                        outData[ix][0] = state.MidLength;
                        outData[ix][1] = state.Temperature;
                        outData[ix][2] = state.TotalPressure;
                        outData[ix][3 + speciesIx] = state.MoleFraction;
                    }
                }

                speciesIx++;
            }

            if (outData.Count == 0) return null;

            double[,] outArray = new double[outData.Count, 5];
            for (int i = 0; i < outData.Count; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    outArray[i, j] = outData[i][j];
                }
            }

            return outArray;
        }

        public string GetStateSignature()
        {
            string outStr = "";

            outStr += GetStateSignature("CO2");
            outStr += GetStateSignature("H2O");

            return outStr;
        }

        public string GetStateSignature(string species)
        {
            if (!StateDefinitions.ContainsKey(species)) return "[]";

            string outStr = "[";

            List<string> outStrList = [];
            foreach (RunEngineState state in StateDefinitions[species])
            {
                string s = "";
                s += "T" + state.Temperature.ToString("0.##");
                s += "P" + state.TotalPressure.ToString("0.##");
                s += "M" + state.MoleFraction.ToString("0.####");
                s += "L" + state.StartLength.ToString("0.####") + ":" + state.EndLength.ToString("0.####");

                outStrList.Add(s);
            }

            outStr += string.Join(";", outStrList);

            outStr += "]";

            return outStr;
        }

        // Initialize state definitions
        public void Initialize(List<string> lineList)
        {
            StateDefinitions = [];
            //Variables = new Dictionary<string, double>();

            double curLength = 0;
            int curStateIndex = 1;

            int lineIx = 0;
            while (lineIx < lineList.Count)
            {
                string line = lineList[lineIx].Trim();

                string[] lList = line.Split('\t');

                string lLower = lList[0].ToLower();

                //switch (lLower)
                //{
                //    case "variable":
                //        {
                //            var sList = lList[1].Split('=');
                //            Variables.Add(sList[0], double.Parse(sList[1]));

                //            break;
                //        }
                //}

                if (lLower.StartsWith("slab"))
                {
                    string stateDefinition = $"Slab{curStateIndex}";
                    curStateIndex++;

                    List<string> speciesList = [];
                    double T = 0, P = 0, l = 0;
                    List<double> mfList = [];

                    while (line != "}")
                    {
                        lineIx++;

                        line = lineList[lineIx].Trim();

                        if (line == "{" || line == "}") continue;

                        lList = line.Split('\t');

                        lLower = lList[0].ToLower();

                        switch (lLower)
                        {
                            case "temperature":
                                {
                                    T = interpretString(lList[1]);

                                    break;
                                }
                            case "pressure":
                            case "totalpressure":
                                {
                                    P = interpretString(lList[1]);

                                    break;
                                }
                            case "length":
                            case "pathlength":
                                {
                                    l = interpretString(lList[1]);

                                    break;
                                }
                        }

                        if (lLower.EndsWith("molefraction"))
                        {
                            string[] cList = lList[0].Split(' ');

                            speciesList.Add(cList[0].ToUpper());
                            mfList.Add(interpretString(lList[1]));
                        }
                    }

                    for (int i = 0; i < speciesList.Count; i++)
                    {
                        string s = speciesList[i];
                        double mf = mfList[i];

                        if (!StateDefinitions.ContainsKey(s))
                        {
                            StateDefinitions.Add(s, []);
                        }

                        RunEngineState state = new()
                        {
                            Temperature = T,
                            TotalPressure = P,
                            MoleFraction = mf,
                            StartLength = curLength,
                            EndLength = curLength + l,
                            StateDefinition = stateDefinition
                        };

                        StateDefinitions[s].Add(state);
                    }

                    curLength += l;
                }

                lineIx++;
            }
        }

        double interpretString(string s)
        {
            //NCalc.Expression _evalExp = new NCalc.Expression(s);
            //foreach (string v in Variables.Keys)
            //    _evalExp.Parameters.Add(v, Variables[v]);

            //var obj = _evalExp.Evaluate();

            return Convert.ToDouble(s);
        }

        // Read homogeneous slab definitions from a disk file
        public void ReadFromFile(string fName)
        {
            FileName = fName;

            using StreamReader sr = new(fName);
            List<string> lineList = [.. sr.ReadToEnd().Split(
                [Environment.NewLine], StringSplitOptions.None)];

            Initialize(lineList);
        }

        public void WriteStateArray(string fileName)
        {
            double[,] stateArray = GetStateArrays();
            int M = stateArray.GetLength(0), N = stateArray.GetLength(1);

            using StreamWriter sw = new(fileName);
            sw.WriteLine("Length\tTemperature\tPressure\tCO2\tH2O");
            sw.WriteLine("double\tdouble\tdouble\tdouble\tdouble");

            for (int i = 0; i < M; i++)
            {
                List<string> lineList = [];
                for (int j = 0; j < N; j++)
                {
                    lineList.Add(stateArray[i, j].ToString());
                }

                sw.WriteLine(string.Join("\t", lineList));
            }
        }

        public void WriteStateSignature(string fileName)
        {
            using StreamWriter sw = new(fileName);
            sw.WriteLine(GetStateSignature());
        }

        public void WriteToStream(StreamWriter sw)
        {
            sw.WriteLine(FileName);
        }

        #endregion Public Methods

        #region Public Properties

        public string FileName { get; set; } = null;

        public RunEngine Parent { get; set; } = null;

        // RunEngineState is the class which defines the homogeneous slab
        // The string key is the state ID, while the value is the state definition
        public Dictionary<string, List<RunEngineState>> StateDefinitions { get; set; } = null;

        #endregion Public Properties

        #region Private Methods

        #endregion Private Methods

        #region Private Properties
        #endregion Private Properties
    }
}
