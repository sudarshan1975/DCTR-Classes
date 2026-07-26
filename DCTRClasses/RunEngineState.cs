namespace DCTRClasses
{
    /// <summary>
    /// Specifies a discretized homogeneous slab within an inhomogeneous gas column
    /// </summary>
    public class RunEngineState
    {
        #region Public Methods
        #endregion Public Methods

        #region Public Properties

        // Path length: ending location of homogeneous slab within inhomogeneous column (cm)
        public double EndLength
        {
            get { return _endLength; }
            set { _endLength = value; }
        }

        // cm
        public double MidLength
        {
            get { return (_endLength + _startLength) / 2; }
        }

        // Dimensionless - 0 to 1
        public double MoleFraction
        {
            get { return _moleFraction; }
            set { _moleFraction = value; }
        }

        // cm
        public double PathLength
        {
            get { return _endLength - _startLength; }
        }

        // Path length: starting location of homogeneous slab within inhomogeneous column (cm)
        public double StartLength
        {
            get { return _startLength; }
            set { _startLength = value; }
        }

        public string StateDefinition
        {
            get { return _stateDefinition; }
            set { _stateDefinition = value; }
        }

        // Temperature (K)
        public double Temperature
        {
            get { return _temperature; }
            set { _temperature = value; }
        }

        // Total pressure (Atmospheres)
        public double TotalPressure
        {
            get { return _totalPressure; }
            set { _totalPressure = value; }
        }

        #endregion Public Properties

        #region Private Methods
        #endregion Private Methods

        #region Private Properties

        // Universal gas constant (J/mol-K)
        static double R = 8.314;

        // Avogadro number (molecules/mol)
        static double NA = 6.02214076e23;

        // Path length: ending location of homogeneous slab within inhomogeneous column (cm)
        double _endLength = 0;

        // Dimensionless - 0 to 1
        double _moleFraction = 0;

        // Path length: starting location of homogeneous slab within inhomogeneous column (cm)
        double _startLength = 0;

        string _stateDefinition = "";

        // Temperature (K)
        double _temperature = 300;

        // Total pressure (Atmospheres)
        double _totalPressure = 1;

        #endregion Private Properties

        // The string key specifies the species (eg.: "HITRAN CO2" or "HITEMP H2O")
        // Within the inner dictionary, the integer key specifies the isotope ID
        // The value is the calculated partition function value for the given species,
        // isotope, and temperature (of the homogeneous slab)
        public Dictionary<string, Dictionary<int, double>> PartitionFunctionValueDictionary { get; set; } =
            new Dictionary<string, Dictionary<int, double>>();

        // Constant value, will be calculated upon initialization
        public double IntensityScalingFactor { get; set; } = 1;

        public void Initialize(string speciesDef, string speciesDatabaseDef)
        {
            speciesDatabaseDef = PartitionFunction.GetDescriptionString(speciesDatabaseDef);
            PartitionFunctionValueDictionary[speciesDef] = PartitionFunction.GetValues(Temperature, speciesDatabaseDef);

            IntensityScalingFactor = 101325 * NA * 1e-6 / R;
            IntensityScalingFactor *= _totalPressure / _temperature;
            IntensityScalingFactor *= _moleFraction * PathLength;
        }

        public string GetDescription()
        {
            List<string> outStrList = [
                "T" + _temperature,
                "P" + _totalPressure,
                "MF" + _moleFraction,
                "L" + PathLength
            ];

            return string.Join(";", outStrList);
        }

        public override string ToString()
        {
            return GetDescription();
        }

        public RunEngineState DeepCopy()
        {
            RunEngineState outState = new()
            {
                Temperature = _temperature,
                TotalPressure = _totalPressure,
                MoleFraction = _moleFraction,
                StartLength = _startLength,
                EndLength = _endLength,
                PartitionFunctionValueDictionary =
                    PartitionFunctionValueDictionary.ToDictionary(item => item.Key, item => item.Value)
            };

            return outState;
        }

        public string GetFileDescription(double StartWaveNumber, double EndWaveNumber)
        {
            string outStr = ((int)StartWaveNumber).ToString("D5") + "_" +
                     ((int)EndWaveNumber).ToString("D5") + "_" +
                     _temperature.ToString() + "_" + _totalPressure.ToString() + "_" +
                     _moleFraction.ToString();

            return outStr;
        }
    }
}
