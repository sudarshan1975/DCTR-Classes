namespace DCTRClasses
{
    /// <summary>
    /// Used to derive the spectral range (start and end wavenumbers in cm-1)
    /// given the name of a binary line data file
    /// </summary>
    public class FileData
    {
        #region Public Methods

        public void SetFileName(string fName)
        {
            FileName = fName;

            string[] splitStrings = StringFunctions.Split(Path.GetFileName(fName), new List<string>() { "_", "-" }, out _);

            StartWavenumber = double.Parse(splitStrings[1]);
            EndWavenumber = double.Parse(splitStrings[2].Replace(".bin", ""));
        }

        public override string ToString()
        {
            return Path.GetFileName(FileName);
        }

        #endregion Public Methods

        #region Public Properties

        public double EndWavenumber
        {
            get { return _endWavenumber; }
            set { _endWavenumber = value; }
        }

        public string FileName
        {
            get { return _fileName; }
            set { _fileName = value; }
        }

        public double StartWavenumber
        {
            get { return _startWavenumber; }
            set { _startWavenumber = value; }
        }

        #endregion Public Properties

        #region Private Methods
        #endregion Private Methods

        #region Private Properties

        double _endWavenumber = 10000;

        string _fileName = "";

        double _startWavenumber = 200;

        #endregion Private Properties
    }
}
