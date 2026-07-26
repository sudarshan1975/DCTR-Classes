using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    /// <summary>
    /// Defines an interface for a radiative heat transfer solver
    /// One implementation of this interface is the RapidSolver class
    /// </summary>
    public interface IRHTSolver
    {
        RunEngine? Parent { get; set; }

        void Initialize();

        void SetWaveNumberRange();

        void Run(double startWNum, double endWNum, bool requiresOffset, StateMachine stateMachine,
            Dictionary<string, DataBuffer> bufferDict);

        void WriteDetails(string folderName);

        void WriteSettings(string folderName);

        void SaveCurrentArray(bool requiresOverlap);

        void FinalizeCalcs(bool requiresOverlap);

        void UpdateReportDictionary(Dictionary<string, string> reportDict);
    }
}
