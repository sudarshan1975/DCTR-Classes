using Serilog;

namespace DCTRClasses
{
    public class Program1
    {
        public static void Main(string[] args)
        {
            if (args.Length < 2)
            {
                Log.Information("No folder or file supplied");
                Log.Information("Terminating run");

                return;
            }

            try
            {
                string fInput = args[0];
                string fOutput = args[1];

                string? inpPath = Path.GetDirectoryName(fInput);
                string? inpPattern = Path.GetFileName(fInput);

                if (inpPath == null || inpPattern == null)
                {
                    Log.Warning($"Could not get directory or file information from string {fInput}");

                    return;
                }

                string[] inpFiles = [.. Directory.EnumerateFiles(inpPath, inpPattern)];

                Log.Information($"Found {inpFiles.Length} files with pattern" +
                    $" \"{inpPattern}\" in directory \"{inpPath}\"");

                LineEnsemble lineEnsemble = new();

                if (!Directory.Exists(fOutput))
                {
                    Log.Information($"Directory \"{fOutput}\" doesn't exist. Attempting to create it.");

                    Directory.CreateDirectory(fOutput);
                }

                foreach (string fName in inpFiles)
                {
                    string outFName = $"{fOutput}\\{Path.GetFileNameWithoutExtension(fName)}.bin";

                    Log.Information($"Reading file \"{fName}\"");

                    lineEnsemble.ReadFromFile(fName);

                    Log.Information($"Writing file \"{outFName}\"");

                    lineEnsemble.WriteBinaryFile(outFName);
                }
            }
            catch (Exception exc)
            {
                Log.Warning($"{exc}{Environment.NewLine}{exc.StackTrace}");
            }
        }
    }
}
