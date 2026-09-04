// See https://aka.ms/new-console-template for more information
using DCTRClasses;
using Serilog;

// Use the code in Program1.cs to convert line data from ASCII to binary format

//string folder = Directory.GetCurrentDirectory(); // Uncomment this line to run from the current directory
string folder = @"D:\RadHeatTransfer\OnlineCode\DCTRPackage"; // Comment this line to run from the current directory

if (args.Length > 0) folder = args[0];

Log.Logger = new LoggerConfiguration()
    .WriteTo.Console()
    .CreateLogger();

try
{
    Log.Information($"Running within base folder: {folder}");

    Log.Information($"Initializing partition function tables");

    // Load partition function data
    PartitionFunction.Initialize($@"{folder}\PartitionFunctionTables.bin");

    // Run engine which implements the DCTR solver
    RunEngine rEngine = new()
    {
        RunType = RunType.DCTR
    };

    Log.Information($@"Reading input file: {folder} \DCTRInput.dat");

    // Reads from the input file "DCTRInput.dat" in the run folder
    rEngine.Initialize(folder);

    string OutputFolder = $@"{folder}\Outputs\{rEngine.RunType}\{rEngine.RunType} Output_{DateTime.Now:yyyy_MM_dd_HH_mm_ss}";

    if (!Directory.Exists(OutputFolder)) Directory.CreateDirectory(OutputFolder);

    rEngine.OutputFolder = OutputFolder;

    // Write run engine parameters to console for verification
    rEngine.Echo();

    Log.Information($"Output folder={OutputFolder}");

    Log.Information($"Performing run");

    // Actual DCTR run
    rEngine.Run();

    Log.Information($"Run successful");
}
catch (Exception exc)
{
    Log.Warning($"{exc}{Environment.NewLine}{exc.StackTrace}");
}

try
{
    Log.Information("Press any key to continue...");

    Console.ReadKey();
}
catch
{

}

