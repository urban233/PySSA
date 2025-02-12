using System.Diagnostics;
using Microsoft.VisualBasic;
using NLog;

namespace WindowsTasks.Util;

public class ConsoleRunner
{
    /// <summary>
    /// Logger for the class.
    /// </summary>
    private static readonly Logger Logger = LogManager.GetCurrentClassLogger();

    public static bool RunCommandInCmd(string theArguments)
    {
        Logger.Debug("Defining cmd.exe process ...");
        string tmpArgs = "/C " + theArguments;
        Process process = new Process
        {
            StartInfo =
            {
                FileName = "cmd.exe",
                UseShellExecute = false,
                CreateNoWindow = true,
                Arguments = tmpArgs
            }
        };
        Logger.Debug($"The arguments string for the process is: {tmpArgs}");
        try
        {
            Logger.Debug("Starting cmd.exe process ...");
            process.Start();
            Logger.Debug("Waiting for cmd.exe process ...");
            process.WaitForExit();
        }
        catch (Exception ex)
        {
            Logger.Error($"Process ended with error: {ex}");
            return false;
        }
        Logger.Info("Process ended.");
        return true;
    }
}