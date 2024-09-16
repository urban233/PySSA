using System.Diagnostics;

namespace PostInstallationRunner.Components;

public class ChimeraXComponent : IComponent
{
    public bool Install()
    {
        try
        {
            Process tmpProcess = new Process
            {
                StartInfo =
                {
                    FileName = "cmd.exe",
                    UseShellExecute = false,
                    CreateNoWindow = true,
                    Arguments = @$"C:\ProgramData\IBCI\PyDD\temp\install.bat"
                }
            };
            tmpProcess.Start();
            tmpProcess.WaitForExit();
        }
        catch (Exception ex)
        {
            return false;
        }
        return true;
    }

    public bool Uninstall()
    {
        try
        {
            Process tmpProcess = new Process
            {
                StartInfo =
                {
                    FileName = "cmd.exe",
                    UseShellExecute = false,
                    CreateNoWindow = true,
                    Arguments = @$"C:\ProgramData\IBCI\PyDD\temp\uninstall.bat"
                }
            };
            tmpProcess.Start();
            tmpProcess.WaitForExit();
        }
        catch (Exception ex)
        {
            return false;
        }
        return true;
    }

    public bool IsInstalled()
    {
        // This is only a simple check, this could be done in a more sophisticated way to avoid false-positive detections
        if (File.Exists(@"C:\ProgramData\IBCI\PyDD\bin\ChimeraX\bin\ChimeraX.exe"))
        {
            return true;
        }
        return false;
    }
}