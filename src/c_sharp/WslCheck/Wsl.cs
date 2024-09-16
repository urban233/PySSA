using System.Diagnostics;

namespace WslCheck;

public class Wsl
{
    public static bool CheckIfWslIsInstalled()
    {
        // PowerShell command to check VirtualMachinePlatform feature
        string powerShellCommand =
            "$vmPlatform = Get-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform | Select-Object State; " +
            "if ($vmPlatform -ne $null) { " +
            "    if ($vmPlatform.State -eq 'Enabled') { " +
            "        Write-Output 'Virtual Machine Platform feature is enabled.'; " +
            "    } elseif ($vmPlatform.State -eq 'Disabled') { " +
            "        Write-Output 'Virtual Machine Platform feature is installed but disabled.'; " +
            "    } " +
            "} else { " +
            "    Write-Output 'Virtual Machine Platform feature is not installed.'; " +
            "} " +
            "$vmPlatform.State";

        // Execute PowerShell command
        string result = ExecutePowerShellCommand(powerShellCommand);

        // Process the result
        if (result.Contains("Virtual Machine Platform feature is enabled."))
        {
            Console.WriteLine("Virtual Machine Platform feature is enabled.");
            return true;
        }
        if (result.Contains("Virtual Machine Platform feature is installed but disabled."))
        {
            Console.WriteLine("Virtual Machine Platform feature is installed but disabled.");
            return false;
        }
        if (result.Contains("Virtual Machine Platform feature is not installed."))
        {
            Console.WriteLine("Virtual Machine Platform feature is not installed.");
            return false;
        }
        Console.WriteLine("Unknown state.");
        return false;
    }

    public static bool Install()
    {
        try
        {
            ExecutePowerShellCommand("wsl --install --no-distribution");
        }
        catch (Exception ex)
        {
            return false;
        }
        return true;
    }

    private static string ExecutePowerShellCommand(string command)
    {
        // Create process start info for PowerShell
        ProcessStartInfo psi = new ProcessStartInfo
        {
            FileName = "powershell.exe",
            Arguments = $"-Command \"{command}\"",
            UseShellExecute = false,
            RedirectStandardOutput = true,
            CreateNoWindow = true,
            Verb = "runas" // This requests the process to run with elevated privileges (admin)
        };

        // Create and start the process
        using (Process process = Process.Start(psi))
        {
            // Read the output (result of PowerShell script)
            using (var reader = process.StandardOutput)
            {
                string result = reader.ReadToEnd();
                return result.Trim(); // Trim whitespace and newlines from the result
            }
        }
    }
}