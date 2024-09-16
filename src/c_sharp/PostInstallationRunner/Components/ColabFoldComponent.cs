using System.Diagnostics;
using PostInstallationRunner.Util;

namespace PostInstallationRunner.Components;

public class ColabFoldComponent : IComponent
{
    /// <summary>
    /// Installs ColabFold.
    /// </summary>
    /// <returns>
    /// True if ColabFold is successfully installed, otherwise false.
    /// </returns>
    public bool Install()
    {
        #region Checks

        if (!File.Exists($"{ConstantPaths.TEMP_DIR}\\alma-colabfold-9-rootfs.tar"))
        {
            return false;
        }

        #endregion

        if (!Directory.Exists(@"C:\ProgramData\localcolabfold\"))
        {
            Directory.CreateDirectory(@"C:\ProgramData\localcolabfold\");
        }
        if (!Directory.Exists(@"C:\ProgramData\localcolabfold\scripts\"))
        {
            Directory.CreateDirectory(@"C:\ProgramData\localcolabfold\scripts\");
        }
        if (!Directory.Exists(@"C:\ProgramData\localcolabfold\storage\"))
        {
            Directory.CreateDirectory(@"C:\ProgramData\localcolabfold\storage\");
        }
        Process process = new Process
        {
            StartInfo =
            {
                FileName = "cmd.exe",
                UseShellExecute = false,
                CreateNoWindow = true,
                Arguments = "/C wsl --import almaColabfold9 C:\\ProgramData\\localcolabfold\\storage C:\\ProgramData\\IBCI\\temp\\alma-colabfold-9-rootfs.tar"
            }
        };
        try
        {
            process.Start();
            process.WaitForExit();
        }
        catch (Exception ex)
        {
            return false;
        }
        return true;
    }

    /// <summary>
    /// Uninstalls ColabFold.
    /// </summary>
    /// <returns>
    /// True if ColabFold is successfully uninstalled, otherwise false.
    /// </returns>
    public bool Uninstall()
    {
        Process process = new Process
        {
            StartInfo =
            {
                FileName = "cmd.exe",
                UseShellExecute = false,
                CreateNoWindow = true,
                Arguments = "/C wsl --unregister almaColabfold9"
            }
        };
        try
        {
            process.Start();
            process.WaitForExit();
        }
        catch (Exception ex)
        {
            return false;
        }

        if (!File.Exists(@"C:\ProgramData\localcolabfold\storage\ext4.vhdx"))
        {
            if (Directory.Exists(@"C:\ProgramData\localcolabfold\"))
            {
                Directory.Delete(@"C:\ProgramData\localcolabfold\", true);
            }
        }
        else
        {
            // Unregistering the WSL2 distro failed therefore return false
            return false;
        }
        return true;
    }

    /// <summary>
    /// Checks if ColabFold is installed or not on the system.
    /// </summary>
    /// <returns>
    /// True if ColabFold is installed, otherwise false.
    /// </returns>
    public bool IsInstalled()
    {
        Process process = new Process
        {
            StartInfo =
            {
                FileName = "cmd.exe",
                UseShellExecute = false,
                CreateNoWindow = true,
                Arguments = "/C wsl -d almaColabfold9 ls /home/rhel_user/localcolabfold/colabfold-conda/bin/colabfold_batch"
            }
        };
        try
        {
            process.Start();
            process.WaitForExit();
        }
        catch (Exception ex)
        {
            return false;
        }
        if (process.ExitCode != 0)
        {
            return false;
        }
        return true;
    }
}