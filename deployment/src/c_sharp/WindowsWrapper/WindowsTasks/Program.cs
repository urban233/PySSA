using WindowsTasks.Util;

namespace WindowsTasks;

class Program
{
    static int Main(string[] args)
    {
        if (!Wsl.CheckIfWslIsInstalled())
        {
            Console.WriteLine("WSL2 is NOT available.");
            return 1;
            //ConsoleRunner.RunCommandInCmd("wsl --install --no-distribution");
        }
        Console.WriteLine("WSL2 is available.");
        return 0;
    }
}
