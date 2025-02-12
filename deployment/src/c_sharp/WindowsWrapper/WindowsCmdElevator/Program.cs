namespace WindowsCmdElevator;

class Program
{
    static void Main(string[] args)
    {
        if (CmdElevator.IsAdministrator())
        {
            // If the current process is already elevated, run the desired command
            CmdElevator.RunCommand(args);
        }
        else
        {
            // If not elevated, restart the process with elevated rights
            //CmdElevator.RestartElevated(args);
            CmdElevator.RestartElevatedWindowsTasks();
        }
    }
}
