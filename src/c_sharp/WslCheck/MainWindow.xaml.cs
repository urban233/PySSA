using System.Windows;

namespace WslCheck;

/// <summary>
/// Interaction logic for MainWindow.xaml
/// </summary>
public partial class MainWindow : Window
{
    public MainWindow()
    {
        InitializeComponent();
        btnEnableWsl.Visibility = Visibility.Hidden;
        if (Wsl.CheckIfWslIsInstalled())
        {
            ProgressTextBlock.Text = "WSL2 is already installed.";
        }
        else
        {
            ProgressTextBlock.Text = "WSL2 is not installed.\nDo you want to enable it now?";
            btnEnableWsl.Visibility = Visibility.Visible;
        }

    }

    private void BtnEnableWsl_OnClick(object sender, RoutedEventArgs e)
    {
        try
        {
            btnEnableWsl.IsEnabled = false;
            if (Wsl.Install())
            {
                ProgressTextBlock.Text = "WSL2 was successfully installed.\nA restart is necessary to apply the changes.";
                btnEnableWsl.Visibility = Visibility.Hidden;
            }
            else
            {
                ProgressTextBlock.Text = "Installation process of WSL2 failed! Please try again.";
                btnEnableWsl.Visibility = Visibility.Visible;
            }
            btnEnableWsl.IsEnabled = true;
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex);
        }
    }
}