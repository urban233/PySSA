==============
 Installation
==============

Step 1 - Installation of PySSA-Installer
========================================
To install the PySSA plugin, you need to first install the PySSA-Installer,
which then installs all needed dependencies.

1. Start pyssa-installer
------------------------
To start the installation process, double click on the
downloaded pyssa-installer.exe in your **Downloads** folder.

.. image:: ./assets/in_pyssa_in_downloads_exe.png
    :width: 60 %
    :align: center

2. Ready to install
-------------------
Click on *Install* to start the installation.

.. image:: ./assets/in_pyssa_in_install.png
    :width: 60 %
    :align: center

3. Installation process
-----------------------
Now the installation process is running. The installation progress is shown in the progress bar.

.. image:: ./assets/in_pyssa_in_installing_process.png
    :width: 60 %
    :align: center

4. Installation is finished
---------------------------
After the installation is finished, click on *Finish* to exit the setup.

.. image:: ./assets/in_pyssa_in_finish.png
    :width: 60 %
    :align: center


Step 2 - Make your computer ready for PySSA
===========================================
.. note::
    To install the PySSA plugin, you need two dependencies: the WSL2 and Local Colabfold.
    Both can be installed through the PySSA-Installer. The PySSA cannot be installed until
    these prerequisites are met.

1. PySSA-Installer
___________________________
Open the PySSA-Installer from your desktop or start menu.

.. image:: ./assets/pyin_in_wsl.png
    :width: 60 %
    :align: center

2. Installation of WSL2
-----------------------
Firstly, click on *Install* button, to start the installation of WSL2.

.. image:: ./assets/pyin_progress_install_wsl.png
    :width: 60 %
    :align: center

The progress of the installation process is shown in the progress bar, together with a
status message below the progress bar.

.. image:: ./assets/pyin_progress_install_wsl.png
    :width: 60 %
    :align: center

.. important::
    After the successful installation of WSL2, the computer needs a restart.
    It is recommended to restart the computer immediately, to be able to proceed with the
    installation of the Local Colabfold.

To finish the installation of the WSL2, you need to restart your computer!

.. image:: ./assets/pyin_in_wsl_finish_restart_yes.png
    :width: 60 %
    :align: center

.. note::
    After the restart, wait until a cmd prompt (black box) opens. It could be necessary
    to create a new user within the cmd prompt. You can choose whatever username and password
    you like. It does not matter for the further installation steps.

    It could also be, that an error message is displayed. If you see this, then just hit enter
    a

3. Installation of Local Colabfold
----------------------------------
.. important::
    The during the installation of the Local Colabfold, 17GB are downloaded.
    Check if you have enough disk space!
    The download can take a long time, depending on your internet connection,
    so please be patient!

Firstly, click on *Install* for the installation of Local Colabfold.

.. image:: ./assets/pyin_in_lc_install.png
    :width: 60 %
    :align: center

Now the installation process is running. The progress of the download is displayed through
the progress bar. The download can take a long time, this is the perfect time to get a
cup of coffee.

.. image:: ./assets/pyin_progress_downl_lc.png
    :width: 60 %
    :align: center

After the download finishes, the local colabfold needs to be imported into the WSL2,
just stay patient and lean back everything will be done automatically.

.. image:: ./assets/pyin_progress_importing_lc_wsl.png
    :width: 60 %
    :align: center

After the installation you can install PySSA.

.. image:: ./assets/pyin_in_lc_ok.png
    :width: 60 %
    :align: center


Step 3 - Installation of PySSA
==============================
Firstly, click on *Install* for the installation of PySSA.

.. image:: ./assets/pyin_in_pyssa.png
    :width: 60 %
    :align: center

Now the download process of PySSA_setup.exe is running. It is shown in the progress bar.

.. image:: ./assets/pyin_progress_downl_pyssa_exe.png
    :width: 60 %
    :align: center

After the download, the setup for installation appears. Click *Next* to continue.

.. image:: ./assets/in_pyssa_exe_next_begin.png
    :width: 60 %
    :align: center

Now you can accept the license agreement with a click on *Next*.
If you don't accept it, the installation can **not** begin.
So you don't come to this step.

.. image:: ./assets/in_pyssa_exe_next_agreement.png
    :width: 60 %
    :align: center

Click on *Install*.

.. image:: ./assets/in_pyssa_exe_install.png
    :width: 60 %
    :align: center

Now the installation process is running. It is shown in the progress bar.
# BOX FOR A "LONG" INSTALLATION TIME!

.. image:: ./assets/in_pyssa_exe_installing_process.png
    :width: 60 %
    :align: center

After the installation is finished, click on *Finish* to exit the setup.

.. image:: ./assets/in_pyssa_exe_finish.png
    :width: 60 %
    :align: center


Step 4 - Launch PySSA for the first time
========================================
- click icon
- open plugin from pymol
- default workspace

18. Start PyMOL-PySSA.
----------------------
To install the PySSA PyMOL plugin, you have to start PyMOL with by clicking on the
*PyMOL-PySSA* desktop icon.

.. image:: assets/images_win/win_pymol_icon.png
    :align: center

19. Navigate to *Plugin*.
-------------------------
After you successfully launched PyMOL, navigate in the menu bar to *Plugin*.

.. image:: assets/images_win/win_pymol_plugin.png
    :align: center

20. Open the *Plugin Manager*.
------------------------------
If you click on *Plugin* in the menu bar, a dropdown menu will occur. There you have to click on
*Plugin Manager*.

.. image:: assets/images_win/win_pymol_plugin_manager_click.png
    :align: center

21. Install New Plugin.
-----------------------
In the Plugin Manager, navigate to *Install New Plugin* and then click on *Choose file ...*.

.. image:: assets/images_win/win_pymol_plugin_manager.png
    :align: center

22. Open the PySSA.zip.
-----------------------
After you clicked on *Choose file ...* a file dialog will open. There you have to click on *.pyssa*
in the top bar of the explorer.

.. image:: assets/images_win/win_plugin_path_info.png
    :align: center

After you clicked on *.pyssa*, the folder will open and there you have to click on *pyssa.zip*.
And then on *Open*.

.. image:: assets/images_win/win_pymol_plugin_zip.png
    :align: center

23. Confirm Installation Path.
------------------------------
Next, a dialog will open, which displays a file path.

Just click on *OK*. And wait a little bit.

Do **not** change the path!

.. image:: assets/images_win/win_pymol_plugin_install_path.png
    :align: center

24. Confirm Installation.
-------------------------
If the installation was successful, a dialog will open which says that
the plugin was installed successfully.

There you have to click on *OK*.

.. image:: assets/images_win/win_pymol_plugin_install_finish.png
    :align: center

25. Open the PySSA Plugin.
--------------------------
To open the PySSA plugin, navigate to *Plugin* and click on *PySSA*.
The plugin should open, after a few seconds.

.. image:: assets/images_win/win_activate_plugin.png
    :align: center

26. Basic PyMOL-PySSA Window Setup.
-----------------------------------
To work with both, PyMOL and PySSA at the same time, you can split the window in two half's,
by dragging the PyMOL window to the right side of your screen and select the PySSA window on the
left side. After that you can resize both windows, if you navigate with your mouse to the center and
move the window to the right side.

.. image:: assets/images_win/win_plugin_screen.png
    :align: center





1. Search for "programs" in Windows search bar.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To remove the PySSA, search in the Windows search bar for "programs".

Then choose "Add or remove programs".

.. image:: assets/images_win/win_search_programs.png
    :align: center

2. Select Uninstall.
~~~~~~~~~~~~~~~~~~~~
Go through the list of programs and find the PySSA program.
After that, you click on *Uninstall*.

.. image:: assets/images_win/win_delete_programs.png
    :align: center

3. Uninstall PySSA.
~~~~~~~~~~~~~~~~~~~
To finally uninstall the PySSA click on *Uninstall* in the new dialog window, which appeared after
the last Uninstall-click. This uninstalls also the open-source version of PyMOL!

.. image:: assets/images_win/win_delete_programs_confirm.png
    :align: center

4. Search for Anaconda.
~~~~~~~~~~~~~~~~~~~~~~~
If you want to also remove Anaconda, search for it in the programs list.

.. image:: assets/images_win/win_anaconda_search.png
    :align: center

5. Select Uninstall.
~~~~~~~~~~~~~~~~~~~~
Click on Anaconda and then on *Uninstall*.

.. image:: assets/images_win/win_anaconda_select_uninstall.png
    :align: center

6. Uninstall Anaconda.
~~~~~~~~~~~~~~~~~~~~~~
To finally uninstall the Anaconda click on *Uninstall* in the new dialog window, which appeared after
the last Uninstall-click.

.. image:: assets/images_win/win_anaconda_uninstall.png
    :align: center

Advanced way
^^^^^^^^^^^^

1. Select unins000.
~~~~~~~~~~~~~~~~~~~
First you select *unins000* with a click on it.

You can find it at: C:\\Users\\your-username\\.pyssa

.. image:: assets/images_win/win_program_dir_uninstall.png
    :align: center

2. Uninstall PySSA.
~~~~~~~~~~~~~~~~~~~
To finish the uninstallation process, click on *Yes*.

.. image:: assets/images_win/win_uninstall_introduction.png
    :align: center

3. Finished Uninstallation.
~~~~~~~~~~~~~~~~~~~~~~~~~~~
If the uninstall was successful, a dialog window appears, which informs you, that the
uninstall has been done.

Then you click on *OK*.

.. image:: assets/images_win/win_finish_uninstall.png
    :align: center




