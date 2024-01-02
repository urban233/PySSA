# # string conversion
#         str_conversion_1 = constants.SETTINGS_DIR.replace("\\", "/")
#         str_conversion_2 = str_conversion_1.replace(":", "")
#         str_conversion_3 = str_conversion_2.replace("C", "c")
#         settings_dir_unix_notation = f"/mnt/{str_conversion_3}"
#         # new way
#         # list of paths of the fasta dir and pdb dir which is needed by colabbatch
#         colabbatch_paths = [self.fasta_path, self.pdb_path]
#
#         # <editor-fold desc="Linux pre-process">
#         # lists with podman commands which need to be used for the prediction process
#
#         podman_machine_start_command = ["podman", "machine", "start"]
#         podman_machine_stop_command = ["podman", "machine", "stop"]
#         podman_container_delete_command = ["podman", "container", "rm", constants.CONTAINER_NAME]
#         podman_container_stop_command = ["podman", "container", "stop", constants.CONTAINER_NAME]
#         podman_container_create_command = ["podman", "run", "-itd", "-v", f"{settings_dir_unix_notation}:/home/ubuntu_colabfold/{os.getlogin()}", "--name", constants.CONTAINER_NAME, constants.IMAGE_NAME]
#         podman_container_start_command = ["podman", "container", "start", constants.CONTAINER_NAME]
#         podman_container_exec_colabbatch_command = ["podman", "container", "exec", constants.CONTAINER_NAME,
#                                                     "/home/ubuntu_colabfold/localcolabfold/colabfold-conda/bin/colabfold_batch"]
#
#         if sys.platform.startswith("linux"):
#             # podman virtual machine gets started to be able to run containers
#             powershell_result = subprocess.run(podman_machine_start_command)
#             if powershell_result.returncode == 0:
#                 logger.info("Podman machine started.")
#             elif powershell_result.returncode == 125:
#                 logger.warning("Podman machine is already running.")
#             else:
#                 logger.error(f"Podman machine could not be started! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("Podman machine could not be started!")
#
#             subprocess.run(podman_container_stop_command)
#             subprocess.run(podman_container_delete_command)
#             powershell_result = subprocess.run(podman_container_create_command)
#             if powershell_result.returncode != 0:
#                 logger.error(f"Podman container could not be created! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("Podman container could not be created!")
#         # </editor-fold>
#
#         # checks cases for args to use in colabbatch command
#         if self.prediction_configuration.templates == "pdb70" and self.prediction_configuration.amber_force_field is True:
#             logger.info("Run prediction with default pdb70 templates and with amber force field correction.")
#             colabbatch_args = ["--amber", "--templates"]
#         elif self.prediction_configuration.templates == "none" and self.prediction_configuration.amber_force_field is True:
#             logger.info("Run prediction with no templates and with amber force field correction.")
#             colabbatch_args = ["--amber"]
#         elif self.prediction_configuration.templates == "pdb70" and self.prediction_configuration.amber_force_field is False:
#             logger.info("Run prediction with default pdb70 templates and with no amber force field correction.")
#             colabbatch_args = ["--templates"]
#         elif self.prediction_configuration.templates == "none" and self.prediction_configuration.amber_force_field is False:
#             logger.info("Run prediction with no templates and with no amber force field correction.")
#             colabbatch_args = []
#         else:
#             logger.error(
#                 f"Invalid prediction configuration. templates: {self.prediction_configuration.templates}, "
#                 f"amber: {self.prediction_configuration.amber_force_field}",
#             )
#             raise ValueError("Invalid prediction configuration.")
#
#         if sys.platform.startswith("linux"):
#             # starts localcolabfold-container, which runs the prediction
#             powershell_result = subprocess.run(podman_container_start_command)
#             if powershell_result.returncode != 0:
#                 logger.error(f"Podman localcolabfold-container could not be started! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("Podman localcolabfold-container could not be started!")
#
#             # concatenating lists to create the complete prediction command run in the container
#             powershell_result = subprocess.run(podman_container_exec_colabbatch_command + colabbatch_paths + colabbatch_args)
#             if powershell_result.returncode != 0:
#                 logger.error(f"Podman exec command which runs colabbatch failed! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("Podman exec command which runs colabbatch failed!")
#
#             powershell_result = subprocess.run(podman_container_stop_command)
#             if powershell_result.returncode != 0:
#                 logger.error(f"Podman container could not be stopped! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("Podman container could not be stopped!")
#             powershell_result = subprocess.run(podman_container_delete_command)
#             # if powershell_result.returncode != 0:
#             #     logger.error(f"Podman container could not be deleted! "
#             #                  f"Command exited with code: {powershell_result.returncode}")
#             #     raise RuntimeError("Podman container could not be deleted")
#
#             # stops podman virtual machine after completing prediction process
#             powershell_result = subprocess.run(podman_machine_stop_command)
#             if powershell_result.returncode != 0:
#                 logger.error(f"Podman machine could not be shutdown! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("Podman machine could not be shutdown!")
#         elif sys.platform.startswith("win32"):
#             if not os.path.exists(constants.WSL_SCRATCH_DIR):
#                 subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch"])
#                 subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch/local_predictions"])
#                 subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch/local_predictions/fasta"])
#                 subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch/local_predictions/pdb"])
#
#             print(settings_dir_unix_notation)
#             subprocess.run(["wsl", "-d", "almaColabfold9", "cp", "-r", f"{settings_dir_unix_notation}/scratch/local_predictions/fasta", "/home/rhel_user/scratch/local_predictions/"])
#             try:
#                 # Run the 'ls' command in WSL to list the directory
#                 command = ["wsl", "-d", "almaColabfold9", "ls", "/home/rhel_user/scratch/local_predictions/fasta"]
#
#                 # Run the command and capture the output
#                 result = subprocess.run(command, capture_output=True, text=True, shell=True)
#
#                 # Check if the command was successful
#                 if result.returncode == 0:
#                     # Print the output
#                     logger.info(result.stdout)
#                     if result.stdout == "":
#                         logger.critical("There are no fasta files in the WSL2 fasta directory!")
#                         raise exception.FastaFilesNotFoundError("There are no fasta files in the WSL2 fasta directory!")
#                 else:
#                     # Print an error message
#                     print(f"Error: {result.stderr}")
#             except Exception as e:
#                 print(f"An error occurred: {e}")
#
#             # concatenating lists to create the complete prediction command run in the container
#             wsl_almacolabfold_exec_colabbatch_command = ["wsl", "-d", "almaColabfold9", "/home/rhel_user/localcolabfold/colabfold-conda/bin/colabfold_batch"]
#             complete_command = wsl_almacolabfold_exec_colabbatch_command + colabbatch_paths + colabbatch_args
#             logger.info(f"complete shell command: {complete_command}")
#             powershell_result = subprocess.run(
#                 complete_command, capture_output=True, text=True, shell=True,
#             )
#             if powershell_result.returncode != 0:
#                 logger.error(f"WSL2 almaColabfold9 exec command which runs colabbatch failed! "
#                              f"Command exited with code: {powershell_result.returncode}")
#                 raise RuntimeError("WSL2 almaColabfold9 exec command which runs colabbatch failed!")
#             elif powershell_result.returncode == 0:
#                 logger.info(f"Shell output: {powershell_result.stdout}")
#                 logger.info("Prediction process was successful, copying results ...")
#                 subprocess.run(
#                     [
#                         "wsl", "-d", "almaColabfold9", "cp", "-r", "/home/rhel_user/scratch/local_predictions/pdb",
#                         f"{settings_dir_unix_notation}/scratch/local_predictions",
#                      ],
#                 )
#                 subprocess.run(["wsl", "-d", "almaColabfold9", "rm", "-rf", "/home/rhel_user/scratch"])
#         else:
#             print("Unsupported operating system detected.")
#         # shuts down the WSL2 environment used by the podman virtual machine
#         powershell_result = subprocess.run(["wsl", "--shutdown"])
#         if powershell_result.returncode != 0:
#             logger.error(f"WSL2 could not be shutdown! "
#                          f"Command exited with code: {powershell_result.returncode}")
#             raise RuntimeError("WSL2 could not be shutdown!")