import pathlib
import os
import zipfile
import subprocess
from setuptools import setup, find_packages
from setuptools import Command

import setup_util


class CreateWinPackage(Command):
    description = 'Create a custom ZIP file containing source and additional files'
    user_options = []

    project_root_path = pathlib.Path(__file__).parent.absolute()
    temp_path: pathlib.Path = pathlib.Path(project_root_path, "_build", "temp")

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            # Clean temporary build directory
            if setup_util.Directory.exists(self.temp_path):
                setup_util.Directory.purge(self.temp_path)
            setup_util.Directory.create_directory(self.temp_path)
            # Copy base windows package from resources to tmp directory

            # Copy new PySSA src files from GitHub repo to tmp/offline_win_package directory

            # Remove unnecessary directories

            # Compress base_win_package to win_package.zip archive

            # Clean directory
        except Exception as e:
            print(e)

        # Run a Windows command using subprocess
        # result = subprocess.run(['powershell', '-Command', 'Write-Host "This is a message from the PS>"'],
        #                         capture_output=True, text=True)
        # print(result.stdout)

        # with zipfile.ZipFile('my-python-project.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        #     # Add source files
        #     for root, dirs, files in os.walk('src'):
        #         for file in files:
        #             filepath = os.path.join(root, file)
        #             zipf.write(filepath, os.path.relpath(filepath, start='src'))
        #
        #     # Add additional files
        #     for root, dirs, files in os.walk('additional-files'):
        #         for file in files:
        #             filepath = os.path.join(root, file)
        #             zipf.write(filepath, os.path.relpath(filepath, start='.'))

        # Run another Windows command using subprocess
        #subprocess.run(['echo', 'ZIP file created successfully.'])


setup(
    name='my-python-project',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    cmdclass={
        'create_win_package': CreateWinPackage,
    },
)
