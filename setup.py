from setuptools import setup
from setuptools import Command
import setup_build_tasks


class CreateWinPackage(Command):
    """Setup command for creating the Windows package that gets deployed."""
    description = 'Create a custom ZIP file containing source and additional files'
    user_options = []

    def initialize_options(self):
        """Override the default initialization options."""
        pass

    def finalize_options(self):
        """Override the default finalization options."""
        pass

    def run(self):
        """Override the default run method."""
        tmp_build_task = setup_build_tasks.CreateWinPackageBuildTask()
        tmp_build_task.execute_task()


class MakeDocs(Command):
    """Build command for creating sphinx documentation."""
    description = 'Activates the .venv and runs the make.bat for building the sphinx docs.'
    user_options = []

    def initialize_options(self):
        """Override the default initialization options."""
        pass

    def finalize_options(self):
        """Override the default finalization options."""
        pass

    def run(self):
        """Override the default run method."""
        tmp_build_task = setup_build_tasks.MakeSphinxDocs()
        tmp_build_task.execute_task()


class MakeVenv(Command):
    """Build command for creating sphinx documentation."""
    description = 'Export the .venv information into a deployable format.'
    user_options = []

    def initialize_options(self):
        """Override the default initialization options."""
        pass

    def finalize_options(self):
        """Override the default finalization options."""
        pass

    def run(self):
        """Override the default run method."""
        tmp_build_task = setup_build_tasks.MakeVenvForDeployment()
        tmp_build_task.execute_task()


setup(
    cmdclass={
        'create_win_package': CreateWinPackage,
        'make_docs': MakeDocs,
        'make_venv': MakeVenv,
    },
)
