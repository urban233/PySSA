from pyssa.internal.thread import tasks
from pyssa.presenter import main_presenter_async
from pyssa.util import constants


class MainTasksManager:
    """A manager class which manages the threads for the main tasks such as prediction and analysis."""

    prediction_task: "tasks.Task"
    distance_analysis_task: "tasks.Task"

    def __init__(self):
        self.prediction_task = None
        self.distance_analysis_task = None

    def start_prediction_task(self, the_prediction_task: "tasks.Task"):
        """Runs a structure prediction"""
        constants.PYSSA_LOGGER.info("Running only a prediction.")
        # No analysis after prediction
        self.prediction_task = the_prediction_task
        self.prediction_task.start()

    def check_if_prediction_task_is_finished(self) -> bool:
        if self.prediction_task is None:
            raise ValueError("Cannot check the state of the prediction task if the task object is None! Please check first if the task object is None!")
        return self.prediction_task.is_finished()

    def start_distance_analysis_task(self, the_distance_analysis_task: "tasks.Task"):
        """Runs a structure distance_analysis"""
        constants.PYSSA_LOGGER.info("Running only a distance_analysis.")
        # No analysis after distance_analysis
        self.distance_analysis_task = the_distance_analysis_task
        self.distance_analysis_task.start()

    def check_if_distance_analysis_task_is_finished(self) -> bool:
        if self.distance_analysis_task is None:
            raise ValueError("Cannot check the state of the distance_analysis task if the task object is None! Please check first if the task object is None!")
        return self.distance_analysis_task.is_finished()
