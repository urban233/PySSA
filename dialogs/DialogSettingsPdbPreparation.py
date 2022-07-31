from pymol import Qt

from utils import constants, tools
from uiForms.auto.auto_DialogSettingsPdbPreparation import Ui_Dialog


class DialogSettingsPdbPreparation(Qt.QtWidgets.QDialog):
    SETTINGS_FULL_FILENAME = constants.SETTINGS_DIR
    PDB_STORAGE_PATH_NODE = constants.PDB_STORAGE_PATH_NODE
    ZIP_STORAGE_PATH_NODE = constants.ZIP_STORAGE_PATH_NODE
    ATTRIBUTE = constants.ATTRIBUTE

    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.ui.txtPdbStorageDir.setEnabled(False)
        self.ui.txtZipStorageDir.setEnabled(False)
        # sets default values

        self.xmlObj = tools.SettingsXml(self.SETTINGS_FULL_FILENAME)
        self.xmlFile = self.xmlObj.load_xml_in_memory()
        self.ui.txtPdbStorageDir.setText(self.xmlObj.get_path(self.xmlFile,
                                                              self.PDB_STORAGE_PATH_NODE,
                                                              self.ATTRIBUTE))
        self.ui.txtZipStorageDir.setText(self.xmlObj.get_path(self.xmlFile,
                                                              self.ZIP_STORAGE_PATH_NODE,
                                                              self.ATTRIBUTE))

        self.ui.btnPdbStorageDir.clicked.connect(self.choosePdbStorageDir)
        self.ui.btnZipStorageDir.clicked.connect(self.chooseZipStorageDir)
        self.ui.btnCancel.clicked.connect(self.cancelDialog)
        self.ui.btnOk.clicked.connect(self.okDialog)

        self.setWindowTitle("PdbPreparation Settings")


    def getDirectoryPath(self):
        tmpDialog = Qt.QtWidgets.QFileDialog()
        tmpDialog.setFileMode(Qt.QtWidgets.QFileDialog.Directory)
        tmpDialog.setOption(Qt.QtWidgets.QFileDialog.ShowDirsOnly)
        tmpDialog.exec_()
        return tmpDialog.directory().path()

    # @SLOT()
    def choosePdbStorageDir(self):
        self.ui.txtPdbStorageDir.setText(self.getDirectoryPath())

    def chooseZipStorageDir(self):
        self.ui.txtZipStorageDir.setText(self.getDirectoryPath())

    def cancelDialog(self):
        self.close()

    def okDialog(self):
        self.xmlObj.set_value(self.xmlFile, self.PDB_STORAGE_PATH_NODE,
                              self.ATTRIBUTE, self.ui.txtPdbStorageDir.text())
        self.xmlObj.set_value(self.xmlFile, self.ZIP_STORAGE_PATH_NODE,
                              self.ATTRIBUTE, self.ui.txtZipStorageDir.text())
        self.xmlObj.save_xml_file(self.xmlFile)
        self.close()
