from xml.dom import minidom


class xml_document():
    root = minidom.Document()

    def __init__(self, pathToXml):
        self.pathToXml = pathToXml

    def create_settings_xml_file(self):
        """This function create the settings xml with the format:

        <?xml version="1.0" ?>
        <root>
            <pdbStorage>
                <pdbpath DEFAULT_ATTRIBUTE=DEFAULT_PDB_PATH>
            </pdbStorage>
            <zipStorage>
                <zippath DEFAULT_ATTRIBUTE=DEFAULT_ZIP_PATH/>
            </zipStorage>
        </root>
        """
        DEFAULT_PDB_PATH = "/home/$USER/Documents"
        DEFAULT_ZIP_PATH = "/home/$USER/Downloads"
        DEFAULT_ATTRIBUTE = "name"

        root_node = self.root.createElement("root")
        self.root.appendChild(root_node)

        pdb_storage_node = self.root.createElement("pdbStorage")
        root_node.appendChild(pdb_storage_node)

        pdb_path_node = self.root.createElement("pdbpath")
        # init node/attribute with default values
        pdb_path_node.setAttribute(DEFAULT_ATTRIBUTE, DEFAULT_PDB_PATH)
        pdb_storage_node.appendChild(pdb_path_node)

        zip_storage_node = self.root.createElement("zipStorage")
        root_node.appendChild(zip_storage_node)

        zip_path_node = self.root.createElement("zippath")
        # init node/attribute with default values
        zip_path_node.setAttribute(DEFAULT_ATTRIBUTE, DEFAULT_ZIP_PATH)
        pdb_storage_node.appendChild(zip_path_node)

    def load_xml_in_memory(self):
        """This function loads a xml file into the memory.

        Note:
            This function should be used once to load the xml file into the
            memory.
        """
        return minidom.parse(self.pathToXml)

    def get_path(self, xml_file, lowerTag, attribute):
        """This functions returns the value of the path node.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            lowerTag (str):
                e.g. pdbpath or zippath node
            attribute (str):
                e.g. name
        """
        path_name = xml_file.getElementsByTagName(lowerTag)
        path = path_name[0].getAttribute(attribute)
        return path

    def set_value(self, xml_file, lowerTag, attribute, value):
        """This function changes a specific value in the xml file

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            lowerTag (str):
                 e.g. pdbpath or zippath node
            attribute (str):
                 e.g. name
            value (str):
                 new value which should be set to the attribute
        """
        pathName = xml_file.getElementsByTagName(lowerTag)
        pathName[0].setAttribute(attribute, value)

    def save_xml_file(self, xml_file):
        """This function saves the opened xml file. The path of the class will
        be used as save path.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
        """
        with open(self.pathToXml, "w") as f:
            f.write(xml_file.toxml())


if __name__ == '__main__':
    xml = xmlDocument("settings.xml")
    xml.save_xml_file()
