import os
import xml.etree.ElementTree as ET


def generate_rst_from_xml(xml_file, template_file, output_file):
  """
  Generate an RST file from an XML file based on a template.

  :param xml_file: Path to the XML file
  :param template_file: Path to the template RST file
  :param output_file: Path to the output RST file
  """
  # Read and parse the XML file
  tree = ET.parse(xml_file)
  root = tree.getroot()

  # Extract data from the XML file
  title = root.find('title').text

  # Extract card data
  card = root.find('card')
  if card is not None:
    card_title = card.find('title').text
    card_content = card.find('content').text
    card_rst = f".. card:: {card_title}\n\n     {card_content}\n"
  else:
    card_rst = ''

  # Extract description
  if root.find('description').text is not None:
    description = root.find('description').text.strip() if root.find('description') is not None else ''
    description_lines = description.split('\n')
    description = '\n'.join(line.lstrip() for line in description_lines)
  else:
    description = ''

  # Extract getting started guides
  getting_started_guides = ''
  for guide in root.find('getting_started_guides'):
    guide_title = guide.find('title').text
    steps = '\n'.join(f"{i + 1}. {step.text}" for i, step in enumerate(guide.find('steps')))
    getting_started_guides += f"\n{guide_title}\n" + "*" * len(guide_title) + f"\n{steps}\n\n"

  # Extract details
  details = '\n'.join(f"- {detail.text}" for detail in root.find('details'))

  # Extract toctree data
  toctree = root.find('toctree')
  if toctree is not None:
    maxdepth = toctree.get('maxdepth')
    items = [item.text for item in toctree.findall('item')]
    toctree_items = '\n'.join(f'    {item}' for item in items)
  else:
    maxdepth = ''
    toctree_items = ''
  if toctree_items.find("None") != -1:
    toctree_items = ""

  # Read the template RST file
  with open(template_file, 'r') as file:
    template = file.read()

  # Replace placeholders with actual data
  output = template.format(
    title=title,
    card=card_rst,
    description=description,
    getting_started_guides=getting_started_guides,
    details=details,
    maxdepth=maxdepth,
    toctree_items=toctree_items
  )

  # Write the output to a new RST file
  with open(output_file, 'w') as file:
    file.write(output)

  print(f"RST file generated successfully and saved to {output_file}")


def process_directory(input_dir, template_file, output_root):
  """
  Process all XML files in a directory and its subdirectories while preserving the directory hierarchy.

  :param input_dir: Path to the directory containing XML files
  :param template_file: Path to the template RST file
  :param output_root: Path to the root directory where output RST files will be saved
  """
  # Iterate over all files and subdirectories in the input directory
  for root, dirs, files in os.walk(input_dir):
    for file in files:
      # Check if the file is an XML file
      if file.endswith(".xml"):
        # Generate the output directory preserving the hierarchy
        relative_path = os.path.relpath(root, input_dir)
        output_dir = os.path.join(output_root, relative_path)

        # Create the output directory if it doesn't exist
        if not os.path.exists(output_dir):
          os.makedirs(output_dir)

        # Generate the output file path
        input_file = os.path.join(root, file)
        output_file = os.path.join(output_dir, os.path.splitext(file)[0] + ".rst")

        # Process the XML file
        generate_rst_from_xml(input_file, template_file, output_file)


if __name__ == "__main__":
  process_directory(
    r"C:\Users\martin\github_repos\PySSA\docs\source\_my_content\help",
    r"C:\Users\martin\github_repos\PySSA\docs\source\_templates\basic_help_template.rst",
    r"C:\Users\martin\github_repos\PySSA\docs\source\help"
  )
