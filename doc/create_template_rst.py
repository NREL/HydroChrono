import os

def write_rst_template(filepath, overwrite=False):
    # Check if the file exists and if user wants to overwrite
    if os.path.exists(filepath) and not overwrite:
        return

    # Extract the title from the filename
    filename = os.path.basename(filepath).replace('.rst', '')  # Remove the .rst extension
    title = filename.replace('_', ' ').title()  # Convert to title case and replace underscores with spaces
    
    # Create a unique label for the document
    label = f"label-{filename}"

    template = f"""
.. _{label}:

{title}
{'=' * len(title)}

Introduction
------------

Provide a brief overview of what this page is about.

Section 1: Topic Name
---------------------

Description or content about this section.

Sub-section 1.1
^^^^^^^^^^^^^^^

Further details or subdivisions of Section 1.
"""

    with open(filepath, 'w') as file:
        file.write(template)

def process_directory(directory_path):
    # Ask user if they want to overwrite existing files
    choice = input("Do you want to overwrite existing .rst files? (y/n): ").lower()
    overwrite = choice == 'y'
    
    # Walk through the directory and its sub-directories
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in filenames:
            if filename.endswith('.rst'):
                filepath = os.path.join(dirpath, filename)
                write_rst_template(filepath, overwrite)

if __name__ == "__main__":
    root_directory = 'C:\\code\\HydroChrono\\doc'
    process_directory(root_directory)
