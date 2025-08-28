# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 20:33:38 2025

@author: Laila
"""

# -*- coding: utf-8 -*-
"""
Phylogenetic Tree Pipeline: Project 17  
Created on: 17.12.24 - 25.01.25  
Author: Laila Rothe  

Description:  
This program maps data onto phylogenetic trees, automating the process of tree creation, alignment and design.
Optional upload to iTOL for advanced visualization and sharing. 
Streamline the analysis of phylogenetic relationships with an emphasis on customization, visual clarity, and integration with online tools.  

Key Features:  
- **Input:** Excel file containing sequence information, metadata, and outgroup specifications.  
- **Output:** A designed phylogenetic tree in multiple formats (PDF, Newick, Nexus, etc.).  
- **Automated steps:**  
  1. Generate dictionaries from the Excel table.  
  2. Perform sequence alignment using MAFFT.  
  3. Construct a phylogenetic tree using IQ-TREE2.  
  4. Design the tree with the following features:  
     - Assign unique colors to each species.  
     - Root the tree at the specified outgroup.  
     - Highlight significant differences (bootstrap > 70%) with an asterisk (*).  
     - Display genus and gene information in the title for context.  
     - Add all boolean metadata (e.g., morphological traits) as visual markers on the tree.  
     - Include a legend for metadata visualization.  
     - Use bold font to emphasize taxa of primary interest.  
  5. Save the designed tree in the desired output path.  
  6. Optionally upload the tree to iTOL for online visualization and collaboration.  

Dependencies:  
To run this program, ensure the following Python packages and software are installed:  
- **Python 3.8 or higher**  
- **Required Python packages (install via pip):**  
  - `pandas`  
  - `ete3`  
  - `numpy`  
  - `matplotlib`  
  - `itolapi`  
- **External software:**  
  - MAFFT (for sequence alignment)  
  - IQ-TREE2 (for phylogenetic tree construction)  

Usage:  
1. Prepare an Excel file with the required data (example table provided; formatting must match the example).  
2. Provide the file paths for the input Excel file, the MAFFT and IQ-TREE executables, and the output folder.  
3. Run the script and follow the prompts.  
4. The pipeline will produce an aligned FASTA file, a phylogenetic tree, and a fully designed tree visualization.  
5. If you have a subscription to iTOL, the pipeline can automatically upload the tree for online access and sharing.  

"""

import subprocess                                                              # For calling the powershell
import pandas as pd                                                            # For working with the input data
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace              # For visuylization directly in python
import os
from pathlib import Path                                                       # For working with file paths
import sys                                                                     # For exiting the program in case of critical errors
from itolapi import Itol                                                       # iTOL API classes for tree handling and exporting
import matplotlib.colors as mcolors                                            # For validating colors

#Function to get all the right paths for the calculation
def get_valid_path(prompt, is_file=True):
    while True:
        path = input(prompt)                                                   # Prompt the user for input
        if os.path.exists(path) and ((is_file and os.path.isfile(path)) or (not is_file and os.path.isdir(path))): # Check if the path exists and whether it's a file or directory based on the is_file flag
            return path                                                        # Return the valid path
        print("Invalid path. Please try again.")


def validate_excel_input(excel_path):
    issues = []  # List to store all warnings and errors

    try:
        xls = pd.ExcelFile(excel_path)
    except Exception as e:
        print(f"Fehler beim Öffnen der Excel-Datei: {e}")
        sys.exit(1)
    
    df = pd.read_excel(xls)

    # Prüfen, ob die erste Zeile alle benötigten Spalten enthält
    required_columns = {'seqid', 'species', 'Gene', 'Genus', 'Outgroup', 'sequence', 'background_color'}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        issues.append(f"🚨 ERROR: Missing required columns in the Excel file! 🚨\n❌ Missing columns: {missing_columns}\n➡ Please update the Excel table to include all required columns and restart the program.")
    
    # Check if 'seqid' values are unique
    if 'seqid' in df.columns and df['seqid'].duplicated().any():
        duplicate_rows = df[df['seqid'].duplicated(keep=False)].sort_values(by='seqid')
        issues.append("\n🚨 Warning: Duplicate 'seqid' values found!")
        for seqid in duplicate_rows['seqid'].unique():
            rows = duplicate_rows[duplicate_rows['seqid'] == seqid].index.tolist()
            issues.append(f"   - Duplicate '{seqid}' in rows: {', '.join(map(str, rows))}")
    
    # Prüfen, ob 'background color' gültige Farben enthält
    if 'background_color' in df.columns:
        invalid_rows = df.loc[~df['background_color'].isin(mcolors.CSS4_COLORS) & df['background_color'].notna()]
    
        if not invalid_rows.empty:
            issues.append("\n🚨 ERROR: Invalid colors detected in 'background_color'! 🚨")
            for index, row in invalid_rows.iterrows():
                issues.append(f"❌ Error: Invalid color '{row['background_color']}' in row {index + 2}")  # Excel row number
        
        if df['background_color'].isna().any():
            missing_rows = df[df['background_color'].isna()].index + 2  # Excel row numbers
            issues.append("\n🚨 ERROR: Missing colors in 'background_color'! 🚨")
            issues.append(f"❌ Missing colors in rows: {list(missing_rows)}")
    
    # If there are any issues, print them all and ask for user decision
    if issues:
        print("\n".join(issues))
        sys.exit("❌ Program terminated due to errors in the Excel file.")
   
    # Prüfen, welche Spalten nur 0 oder 1 enthalten
    binary_columns = [col for col in df.columns if df[col].dropna().isin([0, 1]).all()]
    if binary_columns:
        issues.append(f"Diese Merkmale werden im Baum hervorgehoben: {binary_columns}")
    
    return True

# Function to create dictionaries from Excel
def create_dictionaries_from_excel(excel_path):
    df = pd.read_excel(excel_path)
    data_dict = {}                                                             # Initialize an empty dictionary to store all column-based dictionaries
    
    # Handle Gene, Genus, and Outgroup separately (they should not be tied to seqid)
    data_dict['Gene'] = {'Gene': df['Gene'].iloc[0]}                           # Using the first row as the value for 'Gene'
    data_dict['Genus'] = {'Genus': df['Genus'].iloc[0]}                        # Using the first row as the value for 'Genus'
    data_dict['Outgroup'] = {'Outgroup': df['Outgroup'].iloc[0]}               # Using the first row as the value for 'Outgroup'
    
    for column in df.columns:                                                  # Loop through the other columns to create dictionaries for each seqid
        if column not in ['Gene', 'Genus', 'Outgroup']:                        # Skip the general columns
            column_dict = {}                                                   # Initialize an empty dictionary for the current column
            for index, row in df.iterrows():                                   # Loop through each row to populate the dictionary
                sequence_id = row['seqid']                                     # Ensure we are referencing seqid for each row
                value = row[column]                                            # Check if the value is not NaN
                if pd.notna(value):                                            # Only add to the dictionary if the value is not NaN
                    column_dict[sequence_id] = value
            data_dict[column] = column_dict                                    # Store the column dictionary in the main dictionary
    
    return data_dict

#Function to generate a fasta file from the excel table
def generate_fasta_from_dict(data_dict, output_fasta_path):
    with open(output_fasta_path, 'w') as fasta_file:                           # Open the output FASTA file in write mode
        for seqid, sequence in data_dict['sequence'].items():                  # Iterate over the sequence-related data (sequences and species)
            species = data_dict['species'].get(seqid, 'Unknown')               # Get the species for the current seqid from the 'species' dictionary
            fasta_file.write(f">{species} {seqid}\n")                          # Write the FASTA header in the format: >Species seqid
            fasta_file.write(f"{sequence}\n")                                  # Write the sequence data
    
    return output_fasta_path                                                   # Return the path to use in subsequent steps

# Function to run MAFFT for alignment
def run_mafft(mafft_path, fasta_path, aligned_fasta):
    cmd = [mafft_path, "--auto", fasta_path]
    with open(aligned_fasta, "w") as outfile:
        subprocess.run(cmd, stdout=outfile, stderr=subprocess.DEVNULL, check=True)        # Run the command

# Function to calculate the phylogenetic tree using IQ-TREE
def calculate_phylogenetic_tree(iqtree2_path, output_path, output_folder):
    prefix = os.path.join(output_folder, "iqtree_output")
    cmd = [
        iqtree2_path,
        "-s", output_path,
        "-nt", "AUTO",
        "-m", "TESTMERGE",
        "-B", "1000",
        "-pre", prefix,
    ]
    subprocess.run(cmd, check=True)

    treefile = prefix + ".treefile"
    if not os.path.exists(treefile):
        raise FileNotFoundError(treefile)
    return treefile                          # Return the path to the treefile generated by IQ-TREE

#Function to apply background colors
def apply_node_colors(tree, data_dict):
    for leaf in tree.get_leaves():                                             # Loop through each leaf in the tree   
        normalized_leaf_name = leaf.name.replace("_", "")                      # Remove underscores from leaf names
        for seqid in data_dict['background_color']:                            # Iterate through all seqids in background_color
           if seqid in normalized_leaf_name:                                   # Check if the seqid is present anywhere in the leaf name
                color = data_dict['background_color'][seqid]                   # Get the color for the matched seqid
                species_style = NodeStyle()
                species_style["bgcolor"] = color                               # Assign color to the node style
                leaf.set_style(species_style)                                  # Apply the style to the leaf
                break                                                          # Exit the loop once the color is applied
    
    return tree

#Function to Load the tree and set colors and outgroup
def get_tree_with_colors(tree_path, data_dict):
    tree = Tree(tree_path)                                                     # Load the phylogenetic tree
    default_style = NodeStyle()                                                # Define the default node style
    default_style["bgcolor"] = "white"                                         # Default color for unmatched seqid
    apply_node_colors(tree, data_dict)                                         # Apply colors and styles using seqid from data_dict
    outgroup_dict = data_dict.get('Outgroup', {})                              # Safely get the 'Outgroup' dictionary, default to an empty dictionary
    outgroup_seqids = set(str(seqid) for seqid in outgroup_dict.values() if pd.notna(seqid))  # Ensure all values are strings
    outgroup_leaves = [leaf for leaf in tree.get_leaves() if any(seqid in leaf.name for seqid in outgroup_seqids)]
    if outgroup_leaves:                                                        # Set the outgroup if found
        mrca = tree.get_common_ancestor(outgroup_leaves)                       # Find the MRCA (most recent common ancestor) of all outgroup leaves
        if mrca != tree:                                                       # Check if the MRCA is not already the root
            tree.set_outgroup(mrca)                                            # Root the tree at the MRCA of the outgroup leaves
        tree.ladderize(direction=1)                                            # Ensure the tree is ladderized with the outgroup clade at the bottom
     
    return tree, TreeStyle()                                                   # Return the tree with the applied style

#Function to add boolean values from excel table
def add_boolean_values(tree, data_dict):
    colors = ["red", "green", "blue", "orange", "purple", "cyan", "magenta", "yellow"]  # Define a set of colors for dots
    color_mapping = {}                                                         # Map each characteristic (dictionary) to a color
    available_colors = iter(colors)                                            # Create an iterator for colors

    for characteristic, characteristic_dict in data_dict.items():              # Identify valid dictionaries (containing boolean values) and assign colors
        if all(isinstance(value, int) and value in [0, 1] for value in characteristic_dict.values()):  # Check for boolean values
            if characteristic not in color_mapping:
                try:
                    color_mapping[characteristic] = next(available_colors)     # Assign a unique color
                except StopIteration:
                    raise ValueError("Not enough colors defined for the number of characteristics.")
    
    
    for leaf in tree.get_leaves():                                             # Loop through each leaf in the tree
        column_index = 1                                                       # Start placing dots at column 1 (branch-right)
        for characteristic, characteristic_dict in data_dict.items():          # Check each characteristic dictionary for matches with the leaf name
            if characteristic in color_mapping:                                # Only process valid boolean dictionaries
                for seqid, value in characteristic_dict.items():
                    if seqid in leaf.name and value == 1:                      # Match seqid and check if value is 1
                        dot = CircleFace(radius=10, color=color_mapping[characteristic], style="circle") # Add a colored dot (CircleFace) to the leaf
                        leaf.add_face(dot, column=column_index, position="branch-right")
                        column_index += 1                                      # Move to the next column for the next dot
                            
    ts = TreeStyle()                                                           # Create a TreeStyle object for customizing the visualization
    ts.show_leaf_name = True                                                   # Ensure leaf names are displayed

    for characteristic, color in color_mapping.items():                        # Add a legend for the dots
        
        dot = CircleFace(radius=15, color=color, style="circle")               # Create a legend entry with a colored dot and the characteristic name
        label = TextFace(f" {characteristic}", fsize=16)                       # Formatting the label
        label.margin_left = 50
        ts.legend.add_face(dot, column=0)
        ts.legend.add_face(label, column=1)

    return ts

#Function to add branch support values(bootstraps)
def add_support_values_to_tree(tree, ts):
    for node in tree.traverse():                                               # Traverse all nodes in the tree
        if node.support:                                                       # Check if there's a support value
            support_value = node.support
            if 2 <= support_value <= 100:  
                if support_value >= 70:                                        # Determine color based on the support value
                    color = "black"                                            # Support >= 70 is black
                    support_label = f"{support_value}*"                        # Add * next to values >= 70
                else:
                    color = "red"                                              # Support < 70 is red
                    support_label = f"{support_value}"                         # No * for values < 70
                support_face = TextFace(support_label, fsize=9, fgcolor=color)
                node.add_face(support_face, column=0, position="branch-bottom")# Add the label on the right side of the branch
    return ts  

#Function to create the tree layout with titles and treestyle
def create_tree_with_layout(tree, ts, data_dict, iqtree_command):
    high_interest_dict = data_dict['high_interest']                            # Extract the high_interest dictionary from the data_dict
    genus_dict = data_dict['Genus']                                            # Extract genus dictionary
    gene_dict = data_dict['Gene']                                              # Extract gene dictionary
    
    genus = ", ".join(str(val) for val in genus_dict.values() if pd.notna(val))# Build title and subtitle
    gene = ", ".join(str(val) for val in gene_dict.values() if pd.notna(val))
    title = f"Phylogenetic Tree of {genus} Based on {gene}"
    subtitle = f"IQ-TREE2 parameters: {iqtree_command}"    
    ts.show_scale = True                                                       # Disable the scale bar for this tree.
    ts.show_leaf_name = False                                                  # Disable leaf names from being shown by default.
    ts.margin_left = 100                                                       # Set the left margin for the tree visualization.
    ts.margin_right = 100                                                      # Set the right margin for the tree visualization.
    ts.margin_top = 150                                                        # Set the top margin for the tree visualization.
    ts.margin_bottom = 150                                                     # Set the bottom margin for the tree visualization.
    ts.force_topology = True                                                   # Force the tree topology to remain fixed during visualization.
    ts.title.add_face(TextFace(title, fsize=25, bold=True), column=0)          # Add the title to the tree with specific formatting (font size and bold).
    ts.title.add_face(TextFace(subtitle, fsize=18), column=0)
    ts.scale = 100
    
    for leaf in tree.get_leaves():                                             # Customize leaf names
        seqid = None
        for key in high_interest_dict.keys():                                  # Look for the seqid in the leaf name
            if key in leaf.name:
                seqid = key
                break
        if seqid and high_interest_dict.get(seqid) == 1:                       # Style the leaf name based on the high_interest_dict
            name_face = TextFace(leaf.name, fsize=18, fgcolor="blue", bold=True)
        else:                                                                  # Default styling
            name_face = TextFace(leaf.name, fsize=18, fgcolor="black")
        
        leaf.add_face(name_face, column=0, position="branch-right")
    
    return ts

#Function to save the tree in diff formats
def save_tree(tree, tree_style, output_folder, file_name="phylogenetic_tree"): # define paths for every format compatible to linux and windows
    pdf_path = os.path.join(output_folder, f"{file_name}.pdf")                 #define paths for every format compatible to linux and windows
    jpg_path = os.path.join(output_folder, f"{file_name}.jpg")
    png_path = os.path.join(output_folder, f"{file_name}.png")
    newick_path = os.path.join(output_folder, f"{file_name}.nwk")
    nexus_path = os.path.join(output_folder, f"{file_name}.nex")

    tree.render(pdf_path, tree_style=tree_style)                               # Save as image formats
    tree.render(jpg_path, tree_style=tree_style)
    tree.render(png_path, tree_style=tree_style)
    tree.write(format=1, outfile=newick_path)                                  # Newick format
    tree.write(format=9, outfile=nexus_path)                                   # Nexus format
 
# Function to upload the tree to iTOL
def upload_tree_to_itol(api_key: str, project_name: str, tree_path: str) -> None:     # upload a phylogenetic tree to iTOL.
    itol_instance = Itol()                                                     # Instantiate the iTOL class
    tree_file = Path(tree_path)                                                # Set the tree file to be uploaded
    itol_instance.add_file(tree_file)                                          # Add the tree file to the iTOL instance
    itol_instance.params['APIkey'] = api_key                                   # Set the upload parameters
    itol_instance.params['projectName'] = project_name
    itol_instance.params['treeName'] = 'My Phylogenetic Tree'
    itol_instance.print_variables()                                            # Print the parameters for verification
    
    print("\nUploading the tree. This may take some time depending on the tree size.")     # Upload the tree
    if not itol_instance.upload():                                             # Attempt the upload
        print("There was an error: " + itol_instance.comm.upload_output)       # Print error message
        sys.exit(1)                                                            # Exit the program with an error status
    print("Tree ID: " + itol_instance.comm.tree_id)                            # Print the unique tree ID
    print("iTOL output: " + itol_instance.comm.upload_output)                  # Print the raw server response
    print("Tree Web Page URL: " + itol_instance.get_webpage())                 # Print the URL to view the tree online
    print("Warnings: ", itol_instance.comm.warnings)                           # Print any warnings associated with the upload

# Function to call all the functions
def main(excel_path, show_tree=True):                                                     # Flag to control whether to show the tree: if True then close browser before saving 
    print("Welcome to the Phylogenetic Tree Pipeline ! Please ensure that you provide all required file paths (Excel, FASTA, and tree output paths) correctly.")
 
    print("\n🔍 Prechecking the input data for errors...")
    validate_excel_input(excel_path)  # This function now exits if critical errors are found

    print("\n✅ Input is correct. Starting the calculation...")

 
    data_dict = create_dictionaries_from_excel(excel_path)                     # Step 1: Create dictionaries from the Excel file
    print("\nStep 1: Generating the FASTA file from the data dictionary...")
    fasta_path = generate_fasta_from_dict(data_dict, output_fasta_path)        # Step 2: Generate the FASTA file from the data_dict
    
    print("\nStep 2: Aligning the FASTA sequences using MAFFT...")
    run_mafft(mafft_path, fasta_path, aligned_fasta)                                       # Step 3: Run MAFFT alignment
    
    print("\nStep 3: Calculating the phylogenetic tree using IQ-TREE...")
    treefile_path = calculate_phylogenetic_tree(iqtree2_path, aligned_fasta, output_folder)  # Step 4: Calculate the phylogenetic tree using IQ-TREE

    print("\nStep 4: Tree design process started...")
    tree, ts = get_tree_with_colors(treefile_path, data_dict)                  # Step 5: Design the calculated tree
    tree = apply_node_colors(tree, data_dict)                                  # Apply node colors
    ts = add_boolean_values(tree, data_dict)                                   # Add boolean values (dots) for species with morphology == 1
    ts = add_support_values_to_tree(tree, ts)                                  # Add formatted support values
    ts = create_tree_with_layout(tree, ts, data_dict, iqtree_command)          # Finalize tree layout
    
    if show_tree:
        print("Displaying tree in ETE3 viewer. Please close the viewer to continue...")
        tree.show(tree_style=ts)                                               # Show the tree
        
    print("\nStep 5: Saving the phylogenetic tree in multiple formats...")     # Optional: Save the tree in multiple formats
    save_tree(tree, ts, output_folder)                                         
  
    user_response = input("\nDo you want to upload the tree to iTOL? (yes/no): ").strip().lower()
    if user_response == "yes":       
        api_key = input("Please enter your iTOL API key: ").strip()            # API key and project name are individual
        project_name = input("Please enter the iTOL project name: ").strip()        
        print("\nUploading the tree to iTOL...")
        upload_tree_to_itol(api_key, project_name, treefile_path)                  # Call upload function with user inputs
    elif user_response == "no":
        print("\nTree upload skipped")                                         
                     
    print("\nPipeline completed successfully! Thank you for using the Phylogenetic Tree Pipeline.")
    
# Paths and files
excel_path = get_valid_path("Path to the Excel file: ")
mafft_path = get_valid_path("Path to the MAFFT executable: ")
iqtree2_path = get_valid_path("Path to the IQ-TREE executable: ")
output_folder = get_valid_path("Path to the output folder: ", is_file=False)
output_fasta_path = os.path.join(output_folder, "output_sequences.fasta")      # Construct dynamic file paths based on output folder
aligned_fasta = os.path.join(output_folder, "aligned_sequences.fas")
tree_path = os.path.join(output_folder, "iqtree_output.treefile")
iqtree_command = "-m TESTMERGE -B 1000"

if __name__ == "__main__":
    main(excel_path)
