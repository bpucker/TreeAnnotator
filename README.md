# TreeAnnotator (interactive)
Python script to map functional annotation to a phylogenetic tree

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

# TreeAnnotator for CLI
This script allows the automatic coloration of a tree file based on metadata input in a TSV file. Different output formats are supported (PNG, JPEG, PDF).

## Installing dependencies

```
pip3 install ete3
pip3 install PyQt5
pip3 install matplotlib
```

OR (if externally managed message appears):

```
sudo apt install python3-ete3
sudo apt install python3-matplotlib
sudo apt install python3-pyqt5
```




## Usage

```
Usage:
  python3 TreeAnnotator.py --info <INFO_FILE> --tree <TREE_FILE> --out <DIR>

Mandatory:
  --info     STR      Metadata input file
  --tree     STR      Tree file
  --out      STR      Output directory

Optional:
	--name     STR      Tree name
	--fontsize INT      Font size of leaf labels
	--dotsize  INT      Dot size
	--bgcolors STR      Comma-separated list of colors (in quotation marks)
	--dotcolors STR      Comma-separated list of colors (in quotation marks)

```

`--info` full path to the file containing the metadata.

`--tree` full path to the tree input file.

`--out` full path to the output folder. The folder will be created if it does not exist already.

`--name` specifies the name prefix of the output files.

`--fontsize` specifies the font size of leaf labels. Default: 10.

`--dotsize` specifies the size of dots used to visualize metadata next to leaf labels. Default: 5.

`--bgcolors` specifies the background colors of leaf labels. A list of comma-separated color strings is expected. Quotation marks around the entire string are required to avoid hashtags breaking the inserted command. Default: "#8B5E3C,#4B6F44,#C2B173,#4A5973,#A77651,#71587D,#4A7F94,#925F86,#9FAD6A,#C7A79A,#5B7A74,#B3A3C1,#7C6647,#5C3D3D,#94A397,#6D6A3D,#BBA785,#434F68".

`--dotcolors` specifies the dot colors of leaf labels. A list of comma-separated color strings is expected. Quotation marks around the entire string are required to avoid hashtags breaking the inserted command. Default: "#e6194b,#3cb44b,#ffe119,#4363d8,#f58231,#911eb4,#42d4f4,#f032e6,#bfef45,#fabebe,#469990,#e6beff,#9a6324,#800000,#aaffc3,#808000,#ffd8b1,#000075"



## Example:

```
python3 TreeAnnotator.py \
--info metadata.tsv \
--tree newick_tree.treefile \
--out ./stylish_tree/
```



## References
- Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Mol Biol Evol. 2016 Jun;33(6):1635-8. doi: [10.1093/molbev/msw046](https://doi.org/10.1093/molbev/msw046). Epub 2016 Feb 26. PMID: 26921390; PMCID: PMC4868116.
- J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.
- ChatGPT v4

