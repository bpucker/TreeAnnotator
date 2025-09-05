### Laila Rothe & Boas Pucker ###
### ChatGPT was used to generate code snippets and to optimize code ###

__version__ = """v0.1"""

__usage__ = """
	python3 TreeAnnotator.py
	--tree <INPUT_TREE_FILE>
	--info <INPUT_INFO_FILE>
	--out <OUTPUT_FOLDER>
"""

from ete3.treeview import TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
import subprocess, os, sys, csv                                                     # For interacting with system
from ete3 import Tree              # For visuylization directly in python
import matplotlib.colors as mcolors                                            # For validating colors


# --- end of imports --- #

#important for headless servers (without GUI)
os.environ["QT_QPA_PLATFORM"] = "offscreen"

def load_metadata( info_file ):
	"""! @brief load metadata """
	metadata = {}
	with open( info_file ) as f:
		reader = csv.DictReader(f, delimiter='\t')
		for row in reader:
			metadata[row["Gene"]] = row	
	return metadata



def infer_style_mapping(metadata, max_categorical=20):
	"""Infer style mapping for each metadata column."""
	if not metadata:
		return {}
	sample = next(iter(metadata.values()))
	columns = [col for col in sample if col != "Gene" and col != "name"]
	# Build sets of unique values for each column
	unique_values = {col: set() for col in columns}
	for row in metadata.values():
		for col in columns:
			val = row.get(col, "").strip().lower()
			if val:
				unique_values[col].add(val)
	style_mapping = {}
	used_styles = set()
	for col, values in unique_values.items():
		val_set = set(values)
		# Check for binary columns
		is_binary = val_set <= {"0", "1"} or val_set <= {"yes", "no"} or val_set <= {"true", "false"}
		if is_binary:
			if "bold" not in used_styles:
				style_mapping[col] = "bold"
				used_styles.add("bold")
			elif "font_color" not in used_styles:
				style_mapping[col] = "font_color"
				used_styles.add("font_color")
			else:
				continue  # No more binary styles available
		# Categorical with few unique values
		elif len(val_set) <= max_categorical:
			if "background" not in used_styles:
				style_mapping[col] = "background"
				used_styles.add("background")
			elif "dot" not in used_styles:
				style_mapping[col] = "dot"
				used_styles.add("dot")
			elif "font_color" not in used_styles:
				style_mapping[col] = "font_color"
				used_styles.add("font_color")
			else:
				continue  # All styles used
	return style_mapping


def get_style_value( meta, field, style_type, color_maps ):
	val = meta.get(field, "")
	if style_type in ["background", "dot", "font_color"]:
		return color_maps[field].get(val, "#000000")
	elif style_type == "bold":
		return val.lower() in ("1", "yes", "true")
	return None


def build_color_maps(metadata, style_mapping):
	"""! @brief generate color maps for modifying names and backgrounds """
	
	# PREDEFINED_COLORS = [
												# "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4",
												# "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff",
												# "#9a6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075"
											# ]
	PREDEFINED_COLORS = [	"#8B5E3C", "#4B6F44", "#C2B173", "#4A5973", "#A77651", "#71587D",
							"#4A7F94", "#925F86", "#9FAD6A", "#C7A79A", "#5B7A74", "#B3A3C1",
							"#7C6647", "#5C3D3D", "#94A397", "#6D6A3D", "#BBA785", "#434F68"
						]
	color_maps = {}
	for field, style in style_mapping.items():
		if style in ["background", "dot", "font_color"]:
			unique_values = sorted(set(row.get(field, '') for row in metadata.values()))
			color_maps[field] = {
												val: PREDEFINED_COLORS[i % len(PREDEFINED_COLORS)]
												for i, val in enumerate(unique_values)
											}
	return color_maps


def modify_tree( tree, metadata, style_mapping, color_maps ):
	"""! @brief modify the tree based on metadata """
	
	for leaf in tree.iter_leaves():
		meta = metadata.get(leaf.name)
		if not meta:
			print(f"Skipping leaf '{leaf.name}' — not found in metadata")
			continue
		nstyle = NodeStyle()
		# --- background --- #
		for field, style_type in style_mapping.items():
			if style_type == "background":
				nstyle["bgcolor"] = get_style_value(meta, field, style_type, color_maps)
		leaf.set_style(nstyle)

		# --- text bold --- #
		font_color = "#000000"
		bold = False
		for field, style_type in style_mapping.items():
			if style_type == "font_color":
				font_color = get_style_value(meta, field, style_type, color_maps)
			elif style_type == "bold":
				bold = get_style_value(meta, field, style_type, color_maps)
		text_face = TextFace(leaf.name, fsize=10, fgcolor=font_color, bold=bold)
		leaf.add_face(text_face, column=0, position="aligned")
		print(f"Adding face for: {leaf.name}")
		
		# --- dots --- #
		dot_column = 1
		for field, style_type in style_mapping.items():
			if style_type == "dot":
				dot_color = get_style_value(meta, field, style_type, color_maps)
				dot_face = CircleFace(radius=5, color=dot_color, style="circle")
				dot_face.margin_right = 4
				leaf.add_face(dot_face, column=dot_column, position="aligned")
				dot_column += 1
	
	# --- tree style --- #
	ts = TreeStyle()
	ts.show_leaf_name = False  # show names via faces
	ts.mode = "c"  # circular layout (or use "r" for rectangular)
	ts.scale = 120
	
	return tree, ts


def add_support_values_to_tree(tree, ts):
	"""! @brief function to add branch support values(bootstraps)  """
	
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
	return tree, ts  


def save_tree(tree, tree_style, output_folder, file_name="phylogenetic_tree"): # define paths for every format compatible to linux and windows
	"""! @brief save tree in different file formats """
	
	pdf_path = os.path.join(output_folder, f"{file_name}.pdf")                 #define paths for every format compatible to linux and windows
	jpg_path = os.path.join(output_folder, f"{file_name}.jpg")
	png_path = os.path.join(output_folder, f"{file_name}.png")
	newick_path = os.path.join(output_folder, f"{file_name}.nwk")
	nexus_path = os.path.join(output_folder, f"{file_name}.nex")
	
	tree.render( pdf_path, tree_style=tree_style )                               # Save as image formats
	tree.render (jpg_path, tree_style=tree_style )
	tree.render( png_path, tree_style=tree_style )
	tree.write( format=1, outfile=newick_path )                                  # Newick format
	tree.write( format=9, outfile=nexus_path )                                   # Nexus format


def main( arguments ):
	"""! @brief runs the actual analysis """
	
	tree_file = arguments[ arguments.index('--tree')+1 ]
	info_file = arguments[ arguments.index('--info')+1 ]
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	if output_folder[-1] != "/":
		output_folder += "/"
		
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--mafft' in arguments:
		mafft_path = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft_path = "mafft"
	if '--iqtree' in arguments:
		iqtree2_path = arguments[ arguments.index('--iqtree')+1 ]
	else:
		iqtree2_path = "iqtree"

	file_name="phylogenetic_tree"
	
	# --- load metadata from info file --- #
	metadata = load_metadata( info_file )                     # Generate dictionary based on info file
	style_mapping = infer_style_mapping(metadata, max_categorical=20)
	print("Style mapping:", style_mapping)
	color_maps = build_color_maps(metadata, style_mapping)
	
	# --- load tree from file --- #
	tree = Tree( tree_file, format=1 )
	
	# --- modify tree style based on metadata --- #
	tree, ts = modify_tree( tree, metadata, style_mapping, color_maps )
	
	## --- add support values --- #
	#tree, ts = add_support_values_to_tree( tree, ts )
	
	# --- save tree into file --- #
	save_tree( tree, ts, output_folder, file_name )


if '--tree' in sys.argv and '--info' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
