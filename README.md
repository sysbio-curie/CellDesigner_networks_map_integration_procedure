# CellDesigner_networks_map_integration_procedure
Repository containing a procedure and [scripts](scripts) allowing to add elements on top of [CellDesigner](http://www.celldesigner.org/) maps and integrate them without overlaps with existing entities. It makes usage of Java apps ([BiNoM](https://binom.curie.fr/)) and python [scripts](scripts), included in the folder [scripts](scripts).
This procedure has been used to generate a Metabolic map by adding proteins on top of catalysed reactions. The map used for this work is [ReconMap 2.0](https://www.nature.com/articles/nbt.2488).
Even though this procedure describes how to add catalysing proteins in the vicinity of metabolic reactions, this procedure can be generalised to more other similar tasks.
The procedure consists in several steps:
1) Take an XML map in [CellDesigner](http://www.celldesigner.org/) SBML format.
2) Retrieve reaction glyphs coordinates from the XML map through the [BiNoM](https://binom.curie.fr/) function.
```bash
java -cp binom.jar fr.curie.BiNoM.pathways.navicell.FindReactionPositionsScript Map_file.xml
```
3) Generate multiple file contaning the information on the reactions and proteins to add:
- List of reaction names
- List of protein names
- Boolean matrix with reaction in rows and proteins in columns. The value is '1' if a protein is part of the catalysis of the reaction.
- Correspondance between reaction ids and names in the XML map
- Correspondance between proteins IDs and needed symbols (e.g. Entrez to HUGO symbols)
- List of reaction positions generate through the [BiNoM](https://binom.curie.fr/) function
4) Apply the script `Entities_addition_with_Voronoi_K_means.py` to generate a file in BiNoM Reaction Format (BRF) containing protein coordinates as well as their width and height to be displayed on the map. For information about the needed inputs, use the `--help` option.
```bash
python Entities_addition_with_Voronoi_K_means.py --help
```
5) Apply the script `BRF_to_CellDesigner.py` on the generated file from step 4. to convert the BRF format to CellDesigner map format. For information about the needed inputs, use the `--help` option.
```bash
python BRF_to_CellDesigner.py --help
```
6) Use the [BiNoM](https://binom.curie.fr/) function to merge the original map and the map of proteins together in one unique map. This is possible either from Cytoscape or through command lines. For command lines, a config file has to be created with `width_size` and `height_size` corresponding to the original map. The config file should be as followed:
```
mapsize width_size height_size
First_map_file.xml 1 heigth_size
Second_map_file.xml 1 heigth_size
```
The command for the merging is the following:
```bash
java -cp binom.jar fr.curie.BiNoM.pathways.utils.MergingMapsProcessor --config config_filename --prefixlength 1 --out Merged_filename.xml --mergemaps --mergespecies --verbose
```
