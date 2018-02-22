#!/usr/bin/env python
#coding: utf-8

import xml.etree.cElementTree as ET
import re
import argparse


##########################
### HELP Section START ###
##########################

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
					argparse.ArgumentDefaultsHelpFormatter): pass

parser = argparse.ArgumentParser(prog='BRF_to_CellDesigner',
								formatter_class=CustomFormatter,
								description=
"""

#-------------------------------------------#
# BRF_to_CellDesigner                       #
#                                           #
# Author : SOMPAIRAC Nicolas                #
# Contact : nicolas.sompairac@gmail.com     #
# Version : (February 2018); Institut Curie #
#-------------------------------------------#

This script allows to transform a set of entities in BiNoM Reaction Format (BRF)
with entities names, X and Y coordinates as well as their Width and Height to a 
CellDesigner XML file. The generated CellDesigner map take into account aliases.
The BRF file can be obtained through the script
Entities_addition_with_Voronoi_K_means.py

""")

file_locations = parser.add_argument_group('Location of different files')


file_locations.add_argument('--brf', metavar='BRF_map', type=str,
							default=None,
							help=("Location of the BRF file to convert to "
								"CellDesigner XML format."))
file_locations.add_argument('--result-file', metavar='Res_outfile', type=str,
							default=None,
							help=("Location and name of the resulting file "
								"in CellDesigner XML format."))


options = parser.add_argument_group('Options to use')

options.add_argument('--map-name', metavar='Map_name', type=str,
					help=("Name of the CellDesigner map to put as info in the "
						"XML result file."))
options.add_argument('--map-width', metavar='Map_width', type=str,
					help=("Width of the resulting CellDesigner map."))
options.add_argument('--map-height', metavar='Map_height', type=str,
					help=("Height of the resulting CellDesigner map."))

args = parser.parse_args()


#########################
### HELP Section STOP ###
#########################


def Create_protein_annotation(elem, alias, name, x, y, size):

	celldesigner_speciesAlias = ET.SubElement(elem, "celldesigner:speciesAlias",id=alias,species=name)
	activity = ET.SubElement(celldesigner_speciesAlias, "celldesigner:activity")
	activity.text = "inactive"
	bounds = ET.SubElement(celldesigner_speciesAlias, "celldesigner:bounds",h=str(size),w=str(size),x=str(x),y=str(y))
	font = ET.SubElement(celldesigner_speciesAlias, "celldesigner:font",size="8")
	view = ET.SubElement(celldesigner_speciesAlias, "celldesigner:view",state="usual")
	usualView = ET.SubElement(celldesigner_speciesAlias, "celldesigner:usualView")

	innerPosition = ET.SubElement(usualView, "celldesigner:innerPosition",x="0.0",y="0.0")
	boxSize = ET.SubElement(usualView, "celldesigner:boxSize",height="10.0",width="10.0")
	singleLine = ET.SubElement(usualView, "celldesigner:singleLine",width="1.0")
	paint = ET.SubElement(usualView, "celldesigner:paint",color="ffccffcc",scheme="Color")

	briefView = ET.SubElement(celldesigner_speciesAlias, "celldesigner:briefView")

	innerPosition = ET.SubElement(briefView, "celldesigner:innerPosition",x="0.0",y="0.0")
	boxSize = ET.SubElement(briefView, "celldesigner:boxSize",height="60.0",width="80.0")
	singleLine = ET.SubElement(briefView, "celldesigner:singleLine",width="1.0")
	paint = ET.SubElement(briefView, "celldesigner:paint",color="3fff0000",scheme="Color")

	info = ET.SubElement(celldesigner_speciesAlias, "celldesigner:info",angle="-1.5707963267948966",state="empty")

	return

def Create_species_annotation(elem, ID1, ID2, name, cd_type):

	species = ET.SubElement(elem, "species",compartment="default",id=ID1,initialAmount="0",metaid=ID1,name=name)

	annotation = ET.SubElement(species, "annotation")

	extension = ET.SubElement(annotation, "celldesigner:extension")
	positionToCompartment = ET.SubElement(extension, "celldesigner:positionToCompartment")
	positionToCompartment.text = "inside"
	speciesIdentity = ET.SubElement(extension, "celldesigner:speciesIdentity")
	clas = ET.SubElement(speciesIdentity, "celldesigner:class")
	clas.text = "PROTEIN"
	proteinReference = ET.SubElement(speciesIdentity, "celldesigner:proteinReference")
	proteinReference.text = ID2

	return


def Construct_XML_file(prot_dict, filename, map_name, map_width, map_height):

	sbml = ET.Element('sbml',level="2",version="4",xmlns="http://www.sbml.org/sbml/level2/version4")
	sbml.set('xmlns:celldesigner',"http://www.sbml.org/2001/ns/celldesigner")

	model = ET.SubElement(sbml, 'model',id=map_name,metaid=map_name)

	annotation = ET.SubElement(model, 'annotation')

	celldesigner_extension = ET.SubElement(annotation, 'celldesigner:extension')

	celldesigner_modelVersion = ET.SubElement(celldesigner_extension, 'celldesigner:modelVersion')
	celldesigner_modelVersion.text = "4.0"
	celldesigner_modelDisplay = ET.SubElement(celldesigner_extension, 'celldesigner:modelDisplay', sizeX=map_width, sizeY=map_height)
	celldesigner_listOfCompartmentAliases = ET.SubElement(celldesigner_extension, 'celldesigner:listOfCompartmentAliases')
	celldesigner_listOfComplexSpeciesAliases = ET.SubElement(celldesigner_extension, 'celldesigner:listOfComplexSpeciesAliases')

	celldesigner_listOfSpeciesAliases = ET.SubElement(celldesigner_extension, 'celldesigner:listOfSpeciesAliases')

	sa = 1
	s = 1

	h_w = 10.0

	for prot in prot_dict:

		for prot_alias in prot_dict[prot]:

			Create_protein_annotation(celldesigner_listOfSpeciesAliases, "sa"+str(sa), "s"+str(s), prot_alias[1], prot_alias[2], h_w)

			sa += 1

		s += 1

	celldesigner_listOfGroups = ET.SubElement(celldesigner_extension, 'celldesigner:listOfGroups')
	celldesigner_listOfProteins = ET.SubElement(celldesigner_extension, 'celldesigner:listOfProteins')

	pr = 1

	for prot in prot_dict:

		celldesigner_protein = ET.SubElement(celldesigner_listOfProteins, "celldesigner:protein",id="pr"+str(pr),name=prot,type="GENERIC")

		pr += 1

	celldesigner_listOfGenes = ET.SubElement(celldesigner_extension, 'celldesigner:listOfGenes')
	celldesigner_listOfRNAs = ET.SubElement(celldesigner_extension, 'celldesigner:listOfRNAs')
	celldesigner_listOfAntisenseRNAs = ET.SubElement(celldesigner_extension, 'celldesigner:listOfAntisenseRNAs')
	celldesigner_listOfLayers = ET.SubElement(celldesigner_extension, 'celldesigner:listOfLayers')
	celldesigner_listOfBlockDiagrams = ET.SubElement(celldesigner_extension, 'celldesigner:listOfBlockDiagrams')

	listOfUnitDefinitions = ET.SubElement(model, 'listOfUnitDefinitions')

	unitDefinition = ET.SubElement(listOfUnitDefinitions, 'unitDefinition',id="substance",metaid="substance",name="substance")
	listOfUnits = ET.SubElement(unitDefinition, 'listOfUnits')
	unit = ET.SubElement(listOfUnits, 'unit',kind="mole",metaid="CDMT00001")

	unitDefinition = ET.SubElement(listOfUnitDefinitions, 'unitDefinition',id="volume",metaid="volume",name="volume")
	listOfUnits = ET.SubElement(unitDefinition, 'listOfUnits')
	unit = ET.SubElement(listOfUnits, 'unit',kind="litre",metaid="CDMT00002")

	unitDefinition = ET.SubElement(listOfUnitDefinitions, 'unitDefinition',id="area",metaid="area",name="area")
	listOfUnits = ET.SubElement(unitDefinition, 'listOfUnits')
	unit = ET.SubElement(listOfUnits, 'unit',exponent="2",kind="metre",metaid="CDMT00003")

	unitDefinition = ET.SubElement(listOfUnitDefinitions, 'unitDefinition',id="length",metaid="length",name="length")
	listOfUnits = ET.SubElement(unitDefinition, 'listOfUnits')
	unit = ET.SubElement(listOfUnits, 'unit',kind="metre",metaid="CDMT00004")

	unitDefinition = ET.SubElement(listOfUnitDefinitions, 'unitDefinition',id="time",metaid="time",name="time")
	listOfUnits = ET.SubElement(unitDefinition, 'listOfUnits')
	unit = ET.SubElement(listOfUnits, 'unit',kind="second",metaid="CDMT00005")

	listOfCompartments = ET.SubElement(model, 'listOfCompartments')

	compartment = ET.SubElement(listOfCompartments, 'compartment',id="default",metaid="default",size="1",units="volume")

	listOfSpecies = ET.SubElement(model, 'listOfSpecies')

	pr = 1

	for prot in prot_dict:

		Create_species_annotation(listOfSpecies, "s"+str(pr), "pr"+str(pr), prot, cd_type)

		pr += 1

	# Writting of the XML structure to a file
	tree = ET.ElementTree(sbml)
	tree.write(filename, xml_declaration=True)

	with open(filename, 'rU') as infile:

		oneline_xml = infile.read()

	with open(filename, 'w') as outfile:

		outfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
		out = oneline_xml.replace("><", ">\n<")
		outfile.write(out)

	return


def Read_BRF(filename):

	prot_dict = {}
	# Key = Prot_name
	# Value = Tuple of (Alias, X, Y)

	with open(filename, "rU") as infile:

		for line in infile:

			tmp = re.split("\t|;|:",line)

			if tmp[0].replace("'","") in prot_dict:

				prot_dict[tmp[0].replace("'","")].append((tmp[0], tmp[2], tmp[4]))

			else:

				prot_dict[tmp[0].replace("'","")] = []
				prot_dict[tmp[0].replace("'","")].append((tmp[0], tmp[2], tmp[4]))

	return prot_dict


print "\n"
############
### MAIN ###
############

Prot_dict = Read_BRF(args.brf)

Construct_XML_file(Prot_dict, args.result_file, args.map_name, args.map_width, args.map_height)

print "\n"