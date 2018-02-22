#!/usr/bin/env python
#coding: utf-8

import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import random
from scipy.spatial import Delaunay
from sklearn.cluster import KMeans
import csv
import xml.etree.cElementTree as ET
import argparse


##########################
### HELP Section START ###
##########################

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
					argparse.ArgumentDefaultsHelpFormatter): pass

parser = argparse.ArgumentParser(prog='Entities_addition_with_Voronoi_K_means',
								formatter_class=CustomFormatter,
								description=
"""

#-------------------------------------------#
# Entities_addition_with_Voronoi_K_means    #
#                                           #
# Author : SOMPAIRAC Nicolas                #
# Contact : nicolas.sompairac@gmail.com     #
# Version : (February 2018); Institut Curie #
#-------------------------------------------#

This script allows to add entities in the vicinity of chosen elements in a
CellDesigner map in SBML format. Position of entities to add are chosen by
breaking the map by Voronoi tesselation and use the K-mean algorithm to chosen
the most relevant coordinates to avoid overlapping.

""")

file_locations = parser.add_argument_group('Location of different files')


file_locations.add_argument('--map', metavar='Map', type=str,
							default=None,
							help=("Location of the CellDesigner XML map on top "
								"of which the entities will be added."))
file_locations.add_argument('--matrix', metavar='RxG_mat', type=str,
							default=None,
							help=("Location of the boolean matrix containing "
								"rows as elements (e.g.reactions) around which "
								"entities will be added and columns as entities"
								" (e.g.proteins) which will be added in the "
								"vicinity of elements. "
								"The matrix has to be in CSV format."))
file_locations.add_argument('--ref-names', metavar='Ref_names', type=str,
							default=None,
							help=("Location of the file containg a list of "
								"elements names in rows. These elements are "
								"the references around which entities will be "
								"added. Their order correspond to the matrix."))
file_locations.add_argument('--add-names', metavar='Add_names', type=str,
							default=None,
							help=("Location of the file containg a list of "
								"entities names in rows. These elements will be"
								" added around needed elements. "
								"Their order correspond to the matrix."))
file_locations.add_argument('--ref-corresp', metavar='Ref_symbols', type=str,
							default=None,
							help=("Location of the file containg two columns "
								"of correspondance between elements to add "
								"original name (first columns) and the name of "
								"needed after the addition (second columns). "
								"For example, a correspondance between "
								"CellDesigner reactions ID and the real name."))
file_locations.add_argument('--add-corresp', metavar='Add_symbols', type=str,
							default=None,
							help=("Location of the file containg two columns "
								"of correspondance between entities to add "
								"original name (first columns) and the name of "
								"needed after the addition (second columns). "
								"For example, a correspondance between "
								"differents IDs such as Entrez and HUGO. "
								"If no correspondance exist, leave a '-' and "
								"a '*' will be appended to the name."))
file_locations.add_argument('--ref-pos', metavar='Ref_positions', type=str,
							default=None,
							help=("Location of the file containg the X and Y "
								"coordinates of the reference elements. First "
								"column with names, Second column with X "
								"coordinates and Third columns with Y "
								"coordinates as in the map."))
file_locations.add_argument('--result-file', metavar='Res_outfile', type=str,
							default=None,
							help=("Location and name of the resulting file "
								"containing the coordinates of the entities "
								"to be added in the vicinity needed elements."))


options = parser.add_argument_group('Options to use')

options.add_argument('-p', '--plot', action='store_true',
					help=("Use if wanting to visualize the addition result. "
						"Reference elements will be marked in Green and added "
						"entities will be marked in Red."))
options.add_argument('--height', metavar='Height', type=float, default=5.0,
							help="Height of the added entitities")
options.add_argument('--width', metavar='Width', type=float, default=5.0,
							help="Width of the added entitities")

args = parser.parse_args()


#########################
### HELP Section STOP ###
#########################


def Parse_XML_file(filename):

	tree = ET.parse(filename)
	root = tree.getroot()

	alias_dict = {}
	for species in root.iter('{http://www.sbml.org/2001/ns/celldesigner}speciesAlias'):
		ID = species.get('id')
		for bounds in species.iter('{http://www.sbml.org/2001/ns/celldesigner}bounds'):
			alias_dict[ID] = (float(bounds.get('x')), float(bounds.get('y')))

	return alias_dict


def Parse_reaction_positions(filename):

	react_dict = {}
	react_id_list = []

	with open(filename, "rU") as infile:

		for line in infile:

			tmp = line.split()

			react_dict[tmp[0]] = (float(tmp[1]), float(tmp[2]))
			react_id_list.append(tmp[0])

	return react_dict, react_id_list


def Add_reaction_to_species(spec_dict, react_dict):

	for react in react_dict:

		spec_dict[react] = react_dict[react]

	return spec_dict


def Get_csv_into_matrix(filename):

	infile = csv.reader(open(filename, "rb"), delimiter=",")
	infile_mat = list(infile)
	matrix = np.array(infile_mat).astype("int")

	return matrix


def Get_list_of_data(filename):

	data_list = []

	with open(filename, "rU") as infile:

		for line in infile:

			data_list.append(line[:-1])

	return np.array(data_list)


def Simplify_genes_symbols(genes_list):

	simpler_list = []

	for item in genes_list:

		simpler_list.append(item.split(".")[0])

	return np.array(simpler_list)


def Get_dict_of_corres_genes(filename):

	data_dict = {}

	with open(filename, "rU") as infile:

		for line in infile:

			tmp = line.split()

			if tmp[1] != "-":

				data_dict[tmp[0]] = tmp[1]

			else:

				data_dict[tmp[0]] = tmp[0]+"*"

	return data_dict


def Get_reactions_corresp_genes(mat, react_list, gene_list, corresp_dict):

	react_genes_dict = {} # Key = React_ID | Value = List of implicated genes
	react_index = 0

	for react in mat:

		if 1 in react:

			genes_index = np.where(react==1)

			genes_entrez_list = list(gene_list[genes_index])

			genes_symbol_list = []

			for gene in genes_entrez_list:

				if corresp_dict[gene] not in genes_symbol_list:

					genes_symbol_list.append(corresp_dict[gene])

			react_genes_dict[react_list[react_index]] = genes_symbol_list

		else:

			react_genes_dict[react_list[react_index]] = []

		react_index += 1

	return react_genes_dict


def Extract_existing_only_react(react_genes_dict, react_id_list, stupid_react_dict):

	existing_react_dict = {}

	for react in react_id_list:

		if react in react_genes_dict:

			existing_react_dict[react] = react_genes_dict[react]

		elif stupid_react_dict[react] in react_genes_dict:

			existing_react_dict[react] = react_genes_dict[stupid_react_dict[react]]

		else:

			existing_react_dict[react] = []

	return existing_react_dict


################################################################################


def Get_voronoi(pos_dict, react_id_list):

	Points_list = []
	ID_list = []

	for ID in Species_pos_dict:

		Points_list.append(Species_pos_dict[ID])
		ID_list.append(ID)

	vor = Voronoi(Points_list)

	Voronoi_dict = {} # Key = ID | Value = List of Voronoi vertices coordinates

	i = 0

	for region in vor.point_region:

		Voronoi_dict[ID_list[i]] = []
		
		for vertice in vor.regions[region]:

			Voronoi_dict[ID_list[i]].append(vor.vertices[vertice])

		i+=1

	Voronoi_reactions_dict = {}

	for key in Voronoi_dict:

		if key in react_id_list:

			Voronoi_reactions_dict[key] = Voronoi_dict[key]

	return Voronoi_reactions_dict


def Correct_voronoi_coordinates(spec_pos_dict, vor_react_pos_dict):

	x_list = []
	y_list = []

	for spec in spec_pos_dict:

		x_list.append(spec_pos_dict[spec][0])
		y_list.append(spec_pos_dict[spec][1])

	x_min = min(x_list)
	x_max = max(x_list)
	y_min = min(y_list)
	y_max = max(y_list)

	x_range = x_max/x_min
	y_range = y_max/y_min

	for react in vor_react_pos_dict:

		for vor in vor_react_pos_dict[react]:

			if vor[0] < x_min:

				vor[0] = x_min - x_range/2

			elif vor[0] > x_max:

				vor[0] = x_max + x_range/2

			if vor[1] < y_min:

				vor[1] = y_min - y_range/2

			elif vor[1] > y_max:

				vor[1] = y_max + y_range/2

	return vor_react_pos_dict


def Test_in_shape(pt_list, shape):
	"""
    Test if points in `pt_list` are in `shape`

    `pt_list` should be a `NxK` coordinates of `N` points in `K` dimensions
    `shape` is scipy.spatial.Delaunay object from a `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """

	return shape.find_simplex(pt_list)>=0


def Generate_N_points_in_shape(shape_list, N):

	X_list = [item[0] for item in shape_list]
	Y_list = [item[1] for item in shape_list]

	X_min = min(X_list)
	X_max = max(X_list)
	Y_min = min(Y_list)
	Y_max = max(Y_list)

	N_points_list = []

	shape_delaunay = Delaunay(shape_list)

	while(len(N_points_list)<N):

		test_pt = (random.uniform(X_min, X_max), random.uniform(Y_min, Y_max))

		if Test_in_shape(test_pt, shape_delaunay):

			N_points_list.append(test_pt)

	return np.array(N_points_list)


def Clust_points_with_K_means(pt_list, K):

	kmeans = KMeans(n_clusters=K, random_state=0).fit(pt_list)

	return kmeans.cluster_centers_


def Find_protein_markers(react_voronoi_dict, react_with_genes_dict):

	protein_markers_dict = {}
	# Key = associated Reaction ID
	# Value = list containing 2 lists:
	# Value_list_1 = list of protein/gene names
	# Value_list_2 = list of corresponding proteins/genes coordinates

	total = len(react_voronoi_dict)
	i = 0

	N_points = 100

	for ID in react_voronoi_dict:

		K = len(react_with_genes_dict[ID])

		if K != 0:

			points_inside_voronoi_list = Generate_N_points_in_shape(react_voronoi_dict[ID], N_points)

			centroids_list = Clust_points_with_K_means(points_inside_voronoi_list, K)

			protein_markers_dict[ID] = (react_with_genes_dict[ID], centroids_list)

		i+=1
		print i, "/", total

	return protein_markers_dict


def Plot_result_of_algo_single_shape(shape_list, pt_list, cent_list):

	X_list = [item[0] for item in shape_list]
	Y_list = [item[1] for item in shape_list]
	X_list.append(shape_list[0][0])
	Y_list.append(shape_list[0][1])

	X_points = [item[0] for item in pt_list]
	Y_points = [item[1] for item in pt_list]

	X_centroids = [item[0] for item in cent_list]
	Y_centroids = [item[1] for item in cent_list]

	plt.plot(X_list, Y_list, ms=2)
	plt.plot(X_points, Y_points, "r.", ms=2)
	plt.plot(X_centroids, Y_centroids, "k.", ms=3)
	plt.show()

	return


def Plot_final_result(spec_pos_dict, prot_pos_dict, react_id_list):

	X_spec_list = []
	Y_spec_list = []

	X_react_list = []
	Y_react_list = []

	X_prot_list = []
	Y_prot_list = []

	for spec in spec_pos_dict:

		if spec in react_id_list:

			X_react_list.append(spec_pos_dict[spec][0])
			Y_react_list.append(spec_pos_dict[spec][1])

		else:

			X_spec_list.append(spec_pos_dict[spec][0])
			Y_spec_list.append(spec_pos_dict[spec][1])

	for prot in prot_pos_dict:

		for mark in prot_pos_dict[prot][1]:

			X_prot_list.append(mark[0])
			Y_prot_list.append(mark[1])

	
	plt.plot(X_spec_list, Y_spec_list, "g.", ms=1)
	plt.plot(X_react_list, Y_react_list, "k.", ms=1)
	plt.plot(X_prot_list, Y_prot_list, "r.", ms=1)
	plt.gca().invert_yaxis()
	plt.axis("off")
	plt.show()
	
	return


def Write_proteins_in_binom_reaction_format(prot_pos_dict, filename, width, height):

	occurencies_dict = {}

	with open(filename, "w") as outfile:

		for react in prot_pos_dict:

			if len(prot_pos_dict[react][0]) != 0:

				for pos in xrange(len(prot_pos_dict[react][0])):

					if prot_pos_dict[react][0][pos] not in occurencies_dict:

						outfile.write(prot_pos_dict[react][0][pos])
						outfile.write("\t")
						outfile.write("X:"+str(prot_pos_dict[react][1][pos][0]))
						outfile.write(";Y:"+str(prot_pos_dict[react][1][pos][1]))
						outfile.write(";W:"+str(width)+";H:"+str(height)+"\n")

						occurencies_dict[prot_pos_dict[react][0][pos]] = 1

					else:

						outfile.write(prot_pos_dict[react][0][pos]+"'"*occurencies_dict[prot_pos_dict[react][0][pos]])
						outfile.write("\t")
						outfile.write("X:"+str(prot_pos_dict[react][1][pos][0]))
						outfile.write(";Y:"+str(prot_pos_dict[react][1][pos][1]))
						outfile.write(";W:"+str(width)+";"+str(height)+"\n")

						occurencies_dict[prot_pos_dict[react][0][pos]] = occurencies_dict[prot_pos_dict[react][0][pos]] + 1

	return


print "\n"
############
### MAIN ###
############


################################################################################
# Taking info from files
################################################################################

Reaction_Genes_matrix = Get_csv_into_matrix(args.matrix)

Reactions_list = Get_list_of_data(args.ref_names)
Genes_list = Get_list_of_data(args.add_names)

Genes_list = Simplify_genes_symbols(Genes_list)

Correspondance_genes_dict = Get_dict_of_corres_genes(args.add_corresp)

Stupid_react_ID_Name_corresp_dict = Get_dict_of_corres_genes(args.ref_corresp)

React_with_genes_dict = Get_reactions_corresp_genes(Reaction_Genes_matrix, 
	Reactions_list, Genes_list, Correspondance_genes_dict)

Species_pos_dict = Parse_XML_file(args.map)
Reaction_pos_dict, Reactions_ID_list = Parse_reaction_positions(args.ref_pos)
Species_pos_dict = Add_reaction_to_species(Species_pos_dict, Reaction_pos_dict)

Existing_React_with_genes_dict = Extract_existing_only_react(React_with_genes_dict, 
	Reactions_ID_list, Stupid_react_ID_Name_corresp_dict)


################################################################################
# Actual Algorithm
################################################################################

Voronoi_reactions_pos_dict = Get_voronoi(Species_pos_dict, Reactions_ID_list)

Voronoi_reactions_pos_dict = Correct_voronoi_coordinates(Species_pos_dict, 
	Voronoi_reactions_pos_dict)

Protein_markers_dict = Find_protein_markers(Voronoi_reactions_pos_dict, 
	Existing_React_with_genes_dict)

if args.plot:
	print "\nPlotting the result\n"
	Plot_final_result(Species_pos_dict, Protein_markers_dict, Reactions_ID_list)

Write_proteins_in_binom_reaction_format(Protein_markers_dict, args.result_file,
	args.width, args.height)

print "\n"