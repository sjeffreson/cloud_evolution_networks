import h5py
import numpy as np
import glob, os, re, sys
from typing import List, Dict
from multiprocessing import Process
import pickle

import networkx as nx
import astro_helper as ah

DEFAULT_INPUT_DIR = "/n/holystore01/LABS/itc_lab/Users/sjeffreson/NGC300/"
DEFAULT_OUTPUT_DIR = "/n/holystore01/LABS/itc_lab/Users/sjeffreson/NGC300/MC-iterations/"
class MergertreeProps:
	'''Analysis of cloud evolution network properties. Note that the particles
	and mask attributes are for adding properties to cloud arrays after the fact,
	and so are not analyzed as cloud properties.'''
	def __init__(
		self,
		DiGraph: str,
		min_resolved_time: int = 0,
		min_radius: float = 0.,
		num_mc_iter: int = None,
		num_mc_workers: int = None,
    ):
		with open(DEFAULT_INPUT_DIR+'/'+DiGraph, 'rb') as f:
			self.DiGraph = pickle.load(f)
		first_node = next(iter(self.DiGraph.nodes()))
		self.all_attr_keys = self.DiGraph.nodes[first_node].keys() - set(['particles', 'mask'])
		self.timestep = np.mean(np.diff(sorted(np.unique(list(nx.get_node_attributes(self.DiGraph, 'time').values())))))
		xs_cent, ys_cent, zs_cent = zip(*list(nx.get_node_attributes(self.DiGraph, 'centroid').values()))
		self.width = 2.*np.max(np.sqrt(np.array(xs_cent)**2 + np.array(ys_cent)**2))
		self.min_resolved_time = min_resolved_time
		self.min_radius = min_radius
		self.num_mc_iter = num_mc_iter
		self.num_mc_workers = num_mc_workers

	def get_timestep(self) -> float:
		'''return the timestep of the cloud evolution network'''
		return self.timestep

	def get_width(self) -> float:
		'''return the width of the cloud evolution network'''
		return self.width

	def get_cloud_evol_network(self) -> nx.DiGraph:
		'''return the full cloud evolution network, cut to specifications
		if cut_wcs has been used'''
		return self.DiGraph

	def cut_wcs(self, time_cut: bool=True, radius_cut: bool=True, time_rsln_cut: bool=True):
		'''remove weakly-connected components of the directed graph that touch the
		minimum/maximum time, as these may be incomplete'''

		wcs = nx.weakly_connected_components(self.DiGraph)

		time_dict = nx.get_node_attributes(self.DiGraph, 'time')
		min_time, max_time = np.min(list(time_dict.values())), np.max(list(time_dict.values()))

		wcs_cmplt, wcs_cmplt_lens = [], []
		for wc in wcs:
			G_wc = self.DiGraph.subgraph(wc)

			cnd = True
			if time_cut:
				times_wc = list(nx.get_node_attributes(G_wc, 'time').values())
				min_time_wc, max_time_wc = np.min(times_wc), np.max(times_wc)
				cnd = cnd & (min_time_wc > min_time) & (max_time_wc < max_time)
			if radius_cut:
				xs_cent_wc, ys_cent_wc, zs_cent_wc = zip(*list(nx.get_node_attributes(G_wc, 'centroid').values()))
				xs_cent_wc = np.array(list(xs_cent_wc))
				ys_cent_wc = np.array(list(ys_cent_wc))
				Rs_cent_wc = np.sqrt(xs_cent_wc**2 + ys_cent_wc**2)
				min_radius_wc, max_radius_wc = np.min(Rs_cent_wc), np.max(Rs_cent_wc)
				cnd = cnd & (max_radius_wc < self.width/2.) & (min_radius_wc > self.min_radius)
			if time_rsln_cut:
				times_wc = list(nx.get_node_attributes(G_wc, 'time').values())
				min_time_wc, max_time_wc = np.min(times_wc), np.max(times_wc)
				cnd = cnd & (max_time_wc - min_time_wc > self.min_resolved_time)
	
			if cnd:
				wcs_cmplt.append(list(wc))
				wcs_cmplt_lens.append(len(list(wc)))

		G_cmplt = self.DiGraph.subgraph(set(ah.flatten_list(wcs_cmplt)))
		self.DiGraph = G_cmplt

	def get_cloud_evol_wcs(self, extensive_attr_keys: List[str]) -> Dict[str, np.array]:
		'''Get the cloud evolution of a given attribute for clouds defined
		by the weakly-connected components of the graph. Specify the extensive
		attributes that are to be summed over the weakly-connected components,
		rather than averaged'''

		intensive_attr_keys = list(set(self.all_attr_keys) - set(extensive_attr_keys))
		attrs_wcs = {attr: [] for attr in self.all_attr_keys}

		wcs = nx.weakly_connected_components(self.DiGraph)
		for wc in wcs:
			G_wc = self.DiGraph.subgraph(wc)
			times_wc = np.unique(list(nx.get_node_attributes(G_wc, 'time').values()))
			nodes_wc = list(G_wc.nodes())
			uniquetimes, idcs_uniquetimes = np.unique(times_wc, return_index=True)

			extensive_attrs_wc = {attr: np.array(list(nx.get_node_attributes(G_wc, attr).values())) for attr in extensive_attr_keys}
			intensive_attrs_wc = {attr: np.array(list(nx.get_node_attributes(G_wc, attr).values())) for attr in intensive_attr_keys}
			print(np.shape(intensive_attrs_wc['centroid']))

			# make sure the attributes are ordered by the node key
			extensive_attrs_wc = {attr: np.array([extensive_attrs_wc[attr][nodes_wc.index(node)] for node in nodes_wc]) for attr in extensive_attr_keys}
			intensive_attrs_wc = {attr: np.array([intensive_attrs_wc[attr][nodes_wc.index(node)] for node in nodes_wc]) for attr in intensive_attr_keys}

			# sum over extensive, average over intensive
			attrs_wc = {}
			for attr in extensive_attr_keys:
				attrs_wc[attr] = np.array([np.sum(extensive_attrs_wc[attr][times_wc==time]) for time in uniquetimes])
			for attr in intensive_attr_keys:
				attrs_wc[attr] = np.array([
					np.average(intensive_attrs_wc[attr][times_wc==time],
					weights=extensive_attrs_wc['mass'][times_wc==time], axis=0) for time in uniquetimes
				])

			for attr in self.all_attr_keys:
				attrs_wcs[attr].append(attrs_wc[attr])

		return {attr: np.array(attrs_wcs[attr]) for attr in all_attr_keys}

	def set_mc_rand_nos(self, mc_no: int) -> Dict[str, float]:
		'''get a random number for each node in the graph, with a given random seed'''
		np.random.seed(mc_no)
		r = np.random.uniform(size=len(self.nodes))
		return {node: r for node, r in zip(self.nodes, r)}

	def save_cloud_evol_mc(self):
		'''Save the cloud evolution along Monte Carlo trajectories to disk.
		Saved to disk so the same computation with the same seeds can be used
		for different analyses.'''
		if self.num_mc_iter is None:
			raise ValueError("instance must specify num_mc_iter and num_mc_workers to use this function")

		self.nodes = self.DiGraph.nodes()
		self.formnodes = [
			node for node in self.nodes
			if self.Digraph.in_degree(node) < self.DiGraph.out_degree(node)
		] # formation nodes are those with more children than parents

		mc_iter = 0
		while(mc_iter < self.num_mc_iter):
			for j in range(mc_iter, min([self.num_mc_iter, i + self.num_mc_workers])):
				print("thread "+str(j)+" started")
				p = Process(target=self.mc_iteration, args=(j))
				procs.append(p)
				p.start()

			# block until all the threads finish (i.e. block until all function_x calls finish)    
			for t in procs:
				t.join()

			mc_iter += self.num_mc_workers

	def mc_iteration(self, mc_no: int):
		'''function to do a single MC iteration for all clouds/trajectories
		of the network graph'''

		mc_rand_nos_dict = self.set_mc_rand_nos(mc_no)
		iterations_dict = {node: 0 for node in self.nodes}
		attrs_dict = {attr: [] for attr in self.all_attr_keys}
		lifetimes = []

		edges_visited = []
		formation_events = 0

		for formnode in self.formnodes: # loop through formation nodes of the graph
			N_form = self.DiGraph.out_degree(formnode) - self.DiGraph.in_degree(formnode)
			for j in range(N_form):
				lifetime = 0.
				attrs_dict[attr].append([self.DiGraph.nodes[formnode][attr]])
				formation_events+=1 # book-keeping
				lifetime = walk_trajectory(formnode)
				lifetimes.append(lifetime)

		'''function to walk one trajectory through the graph'''
		def walk_trajectory(node: int) -> int:
			children = self.DiGraph.successors(node)
			parents = self.DiGraph.predecessors(node)
			N_children = self.DiGraph.out_degree(node)
			N_parents = self.DiGraph.in_degree(node)

			N_outcomes = np.max([N_parents, N_children]) # number of MC outcomes at n
			N_term = np.max([0, N_parents-N_children]) # number of MC outcomes resulting in path termination at n

			k = 0 # iterator for the MC option
			while(mc_rand_nos_dict[node] > k/N_outcomes):
				k+=1 # partition into which the random no. falls
			k = (k+iterations_dict[node]) % N_outcomes # if node was accessed before, take the next option
			iterations_dict[node] += 1 # count number of times node has been accessed

			if(k < N_term): # we terminate the path here
				return lifetime
			else: # we continue along the trajectory
				lifetime += self.timestep
				child = list(children)[k-N_term]
				for attr in attrs_dict.keys():
					attrs_dict[attr].append(self.DiGraph.nodes[child][attr])
				edges_visited.append((node, child)) # book-keeping
				return walk_trajectory(child)

		'''algorithm checks: (1) no edge is visited more than once,
		(2) all edges are visited, and (3) number of lifetimes = number
		of formation events'''
		if(len(edges_visited) != len(set(edges_visited))):
			raise ValueError("an edge has been visited more than once!")
			return 0
		elif(len(edges_visited) != graph.number_of_edges()):
			raise ValueError(
				"not all edges have been visited! Edges visited = {:s}, total edges = {:s}".format(
					str(len(edges_visited)), str(graph.number_of_edges())
				))
			return 0
		elif(len(lifetimes) != formation_events):
			raise ValueError(
				"number of lifetimes = {:s}, number of formation events = {:s}".format(
					str(len(lifetimes)), str(formation_events)
				))
			return 0

		# save trajectories with their lifetimes and properties
		attrs_dict['lifetime'] = np.array(lifetimes)

		# save to disk as pickled dict
		filename = DEFAULT_OUTPUT_DIR + "/mc_traj_{:s}.pkl".format(str(mc_no))
		with open(filename, 'wb') as f:
			pickle.dump(attrs_dict, f)
		print("saved MC trajectory to {:s}".format(filename))























