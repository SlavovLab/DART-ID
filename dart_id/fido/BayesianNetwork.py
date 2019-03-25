#!/usr/bin/env python3
# coding: utf-8
# 
# Downloaded from https://noble.gs.washington.edu/proj/fido/
# Fido is described in this paper: http://dx.doi.org/10.1021/pr100594k
# Modified 2018 by Albert Chen, Northeastern

import argparse
import getopt
import io
import itertools
import logging
import networkx as nx
import os
import pandas as pd
import random
import sys

from dart_id.fido.Utilities import *
from dart_id.fido.GraphUtilities import *
from dart_id.helper import *
from pprint import pprint
from time import time

logger = logging.getLogger('root')

class potential_table:
  def __init__(self, domain):
    self.potentials_on_frozensets = {}
    self.domain = domain
  def tuple_to_domain_frozenset(self, t):
    return frozenset( [(var,value) for var,value in t if var in self.domain] )
  def __getitem__(self, item):
    return self.potentials_on_frozensets[self.tuple_to_domain_frozenset(item)]
  def __setitem__(self, item, value):
    self.potentials_on_frozensets[self.tuple_to_domain_frozenset(item)] = value
  def __contains__(self, item):
    return self.tuple_to_domain_frozenset(item) in self.potentials_on_frozensets

class distribution(potential_table):
  def __init__(self, domain_vars_to_outcomes):
    logger.info('\t\tmaking distribution on {} variables'.format(len(domain_vars_to_outcomes)))
    potential_table.__init__(self, list(domain_vars_to_outcomes) )
    self.domain_vars_to_outcomes = domain_vars_to_outcomes
  def get_all_possible(self):
    outcomes = list( itertools.product(*[ single_outcome for var_name,single_outcome in list(self.domain_vars_to_outcomes.items())]) )
    var_outcome_lists = [ list(zip(list(self.domain_vars_to_outcomes), o)) for o in outcomes ]
    return var_outcome_lists
  def display(self):
    for o in self.get_all_possible():
      logger.info('{} --> {}'.format(o, self[o]))
  def marginalized_out(self, var_set):
    remaining = set(self.domain_vars_to_outcomes).difference(var_set)
    marginal_distribution = distribution( dict([ (v, self.domain_vars_to_outcomes[v]) for v in remaining ]) )
    for rc in marginal_distribution.get_all_possible():
      marginal_distribution[rc] = 0.0
    total = 0.0
    for c in self.get_all_possible():
      value_at_c = self[c]
      marginal_distribution[c] += value_at_c
      total += value_at_c
    for rc in marginal_distribution.get_all_possible():
      marginal_distribution[rc] /= total
    return marginal_distribution
  def __mul__(self, other_distribution):
    return distribution.multiply([self, other_distribution])
  def __div__(self, other_distribution):
    ### other_distribution should have a domain that is a subset
    ### of the current distribution
    new_distribution = distribution(self.domain_vars_to_outcomes)
    for o in self.get_all_possible():
    ### note: this will not work if other_distribution contains
    ### any 0 values; but in those cases, any value will work for
    ### belief propagation because it will be multiplied by zero
      if abs(other_distribution[o]) < 1e-5:
        #new_distribution[o] = 1
        new_distribution[o] = 0
      else:
        new_distribution[o] = self[o] / other_distribution[o]
    return new_distribution
  def distribution_constant(self):
    tot = 0.0
    for o in self.all_outcomes():
      tot += self[o]
    return tot
  @staticmethod
  def multiply(distribution_list):
    new_domain_vars_to_outcomes = {}
    for d in distribution_list:
      new_domain_vars_to_outcomes.update(d.domain_vars_to_outcomes)
    joint_product_distribution = distribution(new_domain_vars_to_outcomes)
    for c in joint_product_distribution.get_all_possible():
      prod = 1.0
      for d in distribution_list:
        prod *= d[c]
      joint_product_distribution[c] = prod
    return joint_product_distribution

class identifier(hashable_dict):
  def __str__(self):
    return 'identifier(' + str(self.my_dict) + ')'

def fs_node_ids(fs):
  return [n['id'] for n in fs]

class node:
  def __init__(self, identity, outcomes = (True, False)):
    self.id = identity
    self.outcomes = outcomes
    self.parameter_names_to_nodes = None
  def __str__(self):
    return 'node(' + str(self.id) + ', param_map = ' + str(self.parameter_names_to_nodes) + ')'
  def __repr__(self):
    return self.__str__()
  def potential(self, **Kwargs):
    raise Exception('node::potential is pure virtual')
  def wrapped_potential(self, nodes_to_values):
    ### Use the nodes_to_values dict the parameter_names_to_nodes
    ### multi_dict to make a dict of parameter_names_to_values
    #logger.info('computing potential for node {}'.format(self.id))
    #logger.info('nodes to values {}'.format(nodes_to_values))
    #logger.info('names to nodes {}'.format(self.parameter_names_to_nodes))
    parameter_names_to_values = {}
    for parameter, nodes in list(self.parameter_names_to_nodes.items()):
      parameter_node_values = [ nodes_to_values[n] for n in nodes ]
      if len(parameter_node_values) > 1:
        parameter_names_to_values[ parameter ] = parameter_node_values
      else:
        parameter_names_to_values[ parameter ] = parameter_node_values[0]

    #print(parameter_names_to_values)
    return self.potential(**parameter_names_to_values)
  def __hash__(self):
    return hash(self.id)
  def __eq__(self, rhs):
    return isinstance(rhs, node) and self.id == rhs.id
  def __ne__(self, rhs):
    return not self == rhs

class iid_probability_node(node):
  gamma = 0.5
  def __init__(self, id):
    node.__init__(self, id, (True, False))
  def potential(self, my_value, **Kwargs):
    if my_value == True:
      return iid_probability_node.gamma
    return 1-iid_probability_node.gamma

class noisy_or_node(node):
  alpha = 0.25
  beta = 0.01
  def __init__(self, id):
    node.__init__(self, id, (True, False))
  def potential(self, input, my_value, **Kwargs):
    ### when there is a single input, you should convert it to a
    ### list so that you can measure its cardinality
    if not isinstance(input, list):
      input = [input]
    count = len([ i for i in input if i == True ])
    false_prob = (1-noisy_or_node.beta)* ( (1-noisy_or_node.alpha)**count )
    if my_value == False:
      return false_prob
    return 1-false_prob

class predecessor_table_node(node):
  def __init__(self, id, likelihood_table):
    node.__init__(self, id, (True,))
    self.likelihood_table = likelihood_table
  def potential(self, input, my_value, **Kwargs):
    return self.likelihood_table[my_value][input]

def bayesian_network_union(bn1, bn2):
  if bn1 == None: return bn2
  if bn2 == None: return bn1
  for i in bn2:
    bn1.add_node(bn2.get_node_from_id(i), **bn2.node[i])
    for j in bn2[i]:
      bn1.add_edge(bn2.get_node_from_id(i),bn2.get_node_from_id(j),**bn2[i][j])
  return bn1

def all_bayesian_network_unions(graph_lst):
  ### do not set as Graph or DiGraph so that it works with any
  ### homogeneous list
  result = None
  for g in graph_lst:
    result = bayesian_network_union(result, g)
  return result

### bayesian_network should only accept 'node' types (but can be
### indexed by identifier class)
class bayesian_network(nx.DiGraph):
  def __init__(self, **Kwargs):
    nx.DiGraph.__init__(self, **Kwargs)
    self.id_to_node_map = {}

  def __str__(self):
    sio = io.StringIO()
    for n_id in self:
      sio.write(n_id['id'] + ': ')
      targets = []
      for tgt in self[n_id]:
        targets.append(tgt['id'])

      for t in targets:
        sio.write(t)
    return sio.getvalue()
    #return 'a'

  def get_potential_distribution_for_node(self, n_id):
    n_node = self.get_node_from_id(n_id)
    pred = self.predecessors(n_id)
    for m_id in pred:
      m_node = self.get_node_from_id(m_id)
          

  def subgraph(self, ibunch):
    ### using the self.__class__() constructor will work properly
    ### when inherited from; for instance, when called from a
    ### fido_network, it should allocate sg as a new fido_network
    
    sg = self.__class__()

    identifier_set = set([i for i in ibunch])
    for i in ibunch:
      sg.add_node(self.id_to_node_map[i], **self.node[i])
      for tgt in self[i]:
        if tgt not in identifier_set: continue
        sg.add_edge(self.id_to_node_map[i], self.id_to_node_map[tgt], **self[i][tgt])

    return sg

  def get_node_from_id(self, ident):
    return self.id_to_node_map[ident]

  ### method to maintain the id_to_node_map after a change to the
  ### graph
  def update_id_to_node_map_after_change(self, all_changed_nodes):
    for n in all_changed_nodes:
      if n.id in self and n.id not in self.id_to_node_map:
        self.id_to_node_map[n.id] = n
      elif n.id not in self and n.id in self.id_to_node_map:
        del self.id_to_node_map[n.id]

  ### overridden functions for addition and removal: these functions
  ### maintain the graph and the id : node map
  def add_node(self, node_instance, **Kwargs):
    assert isinstance(node_instance, node)
    nx.DiGraph.add_node(self, node_instance.id, **Kwargs)
    self.update_id_to_node_map_after_change((node_instance,))
  def remove_node(self, node_instance, **Kwargs):
    assert isinstance(node_instance, node)
    nx.DiGraph.remove_node(self, node_instance.id, **Kwargs)
    self.update_id_to_node_map_after_change((node_instance,))
  def add_edge(self, node_a, node_b, **Kwargs):
    assert isinstance(node_a, node)
    assert isinstance(node_b, node)
    nx.DiGraph.add_edge(self, node_a.id, node_b.id, **Kwargs)
    self.update_id_to_node_map_after_change((node_a, node_b))
  def remove_edge(self, node_a, node_b, **Kwargs):
    assert isinstance(node_a, node)
    assert isinstance(node_b, node)
    nx.DiGraph.remove_edge(self, node_a.id, node_b.id, **Kwargs)
    self.update_id_to_node_map_after_change((node_a, node_b))
  def add_nodes_from(self, node_lst, **Kwargs):
    for n in node_lst:
      assert isinstance(n, node)
    id_lst = [n.id for n in node_lst]
    nx.DiGraph.add_nodes_from(self, id_lst, **Kwargs)
    self.update_id_to_node_map_after_change(node_lst)
  def remove_nodes_from(self, node_lst, **Kwargs):
    for n in node_lst:
      assert isinstance(n, node)
    id_lst = [n.id for n in node_lst]
    nx.DiGraph.remove_nodes_from(self, id_lst, **Kwargs)
  def add_edges_from(self, edge_lst, **Kwargs):
    for e in edge_lst:
      n_a, n_b = e
      assert isinstance(n_a, node)
      assert isinstance(n_b, node)
    id_lst = [(n_a.id,n_b.id) for (n_a,n_b) in node_lst]
    nx.DiGraph.add_edges_from(self, id_lst, **Kwargs)
    self.update_id_to_node_map_after_change([n_a for n_a,n_b in edge_lst])
    self.update_id_to_node_map_after_change([n_b for n_a,n_b in edge_lst])
  def remove_edges_from(self, edge_lst, **Kwargs):
    for e in edge_lst:
      n_a, n_b = e
      assert isinstance(n_a, node)
      assert isinstance(n_b, node)
    id_lst = [(n_a.id,n_b.id) for (n_a,n_b) in node_lst]
    nx.DiGraph.remove_edges_from(self, id_lst, **Kwargs)
    self.update_id_to_node_map_after_change([n_a for n_a,n_b in edge_lst])
    self.update_id_to_node_map_after_change([n_b for n_a,n_b in edge_lst])

  def get_parameter_names_to_nodes_map(self, i):
    input_map = multi_dict()
    pred = self.predecessors(i)
    for p in pred:
      input_map.add( self[p][i]['label'], p )
    # each node depends on its own value
    input_map.add('my_value', i)
    return input_map
  def init_parameter_names_to_nodes_maps(self):
    for i in self:
      n = self.get_node_from_id(i)
      n.parameter_names_to_nodes = self.get_parameter_names_to_nodes_map(i)

class fido_network(bayesian_network):
  def __init__(self, fido_connected_protein_threshold = 15, \
    omit_clean_peptide_name = False, \
    proteins_column = 'Proteins', protein_delimiter = ';', \
    leading_protein_column = 'Leading razor protein', decoy_tag = 'REV__', 
    sequence_column = 'Sequence', error_prob_column = 'PEP', **Kwargs):

    bayesian_network.__init__(self, **Kwargs)
    self.fido_connected_protein_threshold = fido_connected_protein_threshold
    self.clean_peptide_id = ~omit_clean_peptide_name

    self.proteins_column = proteins_column
    self.protein_delimiter = protein_delimiter
    self.leading_protein_column = leading_protein_column
    self.decoy_tag = decoy_tag

    self.sequence_column = sequence_column
    self.error_prob_column = error_prob_column

  @staticmethod
  def likelihood_term(connected_fido_network, pep, pep_state, protein_config):
    ### computes the likelihood term from one peptide (proteins
    ### --> peptide and peptide --> spectrum) in a given state
    pep_node = connected_fido_network.get_node_from_id(pep)

    spectrum_id = list(connected_fido_network.successors(pep))[0]
    spectrum_node = connected_fido_network.get_node_from_id(spectrum_id)
    prop_prob_spectrum_given_peptide_states = spectrum_node.likelihood_table[True]
    
    protein_config[pep] = pep_state
    prob_peptide_state_given_proteins = pep_node.wrapped_potential(protein_config)
    prob_at_pep_state = prob_peptide_state_given_proteins * prop_prob_spectrum_given_peptide_states[pep_state]
    del protein_config[pep]

    return prob_at_pep_state

  @staticmethod
  def likelihood_protein_configuration(connected_fido_network, protein_config):
    ### go through the peptide terms and get sum_e Pr(E=e | R=r) and
    ### Pr(D_i | E=e)
    peps = [ pep for pep in connected_fido_network if pep['type'] == 'peptide' ]
    likelihood = 1.0
    for pep in peps:
      likelihood *= sum([ fido_network.likelihood_term(connected_fido_network, pep, pep_state, protein_config) for pep_state in (True,False) ])
      #likelihood *= sum([ prob_peptide_states[i] * prop_prob_spectrum_given_peptide_states[i] for i in [False, True]])
    return likelihood

  @staticmethod
  def prob_protein_configuration(connected_fido_network, protein_config):
    n = len(protein_config)
    r = len([ (prot,state) for prot,state in list(protein_config.items()) if state == True ])
    return pow(iid_probability_node.gamma, r) * pow(1-iid_probability_node.gamma, n-r)


  def FidoMarginalization_Inference(self):
    ### sum over the power set of proteins (for each protein
    ### configuration, compute the likelihood for the peptide
    ### sets)

    all_posteriors = {}

    subgraphs = nx.weakly_connected_component_subgraphs(self, copy=False)
    for sg_number, sg in enumerate(subgraphs):

      #if log_connected_naive_complexity(sg) > 13:
      #logger.info('fido is marginalizing subgraph with {}'.format(log_connected_naive_complexity(sg)))

      prot_nodes = [ p for p in sg if p['type'] == 'protein' ]
      prot_outcomes = itertools.product(*[(True, False) for p in prot_nodes])

      subgraph_posteriors = dict.fromkeys(prot_nodes, 0.0)
      total_prop_sum = 0.0
      for outcome in prot_outcomes:
        outcome_map = dict(list(zip(prot_nodes,outcome)))
        #print(outcome_map)
        likelihood = fido_network.likelihood_protein_configuration(sg, outcome_map)
        #print(likelihood)
        prior = fido_network.prob_protein_configuration(sg, outcome_map)
        prop_posterior = likelihood * prior
        total_prop_sum += prop_posterior
        for prot,state in list(outcome_map.items()):
          if state == True: subgraph_posteriors[prot] += prop_posterior

      for prot in subgraph_posteriors:
        subgraph_posteriors[prot] /= total_prop_sum
      all_posteriors.update(subgraph_posteriors)


    return all_posteriors

  @staticmethod
  def prune_low_scoring_peptides(dg, score_threshold = 5e-2):
    # go through the nodes with peptide type and call prune_node if
    # the probability is below the threshold
    to_prune_list = []
    for n in dg:
      if n['type'] == 'peptide':
        spectrum_id = list(dg.successors(n))[0]
        spectrum_node = dg.get_node_from_id(spectrum_id)
        if spectrum_node.likelihood_table[True][True] < score_threshold:
          to_prune_list.append(n)
    for pep in to_prune_list:
      fido_network.prune_peptide(dg, pep)
    return dg

  @staticmethod
  def prune_peptide(dg, pep):
    #logger.info('pruning peptide {}'.format(pep))
    associated_proteins = dg.predecessors(pep)
    pep_node = dg.get_node_from_id(pep)
    spect_id = list(dg.successors(pep))[0]
    spect_node = dg.get_node_from_id(spect_id)
    for prot_number, prot in enumerate(associated_proteins):
      prot_node = dg.get_node_from_id(prot)

      new_pep_id_dict = pep.my_dict.copy()
      new_pep_id_dict['copy'] = prot_number
      new_pep_id = identifier(new_pep_id_dict)
      cloned_pep_node = noisy_or_node(new_pep_id)

      new_spect_id_dict = spect_id.my_dict.copy()
      new_spect_id_dict['copy'] = prot_number
      new_spect_id = identifier(new_spect_id_dict)
      cloned_spect_node = predecessor_table_node(new_spect_id, spect_node.likelihood_table)

      dg.add_edge(prot_node, cloned_pep_node, **dg[prot][pep])
      dg.add_edge(cloned_pep_node, cloned_spect_node, **dg[pep][spect_id])
    dg.remove_node(pep_node)
    dg.remove_node(spect_node)

  @staticmethod
  def connected_lowest_prune(connected_dg):
    pep_nodes = [ n for n in connected_dg if n['type'] == 'peptide' and len(list(connected_dg.predecessors(n))) > 1 ]
    pep_present_likelihoods = [ connected_dg.get_node_from_id(list(connected_dg.successors(pep))[0]).likelihood_table[True][True] for pep in pep_nodes ]

    #pep_score, pep_to_prune = np.min(list(zip(pep_present_likelihoods, pep_nodes)))
    pep_to_prune = pep_nodes[np.argmin(pep_present_likelihoods)]

    fido_network.prune_peptide(connected_dg, pep_to_prune)

  @staticmethod
  def dynamic_pruned(dg, connected_protein_threshold, **Kwargs):
    connected_protein_threshold = float(connected_protein_threshold)
    dg_lst = nx.weakly_connected_component_subgraphs(dg, copy=False)
    pruned_dg_lst = fido_network.dynamic_pruned_helper(dg_lst, connected_protein_threshold)
    return all_bayesian_network_unions(pruned_dg_lst)
        
  @staticmethod
  def dynamic_pruned_helper(dg_lst, connected_protein_threshold):
    ### may be slow-- linear time per iteration and possibly many
    ### iterations (--> quadratic in number of peptides in a connected component)
    result_lst = []
    for dg in dg_lst:
      if log_connected_naive_complexity(dg) < connected_protein_threshold:
        result_lst.append(dg)
      else:
        #logger.info('pruning connected graph with {} proteins'.format(log_connected_naive_complexity(dg)))
        while True:
          # prune the lowest-scoring peptide until the graph is split
          fido_network.connected_lowest_prune(dg)
          subgraphs = nx.weakly_connected_component_subgraphs(dg, copy=False)
          if len(list(subgraphs)) > 1:
            pruned_subgraphs = fido_network.dynamic_pruned_helper(subgraphs, connected_protein_threshold)
            result_lst.extend(pruned_subgraphs)
            break

    return result_lst

  @staticmethod
  def cluster_proteins_from_lists(peptide_spectrum_proteins_prob_charge):
    logger.info('clustering from lists')
    #         given pep --> prots
    #         get
    #         1. prot --> peps
    #         to get
    #         2. peps --> prot
    #         to get prots with identical peps

    prot_to_peps = {}
    for pep,spect,prots,prob,charge in peptide_spectrum_proteins_prob_charge:
      for prot in prots:
        if prot in prot_to_peps: prot_to_peps[prot].update(set([pep]))
        else:                    prot_to_peps[prot] = set([pep])

    peps_to_prots = {}

    for prot,pep_set in list(prot_to_peps.items()):
      pep_set = frozenset(pep_set)
      if pep_set in peps_to_prots: peps_to_prots[pep_set].update(set([prot]))
      else:                        peps_to_prots[pep_set] = set([prot])

    prot_to_protgroup_and_number = {}
    count = 0
    for pep_set,prot_set in list(peps_to_prots.items()):
      if len(prot_set) > 1:
        for prot in prot_set:
          prot_to_protgroup_and_number[prot] = (prot_set,count)
        count += 1

    new_peptide_spectrum_proteins_prob_charge = []
    for pep,spect,prots,prob,charge in peptide_spectrum_proteins_prob_charge:
      new_prots = []
      for prot in prots:
        if prot in prot_to_protgroup_and_number:
          prot_group_set, prot_group_number = prot_to_protgroup_and_number[prot]
          new_prots.append('group_' + str(prot_group_number))
        else: new_prots.append(prot)
      new_peptide_spectrum_proteins_prob_charge.append((pep,spect,new_prots,prob,charge))

    group_set = set()
    for prot in prot_to_protgroup_and_number:
      group_set.add(frozenset(prot_to_protgroup_and_number[prot][0]))

    return new_peptide_spectrum_proteins_prob_charge

  def cluster_proteins(self):
    logger.info('clustering')
    md = multi_dict()
    ### find the proteins that have identical peptide connectivity

    for n in self.node:
      if n['type'] == 'protein':
        md.add(frozenset(self[n]), n)
    ### remove these proteins from the graph and replace them with
    ### a single group that contains all of them
    removal_protein_set = set()
    group_number = 0
    group_set = set()
    for pep_set in md:
      prot_list = list(md[pep_set])
      if len(prot_list) > 1:
        group_name = 'protein_group_'+str(group_number)
        group_set.add(group_name)
        prot_group_node = iid_probability_node(identifier({'id':group_name, 'type':'protein', 'protein_group':frozenset(prot_list)}))
        for pep_id in pep_set:
          pep_node = self.get_node_from_id(pep_id)
          ### note: this will need to be modified if more
          ### complicated edge labes are used
          self.add_edge(prot_group_node, pep_node, label='input')
        removal_protein_set.update(set(prot_list))
        group_number += 1
    self.remove_nodes_from([ self.get_node_from_id(n_id) for n_id in removal_protein_set ])
    return set([ rp['id'] for rp in removal_protein_set ])

  @staticmethod
  def wrap_peptide_id(peptide, assumed_charge, merge_peptide_charge_states):
    peptide_id = {'id':peptide, 'type':'peptide'}
    if not merge_peptide_charge_states:
      peptide_id['charge'] = assumed_charge
    return identifier(peptide_id)

  @staticmethod
  def change_isoleucine_to_leucine(pepStr):
    ### this will ensure that peptides that are identical aside
    ### from I -> L or L -> I changes will be treated as
    ### degenerate
    return pepStr.replace('I','L')

  @staticmethod
  def remove_bounding_amino_acids(pepStr):
    if len(pepStr) > 1 and pepStr[1] == '.':
      return pepStr[2:-2]
    else: return pepStr

  def load_from_pivdo2(self, filename, merge_peptide_charge_states=True):
    bayesian_network.__init__(self)
    ### load the pivdo2 file
    file = open(filename, 'r')
    pivdo2_lines = file.readlines()
    ### get the charge priors and edge information
    charge_priors = fido_network.charge_priors_from_pivdo2(pivdo2_lines)

    peptide_spectrum_proteins_prob_charge = fido_network.get_graph_data_from_pivdo2(pivdo2_lines)

    self.make_graph_from_charge_priors_and_lists( charge_priors = charge_priors, peptide_spectrum_proteins_prob_charge = peptide_spectrum_proteins_prob_charge, merge_peptide_charge_states = merge_peptide_charge_states )

  @staticmethod
  def charge_priors_from_pivdo2(pivdo_lines):
    ### fixme error: this makes pruning use peptideprophet scores
    ### by making the likelihood equal the probability; fix this
    ### before you do any inference!
    return dict.fromkeys(list(range(-1,10)),0.1)

    ### use -1 as unspecified charge
    #charge_priors = { -1:0.1 }
    for line in pivdo_lines:
      line_list = line.split()
      if line_list[0] == 'd':
        charge = int(line_list[1])
        probability = float(line_list[2])
        charge_priors[charge] = probability
    return charge_priors

  @staticmethod
  def get_graph_data_from_pivdo2(pivdo_lines):
    peptide_spectrum_proteins_prob_charge = []

    spectrum_count = 0
    state = ('e')
    for line in pivdo_lines:
      line_list = line.split()
      line_type = line_list[0]
          
      if line_type == 'd': continue

      if line_type not in state:
        raise Exception('Wanted line_type ' + ''.join(state) + ' and got ' + line_type + ' one line: ' + line)

      if line_type == 'e':
        if state == ('r','p','e'):
          peptide_spectrum_proteins_prob_charge.append( (peptide, 'spectrum_'+str(spectrum_count), protein_list, None, charge) )
          spectrum_count += 1

        peptide = line_list[1]
        ### use -1 as unknown or unspecified charge
        charge = -1
        protein_list = []
        state = ('c', 'r')

      elif line_type == 'c':
        charge = int(line_list[1])
        state = ('r')
      elif line_type == 'r':
        protein = line_list[1]
        protein_list.append(protein)
        state = ('r','p','e')
      elif line_type == 'p':
        probability = float(line_list[1])
        state = ('e')
      else:
        raise Exception('Unrecognized line_type in pivdo2 ' + line_type + ' on line ' + line)
      
      if line_type == 'p':
        peptide_spectrum_proteins_prob_charge.append( (peptide, 'spectrum_'+str(spectrum_count), protein_list, probability, charge) )
        spectrum_count += 1

    return peptide_spectrum_proteins_prob_charge

  def make_graph_from_charge_priors_and_lists(self, charge_priors, peptide_spectrum_proteins_prob_charge, merge_peptide_charge_states):
    for peptide, spectrum, protein_list, peptide_probability, assumed_charge in peptide_spectrum_proteins_prob_charge:

      if self.clean_peptide_id:
        peptide = fido_network.change_isoleucine_to_leucine(peptide)
        peptide = fido_network.remove_bounding_amino_acids(peptide)
      
      peptide_id = fido_network.wrap_peptide_id(peptide, assumed_charge, merge_peptide_charge_states)

      ### add protein --> peptide edges
      pep_node = noisy_or_node( peptide_id )
      for prot in protein_list:
        protein_id = identifier({'id':prot, 'type':'protein'})
        prot_node = iid_probability_node(protein_id)
        if prot_node not in self:
          self.add_node( prot_node )
        self.add_edge( prot_node, pep_node, label='input' )

      if peptide_probability != None:
        ### compute the likelihoods from the prior for that particular charge state
        cp = charge_priors[assumed_charge]

        likelihood_given_present = peptide_probability / cp
        likelihood_given_absent = (1-peptide_probability) / (1-cp)
      
        ### add peptide --> spectrum edges
        spectrum_id = identifier({'id':spectrum, 'type':'spectrum'})

        spectrum_node = predecessor_table_node( spectrum_id,  {True:{True: likelihood_given_present, False: likelihood_given_absent}} )
        self.add_edge(pep_node, spectrum_node, label='input')
      else:
        logger.info('Warning: found search_result {} with no spectrum matches'.format(peptide))
        ### (it may match be a non-best match to a spectrum or
        ### may be found elsewhere in the graph with its own
        ### spectrum match)

  def remove_all_but_maximum_likelihood_spectrum(self):
    pep_ids = [ n for n in self if n['type'] == 'peptide' ]
    for pi in pep_ids:
      spectrum_ids = [ sp for sp in self.successors(pi) ]
      if len(spectrum_ids) == 0:
        logger.error('error: peptide', pi['id'], 'has no spectrum matches')
      spectrum_nodes = [ self.get_node_from_id(spect_id) for spect_id in spectrum_ids ] 
      pep_node = self.get_node_from_id(pi)
      likelihoods_present_peptide = [ sp.likelihood_table[True][True] for sp in spectrum_nodes]
      highest_likelihood_spectrum_node = spectrum_nodes[ index_max(likelihoods_present_peptide) ]
      self.remove_nodes_from([ self.get_node_from_id(n_id) for n_id in self.successors(pi) ])
      self.add_edge(pep_node, highest_likelihood_spectrum_node, label='input')

  def multi_file_load(self, fname_list, filetype = 'pivdo', **Kwargs):
    for fname in fname_list:
      logger.info('loading from {}'.format(fname))
      sub_fn = fido_network(**Kwargs)
      ### use the command line value for filetype
      sub_fn.load(fname, filetype, **Kwargs)
      sub_fn.remove_all_but_maximum_likelihood_spectrum()
      self = bayesian_network_union(self, sub_fn)

    logger.info('prots: {}'.format(len([ n for n in self if n['type'] == 'protein' ])))
    logger.info('peps: {}'.format(len([ n for n in self if n['type'] == 'peptide' ])))

  def load(self, fname, filetype = 'pivdo', **Kwargs):
    if filetype == 'pivdo':
      self.load_from_pivdo2(fname)
    else:
      logger.error('Error: unrecognized filetype in fido_network::load:{}'.format(filetype))
      sys.exit(2)

  def load_dataframe(self, df, merge_peptide_charge_states=True, **Kwargs):
    bayesian_network.__init__(self)

    # get charge info -- dummied for now. see charge_priors_from_pivdo2
    #charge_priors = dict.fromkeys(list(range(-1,10)), 1e-6)
    charge_priors = dict.fromkeys(list(range(-1,10)), 0.1)

    # get edge info
    # each element of this list is a tuple describing the PSM
    peptide_spectrum_proteins_prob_charge = []

    sequence = ''
    prot_list = []
    prob = 0.0
    charge = -1 # dummy charge for now
    spectrum_count = 0

    #prot_cols = ['Leading proteins', 'Proteins']
    #prot_cols = ['Leading proteins']
    prot_cols = [self.proteins_column]

    for i in range(0, df.shape[0]):
      sequence = df[self.sequence_column][i]
      prob = 1 - df[self.error_prob_column][i]

      prot_list = []
      for col in prot_cols:
        prots = df[col].loc[i]
        if pd.isnull(prots): continue
        #prot_list = prot_list + [str.join('_', prot.split('|')[0:2]) for prot in prots.split(';')]
        prot_list = prot_list + [prot for prot in prots.split(self.protein_delimiter)]

      prot_list = np.unique(prot_list)
      if len(prot_list) < 1: continue

      peptide_spectrum_proteins_prob_charge.append( (sequence, 'spectrum_'+str(spectrum_count), prot_list, prob, charge) )
      spectrum_count += 1

    self.make_graph_from_charge_priors_and_lists(\
      peptide_spectrum_proteins_prob_charge=peptide_spectrum_proteins_prob_charge,
      charge_priors=charge_priors,
      merge_peptide_charge_states=merge_peptide_charge_states)

  def load_from_dataframes(self, df_list, **Kwargs):
    for df in df_list:
      sub_fn = fido_network(**Kwargs)
      sub_fn.load_dataframe(df, **Kwargs)
      #sub_fn.remove_all_but_maximum_likelihood_spectrum()
      self = bayesian_network_union(self, sub_fn)
    
    logger.info('prots: {}'.format(len([ n for n in self if n['type'] == 'protein' ])))
    logger.info('peps: {}'.format(len([ n for n in self if n['type'] == 'peptide' ])))


def set_gab(gamma, alpha, beta):
  iid_probability_node.gamma = gamma
  noisy_or_node.alpha = alpha
  noisy_or_node.beta = beta

def comp_roc_fdr(posteriors, decoy_tag='REV__'):
  fps = [0]
  tps = [0]
  est_fdr_list = [0]
  emp_fdr_list = [0]

  fp_count = 0
  tp_count = 0
  total_fdr = 0.0
  est_fdr = 0.0
  emp_fdr = 0.0

  num_decoys = 0
  num_targets = 0

  last_prob = -1
  for i in range(0, posteriors.shape[0]):
      prob = posteriors.index[i]
      prots = posteriors.iloc[i]
      #print(prots)
      fp_change = np.sum(decoy_tag in prot for prot in prots)
      tp_change = len(prots) - fp_change
      #print(fp_change, tp_change)
      
      num_decoys += fp_change
      num_targets += tp_change
      
      if tp_change > 0 and fp_change > 0:
          tp_change = 0
          fp_change = 0
      elif tp_change > 0:
          tp_change = 1
      elif fp_change > 0:
          fp_change = 1
      
      if i != 0 and prob != last_prob:
          fps.append(fp_count)
          tps.append(tp_count)
          
          if est_fdr > est_fdr_list[-1]:
              est_fdr_list.append(est_fdr)
              emp_fdr_list.append(emp_fdr)
      
      fp_count += fp_change
      tp_count += tp_change
      
      total_fdr += ((1-prob) * (fp_change + tp_change))
      est_fdr = (total_fdr / (fp_count + tp_count))
      emp_fdr = (float(fp_count) / (fp_count + tp_count))

      last_prob = prob

  fps.append(fp_count)
  tps.append(tp_count)

  fps.append(num_decoys)
  tps.append(num_targets)

  est_fdr_list = np.array(est_fdr_list)
  emp_fdr_list = np.array(emp_fdr_list)

  return (fps, tps, est_fdr_list, emp_fdr_list)

def antiderivative_at(m, b, xVal):
    return (m * xVal * xVal / 2.0) + (b * xVal)

def squared_antiderivative_at(m, b, xVal):
    u = m*m
    v = 2*m*b
    t = b*b
    
    return (u * xVal * xVal / 3.0) + (v * xVal * xVal / 2.0) + (t * xVal)

def area(x1, y1, x2, y2, max_x):
    m = (y2-y1) / (x2-x1)
    b = y1-(m*x1)
    result = antiderivative_at(m, b, np.min([max_x, x2])) - antiderivative_at(m, b, x1)
    if result < 0.0:
        logger.error('area: {}\n{} {} {} {}'.format(result, m, b, x1, x2))
    return result

def squared_area(x1, y1, x2, y2, max_x):
    if x2 < x1: return 0.0
    
    m = ((y2-y1) / (x2-x1))
    b = (y1-(m*x1))
    
    result = (squared_antiderivative_at(m, b, np.min([max_x, x2])) - \
              squared_antiderivative_at(m, b, x1))
    return result

def roc_N(fps, tps, N=50):
  roc_N = 0.0

  if fps[-1] < N:
    logger.warning('Warning: There are not enough false positives; needed {} and was only given {}. Will proceed using largest available value'.format(N, fps[-1]))
    N = fps[-1]
      
  for i in range(0, len(fps)-1):
    if fps[i] >= N: break
    if fps[i] != fps[i+1]:
      current_area = area(fps[i], tps[i], fps[i+1], tps[i+1], N)
      roc_N += current_area
      
  return roc_N / (N * tps[-1])

def fdr_divergence(est_fdr, emp_fdr, thresh=0.1):
  # FDR divergence
  diff = est_fdr - emp_fdr
  tot = 0.0

  for i in range(0, len(diff)-1):
    if est_fdr[i] >= thresh:
      if i == 0:
        tot = np.inf
      break
    
    tot += squared_area(est_fdr[i], diff[i], est_fdr[i+1], diff[i+1], est_fdr[i+1])

  x_range = np.min([thresh, est_fdr[-2]]) - est_fdr[0]

  if np.isinf(tot):
    logger.error('infinite FDR divergence')
    return tot

  return tot / x_range


def run_internal(df, parameter_map):

  if parameter_map['gamma'] is None or parameter_map['alpha'] is None or \
    parameter_map['beta'] is None:

    logger.info('gamma, alpha, beta, not defined. choosing best parameters...')

    logger.info('estimating parameters with accuracy level {}'.format(\
      parameter_map['parameter_accuracy']))
    if parameter_map['parameter_accuracy'] > 1:
      logger.info('sampling smaller dataframe from inputs...')
      sample_size = -1
      if parameter_map['parameter_accuracy'] == 2:
        sample_size = 300
      elif parameter_map['parameter_accuracy'] == 3:
        sample_size = 100

      dfs = df.sample(n=sample_size, random_state=1).reset_index(drop=True)
    else:
      dfs = df

    logger.info('loading into graph network...')
    fn = fido_network(**parameter_map)
    fn.load_from_dataframes(dfs_list, **parameter_map)

    if parameter_map['omit_clean_peptide_name']:
      logger.info('keeping original peptide sequences...')
    else:
      logger.info('cleaning peptide sequences...')

    if parameter_map['all_psms']:
      logger.info('including all PSMs...')
    else:
      logger.info('trimming down to best spectrum match...')
      fn.remove_all_but_maximum_likelihood_spectrum()

    logger.info('pruning...')
    #fn = fido_network.prune_low_scoring_peptides(fn, 1e-2)
    fn = fido_network.dynamic_pruned(fn, **parameter_map)
    fn.init_parameter_names_to_nodes_maps()

    gamma_search = [0.1, 0.5, 0.9];
    alpha_search = [0.01, 0.04, 0.09, 0.16, 0.25, 0.36];
    beta_search = [0.0, 0.01, 0.025, 0.05];

    logger.info('searching parameters for {} iterations'.format(
      len(gamma_search)*len(alpha_search)*len(beta_search)))

    gamma_best = 0.5
    alpha_best = 0.1
    beta_best = 0.01

    best_objective = -10000000
    iter_number = 1

    for i, gamma in enumerate(gamma_search):
      for j, alpha in enumerate(alpha_search):
        for k, beta in enumerate(beta_search):
          # set gamma, alpha, beta
          set_gab(gamma, alpha, beta)

          prots_to_posteriors = fn.FidoMarginalization_Inference()
          # convert to dataframe
          ids = list(prots_to_posteriors.keys())
          probs = [prots_to_posteriors[_id] for _id in ids]
          ids = [_id['id'] for _id in ids]
          posteriors = pd.DataFrame({'prot': ids, 'prob': probs})\
            .sort_values('prob', ascending=False).reset_index(drop=True)
          # group proteins by posterior probability
          posteriors = posteriors.groupby('prob')['prot'].apply(lambda x: x.values)\
            .sort_index(ascending=False)

          fps, tps, est_fdr, emp_fdr = comp_roc_fdr(posteriors, decoy_tag=parameter_map['decoy_tag'])
          roc50 = roc_N(fps, tps, N=50)
          fdr_mse = fdr_divergence(est_fdr, emp_fdr, thresh=0.1)
          _lambda = 0.15

          current_objective = (_lambda * roc50) - ((1 - _lambda) * fdr_mse)
          if current_objective > best_objective:
            best_objective = current_objective
            gamma_best = gamma
            alpha_best = alpha
            beta_best = beta

          logger.info('iteration {} || gamma: {}, alpha: {}, beta: {}'.format(
            iter_number, gamma, alpha, beta))
          logger.info('                roc50: {:.4f}, fdr_mse: {:.4f}, objective: {:.4f}'.format(\
            roc50, fdr_mse, current_objective))

          iter_number += 1

    logger.info('best objective: {:.4f}, with parameters {}, {}, {}'.format(best_objective, gamma_best, alpha_best, beta_best))

    # set gamma, alpha, beta
    parameter_map['gamma'] = gamma_best
    parameter_map['alpha'] = alpha_best
    parameter_map['beta'] = beta_best

  logger.info('running with parameters gamma: {}, alpha: {}, beta: {}'.format(
    parameter_map['gamma'], parameter_map['alpha'], parameter_map['beta']))

  # run with the best parameters
  set_gab(parameter_map['gamma'], parameter_map['alpha'], parameter_map['beta'])

  fn = fido_network(**parameter_map)
  
  logger.info('loading into graph network...')
  fn.load_from_dataframes([df], **parameter_map)

  if parameter_map['all_psms']:
    logger.info('including all PSMs...')
  else:
    logger.info('trimming down to best spectrum match...')
    fn.remove_all_but_maximum_likelihood_spectrum()

  if parameter_map['group_proteins']:
    logger.info('clustering proteins...')
    fn.cluster_proteins()

  logger.info('pruning...')
  if parameter_map['prune_low_scores']:
    logger.info('first pruning low scoring peptides (confidence < 1e-2)')
    fn = fido_network.prune_low_scoring_peptides(fn, 1e-2)

  fn = fido_network.dynamic_pruned(fn, **parameter_map)
  fn.init_parameter_names_to_nodes_maps()
  
  logger.info('running fido inference...')
  prots_to_posteriors = fn.FidoMarginalization_Inference()
  #pprint(prots_to_posteriors)

  ids = list(prots_to_posteriors.keys())
  probs = [prots_to_posteriors[_id] for _id in ids]
  ids = [_id['id'] for _id in ids]
  
  out_df = pd.DataFrame({'Protein': ids, 'Probability': probs})\
    .sort_values('Probability', ascending=False).reset_index(drop=True)

  #print(np.percentile(out_df['Probability'], np.arange(10, 100, 10)))

  out_df['Error Probability'] = 1 - out_df['Probability']
  out_df['q-value'] = np.cumsum(out_df['Error Probability'].values) / np.arange(1, out_df.shape[0]+1, 1)

  # save output from fido to a file, protein_fdr.txt
  logger.info('Writing Fido protein inference output to file: {}'.format(os.path.join(parameter_map['output'], 'protein_fdr.txt')))
  out_df.to_csv(os.path.join(parameter_map['output'], 'protein_fdr.txt'), 
    sep='\t', index=False)

  # attach PSM razor protein FDR to the input dataframe
  logger.info('Attaching razor protein FDR to output')
  fdr_series = pd.Series(out_df['q-value'].values, index=out_df['Protein'].values)
  df['razor_protein_fdr'] = df[parameter_map['leading_protein_column']].map(fdr_series)

  #print(np.sum(df['razor_protein_fdr'] < 0.1), df.shape[0])

  return df
  


def pd_main():

  # initialize logger
  init_logger(True, 'fido.log', log_to_file=False)

  parser = argparse.ArgumentParser()  
  parser.add_argument('-l', '--connected-protein-threshold', type=int, default=14,
    help='log2 of maximum number of subgraph connected states')
  parser.add_argument('-p', '--omit-clean-peptide-name', action='store_true', default=False,
    help='omit cleaning the peptide names')
  parser.add_argument('-a', '--all-psms', action='store_true', default=False,
    help='use all PSM matches instead the best one')
  parser.add_argument('-g', '--group-proteins', action='store_true', default=False,
    help='use protein group level inference')
  parser.add_argument('-s', '--prune-low-scores', action='store_true', default=False,
    help='prune low scoring PSMs (threshold set at 0.01)')
  parser.add_argument('-c', '--parameter-accuracy', type=int, default=3,
    help='set start parameter\'s accurary level (1-3)\n1 = best (slower)\n2 = relaxed (faster)\n3 = sloppy (very fast)')

  parser.add_argument('--gamma', type=float, default=None,
    help='protein prior probability')
  parser.add_argument('--alpha', type=float, default=None,
    help='peptide emission probability')
  parser.add_argument('--beta', type=float, default=None,
    help='spurious peptide identification probability')

  parser.add_argument('--proteins-column', type=str, default='Proteins',
    help='column name for the proteins of each PSM')
  parser.add_argument('--protein-delimiter', type=str, default=';',
    help='delimiter character separating proteins in the proteins column')
  parser.add_argument('--leading-protein-column', type=str, default='Leading razor protein',
    help='column name for the leading (razor) protein of each PSM')
  parser.add_argument('--decoy-tag', type=str, default='REV__',
    help='expression used to deliminate decoy proteins from target proteins')
  parser.add_argument('--sequence-column', type=str, default='Sequence',
    help='column name for the peptide sequence of each PSM')
  parser.add_argument('--error-prob-column', type=str, default='PEP',
    help='column name for the error probability of each PSM')

  parser.add_argument('-o', '--output', type=str, default=os.getcwd(), help='Folder to output protein_fdr list and intermediate fido files')

  args = parser.parse_args()

  #parameter_map = {
  #  'connected_protein_threshold': 14,
  #  'gamma': 0.5, 'alpha': 0.1, 'beta': 0.01 }
  parameter_map = vars(args)

  logger.info('reading...')
  df = pd.read_csv('/gd/MS/SCoPE/SQC/SQC77/evidence.txt', sep='\t', low_memory=False)
  #df = pd.read_csv('~/git/RTLib/Alignments/SQC2_20180812_1/ev_updated.txt', 
  #  sep='\t', low_memory=False)

  run_internal(df, parameter_map)

def real_main(argv):

  # initialize logger
  init_logger(True, 'fido.log', log_to_file=False)

  try:
    opts, args = getopt.getopt(argv, 'f:gc:m:', ['filetype=', 'group_proteins', 'connected_protein_threshold=', 'marginalization_mode='])
  except getopt.GetoptError:
    print('cmd line error')
    sys.exit(2)

  ### this converts any short option (like -f) into the appropriate
  ### long option (--filetype)
  shorthand_to_longhand = {'-f':'--filetype', '-g':'--group_proteins', '-c':'--connected_protein_threshold', '-m':'--marginalization_mode'}

  longhand_opts = []
  for flag, value in opts:
    if flag in shorthand_to_longhand:
      longhand_opts.append( (shorthand_to_longhand[flag], value) )
    else:
      longhand_opts.append( (flag, value) )

  parameter_map = { \
    'marginalization_mode': 'fido', 
    'filetype': 'pivdo', 
    'connected_protein_threshold': 14 }

  #parameter_map = { }
  for flag, value in longhand_opts:
    ### strip off the leading '--'
    parameter_map[flag[2:]] = value

  print('parameter_map')
  print(parameter_map)

  fn = fido_network(**parameter_map)
  print('args',args)
  print('reading...')

  fn = fido_network(**parameter_map)
  fn.multi_file_load(args, **parameter_map)

  print('trimming down to best spectrum match...')
  fn.remove_all_but_maximum_likelihood_spectrum()

  if 'group_proteins' in parameter_map:
    print('clustering proteins')
    fn.cluster_proteins()
  
  print('pruning...')

  #fn = fido_network.prune_low_scoring_peptides(fn, 1e-2)
  fn = fido_network.dynamic_pruned(fn, **parameter_map)

  fn.init_parameter_names_to_nodes_maps()

  if parameter_map['marginalization_mode'] == 'fido':
    print('fido inference')
    prots_to_posteriors = fn.FidoMarginalization_Inference()
  else:
    print('Error: unrecognized marginalization_mode = ', parameter_map['marginalization_mode'])

  pprint(prots_to_posteriors)

main = real_main

if __name__ == '__main__':
  main(sys.argv[1:])


