#!/usr/bin/env python3
# coding: utf-8
# 
# Downloaded from https://noble.gs.washington.edu/proj/fido/
# Fido is described in this paper: http://dx.doi.org/10.1021/pr100594k
# Modified 2018 by Albert Chen, Northeastern

import getopt
import io
import itertools
import logging
import networkx as nx
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
    for parameter,nodes in list(self.parameter_names_to_nodes.items()):
      parameter_node_values = [ nodes_to_values[n] for n in nodes ]
      if len(parameter_node_values) > 1:
        parameter_names_to_values[ parameter ] = parameter_node_values
      else:
        parameter_names_to_values[ parameter ] = parameter_node_values[0]
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
  alpha = 0.1
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
    #
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
  def __init__(self, fido_connected_protein_threshold = 15, **Kwargs):
    bayesian_network.__init__(self, **Kwargs)
    self.fido_connected_protein_threshold = fido_connected_protein_threshold

  @staticmethod
  def MSE(prot_posteriors_A, prot_posteriors_B):
    ### the key sets of these dict's should be the same
    tot = 0.0
    for prot in prot_posteriors_A:
      err = prot_posteriors_A[prot] - prot_posteriors_B[prot]
      tot += err*err
    return tot / len(prot_posteriors_A)

  @staticmethod
  def largest_absolute_error(prot_posteriors_A, prot_posteriors_B):
    ### the key sets of these dict's should be the same
    absolute_errors = []
    for prot in prot_posteriors_A:
      err = abs(prot_posteriors_A[prot] - prot_posteriors_B[prot])
      absolute_errors.append(err)
    return max(absolute_errors)

  def MetropolisHastings_Inference(self):
    ### propose a new peptide and protein configuration
    
    ### accept or reject the configuration
    pass

  def Gibbs_Inference(self, ground_truth_posteriors, use_connected_block = False, max_time = 60, iters_per_time_log = 500):
    ### propose new protein and peptide blocks
    
    ### over the power set of states that this (protein, peptide)
    ### block can take on, compute the posterior (relative to the
    ### current posterior)

    ### sample a new configuration based on these relative
    ### posteriors
    all_posteriors = {}
    all_times_and_errors = []
    for sg in nx.weakly_connected_component_subgraphs(self):
      prot_nodes = [p for p in sg if p['type'] == 'protein' ]
      pep_nodes = [p for p in sg if p['type'] == 'peptide' ]

      ### start at random protein and peptide configurations
      prot_configuration = dict([ (p,(False, True)[ random.randint(0,1) ]) for p in prot_nodes])

      pep_configuration = dict([ (p,(False, True)[ random.randint(0,1) ]) for p in pep_nodes])
      ### make sure the initial peptide configuration isn't
      ### highly unlikely
      for pep in pep_configuration:
        pep_node = sg.get_node_from_id(pep)
        spectrum_id = sg.successors(pep)[0]
        spectrum_node = sg.get_node_from_id(spectrum_id)
        if spectrum_node.likelihood_table[True][ pep_configuration[pep] ] < 1e-10:
          pep_configuration[pep] = not pep_configuration[pep]

      prot_block_size = min(3, len(prot_nodes))
      pep_block_size = min(3, len(pep_nodes))

      subgraph_counts = dict.fromkeys(prot_nodes, 0.0)

      start_time = time()
      times_and_errors = []
      ### stop when the timer passes the allowed time
      for iter in range(-100,sys.maxsize/2):
        if use_connected_block:
          prot_block = sg.random_connected_protein_block(prot_block_size)
          pep_block = sg.random_connected_peptide_block(prot_block, pep_block_size)
        else:
          prot_block = sg.random_protein_block(prot_nodes, prot_block_size)
          pep_block = sg.random_peptide_block(pep_nodes, pep_block_size)

        prot_changes, pep_changes = sg.sample_new_protein_and_peptide_configuration(prot_configuration, pep_configuration, prot_block, pep_block)

        for prot in prot_changes:
          prot_configuration[prot] = prot_changes[prot]
        for pep in pep_changes:
          pep_configuration[pep] = pep_changes[pep]

        if iter >= 0:
          for prot in prot_nodes:
            if prot_configuration[prot] == True:
              subgraph_counts[prot] += 1

        if iter >= 0:
          subgraph_posteriors = dict([ (n,c/float(iter+1)) for n,c in list(subgraph_counts.items()) ])
          if iter % iters_per_time_log == 0:
            error = fido_network.largest_absolute_error(ground_truth_posteriors, subgraph_posteriors)
            elapsed_time = time() - start_time
            logger.info('{} {}'.format(elapsed_time, error))
            times_and_errors.append( (elapsed_time, error) )

            if elapsed_time > max_time: break

      all_times_and_errors.append(times_and_errors)
      all_posteriors.update(subgraph_posteriors)
    return all_posteriors, all_times_and_errors

  def CollapsedGibbs_Inference(self, ground_truth_posteriors, use_connected_block = False, max_time = 60, iters_per_time_log = 500):
    ### propose new protein blocks

    ### over the power set of states that this protein block can
    ### take on, compute the posteriors (relative to the current
    ### posterior) by marginalizing over all peptide sets

    ### sample a new configuration based on these relative
    ### posteriors

    all_posteriors = {}
    all_times_and_errors = []
    for sg in nx.weakly_connected_component_subgraphs(self):
      ### start at a random protein configuration

      logger.info('collapsed_gibbs_inference on subgraph of size {}'.format(log_connected_naive_complexity(sg)))

      prot_nodes = [p for p in sg if p['type'] == 'protein' ]
      prot_configuration = dict([ (p,(False, True)[ random.randint(0,1) ]) for p in prot_nodes])
      prot_block_size = min(3, len(prot_nodes))
      subgraph_counts = dict.fromkeys(prot_nodes, 0.0)

      start_time = time()
      times_and_errors = []
      for iter in range(-100,sys.maxsize/2):
        if use_connected_block:
          prot_block = sg.random_connected_protein_block(prot_block_size)
        else:
          prot_block = sg.random_protein_block(prot_nodes, prot_block_size)
        prot_changes = sg.sample_new_protein_configuration(prot_configuration, prot_block)

        for prot in prot_changes:
          prot_configuration[prot] = prot_changes[prot]

        if iter >= 0:
          for prot in prot_nodes:
            if prot_configuration[prot] == True:
              subgraph_counts[prot] += 1

        if iter >= 0:
          subgraph_posteriors = dict([ (n,c/float(iter+1)) for n,c in list(subgraph_counts.items()) ])
          if iter % iters_per_time_log == 0:
            elapsed_time = time() - start_time
            times_and_errors.append( ( elapsed_time, fido_network.largest_absolute_error(ground_truth_posteriors, subgraph_posteriors) ) )

            if elapsed_time > max_time: break

      all_times_and_errors.append(times_and_errors)
      all_posteriors.update(subgraph_posteriors)
    return all_posteriors, all_times_and_errors

  def random_protein_block(self, prot_nodes, block_size):
    ### prot_nodes is passed in to save the O(n) time of
    ### listing the proteins
    if block_size == 0:
      return set()

    prot_block = set()
    prot = prot_nodes[random.randint(0,len(prot_nodes)-1)]
    prot_block.add(prot)
    while len(prot_block) < block_size:
      prot = prot_nodes[random.randint(0,len(prot_nodes)-1)]
      ### this will only be fast when the block size is small
      while prot in prot_block:
        prot = prot_nodes[random.randint(0,len(prot_nodes)-1)]
      prot_block.add(prot)
    return prot_block

  def random_peptide_block(self, pep_nodes, block_size):
    ### pep_nodes is passed in to save the O(n) time of
    ### listing the peptides
    if block_size == 0: return set()

    pep_block = set()
    pep = pep_nodes[random.randint(0,len(pep_nodes)-1)]
    pep_block.add(pep)
    while len(pep_block) < block_size:
      pep = pep_nodes[random.randint(0,len(pep_nodes)-1)]
      ### this will only be fast when the block size is small
      while pep in pep_block:
        pep = pep_nodes[random.randint(0,len(pep_nodes)-1)]
      pep_block.add(pep)
    return pep_block

  def random_connected_protein_block(self, prot_nodes, block_size):
    pass

  def sample_new_protein_configuration(self, prot_configuration, prot_block):
    all_outcomes = list(itertools.product(*[(False, True) for p in prot_block]))

    all_prot_changes = []
    prop_posteriors = []
    for prot_outcome in all_outcomes:
      prot_changes = dict(list(zip(prot_block, prot_outcome)))
      new_likelihood = fido_network.likelihood_new_protein_configuration_relative_to_current(self, prot_configuration, prot_changes)
      new_prior = self.prob_new_protein_configuration_relative_to_current(prot_configuration, prot_changes)
      prop_posteriors.append(new_likelihood * new_prior)
      all_prot_changes.append(prot_changes)
    prop_constant = sum(prop_posteriors)

    #prop_constant = min(1e-9, prop_constant)
    return fido_network.sample_outcome_from_posteriors(list(zip(all_prot_changes, [pp/prop_constant for pp in prop_posteriors])))

  def sample_new_protein_and_peptide_configuration(self, prot_configuration, pep_configuration, prot_block, pep_block):
    all_prot_outcomes = list(itertools.product(*[(False, True) for p in prot_block]))
    all_pep_outcomes = list(itertools.product(*[(False, True) for p in pep_block]))

    all_prot_and_pep_changes = []
    prop_posteriors = []
    for prot_outcome in all_prot_outcomes:
      for pep_outcome in all_pep_outcomes:
        prot_changes = dict(list(zip(prot_block, prot_outcome)))
        pep_changes = dict(list(zip(pep_block, pep_outcome)))
        new_likelihood = fido_network.likelihood_new_protein_and_peptide_configuration_relative_to_current(self, prot_configuration, pep_configuration, prot_changes, pep_changes)
        new_prior = self.prob_new_protein_configuration_relative_to_current(prot_configuration, prot_changes)
        prop_posteriors.append(new_likelihood * new_prior)
        all_prot_and_pep_changes.append((prot_changes, pep_changes))

    prop_constant = sum(prop_posteriors)
    return fido_network.sample_outcome_from_posteriors(list(zip(all_prot_and_pep_changes, [pp/prop_constant for pp in prop_posteriors])))

  @staticmethod
  def sample_outcome_from_posteriors(changes_and_posteriors):
    sample_value = random.uniform(0.0, 1.0)
    cumulative = 0.0
    for change, posterior in changes_and_posteriors:
      cumulative += posterior
      if cumulative >= sample_value: return change
    raise Exception('random real value is greater than all cumulative values in a list that should sum to 1.0. sample_value = ' + str(sample_value) + ', highest cumulative = ' + str(cumulative))

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
  def likelihood_protein_and_peptide_configuration(connected_fido_network, protein_config, peptide_config):
    ### go through the peptide terms and get Pr(E=e | R=r) and
    ### Pr(D_i | E=e)
    pep_nodes = [ peps for peps in connected_fido_network if peps['type'] == 'peptide' ]
    likelihood = 1.0
    for pep in pep_nodes:
      likelihood *= fido_network.likelihood_term(connected_fido_network, pep, peptide_config[pep], protein_config)

    return likelihood

  @staticmethod
  def likelihood_new_protein_configuration_relative_to_current(connected_fido_network, old_protein_config, protein_changes):
    changed_pep_ids = set()
    for prot in protein_changes:
      changed_pep_ids.update(set([ pep for pep in connected_fido_network[prot] ]))
    relative_likelihood = 1.0

    ### build the updated protein configuration; this can be done
    ### by only copying the previous state of the proteins that
    ### share peptides with the changed protein set. This is more
    ### efficient than copying the entire old protein
    ### configuration map.
    interesting_proteins = set()
    for pep in changed_pep_ids:
      interesting_proteins.update( set( connected_fido_network.predecessors(pep) ) )
    new_protein_config = {}
    for prot in interesting_proteins:
      if prot in protein_changes:
        new_protein_config[prot] = protein_changes[prot]
      else:
        new_protein_config[prot] = old_protein_config[prot]

    for pep in changed_pep_ids:
      original_likelihood_term = sum([ fido_network.likelihood_term(connected_fido_network, pep, pep_state, old_protein_config) for pep_state in (True,False)])
      new_likelihood_term = sum([ fido_network.likelihood_term(connected_fido_network, pep, pep_state, new_protein_config) for pep_state in (True,False) ])
      relative_likelihood *= new_likelihood_term / original_likelihood_term

    return relative_likelihood

  @staticmethod
  def likelihood_new_protein_and_peptide_configuration_relative_to_current(connected_fido_network, old_protein_config, old_peptide_config, protein_changes, peptide_changes):
    changed_pep_nodes = set()
    ### add the peptides changed by protein changes to the
    ### interesting peptide set
    for prot in protein_changes:
      changed_pep_nodes.update(set([ pep for pep in connected_fido_network[prot] ]))

    ### add the deliberately changed peptides to the interesting
    ### peptide set
    changed_pep_nodes.update(set(peptide_changes))

    relative_likelihood = 1.0

    ### build the updated protein configuration; this can be done
    ### by only copying the previous state of the proteins that
    ### share peptides with the changed protein set. This is more
    ### efficient than copying the entire old protein
    ### configuration map.
    interesting_proteins = set()
    for pep in changed_pep_nodes:
      interesting_proteins.update( set( connected_fido_network.predecessors(pep) ) )
    new_protein_config = {}
    for prot in interesting_proteins:
      if prot in protein_changes:
        new_protein_config[prot] = protein_changes[prot]
      else:
        new_protein_config[prot] = old_protein_config[prot]

    ### do the same for the new peptide configuration
    for pep in changed_pep_nodes:
      original_peptide_state = old_peptide_config[pep]
      ### if this peptide has changed, use its new value
      ### otherwise, use the old value
      if pep in peptide_changes:
        new_peptide_state = peptide_changes[pep]
      else:
        new_peptide_state = original_peptide_state

      original_likelihood_term = fido_network.likelihood_term(connected_fido_network, pep, original_peptide_state, old_protein_config)
      new_likelihood_term = fido_network.likelihood_term(connected_fido_network, pep, new_peptide_state, new_protein_config)

      original_likelihood_term = max(1e-30, original_likelihood_term)
      relative_likelihood *= new_likelihood_term / original_likelihood_term

    return relative_likelihood

  @staticmethod
  def prob_protein_configuration(connected_fido_network, protein_config):
    n = len(protein_config)
    r = len([ (prot,state) for prot,state in list(protein_config.items()) if state == True ])
    return pow(iid_probability_node.gamma, r) * pow(1-iid_probability_node.gamma, n-r)

  def prob_new_protein_configuration_relative_to_current(self, old_protein_config, protein_changes):
    return fido_network.prob_protein_configuration(self, protein_changes)

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
        likelihood = fido_network.likelihood_protein_configuration(sg, outcome_map)
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
  def marginalize_over_set(potential_table_lst, relevant_variables_fs, do_not_marginalize_variables_fs):
    ### marginalize all relevant variables that are not
    ### dissallowed on the distribution composed of the products
    ### of the provided tables in the list. return the resulting
    ### distribution (as a table).
    marg_out_fs = relevant_variables_fs.difference(do_not_marginalize_variables_fs)
    remaining_fs = relevant_variables_fs.difference(marg_out_fs)
    marg_out_lst = list(relevant_variables_fs.difference(do_not_marginalize_variables_fs))
    remaining_lst = list(relevant_variables_fs.difference(marg_out_fs))
    ### build a table of remaining_fs outcomes to the probability with marg_out_fs marginalized out

    marg_out_outcomes = list(itertools.product(*[(False, True) for p in marg_out]))
    remaining_outcomes = list(itertools.product(*[(False, True) for p in remaining]))
    
    remaining_table = {}
    total = 0.0
    for remaining_outcome in remaining_outcomes:
      remaining_vars_and_values = list(zip(remaining_lst, remaining_outcome))
      total_for_remaining_outcome = 0.0
      for marg_out_outcome in marg_out_outcomes:
        marg_out_vars_and_values = list(zip(marg_out_lst, marg_out_outcome))
        outcome_vars_and_values = tuple(remaining_vars_and_values + marg_out_vars_and_values)

        potential_product = 1.0
        for table in potential_table_lst:
          potential_value = table[outcome_vars_and_values]
          potential_product *= potential_value

        total_for_remaining_outcome += potential_product
      total += total_for_remaining_outcome
      remaining_table[remaining_vars_and_values] = total_for_remaining_outcome

    ### normalize so that all probabilities sum to 1.0 (some of
    ### the tables used as inputs may be likelihoods)
    for remaining_vars_and_values in remaining_table:
      remaining_table[remaining_vars_and_values] /= total

    return remaining_table
                
  def connected_junction_tree_to_sufficiently_parameterized_nodes(self, jt):
    jt_node_to_sufficiently_parameterized_node_ids = multi_dict([(fs, []) for fs in jt])
    assigned_nodes = set()
    for fs in jt.nodes():
      ### find the potential functions encompassed by the set fs
      for n in self:
        ### if the incoming edge nodes are a subset of fs (and
        ### the node takes any arguments), then this node
        ### potential belongs to this junction tree node
        node_and_pred = set(self.predecessors(n))
        node_and_pred.add(n)
        if len( node_and_pred.difference(fs) ) == 0:
          if not n in assigned_nodes:
            jt_node_to_sufficiently_parameterized_node_ids.add(fs, n)
            assigned_nodes.add(n)
    return jt_node_to_sufficiently_parameterized_node_ids

  def node_id_to_potential_table(self):
    node_id_to_potential_table = {}
    for n in self:
      n_node = self.get_node_from_id(n)
      pred_nodes = [ self.get_node_from_id(p) for p in self.predecessors(n) ]
      relevant_node_ids_to_possible_outcomes = dict([(p.id, p.outcomes) for p in pred_nodes] + [(n, n_node.outcomes)])
      n_potential_distribution = distribution(relevant_node_ids_to_possible_outcomes)

      for vars_and_outcomes in n_potential_distribution.get_all_possible():
        vars_and_outcomes_dict = dict(vars_and_outcomes)
        outcome_prob = n_node.wrapped_potential(vars_and_outcomes_dict)
        n_potential_distribution[vars_and_outcomes] = outcome_prob
      
      node_id_to_potential_table[n] = n_potential_distribution

    return node_id_to_potential_table

  def connected_junction_tree_belief_propagation(self, jt, node_label_map):
    jt_node_to_sufficiently_parameterized_node_ids = self.connected_junction_tree_to_sufficiently_parameterized_nodes(jt)
    node_id_to_potential_table = self.node_id_to_potential_table()

    jt_node_to_potential_table = {}
    for fs in jt:
      jt_node_to_potential_table[fs] = distribution.multiply( [ node_id_to_potential_table[n] for n in jt_node_to_sufficiently_parameterized_node_ids[fs] ] )

    edges_to_marginals = {}
    logger.info('\tcollecting')
    self.connected_junction_tree_collect_messages(jt, jt_node_to_potential_table, edges_to_marginals, jt.nodes()[0])
    jt_nodes_to_true_marginals = {}
    logger.info('\tdistributing')
    self.connected_junction_tree_distribute_messages(jt, jt_node_to_potential_table, edges_to_marginals, jt_nodes_to_true_marginals, jt.nodes()[0])
    return jt_nodes_to_true_marginals

  def connected_junction_tree_collect_messages(self, jt, jt_node_to_potential_table, edges_to_marginals, fs, parent_fs = None, visited_nodes = set()):
    logger.info('\t\tcollecting on {} with len {}'.format(node_label_map[fs], len(fs)))
    incoming_messages = []
    for m in jt.neighbors(fs):
      if m in visited_nodes: continue
      message_from_m = self.connected_junction_tree_collect_messages(jt, jt_node_to_potential_table, edges_to_marginals, m, fs, visited_nodes.union( set( [fs] ) ) )
      incoming_messages.append(message_from_m)

      edges_to_marginals[(m, fs)] = message_from_m

    if parent_fs != None:
      ### when the parent is not None, calculate the message
      ### from fs to parent_fs
      message_up = jt_node_to_potential_table[fs]
      if len(incoming_messages) > 0:
        incoming_messages_distribution = distribution.multiply(incoming_messages)
        message_up = message_up * incoming_messages_distribution

      ### marginalize out the unecessary variables from the message
      message_up = message_up.marginalized_out(fs.difference(parent_fs))
      return message_up

  def connected_junction_tree_distribute_messages(self, jt, jt_node_to_potential_table, edges_to_marginals, jt_nodes_to_true_marginals, fs, visited_nodes = set()):
    logger.info('\t\tdistributing on {} with len {}'.format(node_label_map[fs], len(fs)))
    ### all incoming edges should be available before this
    ### function is called; therefore, you can find the marginal
    ### distribution of this node
    incident_incoming_edge_set = set([(m, fs) for m in jt.neighbors(fs)])

    fs_density = jt_node_to_potential_table[fs] * distribution.multiply( [ edges_to_marginals[e] for e in incident_incoming_edge_set ] )
    jt_nodes_to_true_marginals[fs] = fs_density

    ### pass messages to all non-visited nodes
    for m in jt.neighbors(fs):
      if m in visited_nodes: continue
      message_to_m = fs_density / edges_to_marginals[(m,fs)]
      message_to_m = message_to_m.marginalized_out(fs.difference(m))
      edges_to_marginals[(fs, m)] = message_to_m
      self.connected_junction_tree_distribute_messages(jt, jt_node_to_potential_table, edges_to_marginals, jt_nodes_to_true_marginals, m, visited_nodes.union(set([fs])))

  def JunctionTreeMarginalization_Inference(self):
    ### build the junction tree
    protein_posteriors = {}
    for sg_number,sg in enumerate(nx.weakly_connected_component_subgraphs(self, copy=False)):
      #if sg_number != 27:
      #  continue

      moralized = to_moralized(sg)
      jt = moralized_to_junction_tree(moralized)
      
      if log_connected_naive_complexity(sg) > 12:
        logger.info('\tJT inference on subgraph {} with {} proteins, treewidth {}'.format(sg_number, log_connected_naive_complexity(sg), treewidth(jt)))
        #continue

      jt_unique_label_map = {}
      for i,fs in enumerate(jt):
        jt_unique_label_map[fs] = i

      ### perform belief propagation
      jt_nodes_to_true_marginals = sg.connected_junction_tree_belief_propagation(jt, node_label_map = jt_unique_label_map)
      
      subgraph_protein_posteriors = dict.fromkeys( [ prot_n for prot_n in sg.nodes() if prot_n['type'] == 'protein' ], 0.0)

      computed_posteriors = set()

      for fs in jt:
        local_dist = jt_nodes_to_true_marginals[fs]
        tot = 0.0
        for o in local_dist.get_all_possible():
          prob = local_dist[o]
          tot += prob
          for n_id, value in o:
            ### some junction tree nodes will share
            ### variables; use the posterior computed the
            ### first time they appear; afterward, ignore.
            if n_id in computed_posteriors: continue
            if n_id['type'] == 'protein':
              if value == True: subgraph_protein_posteriors[n_id] += prob
        ### normalize
        for n_id in local_dist.domain_vars_to_outcomes:
          if n_id in computed_posteriors: continue
          if n_id['type'] == 'protein':
              subgraph_protein_posteriors[n_id] /= tot

        computed_posteriors.update(fs)
      protein_posteriors.update(subgraph_protein_posteriors)

    return protein_posteriors

  @staticmethod
  def prune_low_scoring_peptides(dg, score_threshold = 5e-2):
    # go through the nodes with peptide type and call prune_node if
    # the probability is below the threshold
    to_prune_list = []
    for n in dg:
      if n['type'] == 'peptide':
        spectrum_id = dg.successors(n)[0]
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
    pep_score, pep_to_prune = min(list(zip(pep_present_likelihoods, pep_nodes)))
    fido_network.prune_peptide(connected_dg, pep_to_prune)

  @staticmethod
  def dynamic_treewidth_pruned(dg, treewidth_threshold, **Kwargs):
    treewidth_threshold = float(treewidth_threshold)
    dg_lst = nx.weakly_connected_component_subgraphs(dg)
    pruned_dg_lst = fido_network.dynamic_treewidth_pruned_helper(dg_lst, treewidth_threshold)
    return all_bayesian_network_unions(pruned_dg_lst)
        
  @staticmethod
  def dynamic_treewidth_pruned_helper(dg_lst, treewidth_threshold):
    ### may be quite slow-- tree decomposition is built multiple
    ### times
    result_lst = []
    for dg in dg_lst:
      if log_connected_naive_complexity(dg) < connected_protein_threshold:
        result_lst.append(dg)
      else:
        while True:
          # prune the lowest-scoring peptide until the graph is split
          fido_network.connected_lowest_prune(dg)
          subgraphs = nx.weakly_connected_component_subgraphs(dg)
          if len(subgraphs) > 1:
            pruned_subgraphs = fido_network.dynamic_pruned_helper(subgraphs, connected_protein_threshold)
            result_lst.extend(pruned_subgraphs)
            break
    return result_lst

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
        logger.info('pruning connected graph with {} proteins'.format(log_connected_naive_complexity(dg)))
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
    return dict.fromkeys(list(range(-1,10)),0.9999999)

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

        #likelihood_given_present = peptide_probability / cp
        #likelihood_given_absent = (1-peptide_probability) / (1-cp)
        likelihood_given_present = peptide_probability
        likelihood_given_absent = (1-peptide_probability)
      
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
    charge_priors = dict.fromkeys(list(range(-1,10)),0.9999999)

    # get edge info
    # each element of this list is a tuple describing the PSM
    peptide_spectrum_proteins_prob_charge = []

    sequence = ''
    prot_list = []
    prob = 0.0
    charge = -1 # dummy charge for now
    spectrum_count = 0

    for i in range(0, df.shape[0]):
      prots = df['Proteins'].loc[i]
      if pd.isnull(prots): continue

      sequence = df['Sequence'].loc[i]
      prot_list = [str.join('_', prot.split('|')[0:2]) for prot in prots.split(';')]
      prob = 1 - df['PEP'].loc[i]

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
      sub_fn.remove_all_but_maximum_likelihood_spectrum()
      self = bayesian_network_union(self, sub_fn)
    
    logger.info('prots: {}'.format(len([ n for n in self if n['type'] == 'protein' ])))
    logger.info('peps: {}'.format(len([ n for n in self if n['type'] == 'peptide' ])))


def pd_main():

  # initialize logger
  init_logger(True, 'fido.log', log_to_file=False)

  parameter_map = { \
    'marginalization_mode': 'fido', 
    'filetype': 'pivdo', 
    'connected_protein_threshold': 14 }

  fn = fido_network(**parameter_map)
  
  logger.info('reading...')
  df = pd.read_csv('/gd/MS/SCoPE/SQC/SQC77/evidence.txt', sep='\t', low_memory=False)
  fn.load_from_dataframes([df], **parameter_map)

  logger.info('trimming down to best spectrum match...')
  fn.remove_all_but_maximum_likelihood_spectrum()

  if 'group_proteins' in parameter_map:
    logger.info('clustering proteins')
    fn.cluster_proteins()

  logger.info('pruning...')

  #fn = fido_network.prune_low_scoring_peptides(fn, 1e-2)
  fn = fido_network.dynamic_pruned(fn, **parameter_map)

  fn.init_parameter_names_to_nodes_maps()

  if parameter_map['marginalization_mode'] == 'fido':
    logger.info('fido inference')
    prots_to_posteriors = fn.FidoMarginalization_Inference()
  elif parameter_map['marginalization_mode'] == 'gibbs':
    logger.info('gibbs inference')
    prots_to_posteriors = fn.Gibbs_Inference()
  elif parameter_map['marginalization_mode'] == 'gibbs_collapsed':
    logger.info('gibbs_collapsed inference')
    prots_to_posteriors = fn.Gibbs_Collapsed()
  elif parameter_map['marginalization_mode'] == 'junction_tree':
    logger.info('junction_tree inference')
    prots_to_posteriors = fn.JunctionTreeMarginalization_Inference()
  else:
    logger.info('Error: unrecognized marginalization_mode = {}'.format(parameter_map['marginalization_mode']))

  pprint(prots_to_posteriors)

  pass

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
    'connected_protein_threshold': 100 }

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

  #import pycallgraph
  #pycallgraph.start_trace()

  #prots_to_posteriors = fn.CollapsedGibbs_Inference()
  if parameter_map['marginalization_mode'] == 'fido':
    print('fido inference')
    prots_to_posteriors = fn.FidoMarginalization_Inference()
  elif parameter_map['marginalization_mode'] == 'gibbs':
    print('gibbs inference')
    prots_to_posteriors = fn.Gibbs_Inference()
  elif parameter_map['marginalization_mode'] == 'gibbs_collapsed':
    print('gibbs_collapsed inference')
    prots_to_posteriors = fn.Gibbs_Collapsed()
  elif parameter_map['marginalization_mode'] == 'junction_tree':
    print('junction_tree inference')
    prots_to_posteriors = fn.JunctionTreeMarginalization_Inference()
  else:
    print('Error: unrecognized marginalization_mode = ', parameter_map['marginalization_mode'])

  #print(prots_to_posteriors)
  pprint(prots_to_posteriors)

  #pycallgraph.make_dot_graph('test.jpg', format='jpg', tool='neato')

def profile_main(argv):
  # This is the main function for profiling
  # We've renamed our original main() above to real_main()
  import cProfile, pstats
  prof = cProfile.Profile()
  prof = prof.runctx("real_main(argv)", globals(), locals())
  print("<pre>")
  stats = pstats.Stats(prof)
  stats.sort_stats("time")  # Or cumulative
  stats.print_stats(80)  # 80 = how many to print
  # The rest is optional.
  # stats.print_callees()
  # stats.print_callers()
  print("</pre>")

  # import cProfile, pstats, StringIO

  # prof = cProfile.Profile()

  # prof = prof.runctx("real_main(argv)", globals(), locals())

  # stats = pstats.Stats(prof, stream=stream)
  # stats.sort_stats("time")  # Or cumulative
  # stats.print_stats(80)  # 80 = how many to print
  # # The rest is optional.
  # # stats.print_callees()
  # # stats.print_callers()
  # print "Profile data:\n%s" % stream.getvalue()

main = real_main

if __name__ == '__main__':
  main(sys.argv[1:])


