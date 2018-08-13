# coding: utf-8

import itertools
import networkx as nx
import dart_id.fido.BinomialHeap

from dart_id.fido.Utilities import *

def elimination_edges(g, n):
  return list(itertools.combinations(g.neighbors(n), 2))
def new_elimination_edges(g, n):
  elim_edges = elimination_edges(g,n)
  added_edges = []
  for a,b in elim_edges:
    if b not in g[a]:
      added_edges.append((a,b))
  return added_edges
def number_of_added_elimination_edges(g,n):
  ### note that it is important to make sure to not count edge (a,b)
  ### different from (b,a)
  edge_set = elimination_edges(g, n)
  total_edge_count = len(edge_set)
  existing_edge_count = 0
  for a,b in edge_set:
    if b in g[a]:
      existing_edge_count += 1
  return total_edge_count - existing_edge_count

def to_moralized(g):
  moral_graph = g.to_undirected()
  for n in g.node:
    ### this is not the most efficient way to perform this or
    ### to visualize. all predecessors of n will form a
    ### complete graph. this graph can be represented as a
    ### group containing all of these nodes.
    p_list = list(g.predecessors(n))
    for pair in itertools.combinations( p_list,  2):
      moral_graph.add_edge(*pair)
  return moral_graph

def to_triangulated(g):
  ### this method is O(n^2) and can be made more efficient with a
  ### heap. triangulating a graph consisting of two disjoint graphs is
  ### the same as triangulating each individually, so do that for now.
  return fast_to_triangulated(g)

  triangulated = nx.Graph()
  subgraphs = nx.connected_component_subgraphs(g)
  print('\t\t\ttriangulating subgraphs...')
  triangulated_subgraphs = [connected_to_triangulated(sg) for sg in subgraphs]
  print('\t\t\tassembling subgraphs into graph union...')
  triangulated = all_graph_unions(triangulated_subgraphs)
  return triangulated

def fast_to_triangulated(g):
  ### this method is O(n^2) and can be made more efficient with a
  ### heap. triangulating a graph consisting of two disjoint graphs is
  ### the same as triangulating each individually, so do that for now.
  triangulated = nx.Graph()
  subgraphs = nx.connected_component_subgraphs(g)
  #print('\tlooking into subgraph triangulation:')
  #for sg in subgraphs:
  #  print('\t\tsubgraph size:',len(sg))
  #print('\t\t\ttriangulating subgraphs...')
  triangulated_subgraphs = [fast_connected_to_triangulated(sg) for sg in subgraphs]
  #print('\t\t\tassembling subgraphs into graph union...')
  triangulated = all_graph_unions(triangulated_subgraphs)
  return triangulated

### does minimum neighbors
def fast_connected_to_triangulated(g):
  added_edges = []
  local = g.copy()
  
  nodes_and_number_of_neighbors = [ (len(local[n]),n) for n in local ]
  #print(str(nodes_and_number_of_neighbors))
  nodes_to_heap_references = {}
  bh = BinomialHeap.heap()
  for number_neighbors,n in nodes_and_number_of_neighbors:
    #print 'heap inserting',number_neighbors,n
    nodes_to_heap_references[n] = bh.insert(number_neighbors,n)
  for removal_index, n in enumerate(bh):
    #print('\t\tremoving',n,'with',nodes_to_heap_references[n].ref.key,'neighbors')
     
    # if removal_index == 1:
    #     for n in bh:
    #         print '\t\t\tresult:', n, 'with',nodes_to_heap_references[n].ref.key,'neighbors'
    #     nx.draw_graphviz(local)
    #     P.show()

    edges = new_elimination_edges(local, n)

    ### count the edges actually added incident to each node
    nodes_to_added_edges = counting_dict( [ (nei,0) for nei in local.neighbors(n) ] )

    for a,b in edges:
      nodes_to_added_edges.add(a)
      nodes_to_added_edges.add(b)

    ### decrease edges that would be eliminated for node a by the
    ### number of edges that have actually been added this
    ### iteration (note that it is faster to do this in a batch
    ### manner because each decrease_key operation is log(n))
    for a in nodes_to_added_edges:
      x = nodes_to_heap_references[a]
      #print('\tchanging',a,':',x.ref.key,'+',nodes_to_added_edges[a],'- 1')
      ### store the old number of neighbors
      old_num_neighbors = x.ref.key
      ### remove the node by decreasing its key lower than any other
      x.delete()
      nodes_to_heap_references[a] = bh.insert(old_num_neighbors + nodes_to_added_edges[a] - 1,a)
  
  local.add_edges_from(edges)
  added_edges.extend(edges)
  local.remove_node(n)

  #to_min = [ (len(local[n]), n) for n in local ]
  #print(str(to_min))
  #if len(to_min) > 0:
  #  print('min',min(to_min))

  result = g.copy()
  return result

### does minimum fill-in
def connected_to_triangulated(g):
  added_edges = []
  local = g.copy()
  while len(local) > 0:
    ### count the edges that would be added by each node
    ### elimination; eliminate the node that will result in the
    ### fewest
    number_of_added_edges_for_nodes = [ (number_of_added_elimination_edges(local,n), n) for n in local ]
    number_of_edges, n = min(number_of_added_edges_for_nodes)
    #print('removing',n,'with',number_of_edges,'edges')

    edges = elimination_edges(local, n)
    local.add_edges_from(edges)
    added_edges.extend(edges)
    local.remove_node(n)
  result = g.copy()
  result.add_edges_from(added_edges)
  return result

def to_clique_graph(g):
  cliques = [ c for c in nx.find_cliques(g) ]
  clique_graph = nx.make_max_clique_graph(g)
  ### networkx.clique_graph does not tell which nodes it contains
  
  frozenset_clique_graph = nx.Graph()
  for a in clique_graph:
    a_clique = frozenset(cliques[a-1])
    frozenset_clique_graph.add_node(a_clique)
    for b in clique_graph[a]:
      b_clique = frozenset(cliques[b-1])
      ### weight it with the negative cardinality of the
      ### intersection so that you can use the built in
      ### minimum_spanning_tree function to get the maximum
      ### spanning tree
      overlap_size = len(a_clique.intersection(b_clique))
      frozenset_clique_graph.add_edge( a_clique, b_clique, weight = -overlap_size)

  return frozenset_clique_graph

def moralized_to_junction_tree(moral_graph):
  junction_tree = nx.Graph()
  #print('\t\ttriangulating...')
  triangulated = to_triangulated(moral_graph)
  #print)'\t\tbuilding clique graph...')
  clique_graph = to_clique_graph(triangulated)
  #print('\t\tfinding max spanning tree...')
  junction_tree = nx.minimum_spanning_tree(clique_graph)
  return junction_tree

def verify_junction_forest(junction_forest):
  for jt in nx.connected_component_subgraphs(junction_forest):
      verify_junction_tree(jt)

def verify_junction_tree(junction_tree):
  ### uses my own method for verification

  ### Cut each edge L--R. The tree is only a valid junction tree if the
  ### set of nodes in the left tree does not contain the nodes removed
  ### by moving from R to L (and vice-versa).

  ### Intuitive proof: Each running intersection of any element must
  ### have a boundary (or is ubiquitous throughout the tree). Past the
  ### boundary, the element cannot be featured.
  local_jt = junction_tree.copy()
  for sL, sR in junction_tree.edges():
    local_jt.remove_edge(sL,sR)
    subtrees = nx.connected_component_subgraphs(local_jt)

    ### if this is a proper tree, then it should have only two
    ### connected subgraphs
    if len(subtrees) != 2:
      raise Exception('Error: verify junction tree was not called using a connected tree')
    
    subtreeL = subtrees[0]
    subtreeR = subtrees[1]
    if sR in subtreeL:
      subtreeL, subtreeR = subtreeR, subtreeL

    newR = sR.difference(sL)
    newL = sL.difference(sR)
    if len(newR.intersection(subtreeL)) > 0 or len(newL.intersection(subtreeR)) > 0:
      print('verify_junction_tree: cutting edge',sL,'--',sR,'results in overlapping edges')
      return
    local_jt.add_edge(sL,sR)

def graph_union(g1, g2):
  if g1 == None:
    return g2
  if g2 == None:
    return g1
  for n in g2:
    g1.add_node(n, **g2.node[n])
    for m in g2[n]:
      g1.add_edge(n,m,**g2[n][m])
  return g1

def all_graph_unions(graph_lst):
  ### do not set as Graph or DiGraph so that it works with any
  ### homogeneous list
  result = None
  for g in graph_lst:
    result = graph_union(result, g)
  return result

def treewidth(junction_tree):
  return max([len(fs) for fs in junction_tree])

def type_subgraph(g, t):
  type_nodes = [ n for n in g if n['type'] == t ]
  return g.subgraph(type_nodes)

def log_connected_naive_complexity(dg):
  return len(type_subgraph(dg, 'protein'))

def log_connected_naive_junction_tree_complexity(junction_tree):
  var_sets = [fs for fs in junction_tree]
  res = set()
  for vs in var_sets:
    res.update(vs)
  return len(res)

def log_junction_tree_complexity(junction_tree):
  return log_sum([ len(fs) for fs in junction_tree ])

