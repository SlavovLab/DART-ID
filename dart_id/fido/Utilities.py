# coding: utf-8

import numpy as np

def index_max(lst):
  return lst.index(max(lst))

# used for one to many mapping
class multi_dict(dict):
  def add(self, key, item):
    if key in self: self[key].append(item)
    else:           self[key] = [ item ]

  # performs __getitem__ and if the resulting list has only one
  # element, returns that element not in a list
  def get_flattened(self, key):
    result = self[key]
    # if there is only one item for that key, only return that
    # item. otherwise, return a set containing those items
    if len(result) == 1: return result[0]
    else:                return result

class counting_dict(dict):
  def add(self, item):
    if item in self: self[item] += 1
    else:            self[item] = 0

class hashable_dict:
  def __init__(self, d):
    self.my_dict = d
    self.my_frozenset = frozenset(list(d.items()))
  def __getitem__(self, item):
    return self.my_dict[item]
  def __hash__(self):
    return hash(self.my_frozenset)
  def __eq__(self, rhs):
    return isinstance(rhs, hashable_dict) and self.my_frozenset == rhs.my_frozenset
  def __ne__(self, rhs):
    return not self == rhs
  def __str__(self):
    return 'hashable_dict(' + str(self.my_dict) + ')'
  def __repr__(self):
    return self.__str__()

def log_add(logA,logB):
  if logA == None: return logB
  if logA < logB:  return log_add(logB,logA)
  return np.log2( 1 + pow(2, logB-logA) ) + logA

def log_sum(lst):
  log_total = None
  for v in lst: log_total = log_add(log_total, v)
  return log_total

def prod(lst):
  result = 1.0
  for i in lst: result *= i
  return result

def binomial(n, r):
  ### n! /( r! (n-r)! )
  return prod(list(range(r+1, n+1))) / prod(2, n-r+1)

