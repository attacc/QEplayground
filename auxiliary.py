# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import numpy as np
import os

def float_from_string(x):
  """
  Convert a string in a float 
  """
  y=[]
  for t in x.split():
    try:
      y.append(float(t))
    except ValueError:
      pass
  return y

class Path():
    """ Class that defines a path in the brillouin zone
    """
    def __init__(self,klist,intervals):
        """
        Generation of a path in reciprocal space by specifying a list of k-points
        """
        self.intervals = intervals
        klabels = []
        kpoints = []

        for kline in klist:
            kpoint, klabel = kline
            kpoints.append(kpoint)
            klabels.append(klabel)
        self.kpoints = np.array(kpoints)
        self.klabels = klabels

    def get_klist(self):
        """ 
        Output in the format of quantum espresso == [ [kx, ky, kz, 1], ... ]
        """
        kpoints = self.kpoints
        intervals = self.intervals
        kout  = np.zeros([sum(intervals)+1,4])
        kout[:,3] = 1
        io = 0
        for ik,interval in enumerate(intervals):
          for ip in range(interval):
            kout[io,:3] = kpoints[ik] + float(ip)/interval*(kpoints[ik+1] - kpoints[ik])
            io = io + 1
        kout[io,:3] = kpoints[ik] + float(ip+1)/interval*(kpoints[ik+1] - kpoints[ik])

        return kout

    def get_indexes(self):
        """ get the index of each point of the path
        """

        indexes = []
        index = 0
        for n,label in enumerate(self.intervals):
            indexes.append([index,self.klabels[n]])
            index += self.intervals[n] 
        indexes.append([index,self.klabels[-1]])
        return indexes

