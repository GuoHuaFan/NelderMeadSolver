"""
  Filename    : nedlerMead.py
  Author      : Soumyaroop Roy
  Date        : May 02, 2016
  Description : Nedler-Mead method using numpy
  Algorithm   : https://en.wikipedia.org/wiki/Nelder-Mead_method#One_possible_variation_of_the_NM_algorithm
"""

import sys
import math
import numpy as np

""" Debug level """
_dbgLvl = 0

def setDebugLevel (x_dbgLvl):
  global _dbgLvl
  _dbgLvl = x_dbgLvl

"""
  Function to print an error message

  @param  Error string

  @return None
"""
def _error (x_str):
  print ("ERROR: %s" %x_str)
  sys.exit()

"""
  Function to print a debug message

  @param  Debug string
  @param  Debug level

  @return None
"""
def _debugPrint (x_str, x_dbgLvl=1):
  global _dbgLvl
  if x_dbgLvl <= _dbgLvl:
    print ("Debug (Level %u): %s" %(x_dbgLvl, x_str))

""" Simplex class """

"""
  Member variables key:
    m_dim - Space dimensions
    m_fn  - The function

    m_arr - Simplex points
    m_sorted - Indices of m_arr sorted by f(x) values

    m_x0   - Simplex centroid   (x_0)
    m_xr   - Reflected point    (x_r)
    m_xe   - Expanded point     (x_e)
    m_x1   - Best point         (x_1)
    m_xn   - Second worst point (x_n)
    m_xnp1 - Worst point        (x_{n+1})

    m_coeffR - Reflection coefficient (\alpha)
    m_coeffE - Expansion coefficient (\gamma)
    m_coeffC - Contraction coefficient (\rho)
    m_coeffS - Shrink coefficient (\sigma)
"""

class ConvergenceError(Exception):
  def __init__ (self, x_message):
    self.message = x_message;

class Simplex ():
  def __init__ (self, x_dim, x_fn, x_solverCoeffs, x_size):
    self.m_dim = x_dim
    self.m_fn  = x_fn
    self.setSolverCoeffs(x_solverCoeffs)
    self.m_arr = np.empty([x_dim+1, x_dim])
    self.m_size = x_size

  def setSolverCoeffs (self, x_solverCoeffs):
    self.m_coeffR, self.m_coeffE, self.m_coeffC, self.m_coeffS = x_solverCoeffs.getCoeffs()

  def anchor (self, x_start):
    """ Initialize simplex to x_start """
    for i in xrange(len(self.m_arr)):
      self.m_arr[i, :] = np.array([x_start])
      if i == 0:
        continue
      self.m_arr[i, i - 1] += self.m_size

  def printf (self):
    print "Simplex details:"
    print "  Anchor: %s" %self.m_arr[0]
    print "  Size: %u" %self.m_size
    print "  Points: ["
    for i in self.m_sorted:
      print "    %s, %.2f" %(self.m_arr[i], self.m_fn(self.m_arr[i]))
    print "  ]"

  """ Sort x_i using key f(x_i) """
  """ FIXME: Is there a numpy version of the standard python 'sort()' that
             takes in a key as well? """
  """
    @returns (x1, fx1), (xn, fxn), (xnp1, fxnp1)
             Best, Second worst, Worst so far
  """
  def order (self):
    """ Sort """
    l_fnVal = np.zeros(len(self.m_arr))
    for i in xrange(len(self.m_arr)):
      l_fnVal[i] = self.m_fn(self.m_arr[i])
    self.m_sorted = np.argsort(l_fnVal)

    """ Compute the centroid """
    self.computeX0()

    """ Return best, second worst, and worst points and corresponding
        function values """
    self.m_x1   = self.m_arr[self.m_sorted[0]]
    self.m_xn   = self.m_arr[self.m_sorted[-2]]
    self.m_xnp1 = self.m_arr[self.m_sorted[-1]]

    return ((self.m_x1, self.m_fn(self.m_x1)),
            (self.m_xn, self.m_fn(self.m_xn)),
            (self.m_xnp1, self.m_fn(self.m_xnp1)))

  """ Compute the centroid of points all x in the simplex for which 
      f(x) < max(f(x): for all x in the simplex """
  def computeX0 (self):
    """ Add all the points """
    self.m_x0 = np.sum(self.m_arr, axis=0)
    """ Subtract the point for which f(x) is max """
    self.m_x0 -= self.m_arr[self.m_sorted[-1], :]
    """ Compute average """
    self.m_x0 /= self.m_dim
    _debugPrint("Centroid: %s" %self.m_x0)

  """ Reflected point: x_r = x_0 + \alpha * (x_0 - x_{n+1}) """
  def computeXr (self):
    self.m_xr = np.add(self.m_x0,
                       self.m_coeffR * (np.subtract(self.m_x0, self.m_xnp1)))
    _debugPrint("Reflected point: %s" %self.m_xr)
    return self.m_xr, self.m_fn(self.m_xr)

  """ Expanded point: x_e = x_0 + \gamma * (x_r - x_0) """
  def computeXe (self):
    self.m_xe = np.add(self.m_x0,
                       self.m_coeffE * (np.subtract(self.m_xr, self.m_x0)))
    _debugPrint("Expanded point: %s" %self.m_xe)
    return self.m_xe, self.m_fn(self.m_xe)

  """ Contracted point: x_c = x_0 + \rho * (x_{n+1} - x_0) """
  def computeXc (self):
    self.m_xc = np.add(self.m_x0,
                       self.m_coeffC * (np.subtract(self.m_xnp1, self.m_x0)))
    _debugPrint("Contracted point: %s" %self.m_xe)
    return self.m_xc, self.m_fn(self.m_xc)

  """ Shrink: x_i = x_1 + \sigma * (x_i - x_1) for all 1 < i <= n+1 """
  def shrink (self):
    _debugPrint("Shrinking...")
    for i in self.m_sorted[1:]:
      self.m_arr[i, :] = np.add(self.m_x1,
        self.m_coeffS * (np.subtract(self.m_arr[i], self.m_x1)))

  """ Overwrite the worst point """
  def xnp1 (self, x_new):
    self.m_arr[self.m_sorted[-1], :] = x_new
    self.m_sorted[-1] = -1

class SolverCoeffs ():
  def __init__ (self):
    self.m_alpha = 1.0
    self.m_gamma = 1.0
    self.m_rho   = 0.5
    self.m_sigma = 0.5

  def __init__ (self, x_coeffs):
    self.m_alpha, self.m_gamma, self.m_rho, self.m_sigma = x_coeffs

  def setCoeffs (self, x_coeffs):
    self.m_alpha, self.m_gamma, self.m_rho, self.m_sigma = x_coeffs

  def getCoeffs (self):
    return self.m_alpha, self.m_gamma, self.m_rho, self.m_sigma

class SolverParams ():
  def __init__ (self):
    """ Progress report interval """
    self.m_progressReportInterval = 10000
    
    """ Across iteration delta - used to decide convergence """
    self.m_convergenceDelta = 0.000000005
    
    """ Maximum number of iterations """
    self.m_maxIterations = 20000

    """ Simplex start size """
    self.m_simplexStartSize = 10

  def setProgressReportInterval (self, x_progressReportInterval):
    self.m_progressReportInterval = x_progressReportInterval
  
  def setConvergenceDelta (self, x_convergenceDelta):
    self.m_convergenceDelta = x_convergenceDelta
  
  def setMaxIterations (self, x_maxIterations):
    self.m_maxIterations = x_maxIterations

  def setSimplexStartSize (self, x_simplexStartSize):
    self.m_simplexStartSize = x_simplexStartSize

  def getProgressReportInterval (self):
    return self.m_progressReportInterval
  
  def getConvergenceDelta (self):
    return self.m_convergenceDelta
  
  def getMaxIterations (self):
    return self.m_maxIterations

  def getSimplexStartSize (self):
    return self.m_simplexStartSize

class NedlerMeadSolver ():
  def __init__ (self):
    self.m_solverParams = SolverParams()
    self.m_solverCoeffs = SolverCoeffs()

  def __init__ (self, x_coeffs):
    self.m_solverParams = SolverParams()
    self.m_solverCoeffs = SolverCoeffs(x_coeffs)

  def getSolverParams (self):
    return self.m_solverParams

  def getSolverCoeffs (self):
    return self.m_solverCoeffs

  def solve (self, x_fn, x_start):
    """ Construct the simplex """
    self.m_simplex = Simplex(len(x_start), x_fn, self.m_solverCoeffs,
                             self.m_solverParams.getSimplexStartSize())

    """ Anchor the simplex at a starting point; also orders the points """
    self.m_simplex.anchor(x_start)

    (l_x1, l_fx1), (l_xn, l_fxn), (l_xnp1, l_fxnp1) = self.m_simplex.order()

    #self.m_simplex.printf()

    l_iter = 0
    l_x1Prev = l_x1
    l_fx1Prev = l_fx1

    while l_iter < self.m_solverParams.getMaxIterations():

      """ Perform Reflection """
      l_xr, l_fxr = self.m_simplex.computeXr()
      if l_fxr >= l_fx1:
        if l_fxr < l_fxn:
          """ The reflected point IS better than the second worst point """
          self.m_simplex.xnp1(l_xr)
        else:
          """ The reflected point is NOT better than the second worst point """
          """ Perform Contraction """
          l_xc, l_fxc = self.m_simplex.computeXc()
          if (l_fxc < l_fxnp1):
            """ The contracted point is better than the worst point """
            self.m_simplex.xnp1(l_xc)
          else:
            """ The contracted point is NOT better than the worst point """
            """ Reduction/shrinkage """
            self.m_simplex.shrink()
      else:
        """ The reflected point is the best so far """
        """ Perform Expansion """
        l_xe, l_fxe = self.m_simplex.computeXe()
        if (l_fxe < l_fxr):
          """ The expanded point is better than the reflected point """
          self.m_simplex.xnp1(l_xe)
        else:
          """ The expanded point is NOT better than the reflected point """
          self.m_simplex.xnp1(l_xr)

      """ Simplex ordering; also computes the centroid """
      (l_x1, l_fx1), (l_xn, l_fxn), (l_xnp1, l_fxnp1) = self.m_simplex.order()

      """ Increment the number of iterations """
      l_iter += 1

      _debugPrint("End of iteration %u: Best f(%s) = %f, Worst f(%s) = %f."
                 %(l_iter,
                   l_x1,
                   l_fx1,
                   l_xnp1,
                   l_fxnp1
                   ))

      """ Print progress report """
      if l_iter % self.m_solverParams.getProgressReportInterval() == 0:
        print "Progress status: Iteration %u..." %l_iter

      """ Check terminating condition """
      if math.fabs(l_fxnp1 - l_fx1) <= self.m_solverParams.getConvergenceDelta():
        print "Algorithm converged (%u iterations)" %l_iter
        break;

    """ All done """
    if l_iter == self.m_solverParams.getMaxIterations():
      raise ConvergenceError("Algorithm did not converge. Try increasing the"
      + " number of iterations or starting simplex size")

    return (l_x1, l_fx1)

