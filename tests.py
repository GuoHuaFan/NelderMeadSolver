"""
  Filename    : tests.py
  Author      : Soumyaroop Roy
  Date        : August 23, 2016
  Description : Tests Nedler-Mead Convex Solver
"""

import sys
import nedlerMead
from nedlerMead import NedlerMeadSolver
from nedlerMead import ConvergenceError

class Test:
  def __init__ (self, x_name, x_fn, x_start, x_maxIter):
    self.m_name = x_name
    self.m_fn = x_fn
    self.m_start = x_start
    self.m_maxIter = x_maxIter

  def run (self):
    """ Sample coefficients """
    l_coeffs = (1,   #\alpha
                1,   #\gamma
                0.5, #\rho
                0.5  #\sigma
                )
  
    """ Construct a Nedler Mead optimizer """
    l_nm = NedlerMeadSolver(len(self.m_start), self.m_fn, l_coeffs)
  
    #nedlerMead.setMaxIterations(self.m_maxIter)
    #nedlerMead.setDebugLevel(1)
  
    """ Solve it """
    print "Solving %s function..." %self.m_name
    l_x1, l_fx1 = l_nm.solve(self.m_start)
    return True

def main():

  """ Functions """

  """ Rosenbrock's """
  testRb = Test(x_name="Rosenbrock's",
    x_fn=lambda x: 100 * ((x[1] - x[0] ** 2) ** 2) + (1 - x[0]) ** 2,
    x_start=[2000, -3000],
    x_maxIter=1000)

  """ Wood's """
  testW = Test(x_name="Wood's",
    x_fn=lambda x: ( 100 * (x[1] - x[0] ** 2) ** 2
                     +              (1 - x[0]) ** 2
                     + 90 * (x[3] - x[2] ** 2) ** 2
                     +              (1 - x[2]) ** 2
                     +       10.1 * (x[1] - 1) ** 2
                     +       10.1 * (x[3] - 1) ** 2
                     + 19.8 * (x[1] - 1) * (x[3] - 1)
                     ),
    x_start=[1000, -3000, 1020, -230],
    x_maxIter=6000)

  """ Powell's """
  testP = Test(x_name="Powell's",
    x_fn=lambda x: (  (x[0] + 10 * x[1]) ** 2
                     +  5 * (x[2] - x[3]) ** 2
                     +  (x[1] - 2 * x[2]) ** 4
                     + 10 * (x[0] - x[3]) ** 4
                     ),
    x_start=[1000, -3000, 1020, -230],
    x_maxIter=3000)

  """ Wayburn Seader 1 """
  testWs = Test(x_name="Wayburn Seader 1",
    x_fn=lambda x: ((x[0] ** 6 + x[1] ** 4 - 17) ** 2
                      +  (2 * x[0] + x[1] - 4) ** 2
                      ),
    x_start=[1000, -3000],
    x_maxIter=2000)

  """ Made-up function """
  testCubic = Test(x_name="Made-up Cubic",
    x_fn = lambda x: (x[0] * (x[0] - 1)
                     + (x[1] - 2) * (x[1] + 2)
                     + (x[2]) * (x[2] + 3)
                     ),
    x_start=[2000, -3000, 2],
    x_maxIter=1000)

  l_testsPassed = 0
  l_tests = [testRb, testW, testP, testWs, testCubic]
  print "Runnning test..."
  for i in xrange(len(l_tests)):
    try:
      print "Test %u:" %(i + 1)
      if l_tests[i].run():
        l_testsPassed += 1
    except ConvergenceError as e:
      print e.message
    except:
      print "Unexpected error: %s" %sys.exc_info()[0]
    finally:
      print

  print "Report: %u out of %u tests passed" %(l_testsPassed, len(l_tests))

if __name__ == "__main__":
    main()
