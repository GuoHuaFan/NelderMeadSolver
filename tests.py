"""
  Filename    : tests.py
  Author      : Soumyaroop Roy
  Date        : August 23, 2016
  Description : Tests Nelder-Mead Convex Solver
"""

import sys
import datetime
import traceback
import math

import nelderMead
from nelderMead import NelderMeadSolver
from nelderMead import ConvergenceError

class Test:
  def __init__ (self, x_name, x_fn, x_start, x_maxIter, x_expectedResult):
    self.m_name = x_name
    self.m_fn = x_fn
    self.m_start = x_start
    self.m_maxIter = x_maxIter
    self.m_expectedResult = x_expectedResult

  def run (self):
    """ Sample coefficients """
    l_coeffs = (1,   #\alpha
                1,   #\gamma
                0.5, #\rho
                0.5  #\sigma
                )
  
    """ Set up the solver """
    l_nm = NelderMeadSolver(l_coeffs)
    l_nm.getSolverParams().setMaxIterations(self.m_maxIter)
    l_maxConvergenceDelta = 1.0E-8
    l_maxErrorTolerance = l_maxConvergenceDelta*100
    l_nm.getSolverParams().setConvergenceDelta(l_maxConvergenceDelta)
  
    """ Solve it """
    print "Solving %s function..." %self.m_name
    l_start = datetime.datetime.now()
    l_x1, l_fx1 = l_nm.solve(self.m_fn, self.m_start)
    l_stop = datetime.datetime.now()
    print "Time elapsed %.3f seconds" %(l_stop-l_start).total_seconds()

    if math.fabs(l_fx1 - self.m_expectedResult) < l_maxErrorTolerance:
      print "Acceptable error. Test passed."
      return True

    print "Unacceptable error. Test failed."
    return False

def main():

  """ Functions """

  """ Rosenbrock's """
  testRb = Test(x_name="Rosenbrock's",
    x_fn=lambda x: 100 * ((x[1] - x[0] ** 2) ** 2) + (1 - x[0]) ** 2,
    x_start=[2000, -3000],
    x_maxIter=1000,
    x_expectedResult=0.0)

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
    x_maxIter=6000,
    x_expectedResult=0.0)

  """ Powell's """
  testP = Test(x_name="Powell's",
    x_fn=lambda x: (  (x[0] + 10 * x[1]) ** 2
                     +  5 * (x[2] - x[3]) ** 2
                     +  (x[1] - 2 * x[2]) ** 4
                     + 10 * (x[0] - x[3]) ** 4
                     ),
    x_start=[1000, -3000, 1020, -230],
    x_maxIter=3000,
    x_expectedResult=0.0)

  """ Wayburn Seader 1 """
  testWs = Test(x_name="Wayburn Seader 1",
    x_fn=lambda x: ((x[0] ** 6 + x[1] ** 4 - 17) ** 2
                      +  (2 * x[0] + x[1] - 4) ** 2
                      ),
    x_start=[1000, -3000],
    x_maxIter=2000,
    x_expectedResult=0.0)

  """ Made-up function """
  testCubic = Test(x_name="Made-up Cubic",
    x_fn = lambda x: (x[0] * (x[0] - 1)
                     + (x[1] - 2) * (x[1] + 2)
                     + (x[2]) * (x[2] + 3)
                     ),
    x_start=[2000, -3000, 2],
    x_maxIter=1000,
    x_expectedResult=-6.5)

  l_testsPassed = 0
  l_tests = [testRb, testW, testP, testWs, testCubic]
  print "Runnning tests..."
  for i in xrange(len(l_tests)):
    try:
      print "Test %u:" %(i + 1)
      if l_tests[i].run():
        l_testsPassed += 1
    except ConvergenceError as e:
      print e.message
    except:
      traceback.print_exc()
    finally:
      print

  print "Report: %u out of %u tests passed" %(l_testsPassed, len(l_tests))

if __name__ == "__main__":
    main()
