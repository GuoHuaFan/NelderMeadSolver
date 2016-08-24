"""
  Filename    : tests.py
  Author      : Soumyaroop Roy
  Date        : May 02, 2016
  Description : Tests Nedler-Mead
"""

import sys
import c_NedlerMead as NedLerMead

""" Global debug level """
g_dbgLvl = 0
g_progressReportInterval = 500

"""
  Function to print an error message

  @param  string

  @return None
"""
def error (x_str):
  print ("ERROR: %s" %x_str)
  sys.exit()

"""
  Function to print a debug message

  @param  String
  @param  Debug level

  @return None
"""
def debugPrint (x_str, x_dbgLvl=1):
  if x_dbgLvl <= g_dbgLvl:
    print ("Debug (Level %u): %s" %(x_dbgLvl, x_str))

def main():

  """ Functions """
  l_fn1 = lambda x: (x[0] * (x[0] - 1)
                     + (x[1] - 2) * (x[1] + 2)
                     + (x[2]) * (x[2] + 3)
                     )
  l_start1 = [2000, -3000, 2]

  """ Rosenbrock's """
  l_fnRb = lambda x: 100 * ((x[1] - x[0] ** 2) ** 2) + (1 - x[0]) ** 2
  l_startRb = [2000, -3000]

  """ Wood's """
  l_fnW = lambda x: ( 100 * (x[1] - x[0] ** 2) ** 2
                     +              (1 - x[0]) ** 2
                     + 90 * (x[3] - x[2] ** 2) ** 2
                     +              (1 - x[2]) ** 2
                     +       10.1 * (x[1] - 1) ** 2
                     +       10.1 * (x[3] - 1) ** 2
                     + 19.8 * (x[1] - 1) * (x[3] - 1)
                     )
  l_startW = [1000, -3000, 1020, -230]

  """ Powell's """
  l_fnP = lambda x: (  (x[0] + 10 * x[1]) ** 2
                     +  5 * (x[2] - x[3]) ** 2
                     +  (x[1] - 2 * x[2]) ** 4
                     + 10 * (x[0] - x[3]) ** 4
                     )
  l_startP = [1000, -3000, 1020, -230]

  """ Wayburn Seader 1 """
  l_fnWs = lambda x: ((x[0] ** 6 + x[1] ** 4 - 17) ** 2
                      +  (2 * x[0] + x[1] - 4) ** 2
                      )
  l_startWs = [1000, -3000]

  l_fn = l_fnW
  l_start = l_startW
  l_spaceDim = len(l_start)

  """ Sample coefficients """
  l_coeffs = (1,   #\alpha
              1,   #\gamma
              0.5, #\rho
              0.5  #\sigma
              )

  """ Construct a Nedler Mead optimizer """
  l_nm = c_NedlerMead(l_spaceDim, l_fn, l_coeffs)

  """ Run it """
  l_nm.run(l_start)

if __name__ == "__main__":
    main()
