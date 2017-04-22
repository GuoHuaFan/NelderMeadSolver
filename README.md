# Nelder-Mead Method
This is a Python module for the Nelder-Mead method. The implementation uses numpy. The algorithm used is available [here](https://en.wikipedia.org/wiki/Nelder-Mead_method#One_possible_variation_of_the_NM_algorithm).

## Development Notes
This has been developed and tested on Mac OSX (Darwin Kernel Version 15.6.0) using Python 2.7.

## Usage
The basic usage is as follows:

	import nelderMead
	...
	l_solver = nelderMead.NelderMeadSolver(<number of dimensions>, <function>)
	l_solver.solve(<start point>)

## Tests
To the tests, run `tests.py`

	$ python tests.py
	Runnning tests...
	Test 1:
	Solving Rosenbrock's function...
	Algorithm converged (810 iterations)
	Time elapsed 0.746 seconds
	
	Test 2:
	Solving Wood's function...
	Algorithm converged (5615 iterations)
	Time elapsed 6.206 seconds
	
	Test 3:
	Solving Powell's function...
	Algorithm converged (2563 iterations)
	Time elapsed 2.805 seconds
	
	Test 4:
	Solving Wayburn Seader 1 function...
	Algorithm converged (1657 iterations)
	Time elapsed 1.512 seconds
	
	Test 5:
	Solving Made-up Cubic function...
	Algorithm converged (944 iterations)
	Time elapsed 0.944 seconds
	
	Report: 5 out of 5 tests passed
 
## Code Coverage
Run [Coverage.py](https://coverage.readthedocs.io/en/coverage-4.2/) to generate coverage report

	$ coverage run tests.py
	$ coverage report --omit="/System*"
	Name            Stmts   Miss  Cover
	-----------------------------------
	nelderMead.py     148     25    83%
	tests.py           53      6    89%
	-----------------------------------
	TOTAL             201     31    85%
