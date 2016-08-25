# Nedler-Mead Method
This is a Python module for the Nedler-Mead method. The implementation uses numpy. The algorithm used is available [here](https://en.wikipedia.org/wiki/Nelder-Mead_method#One_possible_variation_of_the_NM_algorithm).

## Usage
The basic usage is as follows:
	import nedlerMead
	...
	l_solver = nedlerMead.NedlerMeadSolver(<number of dimensions>, <function>)
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
