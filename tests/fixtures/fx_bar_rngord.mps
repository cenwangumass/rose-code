NAME fx_bar_rngord
ROWS
 N obj
 G g1
 L l2
 G g3
 E e1
 E e2
 E e3
COLUMNS
 x1 obj 1 g1 100
 x1 e1 1 e2 0.01
 x2 obj 2 g1 1
 x2 l2 300 e3 1
 x3 obj 1 g3 1
 x3 e1 -1 e2 1
 x4 obj 3 l2 1
 x4 g3 50 e3 1
RHS
 rhs g1 200 l2 1000
 rhs g3 100 e1 1
 rhs e2 3 e3 2
RANGES
 rng e3 4000 l2 6000
 rng e1 -800 g1 5000
 rng e2 0 g3 250
BOUNDS
 UP bnd x2 40
 UP bnd x4 500
ENDATA
