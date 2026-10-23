NAME fx_bar_rngGE
ROWS
 N obj
 G g1
 G g2
 G g3
 E e1
 E e2
 E e3
COLUMNS
 x1 obj 1 g1 1
 x1 e1 1 e2 1
 x2 obj 2 g1 1
 x2 g2 1 e3 1
 x3 obj 1 g3 1
 x3 e1 -1 e2 1
 x4 obj 3 g2 1
 x4 g3 1 e3 1
RHS
 rhs g1 2 g2 1
 rhs g3 1 e1 1
 rhs e2 3 e3 2
RANGES
 rng g1 2 g2 20
 rng g3 3 e1 4
 rng e2 -4 e3 0
BOUNDS
 UP bnd x2 4
 UP bnd x4 5
ENDATA
