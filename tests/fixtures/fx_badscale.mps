NAME fx_badscale
ROWS
 N obj
 L c1
 L c2
 G c3
COLUMNS
 x obj -1e-6 c1 1e6
 x c2 1e-6
 y obj -1e6 c1 1e-3
 y c2 1e3
 y c3 1
RHS
 rhs c1 1e7 c2 5e-3
 rhs c3 1e-7
BOUNDS
 UP bnd x 1e4
ENDATA
