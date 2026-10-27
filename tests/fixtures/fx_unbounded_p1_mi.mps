NAME fx_unbounded_p1
ROWS
 N obj
 L c1
 L c2
 L c3
 L c4
 G c5
COLUMNS
 x1 obj -1 c1 1
 x1 c2 -1 c5 1
 x2 obj -1 c1 -1
 x2 c2 1
 x3 obj -3 c3 1
 x3 c5 1
 x4 obj -2 c4 1
 x4 c5 1
RHS
 rhs c1 2 c2 3
 rhs c3 5 c4 7
 rhs c5 4
BOUNDS
 MI bnd x1
 UP bnd x1 5
ENDATA
