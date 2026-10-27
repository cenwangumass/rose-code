NAME fx_bar_nrow
ROWS
 N obj
 N free1
 L r1
 L r2
 G r3
 E r4
COLUMNS
 x1 obj -1 free1 3
 x1 r1 1 r2 1
 x1 r4 1
 x2 obj -2 free1 -1
 x2 r1 1 r3 1
 x2 r4 1
 x3 obj -1 free1 2
 x3 r2 1 r3 1
 x4 obj 1 free1 1
 x4 r3 1 r4 -1
RHS
 rhs r1 10 r2 8
 rhs r3 1 r4 4
BOUNDS
 UP bnd x3 6
ENDATA
