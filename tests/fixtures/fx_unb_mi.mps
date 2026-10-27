NAME fx_unb_mi
ROWS
 N obj
 L c1
 G c2
 E c3
COLUMNS
 x obj 1 c1 2
 x c3 1
 y obj -2 c1 -1
 y c2 3 c3 -1
 t obj -1 c1 1
 t c3 1
RHS
 rhs c1 8 c2 -4
 rhs c3 2
BOUNDS
 MI bnd x
 UP bnd y 5
 MI bnd t
ENDATA
