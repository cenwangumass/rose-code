NAME fx_unb_free
ROWS
 N obj
 G c1
 L c2
COLUMNS
 x obj 1 c1 1
 x c2 2
 y obj 0.5 c1 1
 y c2 -1
 z obj 1 c2 1
RHS
 rhs c1 -5 c2 6
BOUNDS
 FR bnd x
 FR bnd y
 UP bnd z 4
ENDATA
