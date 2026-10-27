* PDHG exit-rule fixture (G rows c2 c3 binding): maximize with objective constant (+5), L/G/E rows, one ranged G row (c4: -3 <= y - z <= 4), a free column, finite upper bounds.
NAME fx_pdhg_max
OBJSENSE
    MAX
ROWS
 N obj
 L c1
 G c2
 G c3
 G c4
 L c5
 E c6
COLUMNS
 x obj 3 c1 1
 x c2 1 c3 1
 x c5 1
 y obj 2 c1 1
 y c2 -1 c4 1
 y c6 1
 z obj -1 c1 1
 z c3 3 c4 -1
 w obj 1 c1 1
 w c5 2 c6 1
RHS
 rhs c1 14 c2 1
 rhs c3 20 c4 -3
 rhs c5 12 c6 3.5
 rhs obj -5
RANGES
 rng c4 7
BOUNDS
 UP bnd x 4
 FR bnd y
 UP bnd z 6
 UP bnd w 5
ENDATA
