NAME fx_bar_max
OBJSENSE
    MAX
ROWS
 N obj
 L r1
 L r2
 G r3
 E r4
COLUMNS
 x1 obj 1 r1 1
 x1 r2 1 r4 1
 x2 obj 2 r1 1
 x2 r3 1 r4 1
 x3 obj 1 r2 1
 x3 r3 1
 x4 obj -1 r3 1
 x4 r4 -1
RHS
 rhs obj 3 r1 10
 rhs r2 8
 rhs r3 1 r4 4
BOUNDS
 UP bnd x3 6
ENDATA
