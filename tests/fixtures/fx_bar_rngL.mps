NAME fx_bar_rngL
ROWS
 N obj
 L r1
 L r2
 L r3
 L r4
 L r5
 L r6
COLUMNS
 x1 obj -1 r1 1
 x1 r2 1 r3 1
 x2 obj -2 r1 1
 x2 r4 1 r5 1
 x3 obj -1 r2 1
 x3 r4 1 r6 1
 x4 obj -1 r3 -1
 x4 r5 1 r6 1
RHS
 rhs r1 10 r2 10
 rhs r3 10 r4 10
 rhs r5 8 r6 6
RANGES
 rng r1 10 r2 15
 rng r3 10 r4 12
 rng r5 8 r6 5.9
BOUNDS
 UP bnd x4 5
ENDATA
