* Large finite bounds around the 1e20 infinity threshold: X1 UP 1e20, X2 UP 5e19, X3 LO -1e20 (+UP 4), X4 LO -5e19 (+UP 6), X5 UP 1e30, rows with rhs 1e20.
NAME          E22
ROWS
 N  OBJ
 L  R1
 G  R2
 L  R3
 L  R4
COLUMNS
    X1        OBJ                1.0   R1                 1.0
    X1        R2                 1.0
    X2        OBJ                1.0   R1                 1.0
    X2        R3                 1.0
    X3        OBJ                1.0   R2                 1.0
    X3        R4                 1.0
    X4        OBJ                1.5   R1                 1.0
    X4        R3                -1.0
    X5        OBJ                1.0   R2                 1.0
    X5        R4                 1.0
RHS
    RHS       R1                10.0   R2                 2.0
    RHS       R3                 1e20  R4                 8.0
RANGES
    RNG       R1                 1e30
BOUNDS
 UP BND       X1                 1e20
 UP BND       X2                 5e19
 MI BND       X3
 UP BND       X3                 4.0
 LO BND       X4                -5e19
 UP BND       X4                 6.0
 UP BND       X5                 1e30
ENDATA
