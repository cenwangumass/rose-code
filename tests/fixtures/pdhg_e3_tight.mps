NAME          E3
ROWS
 N  OBJ
 E  E1
 G  G1
 L  L1
 E  E2
 N  FREEROW
COLUMNS
    X1        OBJ                1.0   E1                 1.0
    X1        G1                 1.0   FREEROW            5.0
    X2        OBJ                2.0   E1                 1.0
    X2        L1                 1.0
    X3        OBJ                3.0   G1                 2.0
    X3        L1                 1.0   E2                 1.0
    X4        OBJ                1.0   E2                 3.0
    X4        E2                 1.0
    X5        OBJ                1.0
RHS
    RHS       E1                 0.0   G1                 0.0
    RHS       L1                 0.0   E2                 0.0
RANGES
    RNG       E1                 3.0   E2                -2.5
    RNG       G1                 2.0
BOUNDS
 UP BND       X1                -1.0
 UP BND       X2                 1e20
 UP BND       X3                 1e30
 LO BND       X4                 1.0
 UP BND       X4                 1.0
ENDATA
