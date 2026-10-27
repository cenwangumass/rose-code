# Standalone reference fixtures

The primal JSON fixtures are compact extracts of the cached reference trace and
run for the named model and profile. Version 1 stores each pivot as
`[entering, leaving]`; a null leaving variable means that the original trace
did not identify one. The result stores the original iteration count,
objective, primal vector, reduced costs and dual vector. The tests check every
pivot and the same result fields as the original trace tests. The scale fixture
contains the eight available release captures (SC50A, BLEND, BORE3D, BOEING2,
DEGEN2, PILOT4, PEROLD and DEGEN3); AFIRO and ADLITTLE had no capture in the
cache. These fixtures require neither the reference runtime nor the capture tools.

The MPS files copied from the analysis corpus include the Netlib models 80BAU3B,
FIT2P, GREENBEA and PEROLD and small generated `fx_` models used by the tests.
Other models were already present. `standalone-sha256.txt` records SHA-256 and
names of the 53 added data files, including all compact extracts and scale data.
The hashes are for these packaged files, not the original unabridged traces.

Decision counts, controlled: 80BAU3B 5022, FIT2P 4588. Release: AFIRO 14,
SC50A 32, BLEND 64, ADLITTLE 104, BORE3D 101, BOEING2 217, DEGEN2 690,
PILOT4 690, PEROLD 1791, DEGEN3 5165, 80BAU3B 13805, FIT2P 7523.
