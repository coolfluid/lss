
coolfluid_add_test( ATEST atest-lss-simple    PYTHON atest-lss-simple.py    LIBS cf3_lss cf3_lss_pardiso )
coolfluid_add_test( ATEST atest-lss-matrices  PYTHON atest-lss-matrices.py  LIBS cf3_lss cf3_lss_pardiso )

