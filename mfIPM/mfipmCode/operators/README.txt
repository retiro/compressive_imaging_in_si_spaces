OPERATORS: MATLAB files for the partial DCT operator and wrapper 
           files for the construction of A_operators.

1. Contents
===========

* Subdirectory "partialDCT": Implements a 1-d partial DCT operator. For details please
                             look in the pdct.m file inside the "PartialDCT" subdirectory.
                             For a quick help, run in the MATLAB prompt:
                               >> help pdct
                           

* Subdirectory "@A_operator": Implements a wrapper for matrix-operators.
                              It can be used to call matrix-operators as they were
                              implicit matrices. 

                          Example:
                          A = A_operator(@(x) pdct(x,1,n,sort(picks(1:m))), @(x) pdct(x,2,n,sort(picks(1:m))));
                          Where pdct is a partial DCT operator, see in the "PartialDCT"
                          subdirectory.
                          Now the pdct operator can be used through A like an
                          implicity matrix, i.e A*x and A'y.

                          The subdirectory "@A_operator" has been copied from the 
                          FPC_AS_v1.21 software package.
                          url: http://www.caam.rice.edu/~optimization/L1/FPC_AS/
 
                          For a  quick help, run in the MATLAB prompt:
                            >> help A_operator
