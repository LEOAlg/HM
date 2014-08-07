HM
==
This is a Matlab and C implementation for the algorithm LEO.
1. Matlab: > EMbayesnet.m 

2. i) C Code compiles as: > g++ -std=c++0x -w -fopenmp combEM.cpp -O1 -larmadillo
   ii) C Code runs as: ./a.out 1
   % Here 1 represents the application which we want to run, it can be any application between {1,25}.
   

Notes
==
This algorithm works well for most cases. But some applications might require some engineering. I have listed below the tweaks that can be tried.

1. Increase the iteration_limit and epsilon variable in Matlab/splitEM to get higher accuracy. Most of the applications should work well with iteration_limit = 3 and epsilon = 10000. 
But some would require higher values for better convergence. 

2. Similarly, the number of samples required for good accuracy would vary for different applications.

3. The algorithm can be made much faster by using diagonal inverse covariance matrix and having a gamma prior on the diagonal elements.

4. The normalization of the application is extremely crucial. Some normalizations might better compared to others.
