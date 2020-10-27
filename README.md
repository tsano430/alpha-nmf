# alpha-nmf

Implementation of "[A damped Newton algorithm for nonnegative matrix factorization based on alpha-divergence][1]"

Abstract
--------

A novel Newton-type algorithm for nonnegative matrix factorization based on $\alpha$-divergence is proposed in this paper. The proposed algorithm is a cyclic coordinate descent algorithm that decreases the objective function value along one coordinate direction at a time by using a damped Newton method for monotone equations. It is proved that the proposed algorithm has the global convergence property in the sense of Zangwill. It is also shown experimentally that the proposed algorithm is fast, independent of the value of $\alpha$ while conventional algorithms become very slow for some values of $\alpha$. 

Install LAPACKE/CBLAS
---------------------

[LAPACKE/CBLAS][2] is used in this repository. 

### Download 

URL: https://github.com/Reference-LAPACK/lapack/archive/v3.9.0.tar.gz

### Install

URL: https://gist.github.com/tsano430/b9ee96c39937876a576da62e7f917f96

[1]: https://ieeexplore.ieee.org/document/9010306
[2]: http://www.netlib.org/lapack/lapacke.html
[3]: http://www.cl.cam.ac.uk/research/dtg/attarchive/facedatabase.html
[4]: http://glaros.dtc.umn.edu/gkhome/cluto/cluto/overview/
