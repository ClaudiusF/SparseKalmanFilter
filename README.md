# SparseKalmanFilter
Kalman Filter algorithm (in RcppArmadillo code) with Smoother and Disturbance based on sparse matrices for time-varying/invariant models.

Two solutions are proposed:
1) kalmanViaListe.cpp: list-based
2) kalmanViaIndici.cpp: index-based

This work was carried out as a thesis project, aiming to develop a Kalman algorithm that would be computationally faster in the presence of a large state vector that develops sparsity within matrices.

To access the thesis work, please feel free to contact me.
