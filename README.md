This repository includes the source codes, test instances, and computational results for the flight test scheduling problem (FTSP) introduced in [1].

We thank Prof. Jordi Pereira for kindly providing us with the source code developed in one of his works[2].

We also use benchmark instances for simple assembly line balancing problem of type 1 (SALBP-1) introduced in [3] and [4].

The proposed algorithms are coded in C++ and run in a single thread. All MIP formulations are solved with CPLEX 20.1.0 implemented in the C++ programming language, using Concert Technology. The CPLEX solver is set to use only one thread, with all other parameters set to their default values.

If you find any problem, please feel free to contact with taohanqiao19@nudt.edu.cn

Below is the BibTex for citing this work.

@article{FTSP2025,
title = {Flight test scheduling: A generic model, lower bounds, and iterated local search},
journal = {European Journal of Operational Research},
year = {2025},
issn = {0377-2217},
doi = {https://doi.org/10.1016/j.ejor.2025.06.001},
author = {Hanqiao Tao and Guopeng Song and Roel Leus and Zhe Liang and Jiang Jiang}
}

**References**

[1] Tao, H., et al. Flight test scheduling: A generic model, lower bounds, and iterated local search (doi: https://doi.org/10.1016/j.ejor.2025.06.001).

[2] Pereira, J. (2016). Procedures for the bin packing problem with precedence constraints. European Journal of Operational Research, 250, 794–806.

[3] Otto, A., Otto, C., & Scholl, A. (2013). Systematic data generation and test design for solution algorithms on the example of SALBPGen for assembly line balancing. European Journal of Operational Research, 228, 33–45.

[4] Álvarez-Miranda, E., Pereira, J., & Vilà, M. (2023b). Analysis of the simple assembly line balancing problem complexity. Computers & Operations Research, 159, 106323.
