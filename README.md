# bode

A low-order adaptive solver for implicit ODE's

The solver has been extracted from:

> Hopkins, T. R., & Wait, R. (1978). A comparison of Galerkin, collocation and the method of lines for partial differential equations. International Journal for Numerical Methods in Engineering, 12(7), 1081-1107. https://doi.org/10.1002/nme.1620120703

The method used appeared originally in

> Prothero, A., & Robinson, A. (1974). On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations. Mathematics of Computation, 28(125), 145-162. https://doi.org/10.1090/S0025-5718-1974-0331793-2, [PDF](https://www.ams.org/journals/mcom/1974-28-125/S0025-5718-1974-0331793-2/S0025-5718-1974-0331793-2.pdf)

---

Hamming Predictor Corrector Method:
- https://dl.acm.org/doi/10.1145/320954.320958
- https://dl.acm.org/doi/pdf/10.1145/321138.321143
- http://www.ccpo.odu.edu/~klinck/Reprints/PDF/hammingJACM1959.pdf
- https://www.tat.physik.uni-tuebingen.de/~kokkotas/Teaching/Num_Methods_files/Comp_Phys6-md.pdf
- https://www.cfm.brown.edu/people/dobrush/am33/Mathematica/hamming.html
- https://www.ams.org/journals/mcom/1973-27-121/S0025-5718-1973-0331803-1/
- https://www.ams.org/journals/mcom/1970-24-109/S0025-5718-1970-0280010-7/S0025-5718-1970-0280010-7.pdf
- https://www.proquest.com/openview/c11b30c59ce7000a0759ac2b81a48171/1?pq-origsite=gscholar&cbl=18750&diss=y
- https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/3j333584h

See also:
- https://github.com/ivan-pi/stiff3
- http://www.unige.ch/~hairer/software.html
