# AG-code-representation

This package tests heuristic assumption claimed in the paper
"On the computation of a representation of an algebraic geometric code",
that takes as input a generator matrix G (k * n) of the WAG code
and computes maximal value t, such that a shortened code G_t is a VSAG code.

Auxiliary routines restoring an equivalent code representation (Y,Q,F)
following [1] are also implemented.

Two examples are provided:
- AG code over F_(2^4) constructed on Elliptic curve;
- AG code over F_(3^9) constructed on Hermitian curve.

[1] Márquez-Corbella, I., Martínez-Moro, E., Pellikaan, R., Ruano, D.: Computational
aspects of retrieving a representation of an algebraic geometry code. Journal of Symbolic
Computation 64, 67–87 (2014) https://doi.org/10.1016/j.jsc.2013.12.007 . Mathematical
and computer algebra techniques in cryptology

AUTHORS:
- Semyon A. Novoselov [Immanuel Kant Baltic Federal University, Kaliningrad, Russia]
- Nikita S. Kolesnikov [Immanuel Kant Baltic Federal University, Kaliningrad, Russia]
- Artyom A. Kuninets [Immanuel Kant Baltic Federal University, Kaliningrad, Russia]
- Ekaterina S. Malygina [MIEM, HSE University, Moscow, Russia]
- Arthur V. Kuleshov [Immanuel Kant Baltic Federal University, Kaliningrad, Russia]
