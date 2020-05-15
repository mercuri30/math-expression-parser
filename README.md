# math-expression-parser
A parser for mathematical expressions.
It is based on Reverse Polish notation for evaluation and neighbour checks for correctness.
Variables are supported.

Input expression examples:
* 2 \* exp(1) ^ ln(7)
* a=1 t=-0.1 x=0.25  sin(2 \* pi \* x) \* exp(-4 \* pi ^ 2 \* a \* t)
* 1 + 2 - 3 \* 4 / 5 ^ 6 \* (2 \* (1 - 5 + (3 \* 7 ^ 9) \* (4 + 6 \* 7 - 3))) + 12  

Much more examples are provided in the unit tests.
