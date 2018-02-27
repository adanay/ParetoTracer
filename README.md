# ParetoTracer

## Introduction
The Pareto Tracer (PT) is a predictor-corrector method for the numerical treatment of sufficiently smooth multi-objective optimization 
problems (MOP). The algorithm performs a continuation along the set of (local) solutions of a given MOP with 𝑘 objectives and can cope 
with equality and box constraints. 
For the detailed explanations of the algorithms behind the PT we refer to: 

**[1]** A. Martín and O. Schütze<br/>
**Pareto Tracer: a predictor–corrector method for multi-objective optimization problems**<br/> 
Engineering Optimization 50 (3): 516-536, 2018<br/>
http://www.tandfonline.com/doi/abs/10.1080/0305215X.2017.1327579?journalCode=geno20

**[2]** A. Martín<br/>
**Pareto Tracer: A Predictor Corrector Method for Multi-objective Optimization Problems**<br/>
MSc Thesis, Cinvestav-IPN, Mexico, 2014<br/>
www.cs.cinvestav.mx/TesisGraduados/2014/TesisAdanayMartin.pdf

## Getting Started 
This implementation of PT was developed using MATLAB 2015. The code is organized in different packages depending on the scope. 
The main package is called **`pt`** which stands for Pareto Tracer. This contains the main entry point to the algorithm: the function 
**`pt.trace`**. 
Below there is a basic example on how to call this function. 

**`[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);`**

Each one of the parameters is explained later in detail. For now, another key function is introduced: **`pt.minimize`**. This function is 
utilized for the corrector phase of the algorithm. For the unconstrained case, it coincides with the Newton method proposed by 
**Fliege et al. in [3]**. A modification was proposed by **Martín et al. in [1]** to handle equality and box constraints. The current 
implementation of the pt.minimize function follows these instructions. A basic example on how to call this function is given below.

**`[result, stats, EXITFLAG] = pt.minimize(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);`**

**[3]** Jörge Fliege, L. M. Graña Drummond, and Benar F. Svaiter.<br/> 
**Newton’s method for multiobjective optimization.** <br/>
SIAM Journal on Optimization, 20(2):602–626, 2009.<br/>

## Ready-to-use Examples
Several ready-to-use examples are provided, i.e., script files containing examples on how to call the PT main functions. The examples
are all grouped in two packages called **`x_trace`** and **`x_min`**. The x stands for examples or experiments, and the rest of the
folder name denotes which function is being tested. Additionally, the current experiments are grouped by the Hessian approximation 
strategy utilized in the experiment and by the function benchmark name. Finally, there is one file per experiment, **which is ready to use 
by just clicking the Run button of the MATLAB interface**. 

![x_trace](img/readme1.png)
  
Both the plotted and printed results of running the script **`x_trace/exact/misc/quad_n100_nobj2.m`** are displayed below. The 
experiments always start on a randomly selected point on the Pareto set or close to it. Note that the starting point is a blue star 
while the last point is a red star. For bi-objective problems, PT always goes left up first, and later it goes right down the optimal 
curve.  

![x_trace](img/readme_quad_n100_nobj2.png)

```
>> x_trace.exact.misc.quad_n100_nobj2

Pareto Tracer
Func: Quad(n=100,nobj=2)
Initial Point: [0.98,0.19,0.22,…]
Initial Fun Value: [33.73,231.47]
Hess: Modif: chol, Approx: bfgs (Note: Approx used only if hess not provided.)
Step in Obj: 10.00

Exit Flag (1): No more solution points found.
Iterations: 65
Solution Points: 66
Correct Stats: Avg Its: 0.03, Avg Lin Search Its: 0.03

Fun Evals: 68
Jac Evals: 68
Hess Evals: 68
Elapsed time is 2.908885 seconds.
```

Analogously, the plotted result (and part of the printed result) of running the script **`x_min/exact/misc/dent.m`** is displayed below. 
The experiments for the pt.minimize function are setup such that the algorithm is executed 10 times starting at different randomly 
selected points.

![x_trace](img/readme_dent.png)
 
```
>> x_min.exact.misc.dent

PT Minimize
Func: Dent(n=2,nobj=2)
Initial Point: [0.73,0.16]
Initial Fun Value: [2.14,1.57]
Hess: Modif: chol, Approx: bfgs (Note: Approx used only if hess not provided.)

Exit (1): First-order optimality measure was less than 1.000000e-06.
Solution Point: [0.29,-0.29]
Solution Fun Value: [1.97,1.40]
Iterations: 4
Optimality: d: -1.548985e-09, ||v||^2: 3.047592e-09
Avg Dir Subprob Its: 6.00, Avg Lin Search Its: 1.00

Fun Evals: 5
Jac Evals: 5
Hess Evals: 5
Elapsed time is 0.232982 seconds.
```

One last example is provided in this section to illustrate the output of PT for a function with more than two objectives and one 
equality constraint.

![x_trace](img/readme_sproblem1_n100_nobj3.png)

```
>> x_trace.exact.eq.sproblem1_n100_nobj3

Pareto Tracer
Func: SProblem with lin eq(n=100,nobj=3,naeq=1)
Initial Point: [0.95,0.38,0.20,…]
Initial Fun Value: [32.03,244.26,143.54]
Hess: Modif: chol, Approx: bfgs (Note: Approx used only if hess not provided.)
Step in Obj: 10.00

Exit Flag (1): No more solution points found.
Iterations: 651
Solution Points: 651
Correct Stats: Avg Its: 1.30, Avg Lin Search Its: 0.93

Fun Evals: 3780
Lin Eq Evals: 2060
Jac Evals: 2039
Hess Evals: 2039
Elapsed time is 123.262057 seconds.
```

**Note**: This section has selected the most basic problems to demonstrate the capabilities of the current implementation of the algorithm. 
PT does not perform well in all benchmark problems and those examples are also included in the experiment set coming with this 
implementation. The WFG benchmark is one example of a very challenging set of functions for PT.

## Parameters
Both **`pt.trace`** and **`pt.minimize`** receive the following set of parameters:
-	**objfun**: It must be either a cell array of function handles or a struct with these fields. I.e., **`objfun = {f, J, H}`** where f, J, H are function handles that represent the objective, Jacobian, and Hessian functions respectively. They can be empty except f. J and H will be approximated if not provided. objfun can also be a function handle. In this case, it will be assumed to be the objective function only. I.e., objfun = f. 
    - **`f`**  is a function handle of the form **`y = f(x)`** where x is a vector of n components and y is a vector of nobj components.
    - **`J`**  is a function handle of the form **`y = J(x)`** where x is a vector of n components and y is a matrix of size (nobj x n).
    - **`H`**  is a function handle of the form **`y = H(x)`** where x is a vector of n components and y is a block matrix of size (n x n x nobj).

-	**x0**: The initial guess. It must be a vector of n dimensions. It must be specified.

-	**funvals0**: The known function values at x0. It can be empty (or any of its fields can be empty). If specified, it must be either a cell array or a struct with these fields:
**`funvals0 = {fx, Jx, Hx, ax, aeqx, cx, ceqx, dcx, dceqx, Jcx, Jceqx}`**. 
    - **`fx = objfun.f(x0)`**  % objective function
    - **`Jx = objfun.J(x0)`**  % Jacobian 
    - **`Hx = objfun.H(x0)`**  % Hessian
    - **`Ax = lincon.A * x – lincon.b`**  % linear inequalities
    - **`Aeqx = lincon.Aeq * x – lincon.beq`**  % linear equalities
    - **`cx = nonlcon.c(x0)`**  % nonlinear inequalities
    - **`ceqx = nonlcon.ceq(x0)`**  % nonlinear equalities
    - **`dcx = norm(ax)^2 + norm(cx)^2`**  % square norm of the inequalities
    - **`dceqx = norm(aeqx)^2 + norm(ceqx)^2`**  % square norm of the equalities
    - **`Jcx = nonlcon.Jc(x0)`**  % Jacobian of the nonlinear inequalities
    - **`Jceqx = nonlcon.Jceq(x0)`**  % Jacobian of the nonlinear equalities

-	**`lb, ub`**: Vectors that represent the box constraints in decision space. They must have n components or be empty.

-	**`lincon`**: It must be either a cell array of matrices or a struct with these fields. I.e., **`lincon = {A, b, Aeq, beq}`** representing the linear inequality and equality constraints. They all can be empty.
A is a matrix of size (na x n) where na is the number of linear inequalities.
    - **`b`** is a vector of na components.
    - **`Aeq`** is a matrix of size (naeq x n) where naeq is the number of linear equalities.
    - **`beq`** is a vector of naeq components.

-	**`nonlcon`**: It must be either a cell array of function handles or a struct with these fields. I.e., **`nonlcon = {c, ceq, Jc, Jceq}`** representing the inequality and equality constraints together with their respective Jacobians. If the Jacobians are not provided, they will be approximated. nonlcon can also be a function handle. In this case, it will be assumed to be the nonlinear inequality constraints function only. I.e., nonlcon = c.
    - **`c`** is a function handle of the form **`y = c(x)`** where x is a vector of n components and y is a vector of nc components.
    - **`ceq`** is a function handle of the form **`y = ceq(x)`** where x is a vector of n components and y is a vector of nceq components.
    - **`Jc`** is a function handle of the form **`y = Jc(x)`** where x is a vector of n components and y is a matrix of size (nc x n).
    - **`Jceq`** is a function handle of the form **`y = Jceq(x)`** where x is a vector of n components and y is a matrix of size (nceq x n).

-	**`multfun`**: It must be either a cell array of function handles or a struct with these fields. I.e., **`multfun = {vH, Hw, Hwv}`** representing the Hessian multiply functions. They all can be empty. If specified, the multiply functions will be utilized instead of the Hessian function.
    - **`vH`** is a function handle of the form **`y = vH(x, v)`** where **`y = [v' * H1; v' * H2; ...; v' * Hnobj]`**. The result y has a size of (nobj x n).
    - **`Hw`** is a function handle of the form **`y = Hw(x, w)`** where **`y = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj`**. The result y has a size of (n x n), i.e., the weighted sum of Hessians.
    - **`Hwv`** is a function handle of the form **`y = Hwv(x, w, v)`** where **`y = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v`**. The result y is a vector of n components.

-	**`opts`**: This is a structure containing all the options that can be passed to the algorithms. There is a separate section dedicated to this. 


