
# MLE and REML

Machine information


```julia
versioninfo()
```

    Julia Version 0.7.0
    Commit a4cb80f3ed (2018-08-08 06:46 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin14.5.0)
      CPU: Intel(R) Core(TM) i5-6267U CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-6.0.0 (ORCJIT, skylake)


## Demo data

For demonstration, we generate a random data set.


```julia
# generate data from a d-variate response variane component model
using Random, LinearAlgebra
Random.seed!(123)
n = 1000   # no. observations
d = 2      # dimension of responses
m = 2      # no. variance components
p = 2      # no. covariates
# n-by-p design matrix
X = randn(n, p)
# p-by-d mean component regression coefficient
B = ones(p, d)  
# a tuple of m covariance matrices
V = ntuple(x -> zeros(n, n), m) 
for i = 1:m-1
  Vi = randn(n, 50)
  copyto!(V[i], Vi * Vi')
end
copyto!(V[m], Matrix(I, n, n)) # last covarianec matrix is idendity
# a tuple of m d-by-d variance component parameters
Σ = ntuple(x -> zeros(d, d), m) 
for i in 1:m
  Σi = randn(d, d)
  copyto!(Σ[i], Σi' * Σi)
end
# form overall nd-by-nd covariance matrix Ω
Ω = zeros(n * d, n * d)
for i = 1:m
  Ω += kron(Σ[i], V[i])
end
Ωchol = cholesky(Ω)
# n-by-d responses
Y = X * B + reshape(Ωchol.L * randn(n*d), n, d);
```

## Maximum likelihood estimation (MLE)

To find the MLE of parameters $(B,\Sigma_1,\ldots,\Sigma_m)$, we take 3 steps:  

**Step 1 (Construct data)**. Construct an instance of `VarianceComponentVariate`, which consists fields  

* `Y`: $n$-by-$d$ responses  
* `X`: $n$-by-$p$ covariate matrix  
* `V=(V[1],...,V[m])`: a tuple of $n$-by-$n$ covariance matrices. The last covariance matrix must be positive definite and usually is the identity matrix. 


```julia
using VarianceComponentModels
vcdata = VarianceComponentVariate(Y, X, V)
fieldnames(typeof(vcdata))
```




    (:Y, :X, :V)



In the absence of covariates $X$, we can simply initialize by `vcdata = VarianceComponentVariate(Y, V)`.

**Step 2 (Construct a model)**. Construct an instance of `VarianceComponentModel`, which consists of fields  

* `B`: $n$-by-$p$ mean regression coefficients  
* `Σ=(Σ[1],...,Σ[m])`: variane component parameters respectively. 

When constructed from a `VarianceComponentVariate` instance, the mean parameters $B$ are initialized to be zero and the tuple of variance component parameters $\Sigma$ to be `(eye(d),...,eye(d))`.


```julia
vcmodel = VarianceComponentModel(vcdata)
fieldnames(typeof(vcmodel))
```




    (:B, :Σ, :A, :sense, :b, :lb, :ub)




```julia
vcmodel
```




    VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}([0.0 0.0; 0.0 0.0], ([1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0]), Array{Float64}(0,4), Char[], Float64[], -Inf, Inf)



The remaining fields `A`, `sense`, `b`, `lb`, `ub` specify (optional) constraints on the mean parameters `B`:

$A * \text{vec}(B) \,\, =(\text{or } \ge \text{or } \le) \,\, b$

$lb \le \text{vec}(B) \le ub$

`A` is an constraint matrix with $pd$ columns, `sense` is a vector of charaters taking values `'<'`, `'='` or `'>'`, and `lb` and `ub` are the lower and upper bounds for `vec(B)`. By default, `A`, `sense`, `b` are empty, `lb` is `-Inf`, and `ub` is `Inf`. If any constraits are non-trivial, final estimates of `B` are enforced to satisfy them.

When a better initial guess is available, we can initialize by calling `vcmodel=VarianceComponentModel(B0, Σ0)` directly.

**Step 3 (Fit model)**. Call optmization routine `fit_mle!`. The keywork `algo` dictates the optimization algorithm: `:MM` (minorization-maximization algorithm) or `:FS` (Fisher scoring algorithm).


```julia
vcmodel_mle = deepcopy(vcmodel)
@time logl, vcmodel_mle, Σse, Σcov, Bse, Bcov = fit_mle!(vcmodel_mle, vcdata; algo = :MM);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -6.253551e+03
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    ******************************************************************************
    
           1  -3.881454e+03
           2  -3.853179e+03
           3  -3.846525e+03
           4  -3.844906e+03
           5  -3.844506e+03
           6  -3.844406e+03
           7  -3.844381e+03
           8  -3.844375e+03
           9  -3.844374e+03
          10  -3.844373e+03
    
      4.546981 seconds (11.23 M allocations: 566.109 MiB, 5.23% gc time)


The output of `fit_mle!` contains  

* final log-likelihood  


```julia
logl
```




    -3844.3731814180883



* fitted model


```julia
fieldnames(typeof(vcmodel_mle))
```




    (:B, :Σ, :A, :sense, :b, :lb, :ub)




```julia
vcmodel_mle
```




    VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}([1.092 1.04727; 0.955346 1.01632], ([0.380637 -0.305465; -0.305465 4.51938], [1.84009 0.265569; 0.265569 2.17275]), Array{Float64}(0,4), Char[], Float64[], -Inf, Inf)



* standard errors of the estimated varianec component parameters


```julia
Σse
```




    ([0.0765136 0.263047; 0.263047 0.904332], [0.0844292 0.0917441; 0.0917441 0.0996927])



* covariance matrix of the variance component parameters estimates


```julia
Σcov
```




    8×8 Array{Float64,2}:
      0.00585433  -0.00467019  -0.00467019  …  -1.07903e-6   -1.557e-7   
     -0.00467019   0.0691937    0.00372555     -1.557e-7     -1.27444e-6 
     -0.00467019   0.00372555   0.0691937      -8.83212e-6   -1.27444e-6 
      0.00372555  -0.055198    -0.055198       -1.27444e-6   -1.04316e-5 
     -7.4779e-6   -1.07903e-6  -1.07903e-6      0.00102878    0.000148477
     -1.07903e-6  -8.83212e-6  -1.557e-7    …   0.000148477   0.00121477 
     -1.07903e-6  -1.557e-7    -8.83212e-6      0.00841698    0.00121477 
     -1.557e-7    -1.27444e-6  -1.27444e-6      0.00121477    0.00993864 



* standard errors of the estimated mean parameters


```julia
Bse
```




    2×2 Array{Float64,2}:
     0.0425562  0.0483834
     0.0430596  0.0492809



* covariance matrix of the mean parameter estimates


```julia
Bcov
```




    4×4 Array{Float64,2}:
      0.00181103   -1.96485e-5    0.000243441  -4.38252e-6 
     -1.96485e-5    0.00185413   -4.38252e-6    0.000246407
      0.000243441  -4.38252e-6    0.00234096   -5.73331e-6 
     -4.38252e-6    0.000246407  -5.73331e-6    0.00242861 



## Restricted maximum likelihood estimation (REML)

[REML (restricted maximum likelihood estimation)](https://en.wikipedia.org/wiki/Restricted_maximum_likelihood) is a popular alternative to the MLE. To find the REML of a variane component model, we replace the above step 3 by  

**Step 3**. Call optmization routine `fit_reml!`.


```julia
vcmodel_reml = deepcopy(vcmodel)
@time logl, vcmodel_reml, Σse, Σcov, Bse, Bcov = fit_reml!(vcmodel_reml, vcdata; algo = :MM);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -4.215053e+03
           1  -3.925799e+03
           2  -3.865114e+03
           3  -3.851105e+03
           4  -3.847732e+03
           5  -3.846903e+03
           6  -3.846698e+03
           7  -3.846647e+03
           8  -3.846634e+03
           9  -3.846631e+03
          10  -3.846630e+03
    
      0.443964 seconds (8.09 k allocations: 62.532 MiB, 2.39% gc time)


The output of `fit_reml!` contains

* the final log-likelihood at REML estimate


```julia
logl
```




    -3844.3777179025046



* REML estimates


```julia
fieldnames(typeof(vcmodel_reml))
```




    (:B, :Σ, :A, :sense, :b, :lb, :ub)




```julia
vcmodel_reml
```




    VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}([1.092 1.04727; 0.955345 1.01632], ([0.380594 -0.305485; -0.305485 4.51994], [1.84285 0.261963; 0.261963 2.17842]), Array{Float64}(0,4), Char[], Float64[], -Inf, Inf)



* standard errors of the estimated varianec component parameters


```julia
Σse
```




    ([0.0765055 0.26305; 0.26305 0.904446], [0.0845559 0.0919325; 0.0919325 0.0999526])



* covariance matrix of the variance component parameters estimates


```julia
Σcov
```




    8×8 Array{Float64,2}:
      0.0058531   -0.00467005  -0.00467005  …  -1.06597e-6   -1.51499e-7 
     -0.00467005   0.0691951    0.00372613     -1.51499e-7   -1.26041e-6 
     -0.00467005   0.00372613   0.0691951      -8.86843e-6   -1.26041e-6 
      0.00372613  -0.0552092   -0.0552092      -1.26041e-6   -1.0486e-5  
     -7.50035e-6  -1.06597e-6  -1.06597e-6      0.00101633    0.000144472
     -1.06597e-6  -8.86843e-6  -1.51499e-7  …   0.000144472   0.0012014  
     -1.06597e-6  -1.51499e-7  -8.86843e-6      0.00845158    0.0012014  
     -1.51499e-7  -1.26041e-6  -1.26041e-6      0.0012014     0.00999052 



* standard errors of the estimated mean parameters


```julia
Bse
```




    2×2 Array{Float64,2}:
     0.0425881  0.0484485
     0.0430919  0.0493475



* covariance matrix of the mean parameter estimates


```julia
Bcov
```




    4×4 Array{Float64,2}:
      0.00181375   -1.96783e-5    0.000239868  -4.34611e-6 
     -1.96783e-5    0.00185691   -4.34611e-6    0.000242745
      0.000239868  -4.34611e-6    0.00234726   -5.73082e-6 
     -4.34611e-6    0.000242745  -5.73082e-6    0.00243518 



## Optimization algorithms

Finding the MLE or REML of variance component models is a non-trivial nonlinear optimization problem. The main complications are the non-convexity of objective function and the positive semi-definiteness constraint of variane component parameters $\Sigma_1,\ldots,\Sigma_m$. In specific applications, users should try different algorithms with different starting points in order to find a better solution. Here are some tips for efficient computation. 

In general the optimization algorithm needs to invert the $nd$ by $nd$ overall covariance matrix $\Omega = \Sigma_1 \otimes V_1 + \cdots + \Sigma_m \otimes V_m$ in each iteration. Inverting a matrix is an expensive operation with $O(n^3 d^3)$ floating operations. When there are only **two** varianec components ($m=2$), this tedious task can be avoided by taking one (generalized) eigendecomposion of $(V_1, V_2)$ and rotating data $(Y, X)$ by the eigen-vectors. 


```julia
vcdatarot = TwoVarCompVariateRotate(vcdata)
fieldnames(typeof(vcdatarot))
```




    (:Yrot, :Xrot, :eigval, :eigvec, :logdetV2)



Two optimization algorithms are implemented: [Fisher scoring](https://books.google.com/books?id=QYqeYTftPNwC&lpg=PP1&pg=PA142#v=onepage&q&f=false) (`mle_fs!`) and the [minorization-maximization (MM) algorithm](http://arxiv.org/abs/1509.07426) (`mle_mm!`). Both take the rotated data as input. These two functions give finer control of the optimization algorithms. Generally speaking, MM algorithm is more stable while Fisher scoring (if it converges) yields more accurate answer.


```julia
vcmodel_mm = deepcopy(vcmodel)
@time mle_mm!(vcmodel_mm, vcdatarot; maxiter=10000, funtol=1e-8, verbose = true);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -6.253551e+03
           1  -3.881454e+03
           2  -3.853179e+03
           3  -3.846525e+03
           4  -3.844906e+03
           5  -3.844506e+03
           6  -3.844406e+03
           7  -3.844381e+03
           8  -3.844375e+03
           9  -3.844374e+03
          10  -3.844373e+03
    
      0.042187 seconds (21.56 k allocations: 1.366 MiB)



```julia
# MM estimates
vcmodel_mm.B
```




    2×2 Array{Float64,2}:
     1.092     1.04727
     0.955346  1.01632




```julia
# MM estimates
vcmodel_mm.Σ
```




    ([0.380637 -0.305465; -0.305465 4.51938], [1.84009 0.265569; 0.265569 2.17275])



Fisher scoring (`mle_fs!`) uses either [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) (keyword `solver=:Ipopt`) or [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl) (keyword `solver=:Knitro`) as the backend solver. Ipopt is open source and installation of [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) package alone is sufficient.


```julia
# Fisher scoring using Ipopt
vcmodel_ipopt = deepcopy(vcmodel)
@time mle_fs!(vcmodel_ipopt, vcdatarot; solver=:Ipopt, maxiter=1000, verbose=true);
```

    This is Ipopt version 3.12.10, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:        0
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       21
    
    Total number of variables............................:        6
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:        0
                         variables with only upper bounds:        0
    Total number of equality constraints.................:        0
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  4.2109423e+03 0.00e+00 1.00e+02   0.0 0.00e+00    -  0.00e+00 0.00e+00   0 
       5  3.8445586e+03 0.00e+00 7.87e-01 -11.0 4.94e-02    -  1.00e+00 1.00e+00f  1 MaxS
      10  3.8443870e+03 0.00e+00 2.25e-01 -11.0 1.38e-02    -  1.00e+00 1.00e+00f  1 MaxS
      15  3.8443742e+03 0.00e+00 6.23e-02 -11.0 3.78e-03    -  1.00e+00 1.00e+00f  1 MaxS
      20  3.8443733e+03 0.00e+00 1.70e-02 -11.0 1.03e-03    -  1.00e+00 1.00e+00f  1 MaxS
      25  3.8443732e+03 0.00e+00 4.61e-03 -11.0 2.79e-04    -  1.00e+00 1.00e+00f  1 MaxS
      30  3.8443732e+03 0.00e+00 1.25e-03 -11.0 7.56e-05    -  1.00e+00 1.00e+00f  1 MaxS
      35  3.8443732e+03 0.00e+00 3.39e-04 -11.0 2.05e-05    -  1.00e+00 1.00e+00f  1 MaxS
      40  3.8443732e+03 0.00e+00 9.19e-05 -11.0 5.55e-06    -  1.00e+00 1.00e+00f  1 MaxS
      45  3.8443732e+03 0.00e+00 2.49e-05 -11.0 1.51e-06    -  1.00e+00 1.00e+00f  1 MaxS
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50  3.8443732e+03 0.00e+00 6.76e-06 -11.0 4.08e-07    -  1.00e+00 1.00e+00f  1 MaxSA
      55  3.8443732e+03 0.00e+00 1.83e-06 -11.0 1.11e-07    -  1.00e+00 1.00e+00f  1 MaxSA
      60  3.8443732e+03 0.00e+00 4.97e-07 -11.0 3.00e-08    -  1.00e+00 1.00e+00f  1 MaxSA
    
    Number of Iterations....: 63
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.4496886481728779e+02    3.8443731733053728e+03
    Dual infeasibility......:   2.2693631701157965e-07    2.5290047251948971e-06
    Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   2.2693631701157965e-07    2.5290047251948971e-06
    
    
    Number of objective function evaluations             = 64
    Number of objective gradient evaluations             = 64
    Number of equality constraint evaluations            = 0
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 0
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 63
    Total CPU secs in IPOPT (w/o function evaluations)   =      1.488
    Total CPU secs in NLP function evaluations           =      0.375
    
    EXIT: Solved To Acceptable Level.
      2.486697 seconds (4.13 M allocations: 201.947 MiB, 3.03% gc time)



```julia
# Ipopt estimates
vcmodel_ipopt.B
```




    2×2 Array{Float64,2}:
     1.092     1.04727
     0.955346  1.01632




```julia
# Ipopt estimates
vcmodel_ipopt.Σ
```




    ([0.380552 -0.305594; -0.305594 4.52106], [1.84008 0.265385; 0.265385 2.17287])



Knitro is a commercial software and users need to follow instructions at [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl) for proper functioning. Following code invokes Knitro as the backend optimization solver.
```julia
using KNITRO

# Fisher scoring using Knitro
vcmodel_knitro = deepcopy(vcmodel)
@time mle_fs!(vcmodel_knitro, vcdatarot; solver=:Knitro, maxiter=1000, verbose=true);

# Knitro estimates
vcmodel_knitro.B

# Knitro estimates
vcmodel_knitro.Σ
```

## Starting point

Here are a few strategies for successful optimization. 

* For $d>1$ (multivariate response), initialize $B, \Sigma$ from univariate estimates.  
* Use REML estimate as starting point for MLE.  
* When there are only $m=2$ variance components, pre-compute `TwoVarCompVariateRotate` and use it for optimization.

## Constrained estimation of `B`


Many applications invoke constraints on the mean parameters `B`. For demonstration, we enforce `B[1,1]=B[1,2]` and all entries of `B` are within [0, 2].


```julia
# set up constraints on B
vcmodel_constr = deepcopy(vcmodel)
vcmodel_constr.A = [1.0 0.0 -1.0 0.0]
vcmodel_constr.sense = '='
vcmodel_constr.b = 0.0
vcmodel_constr.lb = 0.0
vcmodel_constr.ub = 2.0
vcmodel_constr
```




    VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}([0.0 0.0; 0.0 0.0], ([1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0]), [1.0 0.0 -1.0 0.0], '=', 0.0, 0.0, 2.0)



We first try the MM algorithm.


```julia
# MM algorithm for constrained estimation of B
@time mle_mm!(vcmodel_constr, vcdatarot; maxiter=10000, funtol=1e-8, verbose = true);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -6.253551e+03
           1  -3.881820e+03
           2  -3.853477e+03
           3  -3.846807e+03
           4  -3.845184e+03
           5  -3.844783e+03
           6  -3.844683e+03
           7  -3.844658e+03
           8  -3.844652e+03
           9  -3.844650e+03
          10  -3.844650e+03
    
      0.170236 seconds (170.93 k allocations: 8.918 MiB)



```julia
fieldnames(typeof(vcmodel_constr))
```




    (:B, :Σ, :A, :sense, :b, :lb, :ub)




```julia
vcmodel_constr.B
```




    2×2 Array{Float64,2}:
     1.07177   1.07177
     0.955683  1.01591




```julia
vcmodel_constr.Σ
```




    ([0.380624 -0.305498; -0.305498 4.51948], [1.84051 0.265065; 0.265065 2.17336])



Now let's try Fisher scoring.


```julia
# Fisher scoring using Ipopt for constrained estimation of B
vcmodel_constr = deepcopy(vcmodel)
vcmodel_constr.A = [1.0 0.0 -1.0 0.0]
vcmodel_constr.sense = '='
vcmodel_constr.b = 0.0
vcmodel_constr.lb = 0.0
vcmodel_constr.ub = 2.0
vcmodel_constr
@time mle_fs!(vcmodel_constr, vcdatarot; solver=:Ipopt, maxiter=1000, verbose=true);
```

    This is Ipopt version 3.12.10, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:        0
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       21
    
    Total number of variables............................:        6
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:        0
                         variables with only upper bounds:        0
    Total number of equality constraints.................:        0
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  4.2114270e+03 0.00e+00 1.00e+02   0.0 0.00e+00    -  0.00e+00 0.00e+00   0 
       5  3.8448353e+03 0.00e+00 7.87e-01 -11.0 4.94e-02    -  1.00e+00 1.00e+00f  1 MaxS
      10  3.8446636e+03 0.00e+00 2.25e-01 -11.0 1.38e-02    -  1.00e+00 1.00e+00f  1 MaxS
      15  3.8446509e+03 0.00e+00 6.23e-02 -11.0 3.78e-03    -  1.00e+00 1.00e+00f  1 MaxS
      20  3.8446499e+03 0.00e+00 1.70e-02 -11.0 1.03e-03    -  1.00e+00 1.00e+00f  1 MaxS
      25  3.8446498e+03 0.00e+00 4.61e-03 -11.0 2.79e-04    -  1.00e+00 1.00e+00f  1 MaxS
      30  3.8446498e+03 0.00e+00 1.25e-03 -11.0 7.56e-05    -  1.00e+00 1.00e+00f  1 MaxS
      35  3.8446498e+03 0.00e+00 3.39e-04 -11.0 2.05e-05    -  1.00e+00 1.00e+00f  1 MaxS
      40  3.8446498e+03 0.00e+00 9.19e-05 -11.0 5.56e-06    -  1.00e+00 1.00e+00f  1 MaxS
      45  3.8446498e+03 0.00e+00 2.49e-05 -11.0 1.51e-06    -  1.00e+00 1.00e+00f  1 MaxS
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50  3.8446498e+03 0.00e+00 6.76e-06 -11.0 4.08e-07    -  1.00e+00 1.00e+00f  1 MaxSA
      55  3.8446498e+03 0.00e+00 1.83e-06 -11.0 1.11e-07    -  1.00e+00 1.00e+00h  1 MaxSA
      60  3.8446498e+03 0.00e+00 4.97e-07 -11.0 3.00e-08    -  1.00e+00 1.00e+00f  1 MaxSA
    
    Number of Iterations....: 63
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.4484507551949679e+02    3.8446498170293380e+03
    Dual infeasibility......:   2.2694405212011240e-07    2.5301808562731130e-06
    Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   2.2694405212011240e-07    2.5301808562731130e-06
    
    
    Number of objective function evaluations             = 64
    Number of objective gradient evaluations             = 64
    Number of equality constraint evaluations            = 0
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 0
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 63
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.024
    Total CPU secs in NLP function evaluations           =      0.625
    
    EXIT: Solved To Acceptable Level.
      0.729483 seconds (102.83 k allocations: 8.150 MiB, 0.88% gc time)



```julia
vcmodel_constr.B
```




    2×2 Array{Float64,2}:
     1.07177   1.07177
     0.955683  1.01591




```julia
vcmodel_constr.Σ
```




    ([0.380539 -0.305626; -0.305626 4.52116], [1.8405 0.264881; 0.264881 2.17348])


