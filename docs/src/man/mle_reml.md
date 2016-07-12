
# MLE and REML

## Demo data

For demonstration, we generate a random data set.


```julia
# generate data from a d-variate response variane component model
srand(123)
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
  copy!(V[i], Vi * Vi')
end
copy!(V[m], eye(n)) # last covarianec matrix is idendity
# a tuple of m d-by-d variance component parameters
Σ = ntuple(x -> zeros(d, d), m) 
for i in 1:m
  Σi = randn(d, d)
  copy!(Σ[i], Σi' * Σi)
end
# form overall nd-by-nd covariance matrix Ω
Ω = zeros(n * d, n * d)
for i = 1:m
  Ω += kron(Σ[i], V[i])
end
Ωchol = cholfact(Ω)
# n-by-d responses
Y = X * B + reshape(Ωchol[:L] * randn(n*d), n, d);
```

## Maximum likelihood estimation (MLE)

To find the MLE of parameters $(B,\Sigma_1,\ldots,\Sigma_m)$, we take 3 steps:  
**Step 1 (Construct data)**. Construct an instance of `VarianceComponentVariate`, which consists of the responses $Y$, the covariate matrix $X$, and a tuple of covariance matrices $V$. The last covariance matrix must be positive definite and usually is the identity matrix. In the absence of covariates $X$, we can simply initialize by `vcdata = VarianceComponentVariate(Y, V)`.


```julia
using VarianceComponentModels
vcdata = VarianceComponentVariate(Y, X, V)
fieldnames(vcdata)
```




    3-element Array{Symbol,1}:
     :Y
     :X
     :V



**Step 2 (Construct a model)**. Construct an instance of `VarianceComponentModel`. The fields `B` and `Σ` are mean and variane component parameters respectively. When constructed from a `VarianceComponentVariate` instance, the mean parameters $B$ are initialized to be zero and the tuple of variance component parameters $\Sigma$ to be `(eye(d),...,eye(d))`.

The remaining fields `A`, `sense`, `b`, `lb`, `ub` specify (optional) constraints on the mean parameters `B`:

$A * \text{vec}(B) =\ge\le b$

$lb \le \text{vec}(B) \le ub$

`A` is an `m x pd` constraint matrix, `sense` is a `m` vector of charaters taking values `<`, `=`, or `>`, and `lb` and `ub` are the lower and upper bounds for `vec(B)`. By default, `A`, `sense`, `b` are empty, `lb` is `-Inf`, and `ub` is `Inf`. If any constraits are non-trivial, final estimates of `B` are enforced to satisfy them.


```julia
vcmodel = VarianceComponentModel(vcdata)
fieldnames(vcmodel)
```




    7-element Array{Symbol,1}:
     :B    
     :Σ    
     :A    
     :sense
     :b    
     :lb   
     :ub   




```julia
vcmodel
```




    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(2x2 Array{Float64,2}:
     0.0  0.0
     0.0  0.0,(
    2x2 Array{Float64,2}:
     1.0  0.0
     0.0  1.0,
    
    2x2 Array{Float64,2}:
     1.0  0.0
     0.0  1.0),0x4 Array{Float64,2},Char[],Float64[],-Inf,Inf)



When a better initial guess is available, we can initialize by calling `vcmodel=VarianceComponentModel(B0, Σ0)` directly.

**Step 3 (Fit model)**. Call optmization routine `fit_mle!`. The keywork `algo` dictates the optimization algorithm: `:MM` (minorization-maximization algorithm) or `:FS` (Fisher scoring algorithm).


```julia
vcmodel_mle = deepcopy(vcmodel)
@time logl, vcmodel_mle, Σse, Σcov, Bse, Bcov = fit_mle!(vcmodel_mle, vcdata; algo = :MM);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -7.348297e+03
           1  -4.102367e+03
           2  -3.745567e+03
           3  -3.652392e+03
           4  -3.627744e+03
           5  -3.621170e+03
           6  -3.619381e+03
           7  -3.618878e+03
           8  -3.618730e+03
           9  -3.618684e+03
          10  -3.618670e+03
    
      6.187559 seconds (11.40 M allocations: 481.014 MB, 1.48% gc time)


The output of `fit_mle!` contains  

* final log-likelihood  


```julia
logl
```




    -3618.661776037447



* fitted model


```julia
fieldnames(vcmodel_mle)
```




    7-element Array{Symbol,1}:
     :B    
     :Σ    
     :A    
     :sense
     :b    
     :lb   
     :ub   




```julia
vcmodel_mle
```




    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(2x2 Array{Float64,2}:
     1.14929   1.00139
     0.924768  1.02485,(
    2x2 Array{Float64,2}:
      0.302279  -0.478398
     -0.478398   0.803237,
    
    2x2 Array{Float64,2}:
      5.86133   -0.586939
     -0.586939   0.586382),0x4 Array{Float64,2},Char[],Float64[],-Inf,Inf)



* standard errors of the estimated varianec component parameters


```julia
Σse
```




    (
    2x2 Array{Float64,2}:
     0.06167    0.0995727
     0.0995727  0.16077  ,
    
    2x2 Array{Float64,2}:
     0.268915  0.085063
     0.085063  0.026905)



* covariance matrix of the variance component parameters estimates


```julia
Σcov
```




    8x8 Array{Float64,2}:
      0.00380319  -0.00590812  -0.00590812  …   7.45979e-6   -7.52483e-7 
     -0.00590812   0.00991472   0.00917805     -7.52483e-7    7.55745e-7 
     -0.00590812   0.00917805   0.00991472     -7.49213e-6    7.55745e-7 
      0.00917805  -0.0154021   -0.0154021       7.55745e-7   -7.59023e-7 
     -7.39532e-5   7.45979e-6   7.45979e-6     -0.00724213    0.000725236
      7.45979e-6  -7.49213e-6  -7.52483e-7  …   0.000725236  -0.000724568
      7.45979e-6  -7.52483e-7  -7.49213e-6      0.00723572   -0.000724568
     -7.52483e-7   7.55745e-7   7.55745e-7     -0.000724568   0.000723881



* standard errors of the estimated mean parameters


```julia
Bse
```




    2x2 Array{Float64,2}:
     0.0760538  0.0251806
     0.0769259  0.0252831



* covariance matrix of the mean parameter estimates


```julia
Bcov
```




    4x4 Array{Float64,2}:
      0.00578418  -6.1406e-5    -0.00061133    3.61347e-6 
     -6.1406e-5    0.00591759    3.61347e-6   -0.000619951
     -0.00061133   3.61347e-6    0.000634064  -1.76949e-6 
      3.61347e-6  -0.000619951  -1.76949e-6    0.000639237



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
           0  -5.297304e+03
           1  -3.954293e+03
           2  -3.715870e+03
           3  -3.663101e+03
           4  -3.650121e+03
           5  -3.646663e+03
           6  -3.645672e+03
           7  -3.645367e+03
           8  -3.645268e+03
           9  -3.645233e+03
          10  -3.645221e+03
    
      0.462098 seconds (730.53 k allocations: 42.706 MB, 0.77% gc time)


The output of `fit_reml!` contains

* the final log-likelihood at REML estimate


```julia
logl
```




    -3622.050155483128



* REML estimates


```julia
fieldnames(vcmodel_reml)
```




    7-element Array{Symbol,1}:
     :B    
     :Σ    
     :A    
     :sense
     :b    
     :lb   
     :ub   




```julia
vcmodel_reml
```




    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(2x2 Array{Float64,2}:
     1.08201  1.05525
     0.90694  1.00679,(
    2x2 Array{Float64,2}:
      0.301641  -0.477617
     -0.477617   0.802105,
    
    2x2 Array{Float64,2}:
      5.88057   -0.610115
     -0.610115   0.61958 ),0x4 Array{Float64,2},Char[],Float64[],-Inf,Inf)



* standard errors of the estimated varianec component parameters


```julia
Σse
```




    (
    2x2 Array{Float64,2}:
     0.0615463  0.0994049
     0.0994049  0.160551 ,
    
    2x2 Array{Float64,2}:
     0.269798   0.0875812
     0.0875812  0.0284283)



* covariance matrix of the variance component parameters estimates


```julia
Σcov
```




    8x8 Array{Float64,2}:
      0.00378795  -0.00588696  -0.00588696  …   7.78026e-6   -8.13244e-7 
     -0.00588696   0.00988133   0.00914906     -8.13244e-7    8.30149e-7 
     -0.00588696   0.00914906   0.00988133     -7.94198e-6    8.30149e-7 
      0.00914906  -0.0153568   -0.0153568       8.30149e-7   -8.47406e-7 
     -7.44333e-5   7.78026e-6   7.78026e-6     -0.00755281    0.000783641
      7.78026e-6  -7.94198e-6  -8.13244e-7  …   0.000783641  -0.00079582 
      7.78026e-6  -8.13244e-7  -7.94198e-6      0.00767047   -0.00079582 
     -8.13244e-7   8.30149e-7   8.30149e-7     -0.00079582    0.000808168



* standard errors of the estimated mean parameters


```julia
Bse
```




    2x2 Array{Float64,2}:
     0.0615457  0.0994049
     0.0994049  0.160551 



* covariance matrix of the mean parameter estimates


```julia
Bcov
```




    4x4 Array{Float64,2}:
      0.00378788  -0.00588695  -0.00588695   0.00914906
     -0.00588695   0.00988132   0.00914906  -0.0153568 
     -0.00588695   0.00914906   0.00988132  -0.0153568 
      0.00914906  -0.0153568   -0.0153568    0.0257766 



## Optimization algorithms

Finding the MLE or REML of variance component models is a non-trivial nonlinear optimization problem. The main complications are the non-convexity of objective function and the positive semi-definiteness constraint of variane component parameters $\Sigma_1,\ldots,\Sigma_m$. Here are some tips for efficient computation. 

In general the optimization algorithm needs to invert the $nd$ by $nd$ overall covariance matrix $\Omega = \Sigma_1 \otimes V_1 + \cdots + \Sigma_m \otimes V_m$ in each iteration. Inverting a matrix is an expensive operation with $O(n^3 d^3)$ floating operations. When there are only **two** varianec components ($m=2$), this tedious task can be avoided by taking one (generalized) eigendecomposion of $(V_1, V_2)$ and rotating data $(Y, X)$ by the eigen-vectors. 


```julia
vcdatarot = TwoVarCompVariateRotate(vcdata)
fieldnames(vcdatarot)
```




    4-element Array{Symbol,1}:
     :Yrot    
     :Xrot    
     :eigval  
     :logdetV2



Two optimization algorithms are implemented: Fisher scoring (`mle_fs!`) and the [minorization-maximization (MM) algorithm](http://arxiv.org/abs/1509.07426) (`mle_mm!`). Both take the rotated data as input. These two functions give finer control of the optimization algorithms. Generally speaking, MM algorithm is more stable while Fisher scoring (if it converges) yields more accurate answer.


```julia
vcmodel_mm = deepcopy(vcmodel)
@time mle_mm!(vcmodel_mm, vcdatarot; maxiter=10000, funtol=1e-8, verbose = true);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -7.348297e+03
           1  -4.102367e+03
           2  -3.745567e+03
           3  -3.652392e+03
           4  -3.627744e+03
           5  -3.621170e+03
           6  -3.619381e+03
           7  -3.618878e+03
           8  -3.618730e+03
           9  -3.618684e+03
          10  -3.618670e+03
    
      0.166114 seconds (750.53 k allocations: 28.531 MB, 2.44% gc time)



```julia
# MM estimates
vcmodel_mm.B
```




    2x2 Array{Float64,2}:
     1.14929   1.00139
     0.924768  1.02485




```julia
# MM estimates
vcmodel_mm.Σ
```




    (
    2x2 Array{Float64,2}:
      0.302279  -0.478398
     -0.478398   0.803237,
    
    2x2 Array{Float64,2}:
      5.86133   -0.586939
     -0.586939   0.586382)



Fisher scoring (`mle_fs!`) uses either [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) (keyword `solver=:Ipopt`) or [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl) (keyword `solver=:Knitro`) as the backend solver. Ipopt is open source and installation of [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) package alone is sufficient. However Knitro is a commercial software and users need to follow instructions at [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl) for proper functioning.


```julia
# Fisher scoring using Ipopt
vcmodel_ipopt = deepcopy(vcmodel)
@time mle_fs!(vcmodel_ipopt, vcdatarot; solver=:Ipopt, maxiter=1000, verbose=true);
```

    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    ******************************************************************************
    
    This is Ipopt version 3.12.4, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:        0
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       21
    
    Total number of variables............................:        6
                         variables with only lower bounds:        4
                    variables with lower and upper bounds:        0
                         variables with only upper bounds:        0
    Total number of equality constraints.................:        0
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.8183666e+03 0.00e+00 9.89e+00   0.0 0.00e+00    -  0.00e+00 0.00e+00   0 
       5  3.7059581e+03 0.00e+00 5.96e+00  -4.4 7.13e-01    -  1.00e+00 1.00e+00f  1  sigma=8.40e-03MaxS
      10  3.6256399e+03 0.00e+00 2.71e+00 -11.0 1.04e-01    -  1.00e+00 1.00e+00f  1  sigma=3.37e-03MaxS
      15  3.6186620e+03 0.00e+00 2.77e-02 -11.0 1.51e-03    -  1.00e+00 1.00e+00f  1  sigma=5.30e-09MaxS
      20  3.6186618e+03 0.00e+00 5.95e-05 -11.0 3.16e-06    -  1.00e+00 1.00e+00f  1  sigma=5.41e-17MaxS
      25  3.6186618e+03 0.00e+00 1.27e-07 -11.0 6.76e-09    -  1.00e+00 1.00e+00f  1  sigma=5.35e-25MaxSA
    
    Number of Iterations....: 28
    
                                       (scaled)                 (unscaled)
    Objective...............:   7.8367920603164933e+01    3.6186617667669216e+03
    Dual infeasibility......:   3.1921838841006262e-09    1.4740002905497475e-07
    Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
    Complementarity.........:   1.0000000000010052e-11    4.6175293907497364e-10
    Overall NLP error.......:   3.1921838841006262e-09    1.4740002905497475e-07
    
    
    Number of objective function evaluations             = 29
    Number of objective gradient evaluations             = 29
    Number of equality constraint evaluations            = 0
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 0
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 28
    Total CPU secs in IPOPT (w/o function evaluations)   =      2.334
    Total CPU secs in NLP function evaluations           =      0.244
    
    EXIT: Optimal Solution Found.
      2.987685 seconds (8.12 M allocations: 292.023 MB, 2.09% gc time)



```julia
# Ipopt estimates
vcmodel_ipopt.B
```




    2x2 Array{Float64,2}:
     1.14928   1.00139
     0.924761  1.02485




```julia
# Ipopt estimates
vcmodel_ipopt.Σ
```




    (
    2x2 Array{Float64,2}:
      0.30227   -0.478407
     -0.478407   0.803239,
    
    2x2 Array{Float64,2}:
      5.86164   -0.586964
     -0.586964   0.586375)




```julia
# Fisher scoring using Knitro
vcmodel_knitro = deepcopy(vcmodel)
@time mle_fs!(vcmodel_knitro, vcdatarot; solver=:Knitro, maxiter=1000, verbose=true);
```

    
    Knitro 10.1.0 STUDENT LICENSE (problem size limit = 300)
    
    =======================================
                Student License
           (NOT FOR COMMERCIAL USE)
             Artelys Knitro 10.1.0
    =======================================
    
    Knitro presolve eliminated 0 variables and 0 constraints.
    
    algorithm:            1
    The problem is identified as bound constrained only.
    Knitro changing bar_initpt from AUTO to 3.
    Knitro changing bar_murule from AUTO to 4.
    Knitro changing bar_penaltycons from AUTO to 1.
    Knitro changing bar_penaltyrule from AUTO to 2.
    Knitro changing bar_switchrule from AUTO to 1.
    Knitro changing linsolver from AUTO to 2.
    
    Problem Characteristics                    ( Presolved)
    -----------------------
    Objective goal:  Maximize
    Number of variables:                     6 (         6)
        bounded below:                       4 (         4)
        bounded above:                       0 (         0)
        bounded below and above:             0 (         0)
        fixed:                               0 (         0)
        free:                                2 (         2)
    Number of constraints:                   0 (         0)
        linear equalities:                   0 (         0)
        nonlinear equalities:                0 (         0)
        linear inequalities:                 0 (         0)
        nonlinear inequalities:              0 (         0)
        range:                               0 (         0)
    Number of nonzeros in Jacobian:          0 (         0)
    Number of nonzeros in Hessian:          21 (        21)
    
      Iter      Objective      FeasError   OptError    ||Step||    CGits 
    --------  --------------  ----------  ----------  ----------  -------
           0   -5.272317e+03   0.000e+00
          10   -3.618665e+03   0.000e+00   1.050e-01   6.061e-03        0
          20   -3.618662e+03   0.000e+00   4.786e-07   2.620e-08        0
    
    EXIT: Locally optimal solution found.
    
    Final Statistics
    ----------------
    Final objective value               =  -3.61866176676693e+03
    Final feasibility error (abs / rel) =   0.00e+00 / 0.00e+00
    Final optimality error  (abs / rel) =   4.79e-07 / 4.79e-07
    # of iterations                     =         20 
    # of CG iterations                  =          0 
    # of function evaluations           =         22
    # of gradient evaluations           =         22
    # of Hessian evaluations            =         20
    Total program time (secs)           =       0.17739 (     0.177 CPU time)
    Time spent in evaluations (secs)    =       0.11814
    
    ===============================================================================
    
      0.550156 seconds (2.93 M allocations: 71.557 MB, 1.60% gc time)


    ### Could not find a valid license.
        Your machine ID is 1f-aa-f6-5b-46.
        Please contact licensing@artelys.com or your local distributor to obtain a license.
        If you already have a license, please execute `get_machine_ID -v` and send the output to support.



```julia
# Knitro estimates
vcmodel_knitro.B
```




    2x2 Array{Float64,2}:
     1.14928   1.00139
     0.924761  1.02485




```julia
# Knitro estimates
vcmodel_knitro.Σ
```




    (
    2x2 Array{Float64,2}:
      0.30227   -0.478407
     -0.478407   0.803239,
    
    2x2 Array{Float64,2}:
      5.86164   -0.586964
     -0.586964   0.586375)



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




    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(2x2 Array{Float64,2}:
     0.0  0.0
     0.0  0.0,(
    2x2 Array{Float64,2}:
     1.0  0.0
     0.0  1.0,
    
    2x2 Array{Float64,2}:
     1.0  0.0
     0.0  1.0),1x4 Array{Float64,2}:
     1.0  0.0  -1.0  0.0,'=',0.0,0.0,2.0)



We first try the MM algorithm.


```julia
# MM algorithm for constrained estimation of B
@time mle_mm!(vcmodel_constr, vcdatarot; maxiter=10000, funtol=1e-8, verbose = true);
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -7.348297e+03
           1  -4.105736e+03
           2  -3.747828e+03
           3  -3.654068e+03
           4  -3.629236e+03
           5  -3.622605e+03
           6  -3.620799e+03
           7  -3.620290e+03
           8  -3.620140e+03
           9  -3.620094e+03
          10  -3.620079e+03
    
      0.407094 seconds (618.18 k allocations: 20.440 MB)



```julia
fieldnames(vcmodel_constr)
```




    7-element Array{Symbol,1}:
     :B    
     :Σ    
     :A    
     :sense
     :b    
     :lb   
     :ub   




```julia
vcmodel_constr.B
```




    2x2 Array{Float64,2}:
     1.02428   1.02428
     0.925434  1.02474




```julia
vcmodel_constr.Σ
```




    (
    2x2 Array{Float64,2}:
      0.301672  -0.478044
     -0.478044   0.803126,
    
    2x2 Array{Float64,2}:
      5.88053   -0.590251
     -0.590251   0.58695 )



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

    This is Ipopt version 3.12.4, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:        0
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       21
    
    Total number of variables............................:        6
                         variables with only lower bounds:        4
                    variables with lower and upper bounds:        0
                         variables with only upper bounds:        0
    Total number of equality constraints.................:        0
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.8184296e+03 0.00e+00 9.89e+00   0.0 0.00e+00    -  0.00e+00 0.00e+00   0 
       5  3.7077312e+03 0.00e+00 5.97e+00  -4.4 7.13e-01    -  1.00e+00 1.00e+00f  1  sigma=8.41e-03MaxS
      10  3.6271375e+03 0.00e+00 2.72e+00 -11.0 1.04e-01    -  1.00e+00 1.00e+00f  1  sigma=3.39e-03MaxS
      15  3.6200706e+03 0.00e+00 2.76e-02 -11.0 1.51e-03    -  1.00e+00 1.00e+00f  1  sigma=5.40e-09MaxS
      20  3.6200704e+03 0.00e+00 5.77e-05 -11.0 3.08e-06    -  1.00e+00 1.00e+00f  1  sigma=5.11e-17MaxS
      25  3.6200704e+03 0.00e+00 1.21e-07 -11.0 6.43e-09    -  1.00e+00 1.00e+00f  1  sigma=4.70e-25MaxSA
    
    Number of Iterations....: 28
    
                                       (scaled)                 (unscaled)
    Objective...............:   7.8291178537553890e+01    3.6200703544238395e+03
    Dual infeasibility......:   2.9789254301703743e-09    1.3774118411854456e-07
    Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
    Complementarity.........:   1.0000000000009410e-11    4.6238547203511678e-10
    Overall NLP error.......:   2.9789254301703743e-09    1.3774118411854456e-07
    
    
    Number of objective function evaluations             = 29
    Number of objective gradient evaluations             = 29
    Number of equality constraint evaluations            = 0
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 0
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 28
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.463
    
    EXIT: Optimal Solution Found.
      0.522384 seconds (3.67 M allocations: 70.253 MB, 9.65% gc time)



```julia
vcmodel_constr.B
```




    2x2 Array{Float64,2}:
     1.02427   1.02427
     0.925427  1.02474




```julia
vcmodel_constr.Σ
```




    (
    2x2 Array{Float64,2}:
      0.301662  -0.478052
     -0.478052   0.803128,
    
    2x2 Array{Float64,2}:
      5.88085   -0.590275
     -0.590275   0.586943)


