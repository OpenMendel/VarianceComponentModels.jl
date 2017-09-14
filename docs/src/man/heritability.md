
# Heritability Analysis

As an application of the variance component model, this note demonstrates the workflow for heritability analysis in genetics, using a sample data set `cg10k` with **6,670** individuals and **630,860** SNPs. Person IDs and phenotype names are masked for privacy. `cg10k.bed`, `cg10k.bim`, and `cg10k.fam` is a set of Plink files in binary format. `cg10k_traits.txt` contains 13 phenotypes of the 6,670 individuals.


```julia
;ls cg10k.bed cg10k.bim cg10k.fam cg10k_traits.txt
```

    cg10k.bed
    cg10k.bim
    cg10k.fam
    cg10k_traits.txt


Machine information:


```julia
versioninfo()
```

    Julia Version 0.6.0
    Commit 903644385b (2017-06-19 13:05 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin13.4.0)
      CPU: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
      WORD_SIZE: 64
      BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
      LAPACK: libopenblas64_
      LIBM: libopenlibm
      LLVM: libLLVM-3.9.1 (ORCJIT, haswell)


## Read in binary SNP data

We will use the [`SnpArrays.jl`](https://github.com/OpenMendel/SnpArrays.jl) package to read in binary SNP data and compute the empirical kinship matrix. Issue 
```julia
Pkg.clone("https://github.com/OpenMendel/SnpArrays.jl.git")
```
within `Julia` to install the `SnpArrays` package.


```julia
using SnpArrays
```


```julia
# read in genotype data from Plink binary file (~50 secs on my laptop)
@time cg10k = SnpArray("cg10k")
```

     22.902730 seconds (51.62 k allocations: 1005.845 MiB, 0.11% gc time)





    6670×630860 SnpArrays.SnpArray{2}:
     (false, true)   (false, true)   …  (true, true)    (true, true) 
     (true, true)    (true, true)       (false, true)   (true, false)
     (true, true)    (true, true)       (true, true)    (true, true) 
     (true, true)    (true, true)       (false, true)   (true, true) 
     (true, true)    (true, true)       (true, true)    (false, true)
     (false, true)   (false, true)   …  (true, true)    (true, true) 
     (false, false)  (false, false)     (true, true)    (true, true) 
     (true, true)    (true, true)       (true, true)    (false, true)
     (true, true)    (true, true)       (true, true)    (true, true) 
     (true, true)    (true, true)       (false, true)   (true, true) 
     (true, true)    (true, true)    …  (true, true)    (true, true) 
     (false, true)   (false, true)      (true, true)    (false, true)
     (true, true)    (true, true)       (true, true)    (false, true)
     ⋮                               ⋱                               
     (false, true)   (false, true)      (false, true)   (false, true)
     (false, true)   (false, true)      (false, true)   (true, true) 
     (true, true)    (true, true)    …  (false, true)   (true, true) 
     (false, true)   (false, true)      (true, true)    (false, true)
     (true, true)    (true, true)       (false, true)   (true, true) 
     (true, true)    (true, true)       (false, false)  (false, true)
     (true, true)    (true, true)       (true, true)    (false, true)
     (true, true)    (true, true)    …  (true, true)    (true, true) 
     (true, true)    (true, true)       (false, true)   (true, true) 
     (true, true)    (true, true)       (true, true)    (false, true)
     (false, true)   (false, true)      (true, true)    (true, true) 
     (true, true)    (true, true)       (true, true)    (true, true) 



## Summary statistics of SNP data


```julia
people, snps = size(cg10k)
```




    (6670, 630860)




```julia
# summary statistics (~50 secs on my laptop)
@time maf, _, missings_by_snp, = summarize(cg10k);
```

     24


```julia
# 5 number summary and average MAF (minor allele frequencies)
quantile(maf, [0.0 .25 .5 .75 1.0]), mean(maf)
```




    ([0.00841726 0.124063 … 0.364253 0.5], 0.24536516625042462)




```julia
# Pkg.add("Plots")
# Pkg.add("PyPlot")
using Plots
pyplot()

histogram(maf, xlab = "Minor Allele Frequency (MAF)", label = "MAF")
```




<img src="data:image/png;base64,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+/1bb70lh8Oh7Oxsr85z/nEuRnR0tBwOR4nXCy+8YOT4l7JBgwbpkUcecb9/7LHH3A84Lm0Jxeuvv14Oh0P33HNPqcfr37+/HA6HUlNTS5QVFhaW+nN0OBz68MMPJUlvv/22YmNjS73GgOqGgAUYEhwcrOzsbG3YsMFj+969e/Xll1+WCFP33HOPPvvss6psYoVSUlIkSb/88ou+/vprH7fGvqFDh+rLL795Da8zAAAOYUlEQVT0eI0cOdLXzarRvvnmG61Zs0aTJk3y2O50OnXw4EFt3rzZY/svv/yib775psz/acjOzta6deskSUuWLCnzvJMmTSrxs+zWrZskaeTIkTp27Jj79xSozniqMWBIrVq11KtXL6Wmpqp79+7u7UuXLlW7du0UGBjoUb9p06Zq2rRplbXv1KlTuuyyy8osLyoqUlpamrp27apvv/1WKSkpuv7666usfRejcePG6ty5s+36FX0XkF5++WX95S9/UePGjT22165dW127dlVqaqq6du3q3r506VLFxcXp7NmzpR7vvffeU2FhoXr16qW1a9fqyJEjCgsLK1EvKiqqzJ9lUFCQRo8erZdfflkJCQkX8emAykcPFmDQ8OHDtWzZMo8/MkuWLNGIESNK1D1/iPDzzz+Xw+HQJ598ohEjRqh+/fqKiorSc889V2LfDz74QHFxcapdu7YiIiKUlJSk06dPlzjWxx9/rKFDhyokJETDhg0rt+0bNmzQ/v37NXHiRPXv319paWk6d+6c199BQUGBHn30UUVFRcnpdKpt27bl9lgU279/v0aNGqWGDRvqsssu080336x//etfXp//fJ9++qkcDodWr16twYMHq379+ho+fLi7fOHChWrfvr2cTqeaNGmi6dOnl/jcGzduVMeOHVW7dm21b99ea9euVWxsrMdQWNeuXTVw4ECP/b799ls5HA5t2rTJva2oqEjPPvusWrZsKafTqebNm+uVV17x2O+xxx5TgwYN9OOPP+rGG29UnTp11L59e3366aclPt+iRYvcvwuNGjXSgAEDtG/fPh06dEi1atXSokWLSuzzpz/9qdTfyWLHjh3TihUrNHTo0FLLhw8frvfff1+FhYXubampqeUec8mSJWrTpo2ef/55nT17Vu+9916ZdcszbNgw/etf/9LPP/98QfsDVYWABRh02223qaCgwD0UsnXrVm3ZskV33XWX7WNMnDhRrVq10ooVK3Tbbbdp6tSpWrNmjbt81apVGjp0qK6++mp9+OGHmjJliubNm6dRo0aVONb48eMVExOjFStW6KGHHir3vCkpKapTp44GDhyoESNG6NChQ6X+Qa/InXfeqTfeeEMPPvigPvroI/Xt21ejRo3S6tWry9wnNzdXXbt21Q8//KBXX31Vy5cvV926dXXLLbfo0KFDFZ7TsiwVFha6X6UFw3Hjxql169b68MMP9cADD0iSnnvuOU2YMEH9+/fXRx99pMmTJ+ull17S448/7t7v119/Vd++fVW3bl299957evDBBzVhwgRlZWV5/d1I0r333qsnn3xSY8aM0ccff6y//vWvevDBB7VgwQKPegUFBUpISNDYsWP1wQcfyOVyafDgwcrNzXXXmT17tsaMGaNOnTrpgw8+0IIFC9S8eXNlZ2friiuu0O23366FCxd6HHfLli367rvvNHbs2DLbuHnzZp06dUpdunQptfyOO+7Q8ePHtX79evcxt27dWubveUZGhv75z39qxIgRiouLKzd0FxUVlfuzjI2NVUhIiD755JMy2w9UCxaAi/b4449bdevWtSzLskaMGGGNGjXKsizLeuyxx6wbbrjBsizLuuaaa6zRo0eXuo9lWdZnn31mSbImT57s3lZUVGRFR0dbY8eOdW/r2LGj+5jF3njjDUuStWXLFo9jTZw40Vb7CwoKrMsvv9y66667LMuyrNOnT1uhoaFWQkKCR73Ro0db7dq1c79ftGiRJck6fPiwZVmW9b//+7+WJGvt2rUe+8XHx1vXXXddmceZMWOGFRoaav3222/ubadPn7YiIyM9vo/SREVFWZI8XoGBge7yTz75xJJk3XvvvR775ebmWnXq1LGmT5/usf3VV1+16tSpY+Xm5lqWZVkPPvigFRoaauXn57vrrF271pLk8XPp0qWLdccdd3gcKz093ZJkbdy40bIsy9qxY4flcDis5ORkj3oPPvig1aRJE6uoqMiyLMuaNm1aie9x165dliQrNTXVsizLysnJsWrXrm397W9/K/O7WbNmjSXJ2rlzp3vb/fffb0VHR7vPVZonn3zSatCgQYnt06ZNs0JDQy3Lsqw777zTSkxMtCzLsh5++GHrpptusizLstq1a+fxvViWZf3jH/+wJFm7d++2LMuynnrqKcvhcFgZGRnuOmfPni3xc5RktW7dukQ7unTp4v5dBaorerAAw4YPH66VK1fq1KlTWrp0qcdwlB29e/d2/9vhcKht27bav3+/JOn48eP64YcfSgzdxMfHS5LHUJT0+11bdqxevVq5ubnuIR6n06nBgwdrxYoVOnXqlO22r1u3Ti6XS7fccotHL8Stt96q77//vswhx3Xr1qlHjx5yuVzufQIDA9WtWzelp6dXeN4777xT6enp7ldpE/TP/y42b96skydPatiwYR5t7dWrl06ePOkegvr666/Vs2dP1a9f371v7969FRISYvt7KfbJJ5/I4XBo8ODBJc7566+/6sCBA+66gYGBuuWWW9zvW7RooVq1arl/FzZv3qzTp0+X2xN16623Kioqyt2LdebMGaWkpOjuu++Ww+Eoc7+srCw1bNiw3M8yfPhwrVixQqdPn67w93zJkiXq1KmTYmJiJEkjRoyQZVml9mIlJSV5/CyXL19eok7Dhg0vuAcRqCpMcgcM69Onj4KDgzVjxgzt2bNHd955p1f7N2jQwON9rVq1dPToUUnS0aNHZVlWiYnHoaGhcjqdysnJ8dh+fr2ypKSkKDQ0VJ07d3afa8CAAVq0aJFWrVrlDnAVyc7OVk5OjoKDg0stz8rKKnVif3Z2tr766qtS9yv+o1yeRo0a6c9//nO5dc7/LoqXlujQoUOp9fft2+duc2lLSlxxxRUVtut82dnZKioq0uWXX17mOZs0aSJJqlevnoKCPP8THRwc7J5rd+TIEUlSREREmecLCAjQ2LFj9frrr2vWrFlatWqVcnNzPZYLKc3p06fldDrLrdOvXz9J0vTp07V///4y5/ht2bJFP/30k55++mn375bL5dK1116rJUuWeCwDIUnNmjWr8GfpdDq9Cv6ALxCwAMOCg4M1ZMgQvfjii+rZs6ftkGNHgwYN5HA4SsxLysvLU0FBgVwul8f28nopih07dkwfffSRTp06VWpoSElJsR2wXC6XGjVqpP/5n/8ptbysUOJyudS3b1899dRTJcoq+kNv1/nfRfF3tXLlylJDSvPmzSVJ4eHhpc4DO39b7dq1debMGY9tf5wvVXzOgIAAbd68uUR4kqQ2bdrY+CS/K74D78CBA7ryyivLrDdmzBg98cQTWr16tRYuXKhevXopMjKy3GO7XC53GCpLcS/niy++qD59+pTZ41XcSzVt2jRNmzatRPmWLVvKDLllOXr0aKl3IALVCQELqAT33HOPDh06pHHjxhk9br169RQXF6dly5a5J2pLct+R9cfb5u0qHgacN2+eWrdu7VH21ltvacmSJcrJySkR3krTq1cvPffcc6pVq5ZXfzR79eqlxYsXq23btqpbt67Xn+FCdOnSRbVr19avv/6q22+/vcx6nTp1UnJyso4dO+YeJly3bp3y8/M96jVt2lQbN26UZVnuMFd8s0Oxnj17qqioSLm5ue4eoItt/6JFi3TttdeWWa9Jkybq16+fZs+era+//trWGlKtW7fWb7/9VuFyFuPGjVNOTo4mTJhQarllWUpNTVWXLl00a9Ysj7LTp0/rtttuU0pKitcBKyMjQ3/5y1+82geoagQsoBJ06tTJvfq0aTNnztTAgQM1atQojRo1Sjt27NCjjz6qIUOGqH379l4fLyUlRVFRURo/fnypvTxvv/223n///TL/iP7Rrbfeqttuu019+/bVlClT1KFDB504cUI///yzdu/eXeJOuWJJSUlKSUlRt27d9Pe//12RkZE6fPiwvv76a0VERHiESVNcLpdmzpyppKQkZWZmqlu3bgoICNAvv/yiDz/8UKtWrZLT6dQDDzyg119/Xf369dOUKVOUk5OjmTNnlgicQ4cO1dtvv61Jkybptttu06ZNm7RixQqPOldffbUmTpyokSNH6qGHHlKnTp105swZ7dixQxs3btQHH3xgu/2XX365pk+frmnTpqmwsFC33367zp07p/Xr1+uvf/2rOnbs6K47btw43XHHHXK5XCWWkihNly5dVFhYqB9//LHc9cVuuOGGcn/PN23apMzMTM2aNctjbbhi/fr109KlS/XMM89U2KZi+fn52rVrl2666Sbb+wC+wCR3oIa5/fbb9f777+vf//637rjjDj3zzDMaP368Fi9e7PWxDh06pPXr1yshIaHU4cQOHTooLi7Oq5Wzly1bpokTJ2ru3Lnq16+fxo4dq3Xr1rlX4y5NWFiYvvrqK8XFxWnq1Knq3bu3HnjgAWVkZFTqYqdTp07VggUL9Omnn2rQoEEaNmyYFixYoM6dO7vngzVt2lSrV6/W8ePHNWzYMD3//PN6/fXXFR4e7nGsAQMGaPbs2VqxYoUGDhyo7du36/XXXy9xzrlz52rmzJlasmSJ+vfvr4SEBC1btqzUAFKRRx99VPPnz9emTZs0cOBA3X333frll1/UqFEjj3r9+vWT0+nUyJEjbQ25Xn311Wrbtm25S2vYsWTJEtWrV0+DBw8utXz06NHKzMzUxo0bbR9zzZo1qlevnvr06XNRbQMqm8OySnmgFACgXLGxsercuXOZvXLVybp169SnTx/98MMPuuaaa2zt89JLL2nevHnasWNHJbfOO4MGDVKjRo305ptv+ropQLnowQIAP3XgwAF98cUXmjp1qrp162Y7XEnShAkTlJeXV+YNC76we/durV27Vo8++qivmwJUiIAFAH5q7ty56tmzpwIDA73u8alTp47efvttj0cw+dqBAwe0YMECRUdH+7opQIX+H8eyLemiojkbAAAAAElFTkSuQmCC" />




```julia
# proportion of missing genotypes
sum(missings_by_snp) / length(cg10k)
```




    0.0013128198764010824




```julia
# proportion of rare SNPs with maf < 0.05
countnz(maf .< 0.05) / length(maf)
```




    0.07228069619249913



## Empirical kinship matrix

We estimate empirical kinship based on all SNPs by the genetic relation matrix (GRM). Missing genotypes are imputed on the fly by drawing according to the minor allele frequencies.


```julia
# GRM using SNPs with maf > 0.01 (default) (~10 mins on my laptop)
srand(123)
@time Φgrm = grm(cg10k; method = :GRM)
```

    396.943890 seconds (8.43 G allocations: 127.378 GiB, 4.38% gc time)





    6670×6670 Array{Float64,2}:
      0.503024      0.00335505   -0.000120075  …  -5.45185e-5   -0.00278072 
      0.00335505    0.498958     -0.00195952       0.000868471   0.0034285  
     -0.000120075  -0.00195952    0.493828         0.000174648  -0.000381467
      0.000923828  -0.00329169   -0.00194166      -0.00223595   -0.00123508 
     -8.39649e-5   -0.00353358    0.0018709        0.00222858   -0.00171176 
      0.00204208    0.000572952   0.00254025   …   0.000861385   2.99785e-5 
      0.000569323   0.0024786    -0.00185743       0.00117649   -0.00118027 
     -0.000642144   0.00317992   -0.00099777       0.00354182   -0.000260645
     -0.00102913   -0.00123475   -0.00061138       0.00173885    0.00177727 
     -0.00139442    0.00208423    0.000124525     -0.00145156   -0.001011   
     -0.00204555    0.00011055   -0.000419398  …  -0.000198235  -0.00110353 
      0.000947587   0.00167346    0.00184451      -0.000690143  -0.00304087 
      0.000322759  -0.000899805   0.00303981       0.000739331  -0.00118835 
      ⋮                                        ⋱                            
      0.00298012    0.00130003    0.000998861      4.18454e-6    0.00303991 
     -0.00207748    0.00274717   -0.00191741      -0.00107073    0.00368267 
      0.000545569  -0.00244439   -0.00299578   …  -0.000669885   0.00221027 
     -0.00423186   -0.00208514   -0.00108833      -0.000622127  -0.000567483
     -0.00325644   -0.000781353   0.0030423        0.000501423  -0.00010267 
      0.00041055   -0.00200772    0.00274867      -0.00624933   -0.00521365 
      0.00210519    0.000879889  -0.00107817      -0.000797878  -0.000557352
     -0.00230058   -0.000119132   0.000116817  …   0.000867087  -0.00233512 
     -0.0020119     0.00230772   -0.00128837       0.00194798   -0.00048733 
     -0.000944942  -0.000928073  -0.000175096      0.00126911   -0.00303766 
     -5.45185e-5    0.000868471   0.000174648      0.500829      0.000469478
     -0.00278072    0.0034285    -0.000381467      0.000469478   0.500627   



## Phenotypes

Read in the phenotype data and compute descriptive statistics.


```julia
# Pkg.add("DataFrames")
using DataFrames

cg10k_trait = readtable(
    "cg10k_traits.txt"; 
    separator = ' ',
    names = [:FID; :IID; :Trait1; :Trait2; :Trait3; :Trait4; :Trait5; :Trait6; 
             :Trait7; :Trait8; :Trait9; :Trait10; :Trait11; :Trait12; :Trait13],  
    eltypes = [String; String; Float64; Float64; Float64; Float64; Float64; 
               Float64; Float64; Float64; Float64; Float64; Float64; Float64; Float64]
    )
# do not display FID and IID for privacy
cg10k_trait[:, 3:end]
```




<table class="data-frame"><thead><tr><th></th><th>Trait1</th><th>Trait2</th><th>Trait3</th><th>Trait4</th><th>Trait5</th><th>Trait6</th><th>Trait7</th><th>Trait8</th><th>Trait9</th><th>Trait10</th><th>Trait11</th><th>Trait12</th><th>Trait13</th></tr></thead><tbody><tr><th>1</th><td>-1.81573145026234</td><td>-0.94615046147283</td><td>1.11363077580442</td><td>-2.09867121119159</td><td>0.744416614111748</td><td>0.00139171884080131</td><td>0.934732480409667</td><td>-1.22677315418103</td><td>1.1160784277875</td><td>-0.4436280335029</td><td>0.824465656443384</td><td>-1.02852542216546</td><td>-0.394049201727681</td></tr><tr><th>2</th><td>-1.24440094378729</td><td>0.109659992547179</td><td>0.467119394241789</td><td>-1.62131304097589</td><td>1.0566758355683</td><td>0.978946979419181</td><td>1.00014633946047</td><td>0.32487427140228</td><td>1.16232175219696</td><td>2.6922706948705</td><td>3.08263672461047</td><td>1.09064954786013</td><td>0.0256616415357438</td></tr><tr><th>3</th><td>1.45566914502305</td><td>1.53866932923243</td><td>1.09402959376555</td><td>0.586655272226893</td><td>-0.32796454430367</td><td>-0.30337709778827</td><td>-0.0334354881314741</td><td>-0.464463064285437</td><td>-0.3319396273436</td><td>-0.486839089635991</td><td>-1.10648681564373</td><td>-1.42015780427231</td><td>-0.687463456644413</td></tr><tr><th>4</th><td>-0.768809276698548</td><td>0.513490885514249</td><td>0.244263028382142</td><td>-1.31740254475691</td><td>1.19393774326845</td><td>1.17344127734288</td><td>1.08737426675232</td><td>0.536022583732261</td><td>0.802759240762068</td><td>0.234159411749815</td><td>0.394174866891074</td><td>-0.767365892476029</td><td>0.0635385761884935</td></tr><tr><th>5</th><td>-0.264415132547719</td><td>-0.348240421825694</td><td>-0.0239065083413606</td><td>0.00473915802244948</td><td>1.25619191712193</td><td>1.2038883667631</td><td>1.29800739042627</td><td>0.310113660247311</td><td>0.626159861059352</td><td>0.899289129831224</td><td>0.54996783350812</td><td>0.540687809542048</td><td>0.179675416046033</td></tr><tr><th>6</th><td>-1.37617270917293</td><td>-1.47191967744564</td><td>0.291179894254146</td><td>-0.803110740704731</td><td>-0.264239977442213</td><td>-0.260573027836772</td><td>-0.165372266287781</td><td>-0.219257294118362</td><td>1.04702422290318</td><td>-0.0985815534616482</td><td>0.947393438068448</td><td>0.594014812031438</td><td>0.245407436348479</td></tr><tr><th>7</th><td>0.1009416296374</td><td>-0.191615722103455</td><td>-0.567421321596677</td><td>0.378571487240382</td><td>-0.246656179817904</td><td>-0.608810750053858</td><td>0.189081058215596</td><td>-1.27077787326519</td><td>-0.452476199143965</td><td>0.702562877297724</td><td>0.332636218957179</td><td>0.0026916503626181</td><td>0.317117176705358</td></tr><tr><th>8</th><td>-0.319818276367464</td><td>1.35774480657283</td><td>0.818689545938528</td><td>-1.15565531644352</td><td>0.63448368102259</td><td>0.291461908634679</td><td>0.933323714954726</td><td>-0.741083289682492</td><td>0.647477683507572</td><td>-0.970877627077966</td><td>0.220861165411304</td><td>0.852512250237764</td><td>-0.225904624283945</td></tr><tr><th>9</th><td>-0.288334173342032</td><td>0.566082538090752</td><td>0.254958336116175</td><td>-0.652578302869714</td><td>0.668921559277347</td><td>0.978309199170558</td><td>0.122862966041938</td><td>1.4790926378214</td><td>0.0672132424173449</td><td>0.0795903917527827</td><td>0.167532455243232</td><td>0.246915579442139</td><td>0.539932616458363</td></tr><tr><th>10</th><td>-1.15759732583991</td><td>-0.781198583545165</td><td>-0.595807759833517</td><td>-1.00554980260402</td><td>0.789828885933321</td><td>0.571058413379044</td><td>0.951304176233755</td><td>-0.295962982984816</td><td>0.99042002479707</td><td>0.561309366988983</td><td>0.733100030623233</td><td>-1.73467772245684</td><td>-1.35278484330654</td></tr><tr><th>11</th><td>0.740569150459031</td><td>1.40873846755415</td><td>0.734689999440088</td><td>0.0208322841295094</td><td>-0.337440968561619</td><td>-0.458304040611395</td><td>-0.142582512772326</td><td>-0.580392297464107</td><td>-0.684684998101516</td><td>-0.00785381461893456</td><td>-0.712244337518008</td><td>-0.313345561230878</td><td>-0.345419463162219</td></tr><tr><th>12</th><td>-0.675892486454995</td><td>0.279892613829682</td><td>0.267915996308248</td><td>-1.04103665392985</td><td>0.910741715645888</td><td>0.866027618513171</td><td>1.07414431702005</td><td>0.0381751003538302</td><td>0.766355377018601</td><td>-0.340118016143495</td><td>-0.809013958505059</td><td>0.548521663785885</td><td>-0.0201828675962336</td></tr><tr><th>13</th><td>-0.795410435603455</td><td>-0.699989939762738</td><td>0.3991295030063</td><td>-0.510476261900736</td><td>1.51552245416844</td><td>1.28743032939467</td><td>1.53772393250903</td><td>0.133989160117702</td><td>1.02025736886037</td><td>0.499018733899186</td><td>-0.36948273277931</td><td>-1.10153460436318</td><td>-0.598132438886619</td></tr><tr><th>14</th><td>-0.193483122930324</td><td>-0.286021160323518</td><td>-0.691494225262995</td><td>0.0131581678700699</td><td>1.52337470686782</td><td>1.4010638072262</td><td>1.53114620451896</td><td>0.333066483478075</td><td>1.04372480381099</td><td>0.163206783570466</td><td>-0.422883765001728</td><td>-0.383527976713573</td><td>-0.489221907788158</td></tr><tr><th>15</th><td>0.151246203379718</td><td>2.09185108993614</td><td>2.03800472474384</td><td>-1.12474717143531</td><td>1.66557024390713</td><td>1.62535675109576</td><td>1.58751070483655</td><td>0.635852186043776</td><td>0.842577784605979</td><td>0.450761870778952</td><td>-1.39479033623028</td><td>-0.560984107567768</td><td>0.289349776549287</td></tr><tr><th>16</th><td>-0.464608740812712</td><td>0.36127694772303</td><td>1.2327673928287</td><td>-0.826033731086383</td><td>1.43475224709983</td><td>1.74451823818846</td><td>0.211096887484638</td><td>2.64816425140548</td><td>1.02511433146096</td><td>0.11975731603184</td><td>0.0596832073448267</td><td>-0.631231612661616</td><td>-0.207878671782927</td></tr><tr><th>17</th><td>-0.732977488012215</td><td>-0.526223425889779</td><td>0.61657871336593</td><td>-0.55447974332593</td><td>0.947484859025104</td><td>0.936833214138173</td><td>0.972516806335524</td><td>0.290251013865227</td><td>1.01285359725723</td><td>0.516207422283291</td><td>-0.0300689171988194</td><td>0.8787322524583</td><td>0.450254629309513</td></tr><tr><th>18</th><td>-0.167326459622119</td><td>0.175327165487237</td><td>0.287467725892572</td><td>-0.402652532084246</td><td>0.551181509418056</td><td>0.522204743290975</td><td>0.436837660094653</td><td>0.299564933845579</td><td>0.583109520896067</td><td>-0.704415820005353</td><td>-0.730810367994577</td><td>-1.95140580379896</td><td>-0.933504665700164</td></tr><tr><th>19</th><td>1.41159485787418</td><td>1.78722407901017</td><td>0.84397639585364</td><td>0.481278083772991</td><td>-0.0887673728508268</td><td>-0.49957757426858</td><td>0.304195897924847</td><td>-1.23884208383369</td><td>-0.153475724036624</td><td>-0.870486102788329</td><td>0.0955473331150403</td><td>-0.983708050882817</td><td>-0.3563445644514</td></tr><tr><th>20</th><td>-1.42997091652825</td><td>-0.490147045034213</td><td>0.272730237607695</td><td>-1.61029992954153</td><td>0.990787817197748</td><td>0.711687532608184</td><td>1.1885836012715</td><td>-0.371229188075638</td><td>1.24703459239952</td><td>-0.0389162332271516</td><td>0.883495749072872</td><td>2.58988026321017</td><td>3.33539552370368</td></tr><tr><th>21</th><td>-0.147247288176765</td><td>0.12328430415652</td><td>0.617549051912237</td><td>-0.18713077178262</td><td>0.256438107586694</td><td>0.17794983735083</td><td>0.412611806463263</td><td>-0.244809124559737</td><td>0.0947624806136492</td><td>0.723017223849532</td><td>-0.683948354633436</td><td>0.0873751276309269</td><td>-0.262209652750371</td></tr><tr><th>22</th><td>-0.187112676773894</td><td>-0.270777264595619</td><td>-1.01556818551606</td><td>0.0602850568600233</td><td>0.272419757757978</td><td>0.869133161879197</td><td>-0.657519461414234</td><td>2.32388522018189</td><td>-0.999936011525034</td><td>1.44671844178306</td><td>0.971157886040772</td><td>-0.358747904241515</td><td>-0.439657942096136</td></tr><tr><th>23</th><td>-1.82434047163768</td><td>-0.933480446068067</td><td>1.29474003766977</td><td>-1.94545221151036</td><td>0.33584651189654</td><td>0.359201654302844</td><td>0.513652924365886</td><td>-0.073197696696958</td><td>1.57139042812005</td><td>1.53329371326728</td><td>1.82076821859528</td><td>2.22740301867829</td><td>1.50063347195857</td></tr><tr><th>24</th><td>-2.29344084351335</td><td>-2.49161842344418</td><td>0.40383988742336</td><td>-2.36488074752948</td><td>1.4105254831956</td><td>1.42244117147792</td><td>1.17024166272172</td><td>0.84476650176855</td><td>1.79026875432495</td><td>0.648181858970515</td><td>-0.0857231057403538</td><td>-1.02789535292617</td><td>0.491288088952859</td></tr><tr><th>25</th><td>-0.434135932888305</td><td>0.740881989034652</td><td>0.699576357578518</td><td>-1.02405543187775</td><td>0.759529223983713</td><td>0.956656110895288</td><td>0.633299568656589</td><td>0.770733932268516</td><td>0.824988511714526</td><td>1.84287437634769</td><td>1.91045942063443</td><td>-0.502317207869366</td><td>0.132670133448219</td></tr><tr><th>26</th><td>-2.1920969546557</td><td>-2.49465664272271</td><td>0.354854763893431</td><td>-1.93155848635714</td><td>0.941979400289938</td><td>0.978917101414106</td><td>0.894860097289736</td><td>0.463239402831873</td><td>1.12537133317163</td><td>1.70528446191955</td><td>0.717792714479123</td><td>0.645888049108261</td><td>0.783968250169388</td></tr><tr><th>27</th><td>-1.46602269088422</td><td>-1.24921677101897</td><td>0.307977693653039</td><td>-1.55097364660989</td><td>0.618908494474798</td><td>0.662508171662042</td><td>0.475957173906078</td><td>0.484718674597707</td><td>0.401564892028249</td><td>0.55987973254026</td><td>-0.376938143754217</td><td>-0.933982629228218</td><td>0.390013151672955</td></tr><tr><th>28</th><td>-1.83317744236881</td><td>-1.53268787828701</td><td>2.55674262685865</td><td>-1.51827745783835</td><td>0.789409601746455</td><td>0.908747799728588</td><td>0.649971922941479</td><td>0.668373649931667</td><td>1.20058303519903</td><td>0.277963256075637</td><td>1.2504953198275</td><td>3.31370445071638</td><td>2.22035828885342</td></tr><tr><th>29</th><td>-0.784546628243178</td><td>0.276582579543931</td><td>3.01104958800057</td><td>-1.11978843206758</td><td>0.920823858422707</td><td>0.750217689886151</td><td>1.26153730009639</td><td>-0.403363882922417</td><td>0.400667296857811</td><td>-0.217597941303479</td><td>-0.724669537565068</td><td>-0.391945338467193</td><td>-0.650023936358253</td></tr><tr><th>30</th><td>0.464455916345135</td><td>1.3326356122229</td><td>-1.23059563374303</td><td>-0.357975958937414</td><td>1.18249746977104</td><td>1.54315938069757</td><td>-0.60339041154062</td><td>3.38308845958422</td><td>0.823740765148641</td><td>-0.129951318508883</td><td>-0.657979878422938</td><td>-0.499534924074273</td><td>-0.414476569095651</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table>




```julia
describe(cg10k_trait[:, 3:end])
```

    Trait1
    Summary Stats:
    Mean:           0.002211
    Minimum:        -3.204128
    1st Quartile:   -0.645771
    Median:         0.125010
    3rd Quartile:   0.723315
    Maximum:        3.479398
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait2
    Summary Stats:
    Mean:           0.001353
    Minimum:        -3.511659
    1st Quartile:   -0.642621
    Median:         0.033517
    3rd Quartile:   0.657467
    Maximum:        4.913423
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait3
    Summary Stats:
    Mean:           -0.001296
    Minimum:        -3.938436
    1st Quartile:   -0.640907
    Median:         -0.000782
    3rd Quartile:   0.637108
    Maximum:        7.916299
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait4
    Summary Stats:
    Mean:           0.002309
    Minimum:        -3.608403
    1st Quartile:   -0.546086
    Median:         0.228165
    3rd Quartile:   0.715291
    Maximum:        3.127688
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait5
    Summary Stats:
    Mean:           -0.001790
    Minimum:        -4.148749
    1st Quartile:   -0.690765
    Median:         0.031034
    3rd Quartile:   0.734916
    Maximum:        2.717184
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait6
    Summary Stats:
    Mean:           -0.001196
    Minimum:        -3.824792
    1st Quartile:   -0.662796
    Median:         0.036242
    3rd Quartile:   0.741176
    Maximum:        2.589728
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait7
    Summary Stats:
    Mean:           -0.001989
    Minimum:        -4.272455
    1st Quartile:   -0.638923
    Median:         0.069801
    3rd Quartile:   0.710423
    Maximum:        2.653779
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait8
    Summary Stats:
    Mean:           0.000614
    Minimum:        -5.625488
    1st Quartile:   -0.601575
    Median:         -0.038630
    3rd Quartile:   0.527342
    Maximum:        5.805702
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait9
    Summary Stats:
    Mean:           -0.001810
    Minimum:        -5.381968
    1st Quartile:   -0.601429
    Median:         0.106571
    3rd Quartile:   0.698567
    Maximum:        2.571936
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait10
    Summary Stats:
    Mean:           -0.000437
    Minimum:        -3.548506
    1st Quartile:   -0.633641
    Median:         -0.096651
    3rd Quartile:   0.498610
    Maximum:        6.537820
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait11
    Summary Stats:
    Mean:           -0.000616
    Minimum:        -3.264910
    1st Quartile:   -0.673685
    Median:         -0.068044
    3rd Quartile:   0.655486
    Maximum:        4.262410
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait12
    Summary Stats:
    Mean:           -0.000589
    Minimum:        -8.851909
    1st Quartile:   -0.539686
    Median:         -0.141099
    3rd Quartile:   0.350779
    Maximum:        13.211402
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    
    Trait13
    Summary Stats:
    Mean:           -0.000151
    Minimum:        -5.592104
    1st Quartile:   -0.492289
    Median:         -0.141022
    3rd Quartile:   0.324804
    Maximum:        24.174436
    Length:         6670
    Type:           Float64
    Number Missing: 0
    % Missing:      0.000000
    



```julia
Y = convert(Matrix{Float64}, cg10k_trait[:, 3:15])
histogram(Y, layout = 13)
```




<img src="data:image/png;base64,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" />



## Pre-processing data for heritability analysis

To prepare variance component model fitting, we form an instance of `VarianceComponentVariate`. The two variance components are $(2\Phi, I)$.


```julia
using VarianceComponentModels

# form data as VarianceComponentVariate
cg10kdata = VarianceComponentVariate(Y, (2Φgrm, eye(size(Y, 1))))
fieldnames(cg10kdata)
```




    3-element Array{Symbol,1}:
     :Y
     :X
     :V




```julia
cg10kdata
```




    VarianceComponentModels.VarianceComponentVariate{Float64,2,Array{Float64,2},Array{Float64,2},Array{Float64,2}}([-1.81573 -0.94615 … -1.02853 -0.394049; -1.2444 0.10966 … 1.09065 0.0256616; … ; 0.886626 0.487408 … -0.636874 -0.439825; -1.24394 0.213697 … 0.299931 0.392809], Array{Float64}(6670,0), ([1.00605 0.00671009 … -0.000109037 -0.00556144; 0.00671009 0.997916 … 0.00173694 0.00685701; … ; -0.000109037 0.00173694 … 1.00166 0.000938955; -0.00556144 0.00685701 … 0.000938955 1.00125], [1.0 0.0 … 0.0 0.0; 0.0 1.0 … 0.0 0.0; … ; 0.0 0.0 … 1.0 0.0; 0.0 0.0 … 0.0 1.0]))



Before fitting the variance component model, we pre-compute the eigen-decomposition of $2\Phi_{\text{GRM}}$, the rotated responses, and the constant part in log-likelihood, and store them as a `TwoVarCompVariateRotate` instance, which is re-used in various variane component estimation procedures.


```julia
# pre-compute eigen-decomposition (~50 secs on my laptop)
@time cg10kdata_rotated = TwoVarCompVariateRotate(cg10kdata)
fieldnames(cg10kdata_rotated)
```

     48.837361 seconds (39 allocations: 1021.427 MiB, 0.57% gc time)





    5-element Array{Symbol,1}:
     :Yrot    
     :Xrot    
     :eigval  
     :eigvec  
     :logdetV2



## Save intermediate results

We don't want to re-compute SnpArray and empirical kinship matrices again and again for heritibility analysis.


```julia
# # Pkg.add("JLD")
# using JLD
# @save "cg10k.jld"
# whos()
```

To load workspace


```julia
using SnpArrays, JLD, DataFrames, VarianceComponentModels, Plots
pyplot()
@load "cg10k.jld"
whos()
```

                              Base               Module
                           BinDeps  41348 KB     Module
                             Blosc  41202 KB     Module
                        ColorTypes  41457 KB     Module
                            Colors  41480 KB     Module
                            Compat  41196 KB     Module
                             Conda  41205 KB     Module
                              Core               Module
                        DataArrays  41456 KB     Module
                        DataFrames  41684 KB     Module
                    DataStructures  41356 KB     Module
                            FileIO  41310 KB     Module
                 FixedPointNumbers  41695 KB     Module
                              GZip  41181 KB     Module
                              HDF5  41403 KB     Module
                            IJulia 4185781 KB     Module
                             Ipopt  41172 KB     Module
                               JLD  41376 KB     Module
                              JSON  41245 KB     Module
                      LaTeXStrings   4058 bytes  Module
                     LegacyStrings  41212 KB     Module
                        LinearMaps     22 KB     Module
                        MacroTools  41606 KB     Module
                              Main               Module
                      MathProgBase  41353 KB     Module
                           MbedTLS  41269 KB     Module
                          Measures  41175 KB     Module
                           NaNMath  41200 KB     Module
                        PlotThemes  41167 KB     Module
                         PlotUtils  41332 KB     Module
                             Plots  42960 KB     Module
                            PyCall  41711 KB     Module
                            PyPlot  41771 KB     Module
                       RecipesBase  41283 KB     Module
                          Reexport  41160 KB     Module
                          Requires  41172 KB     Module
                               SHA     62 KB     Module
                           Showoff  41163 KB     Module
                         SnpArrays  41218 KB     Module
                 SortingAlgorithms  41178 KB     Module
                  SpecialFunctions  41252 KB     Module
                      StaticArrays  41744 KB     Module
                         StatsBase  41810 KB     Module
                         URIParser  41171 KB     Module
           VarianceComponentModels  41278 KB     Module
                                 Y    677 KB     6670×13 Array{Float64,2}
                               ZMQ  41223 KB     Module
                                 _     77 KB     630860-element BitArray{1}
                             cg10k 1027303 KB     6670×630860 SnpArrays.SnpArray{2}
                       cg10k_trait    978 KB     6670×15 DataFrames.DataFrame
                         cg10kdata 695816 KB     VarianceComponentModels.VarianceCo…
                 cg10kdata_rotated 348299 KB     VarianceComponentModels.TwoVarComp…
                                 h     24 bytes  3-element Array{Float64,1}
                               hST    104 bytes  13-element Array{Float64,1}
                            hST_se    104 bytes  13-element Array{Float64,1}
                               hse     24 bytes  3-element Array{Float64,1}
                               maf   4928 KB     630860-element Array{Float64,1}
                   missings_by_snp   4928 KB     630860-element Array{Int64,1}
                            people      8 bytes  Int64
                              snps      8 bytes  Int64
                      trait57_data 347778 KB     VarianceComponentModels.TwoVarComp…
                     trait57_model    232 bytes  VarianceComponentModels.VarianceCo…
                    traitall_model   2792 bytes  VarianceComponentModels.VarianceCo…
                          traitidx     16 bytes  3-element UnitRange{Int64}
                                Σa   3848 bytes  13×13 Array{Array{Float64,2},2}
                              Σcov   2592 bytes  18×18 Array{Float64,2}
                                Σe   3848 bytes  13×13 Array{Array{Float64,2},2}
                              Φgrm 347569 KB     6670×6670 Array{Float64,2}
                               σ2a    104 bytes  13-element Array{Float64,1}
                               σ2e    104 bytes  13-element Array{Float64,1}


## Heritability of single traits

We use Fisher scoring algorithm to fit variance component model for each single trait.


```julia
# heritability from single trait analysis
hST = zeros(13)
# standard errors of estimated heritability
hST_se = zeros(13)
# additive genetic effects
σ2a = zeros(13)
# enviromental effects
σ2e = zeros(13)

tic()
for trait in 1:13
    println(names(cg10k_trait)[trait + 2])
    # form data set for trait j
    traitj_data = TwoVarCompVariateRotate(cg10kdata_rotated.Yrot[:, trait], cg10kdata_rotated.Xrot, 
        cg10kdata_rotated.eigval, cg10kdata_rotated.eigvec, cg10kdata_rotated.logdetV2)
    # initialize model parameters
    traitj_model = VarianceComponentModel(traitj_data)
    # estimate variance components
    _, _, _, Σcov, _, _ = mle_fs!(traitj_model, traitj_data; solver=:Ipopt, verbose=false)
    σ2a[trait] = traitj_model.Σ[1][1]
    σ2e[trait] = traitj_model.Σ[2][1]
    @show σ2a[trait], σ2e[trait]
    h, hse = heritability(traitj_model.Σ, Σcov)
    hST[trait] = h[1]
    hST_se[trait] = hse[1]
end
toc()
```

    Trait1
    (σ2a[trait], σ2e[trait]) = (0.25978160614793233, 0.7369535197912689)
    Trait2
    (σ2a[trait], σ2e[trait]) = (0.18647130348299173, 0.8129591079735827)
    Trait3
    (σ2a[trait], σ2e[trait]) = (0.3188368159422607, 0.6798809726936244)
    Trait4
    (σ2a[trait], σ2e[trait]) = (0.2651357653143703, 0.7308007669086968)
    Trait5
    (σ2a[trait], σ2e[trait]) = (0.28083388108246, 0.7172036435586534)
    Trait6
    (σ2a[trait], σ2e[trait]) = (0.2824159905728832, 0.7170988773569172)
    Trait7
    (σ2a[trait], σ2e[trait]) = (0.2155274336968625, 0.7815346282986375)
    Trait8
    (σ2a[trait], σ2e[trait]) = (0.194687807263945, 0.8049690651320599)
    Trait9
    (σ2a[trait], σ2e[trait]) = (0.24706855916591713, 0.7512942998567308)
    Trait10
    (σ2a[trait], σ2e[trait]) = (0.098712236297271, 0.9011756660217387)
    Trait11
    (σ2a[trait], σ2e[trait]) = (0.1664264642608195, 0.8322427413046204)
    Trait12
    (σ2a[trait], σ2e[trait]) = (0.0834296761650666, 0.9153609794266608)
    Trait13
    (σ2a[trait], σ2e[trait]) = (0.05893968504298988, 0.940270012443928)
    elapsed time: 0.160999612 seconds





    0.160999612




```julia
# heritability and standard errors
[hST'; hST_se']
```




    2×13 Array{Float64,2}:
     0.260633   0.186578   0.319246   …  0.166648  0.0835307  0.0589863
     0.0799732  0.0869002  0.0741007     0.08862   0.0944407  0.0953238



## Pairwise traits

Joint analysis of multiple traits is subject to intensive research recently. Following code snippet does joint analysis of all pairs of traits, a total of 78 bivariate variane component models.


```julia
# additive genetic effects (2x2 psd matrices) from bavariate trait analysis;
Σa = Array{Matrix{Float64}}(13, 13)
# environmental effects (2x2 psd matrices) from bavariate trait analysis;
Σe = Array{Matrix{Float64}}(13, 13)

tic()
for i in 1:13
    for j in (i+1):13
        println(names(cg10k_trait)[i + 2], names(cg10k_trait)[j + 2])
        # form data set for (trait1, trait2)
        traitij_data = TwoVarCompVariateRotate(cg10kdata_rotated.Yrot[:, [i;j]], cg10kdata_rotated.Xrot, 
            cg10kdata_rotated.eigval, cg10kdata_rotated.eigvec, cg10kdata_rotated.logdetV2)
        # initialize model parameters
        traitij_model = VarianceComponentModel(traitij_data)
        # estimate variance components
        mle_fs!(traitij_model, traitij_data; solver=:Ipopt, verbose=false)
        Σa[i, j] = traitij_model.Σ[1]
        Σe[i, j] = traitij_model.Σ[2]
        @show Σa[i, j], Σe[i, j]
    end
end
toc()
```

    Trait1Trait2
    (Σa[i, j], Σe[i, j]) = ([0.258822 0.174358; 0.174358 0.185108], [0.737892 0.585751; 0.585751 0.814301])
    Trait1Trait3
    (Σa[i, j], Σe[i, j]) = ([0.260236 -0.0144726; -0.0144726 0.319245], [0.736512 -0.11979; -0.11979 0.679488])
    Trait1Trait4
    (Σa[i, j], Σe[i, j]) = ([0.259615 0.222203; 0.222203 0.265149], [0.737116 0.599854; 0.599854 0.730788])
    Trait1Trait5
    (Σa[i, j], Σe[i, j]) = ([0.259574 -0.146827; -0.146827 0.28153], [0.737153 -0.254777; -0.254777 0.71653])
    Trait1Trait6
    (Σa[i, j], Σe[i, j]) = ([0.259476 -0.129115; -0.129115 0.282688], [0.73725 -0.23161; -0.23161 0.716837])
    Trait1Trait7
    (Σa[i, j], Σe[i, j]) = ([0.259115 -0.140455; -0.140455 0.215297], [0.737606 -0.197616; -0.197616 0.781774])
    Trait1Trait8
    (Σa[i, j], Σe[i, j]) = ([0.259778 -0.0327756; -0.0327756 0.194698], [0.736957 -0.127026; -0.127026 0.804959])
    Trait1Trait9
    (Σa[i, j], Σe[i, j]) = ([0.261858 -0.204589; -0.204589 0.246027], [0.734961 -0.307734; -0.307734 0.75232])
    Trait1Trait10
    (Σa[i, j], Σe[i, j]) = ([0.259649 -0.0994858; -0.0994858 0.0956585], [0.737083 -0.303942; -0.303942 0.904218])
    Trait1Trait11
    (Σa[i, j], Σe[i, j]) = ([0.25947 -0.138603; -0.138603 0.164709], [0.737257 -0.359557; -0.359557 0.83395])
    Trait1Trait12
    (Σa[i, j], Σe[i, j]) = ([0.261779 -0.145414; -0.145414 0.0807748], [0.735076 -0.041823; -0.041823 0.9181])
    Trait1Trait13
    (Σa[i, j], Σe[i, j]) = ([0.261125 -0.108774; -0.108774 0.0538214], [0.735674 -0.114123; -0.114123 0.945416])
    Trait2Trait3
    (Σa[i, j], Σe[i, j]) = ([0.186541 0.144056; 0.144056 0.320627], [0.812888 0.0995944; 0.0995944 0.678167])
    Trait2Trait4
    (Σa[i, j], Σe[i, j]) = ([0.186131 0.0746032; 0.0746032 0.265122], [0.813293 0.221109; 0.221109 0.730814])
    Trait2Trait5
    (Σa[i, j], Σe[i, j]) = ([0.186442 -0.0118093; -0.0118093 0.280842], [0.812987 -0.0365191; -0.0365191 0.717195])
    Trait2Trait6
    (Σa[i, j], Σe[i, j]) = ([0.18649 -0.00366533; -0.00366533 0.282471], [0.812941 -0.0206271; -0.0206271 0.717046])
    Trait2Trait7
    (Σa[i, j], Σe[i, j]) = ([0.186104 -0.030665; -0.030665 0.215304], [0.81332 -0.000667009; -0.000667009 0.781755])
    Trait2Trait8
    (Σa[i, j], Σe[i, j]) = ([0.187023 0.0331783; 0.0331783 0.195259], [0.812421 -0.0326343; -0.0326343 0.804415])
    Trait2Trait9
    (Σa[i, j], Σe[i, j]) = ([0.185032 -0.085334; -0.085334 0.245909], [0.814386 -0.0809638; -0.0809638 0.752433])
    Trait2Trait10
    (Σa[i, j], Σe[i, j]) = ([0.186587 -0.123303; -0.123303 0.0987387], [0.812872 -0.273083; -0.273083 0.901229])
    Trait2Trait11
    (Σa[i, j], Σe[i, j]) = ([0.185484 -0.117256; -0.117256 0.167776], [0.81393 -0.296772; -0.296772 0.830934])
    Trait2Trait12
    (Σa[i, j], Σe[i, j]) = ([0.185907 -0.0909104; -0.0909104 0.0827171], [0.813555 0.0457924; 0.0457924 0.916135])
    Trait2Trait13
    (Σa[i, j], Σe[i, j]) = ([0.185979 -0.0720811; -0.0720811 0.0568238], [0.8135 0.0751703; 0.0751703 0.942424])
    Trait3Trait4
    (Σa[i, j], Σe[i, j]) = ([0.3188 -0.154562; -0.154562 0.264323], [0.679917 -0.303223; -0.303223 0.731591])
    Trait3Trait5
    (Σa[i, j], Σe[i, j]) = ([0.319216 0.183527; 0.183527 0.282063], [0.679514 0.33724; 0.33724 0.716008])
    Trait3Trait6
    (Σa[i, j], Σe[i, j]) = ([0.319776 0.165672; 0.165672 0.284448], [0.678972 0.298667; 0.298667 0.715124])
    Trait3Trait7
    (Σa[i, j], Σe[i, j]) = ([0.318838 0.166283; 0.166283 0.215261], [0.67988 0.347706; 0.347706 0.781796])
    Trait3Trait8
    (Σa[i, j], Σe[i, j]) = ([0.320718 0.0566397; 0.0566397 0.197764], [0.678063 0.0451569; 0.0451569 0.801956])
    Trait3Trait9
    (Σa[i, j], Σe[i, j]) = ([0.319001 0.137699; 0.137699 0.246142], [0.679722 0.266704; 0.266704 0.752197])
    Trait3Trait10
    (Σa[i, j], Σe[i, j]) = ([0.31908 -0.076513; -0.076513 0.0996001], [0.679646 -0.142905; -0.142905 0.900298])
    Trait3Trait11
    (Σa[i, j], Σe[i, j]) = ([0.318094 -0.0177494; -0.0177494 0.16629], [0.6806 -0.1144; -0.1144 0.832376])
    Trait3Trait12
    (Σa[i, j], Σe[i, j]) = ([0.321164 0.0843842; 0.0843842 0.0874609], [0.677639 0.0341558; 0.0341558 0.911368])
    Trait3Trait13
    (Σa[i, j], Σe[i, j]) = ([0.323273 0.109443; 0.109443 0.0634295], [0.675635 -0.0060525; -0.0060525 0.935819])
    Trait4Trait5
    (Σa[i, j], Σe[i, j]) = ([0.26525 -0.215125; -0.215125 0.282572], [0.73068 -0.377406; -0.377406 0.715518])
    Trait4Trait6
    (Σa[i, j], Σe[i, j]) = ([0.265715 -0.199714; -0.199714 0.283942], [0.730231 -0.347732; -0.347732 0.715619])
    Trait4Trait7
    (Σa[i, j], Σe[i, j]) = ([0.26407 -0.18238; -0.18238 0.214324], [0.731843 -0.32655; -0.32655 0.782733])
    Trait4Trait8
    (Σa[i, j], Σe[i, j]) = ([0.266229 -0.0965381; -0.0965381 0.196655], [0.729739 -0.151461; -0.151461 0.803044])
    Trait4Trait9
    (Σa[i, j], Σe[i, j]) = ([0.269627 -0.226931; -0.226931 0.247265], [0.726443 -0.416085; -0.416085 0.751086])
    Trait4Trait10
    (Σa[i, j], Σe[i, j]) = ([0.265098 -0.0352926; -0.0352926 0.0981462], [0.730847 -0.226248; -0.226248 0.901736])
    Trait4Trait11
    (Σa[i, j], Σe[i, j]) = ([0.265178 -0.0970634; -0.0970634 0.164885], [0.73076 -0.272291; -0.272291 0.833762])
    Trait4Trait12
    (Σa[i, j], Σe[i, j]) = ([0.267732 -0.140985; -0.140985 0.081029], [0.728323 -0.0834791; -0.0834791 0.917815])
    Trait4Trait13
    (Σa[i, j], Σe[i, j]) = ([0.265695 -0.0970238; -0.0970238 0.0564809], [0.730259 -0.226115; -0.226115 0.942736])
    Trait5Trait6
    (Σa[i, j], Σe[i, j]) = ([0.281198 0.280259; 0.280259 0.281764], [0.716855 0.661013; 0.661013 0.717735])
    Trait5Trait7
    (Σa[i, j], Σe[i, j]) = ([0.280442 0.231918; 0.231918 0.211837], [0.717597 0.674491; 0.674491 0.785172])
    Trait5Trait8
    (Σa[i, j], Σe[i, j]) = ([0.280958 0.163168; 0.163168 0.193315], [0.717089 0.221817; 0.221817 0.806314])
    Trait5Trait9
    (Σa[i, j], Σe[i, j]) = ([0.283544 0.243884; 0.243884 0.240564], [0.714585 0.509072; 0.509072 0.757631])
    Trait5Trait10
    (Σa[i, j], Σe[i, j]) = ([0.281378 -0.0454427; -0.0454427 0.100081], [0.716678 -0.0579778; -0.0579778 0.899822])
    Trait5Trait11
    (Σa[i, j], Σe[i, j]) = ([0.280066 0.0195669; 0.0195669 0.165607], [0.71795 -0.0345589; -0.0345589 0.833047])
    Trait5Trait12
    (Σa[i, j], Σe[i, j]) = ([0.28101 0.0592641; 0.0592641 0.0831831], [0.717036 0.0552788; 0.0552788 0.915608])
    Trait5Trait13
    (Σa[i, j], Σe[i, j]) = ([0.281854 0.0680641; 0.0680641 0.0591899], [0.716223 0.0551992; 0.0551992 0.940027])
    Trait6Trait7
    (Σa[i, j], Σe[i, j]) = ([0.282435 0.220236; 0.220236 0.213997], [0.71708 0.581507; 0.581507 0.783041])
    Trait6Trait8
    (Σa[i, j], Σe[i, j]) = ([0.282435 0.18375; 0.18375 0.192999], [0.717081 0.436932; 0.436932 0.80663])
    Trait6Trait9
    (Σa[i, j], Σe[i, j]) = ([0.284516 0.233768; 0.233768 0.242478], [0.715071 0.477502; 0.477502 0.755765])
    Trait6Trait10
    (Σa[i, j], Σe[i, j]) = ([0.283087 -0.0427658; -0.0427658 0.100634], [0.716449 -0.0599491; -0.0599491 0.899275])
    Trait6Trait11
    (Σa[i, j], Σe[i, j]) = ([0.281046 0.0272144; 0.0272144 0.165044], [0.71843 -0.0516242; -0.0516242 0.833601])
    Trait6Trait12
    (Σa[i, j], Σe[i, j]) = ([0.28256 0.0548537; 0.0548537 0.083133], [0.716961 0.0502064; 0.0502064 0.915658])
    Trait6Trait13
    (Σa[i, j], Σe[i, j]) = ([0.283231 0.0585667; 0.0585667 0.0592752], [0.716314 0.055827; 0.055827 0.939942])
    Trait7Trait8
    (Σa[i, j], Σe[i, j]) = ([0.213998 0.0875641; 0.0875641 0.192993], [0.78304 -0.055939; -0.055939 0.806635])
    Trait7Trait9
    (Σa[i, j], Σe[i, j]) = ([0.219039 0.216925; 0.216925 0.243338], [0.778156 0.463024; 0.463024 0.754935])
    Trait7Trait10
    (Σa[i, j], Σe[i, j]) = ([0.216296 -0.0412106; -0.0412106 0.100663], [0.780785 -0.0868086; -0.0868086 0.899246])
    Trait7Trait11
    (Σa[i, j], Σe[i, j]) = ([0.2142 0.0204227; 0.0204227 0.165077], [0.782833 -0.0478727; -0.0478727 0.833569])
    Trait7Trait12
    (Σa[i, j], Σe[i, j]) = ([0.215054 0.0738562; 0.0738562 0.0814228], [0.782012 0.0366272; 0.0366272 0.917365])
    Trait7Trait13
    (Σa[i, j], Σe[i, j]) = ([0.216093 0.0728515; 0.0728515 0.0570272], [0.781006 0.0409945; 0.0409945 0.942189])
    Trait8Trait9
    (Σa[i, j], Σe[i, j]) = ([0.195154 0.111756; 0.111756 0.246453], [0.804528 0.185842; 0.185842 0.751896])
    Trait8Trait10
    (Σa[i, j], Σe[i, j]) = ([0.195015 -0.015506; -0.015506 0.0990776], [0.804651 0.0118538; 0.0118538 0.900815])
    Trait8Trait11
    (Σa[i, j], Σe[i, j]) = ([0.194421 0.0215044; 0.0215044 0.166226], [0.805231 -0.026247; -0.026247 0.83244])
    Trait8Trait12
    (Σa[i, j], Σe[i, j]) = ([0.194491 -0.00425152; -0.00425152 0.0832711], [0.805162 0.0349872; 0.0349872 0.915518])
    Trait8Trait13
    (Σa[i, j], Σe[i, j]) = ([0.19448 0.00235501; 0.00235501 0.0589351], [0.805173 0.0396048; 0.0396048 0.940275])
    Trait9Trait10
    (Σa[i, j], Σe[i, j]) = ([0.246455 -0.00257997; -0.00257997 0.0984563], [0.751895 0.0743439; 0.0743439 0.901429])
    Trait9Trait11
    (Σa[i, j], Σe[i, j]) = ([0.247001 0.0303415; 0.0303415 0.166421], [0.75136 0.153765; 0.153765 0.832248])
    Trait9Trait12
    (Σa[i, j], Σe[i, j]) = ([0.249421 0.0829968; 0.0829968 0.0890874], [0.749007 0.109331; 0.109331 0.909778])
    Trait9Trait13
    (Σa[i, j], Σe[i, j]) = ([0.24861 0.0916799; 0.0916799 0.0602352], [0.749811 0.100027; 0.100027 0.939032])
    Trait10Trait11
    (Σa[i, j], Σe[i, j]) = ([0.0914658 0.100613; 0.100613 0.166501], [0.908376 0.473847; 0.473847 0.83217])
    Trait10Trait12
    (Σa[i, j], Σe[i, j]) = ([0.0951392 0.0588424; 0.0588424 0.0796744], [0.904735 0.0828862; 0.0828862 0.919115])
    Trait10Trait13
    (Σa[i, j], Σe[i, j]) = ([0.0995192 -0.0257171; -0.0257171 0.0598595], [0.900397 0.163778; 0.163778 0.939368])
    Trait11Trait12
    (Σa[i, j], Σe[i, j]) = ([0.165386 0.0579914; 0.0579914 0.0796005], [0.83327 0.144637; 0.144637 0.919166])
    Trait11Trait13
    (Σa[i, j], Σe[i, j]) = ([0.166417 -0.000985185; -0.000985185 0.0595681], [0.832265 0.200012; 0.200012 0.939646])
    Trait12Trait13
    (Σa[i, j], Σe[i, j]) = ([0.085082 0.0696185; 0.0696185 0.0569655], [0.913729 0.572041; 0.572041 0.942247])
    elapsed time: 3.587337102 seconds





    3.587337102



## 3-trait analysis

Researchers want to jointly analyze traits 5-7. Our strategy is to try both Fisher scoring and MM algorithm with different starting point, and choose the best local optimum. We first form the data set and run Fisher scoring, which yields a final objective value -1.4700991+04.


```julia
traitidx = 5:7
# form data set
trait57_data = TwoVarCompVariateRotate(cg10kdata_rotated.Yrot[:, traitidx], cg10kdata_rotated.Xrot, 
    cg10kdata_rotated.eigval, cg10kdata_rotated.eigvec, cg10kdata_rotated.logdetV2)
# initialize model parameters
trait57_model = VarianceComponentModel(trait57_data)
# estimate variance components
@time mle_fs!(trait57_model, trait57_data; solver=:Ipopt, verbose=true)
trait57_model
```

    This is Ipopt version 3.12.4, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:        0
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       78
    
    Total number of variables............................:       12
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:        0
                         variables with only upper bounds:        0
    Total number of equality constraints.................:        0
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  3.0247565e+04 0.00e+00 1.00e+02   0.0 0.00e+00    -  0.00e+00 0.00e+00   0 
       5  1.6835078e+04 0.00e+00 4.08e+02 -11.0 3.64e-01    -  1.00e+00 1.00e+00f  1 MaxS
      10  1.4742941e+04 0.00e+00 1.10e+02 -11.0 2.35e-01    -  1.00e+00 1.00e+00f  1 MaxS
      15  1.4701394e+04 0.00e+00 1.16e+01 -11.0 7.78e-02  -4.5 1.00e+00 1.00e+00f  1 MaxS
      20  1.4701019e+04 0.00e+00 5.75e-01 -11.0 1.51e-04  -6.9 1.00e+00 1.00e+00f  1 MaxS
      25  1.4701018e+04 0.00e+00 2.40e-02 -11.0 6.38e-06  -9.2 1.00e+00 1.00e+00f  1 MaxS
      30  1.4701018e+04 0.00e+00 9.98e-04 -11.0 2.66e-07 -11.6 1.00e+00 1.00e+00f  1 MaxS
      35  1.4701018e+04 0.00e+00 4.15e-05 -11.0 1.10e-08 -14.0 1.00e+00 1.00e+00h  1 MaxS
      40  1.4701018e+04 0.00e+00 1.72e-06 -11.0 4.59e-10 -16.4 1.00e+00 1.00e+00f  1 MaxSA
      45  1.4701018e+04 0.00e+00 7.17e-08 -11.0 1.91e-11 -18.8 1.00e+00 1.00e+00h  1 MaxSA
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    
    Number of Iterations....: 49
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4720359684330265e+02    1.4701017692082147e+04
    Dual infeasibility......:   5.6081357364421780e-09    1.8435742302386474e-07
    Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   5.6081357364421780e-09    1.8435742302386474e-07
    
    
    Number of objective function evaluations             = 50
    Number of objective gradient evaluations             = 50
    Number of equality constraint evaluations            = 0
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 0
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 49
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.014
    Total CPU secs in NLP function evaluations           =      0.076
    
    EXIT: Optimal Solution Found.
      0.097955 seconds (55.15 k allocations: 5.632 MiB)





    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(Array{Float64}(0,3), ([0.280777 0.279441 0.232208; 0.279441 0.28422 0.219831; 0.232208 0.219831 0.212832], [0.717266 0.66183 0.674206; 0.66183 0.715287 0.581891; 0.674206 0.581891 0.784183]), Array{Float64}(0,0), Char[], Float64[], -Inf, Inf)



We then run the MM algorithm, starting from the Fisher scoring answer. MM finds an improved solution with objective value 8.955397e+03.


```julia
# trait59_model contains the fitted model by Fisher scoring now
@time mle_mm!(trait57_model, trait57_data; verbose=true)
trait57_model
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -1.470102e+04
           1  -1.470102e+04
    
      0.003006 seconds (21.01 k allocations: 1.551 MiB)





    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(Array{Float64}(0,3), ([0.280777 0.279441 0.232208; 0.279441 0.28422 0.219831; 0.232208 0.219831 0.212832], [0.717266 0.66183 0.674206; 0.66183 0.715287 0.581891; 0.674206 0.581891 0.784183]), Array{Float64}(0,0), Char[], Float64[], -Inf, Inf)



Do another run of MM algorithm from default starting point. It leads to a slightly better local optimum -1.470104e+04, slighly worse than the Fisher scoring result. Follow up anlaysis should use the Fisher scoring result.


```julia
# default starting point
trait57_model = VarianceComponentModel(trait57_data)
@time _, _, _, Σcov, = mle_mm!(trait57_model, trait57_data; verbose=true)
trait57_model
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -3.024757e+04
           1  -2.040300e+04
           2  -1.656070e+04
           3  -1.528529e+04
           4  -1.490986e+04
           5  -1.480638e+04
           6  -1.477811e+04
           7  -1.476968e+04
           8  -1.476639e+04
           9  -1.476444e+04
          10  -1.476286e+04
          20  -1.475000e+04
          30  -1.474011e+04
          40  -1.473248e+04
          50  -1.472658e+04
          60  -1.472200e+04
          70  -1.471840e+04
          80  -1.471555e+04
          90  -1.471328e+04
         100  -1.471145e+04
         110  -1.470997e+04
         120  -1.470875e+04
         130  -1.470775e+04
         140  -1.470691e+04
         150  -1.470621e+04
         160  -1.470562e+04
         170  -1.470511e+04
         180  -1.470469e+04
         190  -1.470432e+04
         200  -1.470400e+04
         210  -1.470372e+04
         220  -1.470348e+04
         230  -1.470326e+04
         240  -1.470308e+04
         250  -1.470291e+04
         260  -1.470276e+04
         270  -1.470263e+04
         280  -1.470251e+04
         290  -1.470241e+04
         300  -1.470231e+04
         310  -1.470223e+04
         320  -1.470215e+04
         330  -1.470208e+04
         340  -1.470201e+04
         350  -1.470195e+04
         360  -1.470190e+04
         370  -1.470185e+04
         380  -1.470180e+04
         390  -1.470176e+04
         400  -1.470172e+04
         410  -1.470168e+04
         420  -1.470165e+04
         430  -1.470162e+04
         440  -1.470159e+04
         450  -1.470156e+04
         460  -1.470153e+04
         470  -1.470151e+04
         480  -1.470149e+04
         490  -1.470147e+04
         500  -1.470145e+04
         510  -1.470143e+04
         520  -1.470141e+04
         530  -1.470139e+04
         540  -1.470138e+04
         550  -1.470136e+04
         560  -1.470135e+04
         570  -1.470133e+04
         580  -1.470132e+04
         590  -1.470131e+04
         600  -1.470130e+04
         610  -1.470129e+04
         620  -1.470128e+04
         630  -1.470127e+04
         640  -1.470126e+04
         650  -1.470125e+04
         660  -1.470124e+04
         670  -1.470123e+04
         680  -1.470122e+04
         690  -1.470122e+04
         700  -1.470121e+04
         710  -1.470120e+04
         720  -1.470120e+04
         730  -1.470119e+04
         740  -1.470118e+04
         750  -1.470118e+04
         760  -1.470117e+04
         770  -1.470117e+04
         780  -1.470116e+04
         790  -1.470116e+04
         800  -1.470115e+04
         810  -1.470115e+04
         820  -1.470114e+04
         830  -1.470114e+04
         840  -1.470114e+04
         850  -1.470113e+04
         860  -1.470113e+04
         870  -1.470112e+04
         880  -1.470112e+04
         890  -1.470112e+04
         900  -1.470111e+04
         910  -1.470111e+04
         920  -1.470111e+04
         930  -1.470111e+04
         940  -1.470110e+04
         950  -1.470110e+04
         960  -1.470110e+04
         970  -1.470109e+04
         980  -1.470109e+04
         990  -1.470109e+04
        1000  -1.470109e+04
        1010  -1.470109e+04
        1020  -1.470108e+04
        1030  -1.470108e+04
        1040  -1.470108e+04
        1050  -1.470108e+04
        1060  -1.470108e+04
        1070  -1.470107e+04
        1080  -1.470107e+04
        1090  -1.470107e+04
        1100  -1.470107e+04
        1110  -1.470107e+04
        1120  -1.470107e+04
    
      0.794377 seconds (168.12 k allocations: 15.640 MiB, 0.80% gc time)





    VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(Array{Float64}(0,3), ([0.2808 0.279454 0.232256; 0.279454 0.284312 0.219977; 0.232256 0.219977 0.213052], [0.717243 0.661816 0.674158; 0.661816 0.715193 0.581746; 0.674158 0.581746 0.783965]), Array{Float64}(0,0), Char[], Float64[], -Inf, Inf)



Heritability from 3-variate estimate and their standard errors.


```julia
h, hse = heritability(trait57_model.Σ, Σcov)
[h'; hse']
```




    2×3 Array{Float64,2}:
     0.281351   0.284453  0.213689
     0.0778252  0.077378  0.084084



## 13-trait joint analysis

In some situations, such as studying the genetic covariance, we need to jointly analyze 13 traits. We first try the **Fisher scoring algorithm**.


```julia
# initialize model parameters
traitall_model = VarianceComponentModel(cg10kdata_rotated)
# estimate variance components using Fisher scoring algorithm
@time mle_fs!(traitall_model, cg10kdata_rotated; solver=:Ipopt, verbose=true)
```

    This is Ipopt version 3.12.4, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:        0
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:    16653
    
    Total number of variables............................:      182
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:        0
                         variables with only upper bounds:        0
    Total number of equality constraints.................:        0
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  1.3113371e+05 0.00e+00 1.00e+02   0.0 0.00e+00    -  0.00e+00 0.00e+00   0 
       5  8.2233766e+04 0.00e+00 6.03e+02 -11.0 2.32e+00    -  1.00e+00 1.00e+00f  1 MaxS
      10  1.1960260e+05 0.00e+00 8.76e+02 -11.0 6.20e+01  -5.4 1.00e+00 1.00e+00h  1 MaxS
      15  2.4416551e+05 0.00e+00 2.50e+02 -11.0 8.69e+02  -7.8 1.00e+00 1.00e+00f  1 MaxS



    DomainError:
    log will only return a complex result if called with a complex argument. Try log(complex(x)).

    

    Stacktrace:

     [1] nan_dom_err at ./math.jl:300 [inlined]

     [2] log at ./math.jl:419 [inlined]

     [3] logdet(::Array{Float64,2}) at ./linalg/generic.jl:1244

     [4] VarianceComponentModels.TwoVarCompModelRotate(::VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}) at /Users/huazhou/.julia/v0.6/VarianceComponentModels/src/VarianceComponentModels.jl:127

     [5] eval_f(::VarianceComponentModels.TwoVarCompOptProb{VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}},VarianceComponentModels.TwoVarCompVariateRotate{Float64,Array{Float64,2},Array{Float64,2}},Array{Float64,2},Array{Float64,1},VarianceComponentModels.VarianceComponentAuxData{Array{Float64,2},Array{Float64,1}}}, ::Array{Float64,1}) at /Users/huazhou/.julia/v0.6/VarianceComponentModels/src/two_variance_component.jl:683

     [6] (::Ipopt.#eval_f_cb#4{VarianceComponentModels.TwoVarCompOptProb{VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}},VarianceComponentModels.TwoVarCompVariateRotate{Float64,Array{Float64,2},Array{Float64,2}},Array{Float64,2},Array{Float64,1},VarianceComponentModels.VarianceComponentAuxData{Array{Float64,2},Array{Float64,1}}}})(::Array{Float64,1}) at /Users/huazhou/.julia/v0.6/Ipopt/src/IpoptSolverInterface.jl:53

     [7] eval_f_wrapper(::Int32, ::Ptr{Float64}, ::Int32, ::Ptr{Float64}, ::Ptr{Void}) at /Users/huazhou/.julia/v0.6/Ipopt/src/Ipopt.jl:89

     [8] solveProblem(::Ipopt.IpoptProblem) at /Users/huazhou/.julia/v0.6/Ipopt/src/Ipopt.jl:304

     [9] optimize!(::Ipopt.IpoptMathProgModel) at /Users/huazhou/.julia/v0.6/Ipopt/src/IpoptSolverInterface.jl:120

     [10] #mle_fs!#29(::Int64, ::Symbol, ::Symbol, ::Bool, ::Function, ::VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}, ::VarianceComponentModels.TwoVarCompVariateRotate{Float64,Array{Float64,2},Array{Float64,2}}) at /Users/huazhou/.julia/v0.6/VarianceComponentModels/src/two_variance_component.jl:893

     [11] (::VarianceComponentModels.#kw##mle_fs!)(::Array{Any,1}, ::VarianceComponentModels.#mle_fs!, ::VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}, ::VarianceComponentModels.TwoVarCompVariateRotate{Float64,Array{Float64,2},Array{Float64,2}}) at ./<missing>:0

     [12] include_string(::String, ::String) at ./loading.jl:515


From the output we can see the Fisher scoring algorithm ran into some numerical issues. Let's try the **MM algorithm**.


```julia
# reset model parameters
traitall_model = VarianceComponentModel(cg10kdata_rotated)
# estimate variance components using Fisher scoring algorithm
@time mle_mm!(traitall_model, cg10kdata_rotated; verbose=true)
```

    
         MM Algorithm
      Iter      Objective  
    --------  -------------
           0  -1.311337e+05
           1  -8.002108e+04
           2  -5.806935e+04
           3  -4.926111e+04
           4  -4.611059e+04
           5  -4.511606e+04
           6  -4.482679e+04
           7  -4.474294e+04
           8  -4.471496e+04
           9  -4.470174e+04
          10  -4.469246e+04
          20  -4.462243e+04
          30  -4.456888e+04
          40  -4.452774e+04
          50  -4.449601e+04
          60  -4.447134e+04
          70  -4.445199e+04
          80  -4.443665e+04
          90  -4.442436e+04
         100  -4.441442e+04
         110  -4.440630e+04
         120  -4.439961e+04
         130  -4.439405e+04
         140  -4.438938e+04
         150  -4.438544e+04
         160  -4.438210e+04
         170  -4.437923e+04
         180  -4.437676e+04
         190  -4.437463e+04
         200  -4.437277e+04
         210  -4.437115e+04
         220  -4.436972e+04
         230  -4.436846e+04
         240  -4.436735e+04
         250  -4.436636e+04
         260  -4.436548e+04
         270  -4.436469e+04
         280  -4.436399e+04
         290  -4.436335e+04
         300  -4.436278e+04
         310  -4.436226e+04
         320  -4.436179e+04
         330  -4.436137e+04
         340  -4.436098e+04
         350  -4.436063e+04
         360  -4.436030e+04
         370  -4.436001e+04
         380  -4.435974e+04
         390  -4.435949e+04
         400  -4.435926e+04
         410  -4.435905e+04
         420  -4.435886e+04
         430  -4.435868e+04
         440  -4.435851e+04
         450  -4.435836e+04
         460  -4.435822e+04
         470  -4.435809e+04
         480  -4.435797e+04
         490  -4.435785e+04
         500  -4.435775e+04
         510  -4.435765e+04
         520  -4.435756e+04
         530  -4.435747e+04
         540  -4.435739e+04
         550  -4.435732e+04
         560  -4.435725e+04
         570  -4.435718e+04
         580  -4.435712e+04
         590  -4.435706e+04
         600  -4.435701e+04
         610  -4.435696e+04
         620  -4.435691e+04
         630  -4.435687e+04
         640  -4.435683e+04
         650  -4.435679e+04
         660  -4.435675e+04
         670  -4.435671e+04
         680  -4.435668e+04
         690  -4.435665e+04
         700  -4.435662e+04
         710  -4.435659e+04
         720  -4.435657e+04
         730  -4.435654e+04
         740  -4.435652e+04
         750  -4.435649e+04
         760  -4.435647e+04
         770  -4.435645e+04
         780  -4.435643e+04
         790  -4.435642e+04
         800  -4.435640e+04
         810  -4.435638e+04
         820  -4.435637e+04
         830  -4.435635e+04
         840  -4.435634e+04
         850  -4.435633e+04
         860  -4.435631e+04
         870  -4.435630e+04
         880  -4.435629e+04
         890  -4.435628e+04
         900  -4.435627e+04
         910  -4.435626e+04
         920  -4.435625e+04
         930  -4.435624e+04
         940  -4.435623e+04
         950  -4.435622e+04
         960  -4.435621e+04
         970  -4.435621e+04
         980  -4.435620e+04
         990  -4.435619e+04
        1000  -4.435619e+04
        1010  -4.435618e+04
        1020  -4.435617e+04
        1030  -4.435617e+04
        1040  -4.435616e+04
        1050  -4.435616e+04
        1060  -4.435615e+04
        1070  -4.435615e+04
        1080  -4.435614e+04
        1090  -4.435614e+04
    
      3.551301 seconds (178.42 k allocations: 70.115 MiB, 0.42% gc time)





    (-44356.138529861186, VarianceComponentModels.VarianceComponentModel{Float64,2,Array{Float64,2},Array{Float64,2}}(Array{Float64}(0,13), ([0.272384 0.190358 … -0.128222 -0.0980655; 0.190358 0.21692 … -0.0689912 -0.0444349; … ; -0.128222 -0.0689912 … 0.118227 0.0909188; -0.0980655 -0.0444349 … 0.0909188 0.107456], [0.724562 0.56992 … -0.0590518 -0.124939; 0.56992 0.782639 … 0.0238629 0.0475408; … ; -0.0590518 0.0238629 … 0.880671 0.550889; -0.124939 0.0475408 … 0.550889 0.891929]), Array{Float64}(0,0), Char[], Float64[], -Inf, Inf), ([0.0111619 0.0131088 … 0.0128956 0.0127641; 0.0131091 0.0151759 … 0.017162 0.0171466; … ; 0.0128956 0.017162 … 0.0173994 0.0182002; 0.0127643 0.0171461 … 0.0182003 0.0187848], [0.0112235 0.0133094 … 0.0130111 0.0127861; 0.01331 0.0158262 … 0.017867 0.017798; … ; 0.013011 0.0178666 … 0.0179487 0.0187579; 0.012786 0.0177975 … 0.0187578 0.0193328]), [0.000124587 7.24074e-5 … -3.35716e-7 -1.40982e-5; 7.24411e-5 0.000171849 … -2.05381e-5 -3.17975e-6; … ; -3.60221e-7 -2.05683e-5 … 0.000351859 -1.5168e-5; -1.40799e-5 -3.16738e-6 … -1.51641e-5 0.000373756], Array{Float64}(0,13), Array{Float64}(0,0))



It converges after ~1000 iterations.

## Save analysis results


```julia
#using JLD
#@save "copd.jld"
#whos()
```
