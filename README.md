# VarianceComponentModels

[![Build Status](https://travis-ci.org/OpenMendel/VarianceComponentModels.jl.svg?branch=master)](https://travis-ci.org/OpenMendel/VarianceComponentModels.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/5yyf2m4y8p68glbh/branch/master?svg=true)](https://ci.appveyor.com/project/Hua-Zhou/variancecomponentmodels-jl-cw40h/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/OpenMendel/VarianceComponentModels.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/VarianceComponentModels.jl?branch=master)
[![codecov](https://codecov.io/gh/OpenMendel/VarianceComponentModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMendel/VarianceComponentModels.jl)

This [Julia](http://julialang.org/) package provides computational routines for fitting and testing variance component models. It is one component of the umbrella [OpenMendel](https://openmendel.github.io) project.

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/VarianceComponentModels.jl/latest)

## Installation

*Note: Three OpenMendel packages - [SnpArrays](https://github.com/OpenMendel/SnpArrays.jl), [Search](https://github.com/OpenMendel/Search.jl), and [MendelBase](https://github.com/OpenMendel/MendelBase.jl) must be installed before any Mendel analysis packages will run.*

Within Julia, use the package manager to install MendelVarianceComponentModels:

    Pkg.clone("git@github.com:OpenMendel/MendelVarianceComponentModels.jl.git")

This package supports Julia v0.4.

## Data Files

To run this analysis package you will need to prepare a Control file and have your data files available. The Control file holds the names of your data files and any optional parameters for the analysis. Details on the general format and contents of the Control and data files can be found on the OpenMendel [documentation page](https://openmendel.github.io/). Descriptions of the specific options available within the MendelVarianceComponentModels analysis package are in its [documentation page](https://openmendel.github.io/MendelVarianceComponentModels.jl).

There are example data files in the "docs" subfolder of the MendelVarianceComponentModels package, for example, ~/.julia/v0.4/MendelVarianceComponentModels/docs.

## Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:     julia> using MendelVarianceComponentModels

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:     julia> VarianceComponentModels("Control_file.txt")

*Note: The package is called* MendelVarianceComponentModels *but the analysis function is called simply* VarianceComponentModels.

## Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.

<!--- ## Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.