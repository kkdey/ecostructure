# ecostructure
R package for Clustering and Visualization of structure in ecological data

This package is targeted at fitting [STRUCTURE](http://www.genetics.org/content/155/2/945)
type model on ecological species abundance and presence-absence data both at local
and large geographic scales - together with high quality visualizations of the fitted
models.

### Authors 

- [Kushal K Dey](https://kkdey.github.io) 
- [Alexander E White](http://www.alexwhitebiology.com/)

Advised by 

- [Matthew Stephens](http://stephenslab.uchicago.edu/)
- [Trevor Price](https://pondside.uchicago.edu/ecol-evol/people/price.html)


### Citation

If you find **ecostructure**, please cite our papers

White, Alexander E. and Dey, Kushal K. and Mohan, Dhananjai and Stephens, Matthew and Price, Trevor D. *Regional influences on community structure across the tropical-temperate divide*. Nature Communications. 2019. 10 (1). 2646. 10.1038/s41467-019-10253-6.

White, Alexander E. and Dey, Kushal K. and and Stephens, Matthew and Price, Trevor D. *Dispersal syndromes drive the formation of biogeographical regions, illustrated by the case of Wallace’s Line*. Global Ecology and Biogeography. 2021.  https://doi.org/10.1111/geb.13250


### Installation

Install **ecostructure** following the instructions below.

```R
install.packages(devtools)
install.packages("sf")
devtools::install_github("kkdey/methClust")
devtools::install_github("kkdey/maptpx"). # this is an updated version of CRAN package maptpx
devtools::install_github("kkdey/CountClust")
devtools::install_github("kkdey/ecostructure")
```
Then load **ecostructure**

```R
library(ecostructure)
```

### Demo

Some examples of produced using our **ecostructure** package

<img src="bin/ecostructure.2.001.jpeg" alt="misc" height="600" width="1000" align = "middle">

If you want to try **ecostructure** and replicate figures like this, please check our tutorial [here](https://kkdey.github.io/ecostructure/).


### Questions?

For any queries or concerns related to the software, you can open an issue [here](https://github.com/kkdey/ecostructure/issues). Or you can contact 
us. Kushal Dey - *kshldey@gmail.com* or Alex White -
*aewhite100@gmail.com*.

We thank Prof. Trevor D. Price and Dhananjai Mohan for data collection for our
accompanying paper. We thank our mentors Prof. Matthew Stephens and 
Prof. Trevor Price for helpful suggestions and discussions. 

You are also most welcome to contribute to **ecostructure** !!









