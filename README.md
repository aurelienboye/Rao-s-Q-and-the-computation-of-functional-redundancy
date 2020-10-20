Exploration of the different implementations of RaoQ: which one to use
to compute functional redundancy indices ?
================

A. Boyé
12/03/2018, updated the 20 octobre, 2020

-----

 A clean version of what's below, with the equation properly printed, can be found [here](https://github.com/aurelienboye/Rao-s-Q-and-the-computation-of-functional-redundancy/Rao_redundancy.html). 

-----

# The different implementations of Rao’s quadratic entropy in R

## Package SYNCSA

`rao.diversity` of the `SYNCSA` package uses the formula found in Rao
([1982](#ref-rao1982)), which is the most commonly cited version of this
index (e.g. Botta-Dukát [2005](#ref-botta2005rao); Ricotta
[2005](#ref-ricotta2005); Lepš *et al.* [2006](#ref-leps2006))
:

\(RaoQ = \sum \limits_{i=1}^{S-1} \sum \limits_{j=i+1}^S d_{ij} p_i p_j\)

with \(d_{ij}\) the trait dissimilarity between species \(i\) and \(j\),
and with \(p_i\) and \(p_j\) their relative abundance

This `rao.diversity` function only uses relative abundances (**contrary
to the `divc` function of ADE4** that can also use raw abundance; [see
the code of the different functions](#code_details) and the help of
[rao.diversity](https://www.rdocumentation.org/packages/SYNCSA/versions/1.3.2/topics/rao.diversity)
and [divc](https://pbil.univ-lyon1.fr/ade4/ade4-html/divc.html)).

Also, `rao.diversity` the square-root of Gower’s distance and therefore
computes
:

\(RaoQ = \sum \limits_{i=1}^{S-1} \sum \limits_{j=i+1}^S \sqrt{d^{gower}_{ij}} p_i p_j\)

with \(\sum \limits_{i=1}^S p_i = 1\)

> This transformation makes Gower’s distance *metric*, and likely
> *euclidean* (follwoing Legendre & Legendre
> [2012](#ref-Legendre_Legendre2012), p . 297, S15)

## Package FD

The `dbFD` function from the `FD` package is based on the `divc` from
the `ADE4` package that computes Rao’s Q based on Champely & Chessel
([2002](#ref-Champely_Chessel_2002)) formula :
\(RaoQ = \sum \limits_{i=1}^n \sum \limits_{j=1}^n p_i p_j \frac{(\delta_{ij})^2}{2}\)

The exact formula used in `divc` is :
\(RaoQ = \sum \limits_{i=1}^n \sum \limits_{j=1}^n p_i p_j \delta_{ij}^2 \times \frac{1}{2} \times \frac{1}{\sum p^2}\)

with \(\sum p^2 = 1\) in the case of relative abundances

> It is also possible in `divc` and `dbFD` to use raw abundances in
> which case **the influence of divinding by \(\sum p^2\) is not clear
> to me**
> 
> In the book *Quantifying Functional Biodiversity* (Pla *et al.*
> [2011](#ref-pla2011)) it is said p.40 that an unbiasied estimator of
> Rao’s Q when using raw abundances would be
> :
> 
> \(\hat{Rao} = \frac{n}{n-1} 2 \sum \limits_{i>j}^S d_{ij} \frac{n_i n_j}{n^2}\)
> 
> with \(n_i\) the number of individuals of species \(i\) and
> \(n = \sum \limits_{i=1}^S n_i\) the total number of individuals
> 
> This formula is close to the one found in Schleuter *et al.*
> ([2010](#ref-schleuter2010) Table 1 Eq. 3.5) that defines Rao’s Q
> according to Rao ([1982](#ref-rao1982)) and Champely & Chessel
> ([2002](#ref-Champely_Chessel_2002))
> as:
> 
> \(FD_Q = \sum \limits_{S \in S_c} \sum \limits_{S' \in S_c} \frac{A_sA_{s'}}{A^2} dist(s,s')\)

Champely & Chessel ([2002](#ref-Champely_Chessel_2002)) and after them
Pavoine (Pavoine [2005](#ref-pavoine2005phd), [2012](#ref-pavoine2012);
Pavoine *et al.* [2005](#ref-pavoine2005_rao_dissim); Pavoine & Bonsall
[2009](#ref-pavoine_bonsall2009)) have adapted the original formula of
Rao ([1982](#ref-rao1982))
as:

\(RaoQ = \sum \limits_{i=1}^n \sum \limits_{j=1}^n p_i p_j \frac{(\delta_{ij})^2}{2}\)
for

1)  maximizing this function *i.e.* estimate the maximum value of
    \(RaoQ\) possible for the set
    \(p=(p_1...p_n), p_i \geqslant 0, \sum \limits_{ i=1}^n p_i=1\)
    (*e.g.* Champely & Chessel [2002](#ref-Champely_Chessel_2002);
    Pavoine *et al.* [2005](#ref-pavoine2005_rao_dissim); Pavoine &
    Bonsall [2009](#ref-pavoine_bonsall2009) ; Pavoine
    [2012](#ref-pavoine2012))

2)  Generalize Rao’s Q and link its formulation to the one of the
    variance (Pavoine [2005](#ref-pavoine2005phd) *p. *51-54; Pavoine
    *et al.* [2005](#ref-pavoine2005_rao_dissim))

de Bello *et al.* ([2016](#ref-deBello2016rao)) noted that `divc` et
`dbFD` express Rao’s Q in terms of variance and that therefore
**Simpson’s index is not the maximum value possible for this version
of Rao’s Q** (which is expected from Rao [1982](#ref-rao1982) original
formula when the dissimilarity all species pair is 1; Botta-Dukát
[2005](#ref-botta2005rao)). This is critical when measuring functional
redundancy (see [Methods to compute functional redundancy](#redond)).

> **de Bello *et al.* ([2016](#ref-deBello2016rao))**
> 
> Both the Simpson index and the Rao index can be computed with the
> “divc” function in the “ade4” package in R. The same function is
> also implemented in the widely used “dbFD” function in the “FD”
> package. However, this function computes trait variance by squaring
> dissimilar- ity values and dividing by two (de Bello et al. 2011).
> Consequently, it does not produce Rao as such, but a speci c case of
> Rao, expressed in terms of variance. We recall here that one of the
> properties expected from the Rao index that it is a generalization of
> Simpson, so that when all species are completely dissimilar (i.e., dij
> = 1 for every i and j), an ideal algorithm should provide Simpson
> diversity. Simpson diversity is obtained with “divc” only when dij =
> sqrt(2) for all species pairs, which does not seem intuitive and it is
> certainly not explained in the help function. In particular, the upper
> bound for Rao diversity (when the dissimilarity between all species
> pairs is equal to 1, i.e., dij = 1) with the “divc” function is half
> of the Simpson diversity, not the Simpson diversity obtained with Eqs.
> 2 and 3.

However, the formula used in `divc` et `dbFD` allows to estimate the
maximum value possible of Rao’s Q (for a fixed dissimilarity matrix
among species) and therefore allows to express Rao’s Q as a relative
value which might prove usefull when to compare samples that have
important richness differences (Pla *et al.* [2011](#ref-pla2011))

> **Pla *et al.* ([2011](#ref-pla2011)) p. 42-43**
> 
> The expression of quadratic entropy as an **absolute value is not
> useful when the comparisons have to be done between communities with
> very different numbers of species or when different sets of traits
> were used to define the distance matrix. To get a relative expression
> the maximum has to be estimated from the data. The distance matrix
> does not depend on the abundance of species and is fixed for a given
> set of species, but changes in the relative abundance of these species
> may lead to the maximum diversity index (\(Rao_{max}\))**. There are
> two types of abundance vectors that define two subclasses of maximum:
> (a) weak maximization, when some of the \(w_i\) abundances that
> maximize \(Rao_{max}\) are zero; and (b) strong maximization, when all
> the \(w_i\) values that maximize \(Rao_{max}\) are positives.
> 
> The maximization process relies on the dissimilarity matrix and on any
> ultra-metric matrix that belongs to the strong subclass (Pavoine et
> al. 2005). The drawback arising from having only some species to
> maximize the Rao’s quadratic index when dissimilarity between species
> are based on functional traits is the absence of distance measures
> that guaranty the ultrametric condition and then ecological meaningful
> expression of the functional diversity using relative Rao index.
> Taxonomic or phylogenic dissimilarity trees may have ultrametric
> distances and give a maximum value of Rao that relies on total
> abundance distributed among all the species presents.

## Other versions of Rao’s Q in R

### Package `picante`

The `raoD` in the `picante` package seems to compute Rao’s Q only with
phylogenetic diversity ([see the package
description](https://cran.r-project.org/web/packages/picante/picante.pdf)).

It uses relative abundances but with a new and different formula from
those introduced here before. This function expresses Rao’s Q as an
abundance-weighted version of the MPD (Mean Pairwise Distance). There is
however a difference between the original formula of Rao’s Q and this
one in the way the dissimilarity of a species with itself (which could
represent intra-specific trait variability) is handled (c.f. de Bello
*et al.* [2016](#ref-deBello2016rao)).

> **Details** :
> 
> \(D_{kl}=\sum t_{ij} x_{ki} x_{lj}\) with \(x_{ki}\) the relative
> abundance of species \(i\) in the community \(k\) and \(t_{ij}\) the
> phylogenetic distance matrix ([see the code of the
> function](#code_details_raoD))

### `melodic` from de Bello *et al.* ([2016](#ref-deBello2016rao))

de Bello *et al.* ([2016](#ref-deBello2016rao)) gives in the
Supplementary a function (`melodic`) calculating both MPD and Rao’s Q

-----

# Behaviour of the `divc` and `rao.diversity` functions

``` r
# Packages 
library(SYNCSA)
```

    ## Warning: package 'SYNCSA' was built under R version 3.5.2

``` r
library(FD)
```

    ## Warning: package 'ape' was built under R version 3.5.2

    ## Warning: package 'geometry' was built under R version 3.5.2

    ## Warning: package 'vegan' was built under R version 3.5.2

    ## Warning: package 'permute' was built under R version 3.5.2

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.5.2

``` r
library(magrittr)
```

## Difference between the two functions

We use the `dummy` dataset from the examples of `dbFD`

Rao’s Q computed with `dbFD` on raw abundances

``` r
# dbFD
FD_raw <- dbFD(dummy$trait, dummy$abun, calc.FRic = FALSE)
```

Rao’s Q computed with `dbFD` on relative abundances

``` r
# dbFD
abun_rel <- dummy$abun / rowSums(dummy$abun)

FD_rel <- dbFD(dummy$trait, abun_rel, calc.FRic = FALSE)
```

Rao’s Q computed with `divc` with both raw and relative abundances

``` r
# On raw abundances
divc_raw <- divc(as.data.frame(t(dummy$abun)),gowdis(dummy$trait))
```

    ## Warning in divc(as.data.frame(t(dummy$abun)), gowdis(dummy$trait)):
    ## Euclidean property is expected for distance

``` r
# On relative abundances
divc_rel <- divc(as.data.frame(t(abun_rel)),gowdis(dummy$trait))
```

    ## Warning in divc(as.data.frame(t(abun_rel)), gowdis(dummy$trait)): Euclidean
    ## property is expected for distance

Rao’s Q computed with `rao.diversity` with both raw and relative
abundances

``` r
# On raw abundances
SYN_raw <- rao.diversity(dummy$abun, dummy$trait)
```

    ## Warning: Warning: NA in traits matrix

``` r
# On relative abundances
SYN_rel <- rao.diversity(abun_rel, dummy$trait)
```

    ## Warning: Warning: NA in traits matrix

Comparison

``` r
# dbFD raw abundances
FD_raw$RaoQ
```

    ##       com1       com2       com3       com4       com5       com6 
    ## 0.12835440 0.04063622 0.06497224 0.03003235 0.11858388 0.11738479 
    ##       com7       com8       com9      com10 
    ## 0.07868047 0.09308505 0.12431265 0.12490783

``` r
#dbFD relative abundances
FD_rel$RaoQ
```

    ##       com1       com2       com3       com4       com5       com6 
    ## 0.12835440 0.04063622 0.06497224 0.03003235 0.11858388 0.11738479 
    ##       com7       com8       com9      com10 
    ## 0.07868047 0.09308505 0.12431265 0.12490783

``` r
# divc raw abundances
divc_raw
```

    ##        diversity
    ## com1  0.12835440
    ## com2  0.04063622
    ## com3  0.06497224
    ## com4  0.03003235
    ## com5  0.11858388
    ## com6  0.11738479
    ## com7  0.07868047
    ## com8  0.09308505
    ## com9  0.12431265
    ## com10 0.12490783

``` r
# divc relative abundances
divc_rel
```

    ##        diversity
    ## com1  0.12835440
    ## com2  0.04063622
    ## com3  0.06497224
    ## com4  0.03003235
    ## com5  0.11858388
    ## com6  0.11738479
    ## com7  0.07868047
    ## com8  0.09308505
    ## com9  0.12431265
    ## com10 0.12490783

``` r
# SYNCSA (rao.diversity) raw abundances
SYN_raw 
```

    ## $call
    ## rao.diversity(comm = dummy$abun, traits = dummy$trait)
    ## 
    ## $Simpson
    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.6562500 0.5312500 0.6111111 0.2187500 0.6562500 0.7573696 0.6446281 
    ##      com8      com9     com10 
    ## 0.6280992 0.7438017 0.6250000 
    ## 
    ## $FunRao
    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.5150748 0.3202589 0.4282582 0.1628977 0.5038237 0.5580249 0.4305683 
    ##      com8      com9     com10 
    ## 0.4663877 0.5580015 0.4947682 
    ## 
    ## $FunRedundancy
    ##       com1       com2       com3       com4       com5       com6 
    ## 0.14117516 0.21099107 0.18285296 0.05585227 0.15242633 0.19934470 
    ##       com7       com8       com9      com10 
    ## 0.21405980 0.16171148 0.18580012 0.13023176

``` r
# SYNCSA (rao.diversity) relative abundances
SYN_rel
```

    ## $call
    ## rao.diversity(comm = abun_rel, traits = dummy$trait)
    ## 
    ## $Simpson
    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.6562500 0.5312500 0.6111111 0.2187500 0.6562500 0.7573696 0.6446281 
    ##      com8      com9     com10 
    ## 0.6280992 0.7438017 0.6250000 
    ## 
    ## $FunRao
    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.5150748 0.3202589 0.4282582 0.1628977 0.5038237 0.5580249 0.4305683 
    ##      com8      com9     com10 
    ## 0.4663877 0.5580015 0.4947682 
    ## 
    ## $FunRedundancy
    ##       com1       com2       com3       com4       com5       com6 
    ## 0.14117516 0.21099107 0.18285296 0.05585227 0.15242633 0.19934470 
    ##       com7       com8       com9      com10 
    ## 0.21405980 0.16171148 0.18580012 0.13023176

> **Summary**
> 
>   - `divc` and `rao.diversity` give different results
>   - `divc`, `dbFD` and `rao.diversity` give the same results whether
>     raw or relative abundances are used as input
>   - `divc` send back a warning: it needs the distance matrix to be
>     euclidean
>       - In contrast, the distance matrix is not always euclidean in
>         `dbFD` because it uses the uncorrected species-species
>         distance matrix [see the help of
>         `dbFD`](https://www.rdocumentation.org/packages/FD/versions/1.0-12/topics/dbFD)).

## Equality of RaoQ & Simpson when \(d_{ij}=1\)

Rao’s Q corresponds to the expected dissimilarity between two
individuals drawn at random from the community and is the counterpart of
**Simpson’s diversity index** which is the probability for two
individuals drawn at random to belong to to the same species. The
Simpson index is thus equivalent to Rao’s Q when all species are
functionally completely different *i.e.* unique. Indeed if \(d_{ij}=1\)
for all \(i \neq j\)
then

\(RaoQ = \sum \limits_{i=1}^{S-1} \sum \limits_{j=i+1}^S d_{ij} p_i p_j\)

simplfies to Simpson’s diversity index
\(1-\sum \limits_{i=1}^{S} p_i^2\) (Botta-Dukát
[2005](#ref-botta2005rao))

de Bello *et al.* ([2016](#ref-deBello2016rao)) noted that this is not
true with the “variance”-formula used in `divc` et `dbFD`.

<u> Verification: </u>

``` r
# Create a distance matrix with all species-species distance equal to 1
dist_sp<- 1-diag(8) %>%
  # Add species names
  set_rownames(rownames(dummy$trait)) %>%
  # Transform into a distance amtrix
  as.dist(.)

dist_sp
```

    ##     sp1 sp2 sp3 sp4 sp5 sp6 sp7
    ## sp2   1                        
    ## sp3   1   1                    
    ## sp4   1   1   1                
    ## sp5   1   1   1   1            
    ## sp6   1   1   1   1   1        
    ## sp7   1   1   1   1   1   1    
    ## sp8   1   1   1   1   1   1   1

``` r
# Computation with dbFD
ex2_FD <- dbFD(dist_sp, dummy$abun, calc.FRic = FALSE)
```

``` r
ex2_FD$RaoQ
```

    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.3281250 0.2656250 0.3055556 0.1093750 0.3281250 0.3786848 0.3223140 
    ##      com8      com9     com10 
    ## 0.3140496 0.3719008 0.3125000

``` r
# Test with rao.diversity
#------------------

# Create a trait matrix so that all Gower distances between species paires are equal to 1
trait_sp <- as.data.frame(matrix(LETTERS[1:8],ncol=8,nrow=8)) %>%
  # Add species names
  set_rownames(rownames(dummy$trait))

# Verification
gowdis(trait_sp)
```

    ##     sp1 sp2 sp3 sp4 sp5 sp6 sp7
    ## sp2   1                        
    ## sp3   1   1                    
    ## sp4   1   1   1                
    ## sp5   1   1   1   1            
    ## sp6   1   1   1   1   1        
    ## sp7   1   1   1   1   1   1    
    ## sp8   1   1   1   1   1   1   1

``` r
# Computation with rao.diversity
ex2_SYN <- rao.diversity(dummy$abun, trait_sp)
ex2_SYN
```

    ## $call
    ## rao.diversity(comm = dummy$abun, traits = trait_sp)
    ## 
    ## $Simpson
    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.6562500 0.5312500 0.6111111 0.2187500 0.6562500 0.7573696 0.6446281 
    ##      com8      com9     com10 
    ## 0.6280992 0.7438017 0.6250000 
    ## 
    ## $FunRao
    ##      com1      com2      com3      com4      com5      com6      com7 
    ## 0.6562500 0.5312500 0.6111111 0.2187500 0.6562500 0.7573696 0.6446281 
    ##      com8      com9     com10 
    ## 0.6280992 0.7438017 0.6250000 
    ## 
    ## $FunRedundancy
    ##  com1  com2  com3  com4  com5  com6  com7  com8  com9 com10 
    ##     0     0     0     0     0     0     0     0     0     0

> **Summary**
> 
>   - As expected, only the `rao.diversity` function simplifies to
>     Simpson’s index when \(d_{ij}=1\) while `dbFD` gives
>     \(Simpson~D \div 2\)
> 
> This specificity of `divc` was highlighted in de Bello *et al.*
> ([2016](#ref-deBello2016rao))

-----

# Summary

## Methods to compute functional redundancy

Ricotta *et al.* ([2016](#ref-ricotta2016)) has summarized and analyzed
the behavior of the different existing methods to measure redundancy.
Some of them are based on the relationship between Rao’s quadratic
entropy and the Simpson diversity index such as:

The functional redundancy of de Bello *et al.*
([2007](#ref-deBello2007)) :

> \(FRED = SIMD - RaoQ\) with SIMD Simpson’s index \((1-D)\)

de Bello *et al.* ([2016](#ref-deBello2016rao)) another index based on
the *taxonomic distinctness* (\(\Delta^*\)) proposed by Warwick & Clarke
([1995](#ref-warwick1995)) express here in terms of Rao’s Q and
Simpson’s index
> :

> \(\Delta^* = \frac{\sum \limits_{i \neq j} p_i p_j d_{ij}}{\sum \limits_{i \neq j} p_i p_j} = \frac{RaoQ}{SimD}\)

(<span class="citeproc-not-found" data-reference-id="vanderlinden2016">**???**</span>)
proposed another measure with these two indices :

> \(FRED = 1 - (\frac{RaoQ}{SIMD})\)

Finally, one new article of interest if you want to use functional
redundancy indices is Galland *et al.*
([2020](#ref-GALLAND2020106488))

## Most appropriated function to compute Rao’s Q in the context of functional redundancy ?

The `rao.diversity` function appears as the most appropriate **to
measure functional redundancy** as it allows for the equivalence between
\(RaoQ\) and \(Simpson~D\) when species are all functionally unique
(*c.f* Botta-Dukát [2005](#ref-botta2005rao)), which is essential to
compare both indices when measuring redundancy. In this perspective,
note that the formula mentioned in de Bello *et al.*
([2007](#ref-deBello2007)) is indeed
\(RaoQ = \sum \limits_{i=1}^{S-1} \sum \limits_{j=i+1}^S d_{ij} p_i p_j\).

## Additional comments

Rao’s Q could account for intra-specific variability by giving non-zero
distances in the diagonal of species-species distance matrix (de Bello
*et al.* [2016](#ref-deBello2016rao))

> This index can be understood as the average of all dissimilarities in
> the dissimilarity matrix, including the diagonal (Online resource 1 in
> the ESM), which represents the dissimilarity between one species and
> itself (generally assumed to be zero in most existing approaches)

It is also possible to **partition** Rao’s Q to test for the effect of
multiple factors (space, time, habitat…) in an ANOVA-like approach or an
AMOVA-like approach for hierarchical structures (c.f. Rao
[2010](#ref-rao2010); et Pavoine [2012](#ref-pavoine2012))

> ANOQE has more potential than AMOVA. It can be employed to analyse
> nested, crossed, fixed, random or mixed factors, and has sev- eral
> ramifications allowing multivariate factorial analysis (Pavoine,
> Dufour & Chessel 2004; Pavoine & Bailly 2007) as well as decomposition
> of diversity across a taxonomic (Ricotta 2005) or phylogenetic
> (Pavoine, Love & Bonsall 2009; Pavoine, Baguette & Bonsall 2010) tree,
> all of which render the approach practical for use by biologists.
> Recently, Rao (2010) re-emphasised the importance that ANOQE is likely
> to have for future research on diversity. \*\*(in Pavoine
> [2012](#ref-pavoine2012))

-----

# Code of the functions

## `dbFD` et `divc`

``` r
library(FD)

divc
```

    ## function (df, dis = NULL, scale = FALSE) 
    ## {
    ##     if (!inherits(df, "data.frame")) 
    ##         stop("Non convenient df")
    ##     if (any(df < 0)) 
    ##         stop("Negative value in df")
    ##     if (!is.null(dis)) {
    ##         if (!inherits(dis, "dist")) 
    ##             stop("Object of class 'dist' expected for distance")
    ##         if (!is.euclid(dis)) 
    ##             warning("Euclidean property is expected for distance")
    ##         dis <- as.matrix(dis)
    ##         if (nrow(df) != nrow(dis)) 
    ##             stop("Non convenient df")
    ##         dis <- as.dist(dis)
    ##     }
    ##     if (is.null(dis)) 
    ##         dis <- as.dist((matrix(1, nrow(df), nrow(df)) - diag(rep(1, 
    ##             nrow(df)))) * sqrt(2))
    ##     div <- as.data.frame(rep(0, ncol(df)))
    ##     names(div) <- "diversity"
    ##     rownames(div) <- names(df)
    ##     for (i in 1:ncol(df)) {
    ##         if (sum(df[, i]) < 1e-16) 
    ##             div[i, ] <- 0
    ##         else div[i, ] <- (t(df[, i]) %*% (as.matrix(dis)^2) %*% 
    ##             df[, i])/2/(sum(df[, i])^2)
    ##     }
    ##     if (scale == TRUE) {
    ##         divmax <- divcmax(dis)$value
    ##         div <- div/divmax
    ##     }
    ##     return(div)
    ## }
    ## <bytecode: 0x7fa324f92580>
    ## <environment: namespace:ade4>

``` r
dbFD
```

    ## function (x, a, w, w.abun = TRUE, stand.x = TRUE, ord = c("podani", 
    ##     "metric"), asym.bin = NULL, corr = c("sqrt", "cailliez", 
    ##     "lingoes", "none"), calc.FRic = TRUE, m = "max", stand.FRic = FALSE, 
    ##     scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward", 
    ##     km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, km.crit = c("calinski", 
    ##         "ssi"), calc.CWM = TRUE, CWM.type = c("dom", "all"), 
    ##     calc.FDiv = TRUE, dist.bin = 2, print.pco = FALSE, messages = TRUE) 
    ## {
    ##     tol <- .Machine$double.eps
    ##     corr <- match.arg(corr)
    ##     ord <- match.arg(ord)
    ##     CWM.type <- match.arg(CWM.type)
    ##     km.crit <- match.arg(km.crit)
    ##     if (!is.logical(messages)) 
    ##         stop("'messages' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(stand.FRic)) 
    ##         stop("'stand.FRic' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(stand.x)) 
    ##         stop("'stand.x' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(w.abun)) 
    ##         stop("'w.abun' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(calc.FRic)) 
    ##         stop("'calc.FRic' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(calc.FDiv)) 
    ##         stop("'calc.FDiv' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(calc.FGR)) 
    ##         stop("'calc.FGR' musts be TRUE or FALSE.", "\n")
    ##     if (!is.logical(calc.CWM)) 
    ##         stop("'calc.CWM' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(scale.RaoQ)) 
    ##         stop("'scale.RaoQ' must be TRUE or FALSE.", "\n")
    ##     if (!is.logical(print.pco)) 
    ##         stop("'print.pco' must be TRUE or FALSE.", "\n")
    ##     if (is.matrix(x) | is.data.frame(x)) {
    ##         is.dist.x <- FALSE
    ##         s.x <- dim(x)[1]
    ##         t.x <- dim(x)[2]
    ##         if (is.null(row.names(x))) 
    ##             stop("'x' must have row names.", "\n")
    ##         else x.rn <- row.names(x)
    ##     }
    ##     if (is.vector(x) | is.factor(x)) {
    ##         is.dist.x <- FALSE
    ##         s.x <- length(x)
    ##         t.x <- 1
    ##         if (is.null(names(x))) 
    ##             stop("'x' must have names.", "\n")
    ##         else x.rn <- names(x)
    ##     }
    ##     if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    ##         is.dist.x <- TRUE
    ##         s.x <- attr(x, "Size")
    ##         t.x <- 1
    ##         if (is.null(attr(x, "Labels"))) 
    ##             stop("'x' must have labels.", "\n")
    ##         else x.rn <- attr(x, "Labels")
    ##     }
    ##     if (missing(a)) {
    ##         ab.names <- list("Community1", x.rn)
    ##         a <- matrix(1, 1, s.x, dimnames = ab.names)
    ##     }
    ##     else {
    ##         if (is.matrix(a) | is.data.frame(a)) {
    ##             s.a <- dim(a)[2]
    ##             ab.t <- t(a)
    ##             if (is.null(row.names(ab.t))) 
    ##                 stop("'a' must have column names.", "\n")
    ##             else ab.t.row <- row.names(ab.t)
    ##             a <- as.matrix(a)
    ##         }
    ##         if (is.vector(a)) {
    ##             s.a <- length(a)
    ##             if (is.null(names(a))) 
    ##                 stop("'a' must have names.", "\n")
    ##             else ab.t.row <- names(a)
    ##             ab.names <- list("Community1", ab.t.row)
    ##             a <- matrix(a, 1, s.a, dimnames = ab.names)
    ##         }
    ##         if (s.x != s.a) 
    ##             stop("Different number of species in 'x' and 'a'.", 
    ##                 "\n")
    ##         if (any(ab.t.row != x.rn)) 
    ##             stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
    ##                 "\n")
    ##     }
    ##     a <- as.matrix(a)
    ##     a[which(is.na(a))] <- 0
    ##     abun.sum <- apply(a, 1, sum)
    ##     if (any(abun.sum == 0)) 
    ##         stop("At least one community has zero-sum abundances (no species).", 
    ##             "\n")
    ##     abun.sum2 <- apply(a, 2, sum)
    ##     if (any(abun.sum2 == 0)) 
    ##         stop("At least one species does not occur in any community (zero total abundance across all communities).", 
    ##             "\n")
    ##     if (!missing(w) & is.dist.x) 
    ##         stop("When 'x' is a distance matrix, 'w' should be left missing.", 
    ##             "\n")
    ##     if (!missing(w) & !is.dist.x) {
    ##         if (!is.numeric(w) | length(w) != t.x) 
    ##             stop("'w' should be a numeric vector of length = number of traits.", 
    ##                 "\n")
    ##         else w <- w/sum(w)
    ##     }
    ##     if (missing(w)) 
    ##         w <- rep(1, t.x)/sum(rep(1, t.x))
    ##     if (is.matrix(x) | is.data.frame(x)) {
    ##         x <- data.frame(x)
    ##         if (t.x >= 2) {
    ##             x.class <- sapply(x, data.class)
    ##             if (any(x.class == "character")) 
    ##                 x[, x.class == "character"] <- as.factor(x[, 
    ##                   x.class == "character"])
    ##             else x <- x
    ##             if (all(x.class == "numeric") & all(!is.na(x))) {
    ##                 if (length(unique(w)) == 1) {
    ##                   x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
    ##                   x.dist <- dist(x.s)
    ##                 }
    ##                 else {
    ##                   x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
    ##                 }
    ##             }
    ##             else {
    ##                 x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
    ##             }
    ##         }
    ##         if (t.x == 1) {
    ##             if (is.numeric(x[, 1])) {
    ##                 if (all(!is.na(x))) {
    ##                   x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
    ##                   x.dist <- dist(x.s)
    ##                 }
    ##                 if (any(is.na(x))) {
    ##                   pos.NA <- which(is.na(x), arr.ind = TRUE)
    ##                   x <- na.omit(x)
    ##                   x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
    ##                   x.dist <- dist(x.s)
    ##                   row.excl.ab <- pos.NA[, 1]
    ##                   a <- a[, -row.excl.ab]
    ##                   if (messages) 
    ##                     cat("Warning: Species with missing trait values have been excluded.", 
    ##                       "\n")
    ##                 }
    ##             }
    ##             if (is.factor(x[, 1]) | is.character(x[, 1])) {
    ##                 if (is.ordered(x[, 1])) 
    ##                   x <- x
    ##                 else x[, 1] <- as.factor(x[, 1])
    ##                 if (any(is.na(x))) {
    ##                   pos.NA <- which(is.na(x), arr.ind = TRUE)
    ##                   x <- na.omit(x)
    ##                   row.excl.ab <- pos.NA[, 1]
    ##                   a <- a[, -row.excl.ab]
    ##                   x.rn <- x.rn[-pos.NA]
    ##                   if (messages) 
    ##                     cat("Warning: Species with missing trait values have been excluded.", 
    ##                       "\n")
    ##                 }
    ##                 if (is.ordered(x[, 1])) {
    ##                   x.s <- data.frame(rank(x[, 1]))
    ##                   names(x.s) <- x.rn
    ##                   x.dist <- dist(x.s)
    ##                 }
    ##                 else {
    ##                   x.f <- as.factor(x[, 1])
    ##                   x.dummy <- diag(nlevels(x.f))[x.f, ]
    ##                   x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    ##                   sequence <- 1:10
    ##                   if (all(dist.bin != sequence[any(sequence)])) 
    ##                     stop("'dist.bin' must be an integer between 1 and 10.", 
    ##                       "\n")
    ##                   x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    ##                 }
    ##             }
    ##         }
    ##     }
    ##     if (is.vector(x) & is.numeric(x)) {
    ##         if (any(is.na(x))) {
    ##             pos.NA <- which(is.na(x))
    ##             x <- na.omit(x)
    ##             a <- a[, -pos.NA]
    ##             x.rn <- x.rn[-pos.NA]
    ##             if (messages) 
    ##                 cat("Warning: Species with missing trait values have been excluded.", 
    ##                   "\n")
    ##         }
    ##         else x <- x
    ##         x.s <- scale(x, center = T, scale = stand.x)
    ##         x.dist <- dist(x.s)
    ##         x <- data.frame(x)
    ##         dimnames(x) <- list(x.rn, "Trait")
    ##     }
    ##     if (is.vector(x) & is.character(x)) {
    ##         x <- as.factor(x)
    ##         if (any(is.na(x))) {
    ##             pos.NA <- which(is.na(x))
    ##             x <- na.omit(x)
    ##             a <- a[, -pos.NA]
    ##             x.rn <- x.rn[-pos.NA]
    ##             if (messages) 
    ##                 cat("Warning: Species with missing trait values have been excluded.", 
    ##                   "\n")
    ##         }
    ##         else x <- x
    ##         dimnames(x) <- list(x.rn, "Trait")
    ##         x.dummy <- diag(nlevels(x))[x, ]
    ##         x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    ##         sequence <- 1:10
    ##         if (all(dist.bin != sequence[any(sequence)])) 
    ##             stop("'dist.bin' must be an integer between 1 and 10.", 
    ##                 "\n")
    ##         x <- data.frame(x)
    ##         x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    ##     }
    ##     if (is.ordered(x)) {
    ##         if (any(is.na(x))) {
    ##             pos.NA <- which(is.na(x))
    ##             x <- na.omit(x)
    ##             a <- a[, -pos.NA]
    ##             x.rn <- x.rn[-pos.NA]
    ##             cat("Warning: Species with missing trait values have been excluded.", 
    ##                 "\n")
    ##         }
    ##         else x <- x
    ##         x <- data.frame(x)
    ##         dimnames(x) <- list(x.rn, "Trait")
    ##         x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
    ##     }
    ##     if (is.factor(x) & !is.ordered(x)) {
    ##         if (any(is.na(x))) {
    ##             pos.NA <- which(is.na(x))
    ##             x <- na.omit(x)
    ##             a <- a[, -pos.NA]
    ##             x.rn <- x.rn[-pos.NA]
    ##             if (messages) 
    ##                 cat("Warning: Species with missing trait values have been excluded.", 
    ##                   "\n")
    ##         }
    ##         else x <- x
    ##         x.dummy <- diag(nlevels(x))[x, ]
    ##         x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    ##         sequence <- 1:10
    ##         if (all(dist.bin != sequence[any(sequence)])) 
    ##             stop("'dist.bin' must be an integer between 1 and 10.", 
    ##                 "\n")
    ##         x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    ##         x <- data.frame(x)
    ##         dimnames(x) <- list(x.rn, "Trait")
    ##     }
    ##     if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    ##         if (any(is.na(x))) 
    ##             stop("When 'x' is a distance matrix, it cannot have missing values (NA).", 
    ##                 "\n")
    ##         x.dist <- x
    ##     }
    ##     if (any(is.na(x.dist))) 
    ##         stop("NA's in the distance matrix.", "\n")
    ##     if (!is.dist.x) {
    ##         no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    ##         if (any(no.traits == 0)) 
    ##             stop("At least one species has no trait data.", "\n")
    ##     }
    ##     c <- dim(a)[1]
    ##     if (!w.abun) 
    ##         for (h in 1:c) {
    ##             abpos <- which(a[h, ] > 0)
    ##             a[h, abpos] <- 1
    ##         }
    ##     attr(x.dist, "Labels") <- x.rn
    ##     if (is.euclid(x.dist)) 
    ##         x.dist2 <- x.dist
    ##     if (!is.euclid(x.dist)) {
    ##         if (corr == "lingoes") {
    ##             x.dist2 <- lingoes(x.dist)
    ##             if (messages) 
    ##                 cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.", 
    ##                   "\n")
    ##         }
    ##         if (corr == "cailliez") {
    ##             x.dist2 <- cailliez(x.dist)
    ##             if (messages) 
    ##                 cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.", 
    ##                   "\n")
    ##         }
    ##         if (corr == "sqrt") {
    ##             x.dist2 <- sqrt(x.dist)
    ##             if (!is.euclid(x.dist2)) 
    ##                 stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.", 
    ##                   "\n")
    ##             if (is.euclid(x.dist2)) 
    ##                 if (messages) 
    ##                   cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.", 
    ##                     "\n")
    ##         }
    ##         if (corr == "none") {
    ##             x.dist2 <- quasieuclid(x.dist)
    ##             if (messages) 
    ##                 cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.", 
    ##                   "\n")
    ##         }
    ##     }
    ##     x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
    ##     traits <- round(x.pco$li, .Machine$double.exponent)
    ##     nb.sp <- numeric(c)
    ##     for (i in 1:c) {
    ##         sp.pres <- which(a[i, ] > 0)
    ##         traits.sp.pres <- traits[sp.pres, , drop = F]
    ##         traits.sp.pres[traits.sp.pres != 0 & abs(traits.sp.pres) < 
    ##             tol] <- 0
    ##         nb.sp[i] <- nrow(unique(traits.sp.pres))
    ##     }
    ##     names(nb.sp) <- row.names(a)
    ##     min.nb.sp <- min(nb.sp)
    ##     if (min.nb.sp < 3) 
    ##         if (messages) 
    ##             cat("FEVe: Could not be calculated for communities with <3 functionally singular species.", 
    ##                 "\n")
    ##     if (min.nb.sp < 2) 
    ##         if (messages) 
    ##             cat("FDis: Equals 0 in communities with only one functionally singular species.", 
    ##                 "\n")
    ##     if (calc.FRic) {
    ##         x.class2 <- sapply(x, data.class)
    ##         if (all(x.class2 == "factor" | x.class2 == "ordered")) {
    ##             if (length(x.class2) == 1 & x.class2[1] == "ordered") {
    ##                 traits.FRic1 <- rank(x[, 1])
    ##                 names(traits.FRic1) <- x.rn
    ##                 traits.FRic <- data.frame(traits.FRic1)
    ##                 qual.FRic = 1
    ##                 if (messages) 
    ##                   cat("FRic: Only one ordinal trait present in 'x'. FRic was measured as the range of the ranks, NOT as the convex hull volume.", 
    ##                     "\n")
    ##                 if (calc.FDiv) {
    ##                   calc.FDiv <- FALSE
    ##                   if (messages) 
    ##                     cat("FDiv: Cannot be computed when 'x' is a single ordinal trait.", 
    ##                       "\n")
    ##                 }
    ##                 if (stand.FRic) {
    ##                   traits.range <- range(traits.FRic[, 1])
    ##                   FRic.all <- traits.range[2] - traits.range[1]
    ##                 }
    ##             }
    ##             else {
    ##                 traits.FRic <- x
    ##                 qual.FRic = 1
    ##                 if (messages) 
    ##                   cat("FRic: Only categorical and/or ordinal trait(s) present in 'x'. FRic was measured as the number of unique trait combinations, NOT as the convex hull volume.", 
    ##                     "\n")
    ##                 if (stand.FRic) 
    ##                   FRic.all <- nrow((unique(traits.FRic)))
    ##                 if (calc.FDiv) {
    ##                   calc.FDiv <- FALSE
    ##                   if (messages) 
    ##                     cat("FDiv: Cannot be computed when only categorical and/or ordinal trait(s) present in 'x'.", 
    ##                       "\n")
    ##                 }
    ##             }
    ##         }
    ##         else {
    ##             if (x.pco$nf == 1) {
    ##                 traits.FRic <- x.pco$li
    ##                 qual.FRic = 1
    ##                 if (messages) 
    ##                   cat("FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume.", 
    ##                     "\n")
    ##                 if (calc.FDiv) {
    ##                   calc.FDiv <- FALSE
    ##                   if (messages) 
    ##                     cat("FDiv: Cannot not be computed when 'x' contains one single continuous trait or dimension.", 
    ##                       "\n")
    ##                 }
    ##                 if (stand.FRic) {
    ##                   traits.range <- range(traits.FRic[, 1])
    ##                   FRic.all <- traits.range[2] - traits.range[1]
    ##                 }
    ##             }
    ##             if (x.pco$nf > 1) {
    ##                 warning <- FALSE
    ##                 m.max <- min.nb.sp - 1
    ##                 if (m == "min") {
    ##                   warning <- TRUE
    ##                   if (min.nb.sp < 4) {
    ##                     nb.sp2 <- nb.sp[nb.sp > 3]
    ##                     m.min <- floor(log2(min(nb.sp2)))
    ##                     if (messages) 
    ##                       cat("FRic: To respect s >= 2^t, FRic could not be calculated for communities with <4 functionally singular species.", 
    ##                         "\n")
    ##                   }
    ##                   else m.min <- floor(log2(min.nb.sp))
    ##                 }
    ##                 else {
    ##                   if (min.nb.sp < 3) {
    ##                     nb.sp2 <- nb.sp[nb.sp > 2]
    ##                     m.max <- min(nb.sp2) - 1
    ##                     if (messages) 
    ##                       cat("FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species.", 
    ##                         "\n")
    ##                   }
    ##                   else m.max <- m.max
    ##                 }
    ##                 if (is.numeric(m) & m <= 1) 
    ##                   stop("When 'm' is an integer, it must be >1.", 
    ##                     "\n")
    ##                 if (is.numeric(m) & m > m.max) 
    ##                   m <- m.max
    ##                 if (m == "min") 
    ##                   m <- m.min
    ##                 if (m == "max") 
    ##                   m <- m.max
    ##                 if (!is.numeric(m) & m != "min" & m != "max") 
    ##                   stop("'m' must be an integer >1, 'min', or 'max'.", 
    ##                     "\n")
    ##                 if (m < x.pco$nf) {
    ##                   traits.FRic <- x.pco$li[, 1:m]
    ##                   if (x.pco$nf - m == 1) 
    ##                     if (messages) 
    ##                       cat("FRic: Dimensionality reduction was required. The last PCoA axis (out of", 
    ##                         x.pco$nf, "in total) was removed.", "\n")
    ##                   if (x.pco$nf - m > 1) 
    ##                     if (messages) 
    ##                       cat("FRic: Dimensionality reduction was required. The last", 
    ##                         x.pco$nf - m, "PCoA axes (out of", x.pco$nf, 
    ##                         "in total) were removed.", "\n")
    ##                   if (is.euclid(x.dist)) {
    ##                     qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
    ##                     if (messages) 
    ##                       cat("FRic: Quality of the reduced-space representation =", 
    ##                         qual.FRic, "\n")
    ##                   }
    ##                   if (!is.euclid(x.dist) & corr != "none") {
    ##                     qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
    ##                     if (messages) 
    ##                       cat("FRic: Quality of the reduced-space representation (based on corrected distance matrix) =", 
    ##                         qual.FRic, "\n")
    ##                   }
    ##                   if (!is.euclid(x.dist) & corr == "none") {
    ##                     delta <- -0.5 * bicenter.wt(x.dist * x.dist)
    ##                     lambda <- eigen(delta, symmetric = TRUE, 
    ##                       only.values = TRUE)$values
    ##                     sum.m <- sum(lambda[1:m])
    ##                     sum.n <- sum(lambda)
    ##                     lambda.neg <- c(lambda[lambda < 0])
    ##                     max.neg <- abs(min(lambda.neg))
    ##                     qual.FRic <- (sum.m + (length(lambda[1:m]) * 
    ##                       max.neg))/(sum.n + ((length(lambda) - 1) * 
    ##                       max.neg))
    ##                     if (messages) 
    ##                       cat("FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) =", 
    ##                         qual.FRic, "\n")
    ##                   }
    ##                 }
    ##                 if (m >= x.pco$nf) {
    ##                   qual.FRic = 1
    ##                   traits.FRic <- x.pco$li
    ##                   if (x.pco$nf == 2) 
    ##                     if (messages) 
    ##                       cat("FRic: No dimensionality reduction was required. The 2 PCoA axes were kept as 'traits'.", 
    ##                         "\n")
    ##                   if (x.pco$nf > 2) 
    ##                     if (messages) 
    ##                       cat("FRic: No dimensionality reduction was required. All", 
    ##                         x.pco$nf, "PCoA axes were kept as 'traits'.", 
    ##                         "\n")
    ##                 }
    ##                 if (stand.FRic) {
    ##                   hull.all <- convhulln(traits.FRic, "FA")
    ##                   FRic.all <- hull.all$vol
    ##                 }
    ##             }
    ##         }
    ##     }
    ##     if (!calc.FRic & calc.FDiv) 
    ##         cat("FDiv: Cannot be computed when 'calc.FRic' is FALSE.", 
    ##             "\n")
    ##     if (calc.FRic & calc.FDiv) 
    ##         if (min.nb.sp < 3) 
    ##             if (messages) 
    ##                 cat("FDiv: Could not be calculated for communities with <3 functionally singular species.", 
    ##                   "\n")
    ##     if (calc.FGR) {
    ##         if (clust.type == "kmeans") {
    ##             tr.clust <- cascadeKM(traits, km.inf.gr, km.sup.gr, 
    ##                 km.iter, km.crit)
    ##             cat("FGR: Summary of kmeans clustering\n")
    ##             cat("\nPartition\n")
    ##             print(tr.clust$partition)
    ##             cat("\nResults\n")
    ##             print(tr.clust$results)
    ##             cat("\nSize\n")
    ##             print(tr.clust$size)
    ##             plot(tr.clust)
    ##             part.names <- colnames(tr.clust$partition)
    ##             part.names <- as.numeric(substr(part.names, 1, 1))
    ##             cat("\nFGR: How many groups?", "\n")
    ##             cut.g <- toupper(scan(file = "", what = "character", 
    ##                 nlines = 1, quiet = T))
    ##             cut.gr <- as.integer(cut.g)
    ##             if (cut.gr < km.inf.gr | cut.gr > km.sup.gr) 
    ##                 stop("You must type an integer between 'km.ing.gr' and 'km.sup.gr'.", 
    ##                   "\n")
    ##             spfgr.all <- tr.clust$partition[, part.names == cut.gr]
    ##             names(spfgr.all) <- x.rn
    ##         }
    ##         else {
    ##             tr.clust <- hclust(x.dist, method = clust.type)
    ##             plot(tr.clust, main = "Cluster dengrogram of species based on functional traits")
    ##             cat("FGR: Do you want to cut the dendrogram from height or from the number of groups? Type 'h' for height, 'g' for groups.", 
    ##                 "\n")
    ##             cut <- toupper(scan(file = "", what = "character", 
    ##                 nlines = 1, quiet = T))
    ##             if (cut == "H") {
    ##                 cat("FGR: At what height do you want the dendrogram to be cut?", 
    ##                   "\n")
    ##                 cut.d <- toupper(scan(file = "", what = "character", 
    ##                   nlines = 1, quiet = T))
    ##                 cut.dist <- as.numeric(cut.d)
    ##                 spfgr.all <- cutree(tr.clust, h = cut.dist)
    ##             }
    ##             if (cut == "G") {
    ##                 cat("FGR: How many groups?", "\n")
    ##                 cut.g <- toupper(scan(file = "", what = "character", 
    ##                   nlines = 1, quiet = T))
    ##                 cut.gr <- as.integer(cut.g)
    ##                 spfgr.all <- cutree(tr.clust, k = cut.gr)
    ##             }
    ##             if (cut != "H" & cut != "G") 
    ##                 stop("You must type 'h' or 'g'", "\n")
    ##         }
    ##         a.t <- t(a)
    ##         by.gr <- list(spfgr.all)
    ##         gr.abun <- aggregate(a.t, by.gr, sum)
    ##         lab <- paste("group", gr.abun[, 1], sep = "")
    ##         gr.abun <- data.frame(t(gr.abun[, -1]))
    ##         colnames(gr.abun) <- lab
    ##         rownames(gr.abun) <- rownames(a)
    ##     }
    ##     if (is.matrix(x) | is.data.frame(x) & calc.CWM) {
    ##         CWM <- functcomp(x, a, CWM.type = CWM.type)
    ##     }
    ##     if (calc.CWM & class(x)[1] == "dist" | class(x)[1] == "dissimilarity") 
    ##         if (messages) 
    ##             cat("CWM: When 'x' is a distance matrix, CWM cannot be calculated.", 
    ##                 "\n")
    ##     divc <- function(df, dis = NULL, scale = FALSE) {
    ##         if (!inherits(df, "data.frame")) 
    ##             stop("Non convenient df")
    ##         if (any(df < 0)) 
    ##             stop("Negative value in df")
    ##         if (!is.null(dis)) {
    ##             if (!inherits(dis, "dist")) 
    ##                 stop("Object of class 'dist' expected for distance")
    ##             dis <- as.matrix(dis)
    ##             if (nrow(df) != nrow(dis)) 
    ##                 stop("Non convenient df")
    ##             dis <- as.dist(dis)
    ##         }
    ##         if (is.null(dis)) 
    ##             dis <- as.dist((matrix(1, nrow(df), nrow(df)) - diag(rep(1, 
    ##                 nrow(df)))) * sqrt(2))
    ##         div <- as.data.frame(rep(0, ncol(df)))
    ##         names(div) <- "diversity"
    ##         rownames(div) <- names(df)
    ##         for (i in 1:ncol(df)) {
    ##             if (sum(df[, i]) < 1e-16) 
    ##                 div[i, ] <- 0
    ##             else div[i, ] <- (t(df[, i]) %*% (as.matrix(dis)^2) %*% 
    ##                 df[, i])/2/(sum(df[, i])^2)
    ##         }
    ##         if (scale == TRUE) {
    ##             divmax <- divcmax(dis)$value
    ##             div <- div/divmax
    ##         }
    ##         return(div)
    ##     }
    ##     RaoQ <- divc(data.frame(t(a)), x.dist, scale = scale.RaoQ)
    ##     RaoQ <- RaoQ[, 1]
    ##     names(RaoQ) <- rownames(a)
    ##     disp <- fdisp(x.dist, a)
    ##     FDis <- disp$FDis
    ##     nbsp <- rep(NA, c)
    ##     names(nbsp) <- row.names(a)
    ##     FRic <- rep(NA, c)
    ##     names(FRic) <- row.names(a)
    ##     FEve <- rep(NA, c)
    ##     names(FEve) <- row.names(a)
    ##     FGR <- rep(NA, c)
    ##     names(FGR) <- row.names(a)
    ##     FDiv <- rep(NA, c)
    ##     names(FDiv) <- row.names(a)
    ##     for (i in 1:c) {
    ##         sppres <- which(a[i, ] > 0)
    ##         S <- length(sppres)
    ##         nbsp[i] <- S
    ##         tr <- data.frame(traits[sppres, ])
    ##         if (calc.FRic) 
    ##             tr.FRic <- data.frame(traits.FRic[sppres, ])
    ##         ab <- as.matrix(a[i, sppres])
    ##         abundrel <- ab/sum(ab)
    ##         if (calc.FRic) {
    ##             if (all(x.class2 == "factor" | x.class2 == "ordered")) {
    ##                 if (length(x.class2) == 1 & x.class2[1] == "ordered") {
    ##                   tr.range <- range(tr.FRic[, 1])
    ##                   t.range <- tr.range[2] - tr.range[1]
    ##                   if (!stand.FRic) 
    ##                     FRic[i] <- t.range
    ##                   if (stand.FRic) 
    ##                     FRic[i] <- t.range/FRic.all
    ##                 }
    ##                 else {
    ##                   if (!stand.FRic) 
    ##                     FRic[i] <- nrow((unique(tr.FRic)))
    ##                   if (stand.FRic) 
    ##                     FRic[i] <- nrow((unique(tr.FRic)))/FRic.all
    ##                 }
    ##             }
    ##             else {
    ##                 if (dim(tr.FRic)[2] > 1 & nb.sp[i] >= 3) {
    ##                   if (warning) 
    ##                     thresh <- 4
    ##                   if (!warning) 
    ##                     thresh <- 3
    ##                   if (nb.sp[i] >= thresh) {
    ##                     convhull <- convhulln(tr.FRic, "FA")
    ##                     if (!stand.FRic) 
    ##                       FRic[i] <- convhull$vol
    ##                     if (stand.FRic) 
    ##                       FRic[i] <- convhull$vol/FRic.all
    ##                   }
    ##                   else {
    ##                   }
    ##                 }
    ##                 if (dim(tr.FRic)[2] == 1) {
    ##                   tr.range <- range(tr.FRic[, 1])
    ##                   t.range <- tr.range[2] - tr.range[1]
    ##                   if (!stand.FRic) 
    ##                     FRic[i] <- t.range
    ##                   if (stand.FRic) 
    ##                     FRic[i] <- t.range/FRic.all
    ##                 }
    ##             }
    ##         }
    ##         if (nb.sp[i] >= 3) {
    ##             tr.dist <- dist(tr)
    ##             linkmst <- mst(tr.dist)
    ##             mstvect <- as.dist(linkmst)
    ##             abund2 <- matrix(0, nrow = S, ncol = S)
    ##             for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] + 
    ##                 abundrel[r]
    ##             abund2vect <- as.dist(abund2)
    ##             EW <- rep(0, S - 1)
    ##             flag <- 1
    ##             for (m in 1:((S - 1) * S/2)) {
    ##                 if (mstvect[m] != 0) {
    ##                   EW[flag] <- tr.dist[m]/(abund2vect[m])
    ##                   flag <- flag + 1
    ##                 }
    ##             }
    ##             minPEW <- rep(0, S - 1)
    ##             OdSmO <- 1/(S - 1)
    ##             for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), 
    ##                 OdSmO)
    ##             FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    ##         }
    ##         if (calc.FDiv & calc.FRic) {
    ##             if (any(x.class2 == "numeric") & dim(tr.FRic)[2] > 
    ##                 1 & nb.sp[i] >= 3) {
    ##                 vert0 <- convhulln(tr.FRic, "Fx TO 'vert.txt'")
    ##                 vert1 <- scan("vert.txt", quiet = T)
    ##                 vert2 <- vert1 + 1
    ##                 vertices <- vert2[-1]
    ##                 trvertices <- tr.FRic[vertices, ]
    ##                 baryv <- apply(trvertices, 2, mean)
    ##                 distbaryv <- rep(0, S)
    ##                 for (j in 1:S) distbaryv[j] <- (sum((tr.FRic[j, 
    ##                   ] - baryv)^2))^0.5
    ##                 meandB <- mean(distbaryv)
    ##                 devdB <- distbaryv - meandB
    ##                 abdev2 <- abundrel * devdB
    ##                 ababsdev2 <- abundrel * abs(devdB)
    ##                 FDiv[i] <- (sum(abdev2) + meandB)/(sum(ababsdev2) + 
    ##                   meandB)
    ##             }
    ##         }
    ##         if (calc.FGR) 
    ##             FGR[i] <- length(unique(spfgr.all[sppres]))
    ##     }
    ##     res <- list()
    ##     res$nbsp <- nbsp
    ##     res$sing.sp <- nb.sp
    ##     if (calc.FRic) 
    ##         res$FRic <- FRic
    ##     if (calc.FRic) 
    ##         res$qual.FRic <- qual.FRic
    ##     res$FEve <- FEve
    ##     if (calc.FDiv) 
    ##         res$FDiv <- FDiv
    ##     res$FDis <- FDis
    ##     res$RaoQ <- RaoQ
    ##     if (calc.FGR) {
    ##         res$FGR <- FGR
    ##         res$spfgr <- spfgr.all
    ##         res$gr.abun <- gr.abun
    ##     }
    ##     if (is.matrix(x) | is.data.frame(x) & calc.CWM) 
    ##         res$CWM <- CWM
    ##     if (print.pco) {
    ##         res$x.values <- x.pco$eig
    ##         res$x.axes <- x.pco$li
    ##     }
    ##     invisible(res)
    ## }
    ## <bytecode: 0x7fa324b90708>
    ## <environment: namespace:FD>

## `rao.diversity`

``` r
library(SYNCSA)

rao.diversity
```

    ## function (comm, traits = NULL, phylodist = NULL, checkdata = TRUE, 
    ##     ord = "metric", put.together = NULL, standardize = TRUE, 
    ##     ...) 
    ## {
    ##     diver.internal <- function(community, distance) {
    ##         if (any(is.na(distance))) {
    ##             distance.na <- ifelse(is.na(distance), 0, 1)
    ##             inter.na <- community %*% distance.na
    ##             adjustment <- rowSums(sweep(community, 1, inter.na, 
    ##                 "*", check.margin = FALSE))
    ##             distance[is.na(distance)] <- 0
    ##             inter <- community %*% distance
    ##             res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
    ##             res <- ifelse(adjustment > 0, res/adjustment, res)
    ##         }
    ##         else {
    ##             inter <- community %*% distance
    ##             res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
    ##         }
    ##         return(res)
    ##     }
    ##     res <- list(call = match.call())
    ##     if (inherits(comm, "metacommunity.data")) {
    ##         if (!is.null(traits) | !is.null(phylodist) | !is.null(put.together)) {
    ##             stop("\n When you use an object of class metacommunity.data the arguments traits, phylodist and put.together must be null. \n")
    ##         }
    ##         traits <- comm$traits
    ##         phylodist <- comm$phylodist
    ##         put.together <- comm$put.together
    ##         comm <- comm$community
    ##     }
    ##     list.warning <- list()
    ##     if (checkdata) {
    ##         organize.temp <- organize.syncsa(comm, traits = traits, 
    ##             phylodist = phylodist, check.comm = TRUE)
    ##         if (!is.null(organize.temp$stop)) {
    ##             organize.temp$call <- match.call()
    ##             return(organize.temp)
    ##         }
    ##         list.warning <- organize.temp$list.warning
    ##         comm <- organize.temp$community
    ##         traits <- organize.temp$traits
    ##         phylodist <- organize.temp$phylodist
    ##     }
    ##     if (length(list.warning) > 0) {
    ##         res$list.warning <- list.warning
    ##     }
    ##     if (any(is.na(comm))) {
    ##         stop("\n community data with NA\n")
    ##     }
    ##     comm <- as.matrix(comm)
    ##     N <- nrow(comm)
    ##     S <- ncol(comm)
    ##     dist.1 <- 1 - diag(x = rep(1, S))
    ##     if (!is.null(traits)) {
    ##         traits <- as.data.frame(traits)
    ##         m <- ncol(traits)
    ##         weights <- rep(1, m)
    ##         make.names <- is.null(colnames(traits))
    ##         colnames(traits) <- colnames(traits, do.NULL = FALSE, 
    ##             prefix = "T")
    ##         names(weights) <- colnames(traits)
    ##         if (!is.null(put.together)) {
    ##             if (!inherits(put.together, "list")) {
    ##                 stop("\n put.together must be a object of class list\n")
    ##             }
    ##             if (make.names) {
    ##                 for (k in 1:length(put.together)) {
    ##                   put.together[[k]] <- paste("T", put.together[[k]], 
    ##                     sep = "")
    ##                 }
    ##             }
    ##             if (max(table(unlist(put.together))) > 1) {
    ##                 stop("\n The same trait appears more than once in put.together\n")
    ##             }
    ##             if (length(setdiff(unlist(put.together), colnames(traits))) > 
    ##                 0) {
    ##                 stop("\n Check traits names in put.together\n")
    ##             }
    ##             for (k in 1:length(put.together)) {
    ##                 weights[put.together[[k]]] <- 1/length(put.together[[k]])
    ##             }
    ##         }
    ##         dist.functional <- sqrt(as.matrix(FD::gowdis(x = traits, 
    ##             asym.bin = NULL, ord = ord, w = weights, ...)))
    ##         if (checkdata) {
    ##             if (any(is.na(dist.functional))) {
    ##                 warning("Warning: NA in distance between species", 
    ##                   call. = FALSE)
    ##             }
    ##         }
    ##     }
    ##     if (!is.null(phylodist)) {
    ##         dist.phylogenetic <- as.matrix(phylodist)
    ##         if (checkdata) {
    ##             if (any(is.na(dist.phylogenetic))) {
    ##                 warning("Warning: NA in phylodist", call. = FALSE)
    ##             }
    ##         }
    ##         if (standardize) {
    ##             dist.phylogenetic <- dist.phylogenetic/max(dist.phylogenetic, 
    ##                 na.rm = TRUE)
    ##         }
    ##     }
    ##     comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
    ##     SD <- diver.internal(comm, dist.1)
    ##     res$Simpson <- SD
    ##     if (!is.null(traits)) {
    ##         FD <- diver.internal(comm, dist.functional)
    ##         res$FunRao <- FD
    ##         res$FunRedundancy <- SD - FD
    ##     }
    ##     if (!is.null(phylodist)) {
    ##         PD <- diver.internal(comm, dist.phylogenetic)
    ##         res$PhyRao <- PD
    ##         res$PhyRedundancy <- SD - PD
    ##     }
    ##     return(res)
    ## }
    ## <bytecode: 0x7fa326a0c5f0>
    ## <environment: namespace:SYNCSA>

## `raoD`

``` r
library(picante)
```

    ## Warning: package 'picante' was built under R version 3.5.2

``` r
raoD
```

    ## function (comm, phy = NULL) 
    ## {
    ##     res <- list()
    ##     if (is.null(phy)) {
    ##         tij <- 1 - diag(x = rep(1, length(comm[1, ])))
    ##     }
    ##     else {
    ##         if (!is.ultrametric(phy)) {
    ##             stop("Phylogeny must be ultrametric")
    ##         }
    ##         dat <- match.phylo.comm(phy, comm)
    ##         comm <- dat$comm
    ##         phy <- dat$phy
    ##         tij <- cophenetic(phy)/2
    ##     }
    ##     x <- as.matrix(comm)
    ##     S <- length(x[1, ])
    ##     N <- length(x[, 1])
    ##     total <- apply(x, 1, sum)
    ##     samp.relabund <- total/sum(x)
    ##     x.combined <- matrix(apply(x, 2, sum), nrow = 1)/sum(x)
    ##     x <- sweep(x, 1, total, "/")
    ##     D <- vector(length = N)
    ##     names(D) <- rownames(x)
    ##     for (k in 1:N) D[k] <- sum(tij * outer(as.vector(t(x[k, ])), 
    ##         as.vector(t(x[k, ]))))
    ##     res$Dkk <- D
    ##     Dkl <- matrix(nrow = N, ncol = N)
    ##     for (k in 1:N) {
    ##         for (l in 1:N) {
    ##             Dkl[k, l] <- sum(tij * outer(as.vector(t(x[k, ])), 
    ##                 as.vector(t(x[l, ]))))
    ##         }
    ##     }
    ##     row.names(Dkl) <- row.names(x)
    ##     colnames(Dkl) <- row.names(x)
    ##     H <- Dkl
    ##     res$Dkl <- Dkl
    ##     for (k in 1:N) {
    ##         for (l in 1:N) {
    ##             H[k, l] <- Dkl[k, l] - (Dkl[k, k] + Dkl[l, l])/2
    ##         }
    ##     }
    ##     res$H <- H
    ##     res$total <- sum(tij * outer(as.vector(t(x.combined)), as.vector(t(x.combined))))
    ##     res$alpha <- sum(res$Dkk * samp.relabund)
    ##     res$beta <- res$total - res$alpha
    ##     res$Fst <- res$beta/res$total
    ##     return(res)
    ## }
    ## <bytecode: 0x7fa324e2b348>
    ## <environment: namespace:picante>

# References

<div id="refs" class="references">

<div id="ref-botta2005rao">

Botta-Dukát, Z. (2005). Rao’s quadratic entropy as a measure of
functional diversity based on multiple traits. *Journal of vegetation
science*, 16, 533–540.

</div>

<div id="ref-Champely_Chessel_2002">

Champely, S. & Chessel, D. (2002). Measuring biological diversity using
euclidean metrics. *Environmental and Ecological Statistics*, 9,
167–177.

</div>

<div id="ref-deBello2016rao">

de Bello, F., Carmona, C.P., Lepš, J., Szava-Kovats, R. & Pärtel, M.
(2016). Functional diversity through the mean trait dissimilarity:
Resolving shortcomings with existing paradigms and algorithms.
*Oecologia*, 180, 933–940.

</div>

<div id="ref-deBello2007">

de Bello, F., Lepš, J., Lavorel, S. & Moretti, M. (2007). Importance of
species abundance for assessment of trait composition: An example based
on pollinator communities. *Community Ecology*, 8, 163–170.

</div>

<div id="ref-GALLAND2020106488">

Galland, T., Carmona\], C. \[Pérez, Götzenberger, L., Valencia, E. &
Bello\], F. \[de. (2020). Are redundancy indices redundant? An
evaluation based on parameterized simulations. *Ecological Indicators*,
116, 106488.

</div>

<div id="ref-Legendre_Legendre2012">

Legendre, P. & Legendre, L. (2012). *Numerical ecology, 3rd english
edition*. Elsevier.

</div>

<div id="ref-leps2006">

Lepš, J., de Bello, F., Lavorel, S. & Berman, S. (2006). Quantifying and
interpreting functional diversity of natural communities: Practical
considerations matter. *Preslia*, 78, 481–501.

</div>

<div id="ref-pavoine2005phd">

Pavoine, S. (2005). Méthodes statistiques pour la mesure de la
biodiversité. PhD thesis. Lyon 1.

</div>

<div id="ref-pavoine2012">

Pavoine, S. (2012). Clarifying and developing analyses of biodiversity:
Towards a generalisation of current approaches. *Methods in Ecology and
Evolution*, 3, 509–518.

</div>

<div id="ref-pavoine_bonsall2009">

Pavoine, S. & Bonsall, M.B. (2009). Biological diversity: Distinct
distributions can lead to the maximization of rao’s quadratic entropy.
*Theoretical population biology*, 75, 153–163.

</div>

<div id="ref-pavoine2005_rao_dissim">

Pavoine, S., Ollier, S. & Pontier, D. (2005). Measuring diversity from
dissimilarities with rao’s quadratic entropy: Are any dissimilarities
suitable? *Theoretical population biology*, 67, 231–239.

</div>

<div id="ref-pla2011">

Pla, L., Casanoves, F. & Di Rienzo, J. (2011). *Quantifying functional
biodiversity*. Springer Science.

</div>

<div id="ref-rao1982">

Rao, C.R. (1982). Diversity and dissimilarity coefficients: A unified
approach. *Theoretical population biology*, 21, 24–43.

</div>

<div id="ref-rao2010">

Rao, C.R. (2010). Quadratic entropy and analysis of diversity. *Sankhya
A*, 72, 70–80.

</div>

<div id="ref-ricotta2005">

Ricotta, C. (2005). A note on functional diversity measures. *Basic and
Applied Ecology*, 6, 479–486.

</div>

<div id="ref-ricotta2016">

Ricotta, C., de Bello, F., Moretti, M., Caccianiga, M., Cerabolini, B.E.
& Pavoine, S. (2016). Measuring the functional redundancy of biological
communities: A quantitative guide. *Methods in Ecology and Evolution*,
7, 1386–1395.

</div>

<div id="ref-schleuter2010">

Schleuter, D., Daufresne, M., Massol, F. & Argillier, C. (2010). A
user’s guide to functional diversity indices. *Ecological Monographs*,
80, 469–484.

</div>

<div id="ref-warwick1995">

Warwick, R.M. & Clarke, K. (1995). New ’biodiversity’ measures reveal a
decrease in taxonomic distinctness with increasing stress. *Marine
Ecology Progress Series*, 301–305.

</div>

</div>
