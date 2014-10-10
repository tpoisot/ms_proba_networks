% The structure of probabilistic networks
% The Magical Mistery Stouffer Group
% Working paper -- Oct. 2014

# Introduction

Ecological networks are an efficient way to represent the interactions between
individual, populations, or species. Historically, their study has focused on
(i) linking their structure to community or ecosystem-level properties such as
stability [@mcca14], the maintenance of species richness [@bast09;@haer14],
ecosystem functioning [@theb03;@duff02], and (ii) describing the overall
structure of networks, with a particular attention on food webs [@dunn06]
and plant-pollinator interactions [@basc03;@jord87]. To a large extent,
the description of network structure enabled questions about how it ties
into functional properties, and it is no surprise that the methodology to
describe networks is large.

Most measures of network structure function in the following way. Given a
network as input, they return a *property* based on one or several *units*
within this network. Some of the properties are *direct* properties (they only
require knowledge of the unit on which they are applied), and some others are
*emerging* properties (they require knowledge of higher-order structures). For
example, connectance, the proportion of realized interactions, is a direct
property of a network. The degree of a node (how many interactions it is
involved in) is a direct property of the node, whereas the degree distribution
is an emerging property of all nodes. Establishing a difference between
direct and emerging properties is important when interpreting their values:
direct properties are conceptually equivalent to means, whereas emerging
properties are conceptually equivalent to variances.

In the recent years, the interpretation of the values of network structure
(as indicators of the action of ecological or evolutionary processes) has
been somewhat complicated by the observation that network structure varies
through space, and time; species from the same pool do not interact in a
consistent way [@pois12c]. Empirical and theoretical studies suggest that the
network is not the right unit to understand this variation; rather, network
variation is an emerging property of the response of ecological interactions
to environmental factors and chance events [@pois14]. Interactions can
vary because of local mis-matching in phenology [@oles11a], populations
fluctuations preventing the interaction [@cana14], or a combination of both
[@olit14;@cham14]. @olit14 show that accounting for neutral (population-size
driven) and trait-based effects allows predicting the cumulative change in
network structure, but not the change at the level of individual interactions.

Taken together, these considerations highlight the need to amend our current
methodology on ecological network to give more importance to the variation
at the interaction level. Because the methodology to describe networks has
first been crafted at a time when assuming that interactions did not vary,
it is unsuited to address the questions that probabilistic networks allows
asking. In this paper, we show that several direct and emerging core properties
of ecological networks (both bipartite and unipartite) can be re-formulated
in a probabilistic context; we conclude by showing how this methodology can be
applied to exploit the information contained in the variability and networks,
and reduce the computational burden of current methods in network analysis.

# Metrics

Throughout this section, we will assume the following notation. $\mathbf{A}$
is a matrix wherein $A_{ij}$ is $\text{P}(ij)$, *i.e.* the probability that
species $i$ establishes an interaction with species $j$. If $\mathbf{A}$
represents a unipartite network (*e.g.* a food web), it is a square matrix and
the probabilities of each species interacting with itself. If $\mathbf{A}$
represents a bipartite network (*e.g.* a pollination network), it will
most likely not be square. We call $S$ the number of species, and $R$ and
$C$ respectively the number of rows and columns. $S = R + C$ in unipartite
networks, and $S = R+C$ in bipartite networks.

Note that all of the measures defined below can be applied on a bipartite
network that has been made unipartite; the unipartite transformation of a
bipartite matrix $\mathbf{A}$ is the block matrix

\begin{equation}
\mathbf{B} = 
\begin{pmatrix}
0_{(R,R)} & \mathbf{A}\\
0_{(C,R)} & 0_{(C,C)}
\end{pmatrix},
\end{equation}

where $0_{(C,R)}$ is a matrix of $C$ rows and $R$ columns filled with $0$s,
etc.

We assume that all interactions are independent (so that $\text{P}(ij|kl)
= \text{P}(ij)\text{P}(kl)$ for any species), and can be represented as
Bernoulli trials (so that $0 \leq \text{P}(ij) \leq 1$). The later condition
allows to derive estimates for the *variance* of the measures, since (i)
the variance of a single event $X_i$ of probability $p$ is $\text{var}(X)
= p(1-p)$, its expected value is $\text{E}(X)=p$, (ii) the variance of
additive independent events is the sum of their individual variances, and
(iii) the variance of multiplicative independent events is

\begin{equation}
\text{var}(X_1 X_2 ... X_n) = \prod_i \left(\text{var}(X_i) + [\text{E}(X_i)]^2\right) - \prod_i [\text{E}(X_i)]^2
\end{equation}

As a final note, all of the measures described below can be applied on the
binary (0/1) versions of the networks, and will give the exact value of the
non-probabilistic measure. And ain't that nice?

## Direct properties

### Connectance and number of interactions

Connectance is the proportion of realized upon possible interactions, defined
as $Co = L/(R\times C)$, where $L$ is the total number of interactions. As
all interactions in a probabilistic network are assumed to be indpendent,
the expected value of $L$, is

\begin{equation}
\hat L = \sum A_{ij}, 
\end{equation}

and ${\hat{Co}} = \hat L / (R\times C)$.

The variance of the number of interactions is $\text{var}(\hat L) = \sum
(A_{ij}(1-A_{ij}))$.

### Node degree

The degree distribution of a network is the distribution of the number of
interactions established and received by each node. The expected degree of
species $i$ is

\begin{equation}
\hat k_i = \sum_j(A_{ij} + A_{ji})
\end{equation}

The variance of the degree of each species is $\text{var}(\hat k_i) =
\sum_j(A_{ij}(1-A_{ij})+A_{ji}(1-A_{ji}))$. Note also that as expected,
$\sum \hat k_i = 2\hat L$.

### Average generality and vulnerability

By simplification of the above, generality $\hat g_i$ and vulnerability
$\hat v_i$ are given by, respectively, $\sum_j A_{ij}$ and $\sum_j A_{ji}$,
with their variances $\sum_j A_{ij}(1-A_{ij})$ and $\sum_j A_{ji}(1-A_{ji})$.

## Emerging properties

### Path length

Networks can be used to describe indirect interactions between species,
through the use of paths. The existence of a path of length 2 between species
$i$ and $j$ mean that they are connected through at least one additional
species $k$. In a probabilistic network, unless some elements are $0$,
all pairs of species $i$ and $j$ are connected through a path of length
1, with probability $A_{ij}$. The expected number of paths of length $k$
between species $i$ and $j$ is given by

\begin{equation}
\hat{n^{(2)}_{ij}} = \left(\mathbf{A}^k\right)_{ij},
\end{equation}

where $\mathbf{A}^k$ is the matrix multiplied by itself $k$ times.

It is possible to calculate the probability of having at least one path
between the two species: this can be done by calculating the probability of
having 0 paths, then multiplying the resulting array of probabilities. For
the example of length 2, species $i$ and $j$ are connected through $k$ with
probability $A_{ik}A_{kj}$, and so this path does not exist with probability
$1-A_{ik}A_{kj}$. For any pair $i$, $j$, let $\mathbf{m}$ be the vector
such as $m_{k} = A_{ik}A_{kj}$ for all $k \notin (i,j)$. The probability
of not having any path of length 2 is $\prod (1-\mathbf{m})$. Therefore,
the probability of having a path of length 2 between $i$ and $j$ is

\begin{equation}
\hat{p}^{(2)}_{ij} = 1 - \prod (1-\mathbf{m}) .
\end{equation}

In most situations, one would be interested in knowing the probability of
having a path of length 2 *without* having a path of length 1; this is simply
expressed as $(1-A_{ij})\hat{p}^{(2)}_{ij}$. One can, by the same logic,
generate the expression for having at least one path of length 3:

\begin{equation}
\hat{p}^{(3)}_{ij} = (1-A_{ij})(1-\hat{p}^{(2)}_{ij})\left(1 - \prod (1-\mathbf{m})\right)\prod_{x,y}\left((1-A_{iy})(1-A_{xj})\right),
\end{equation}

where $\mathbf{m}$ is the vector of all $A_{ix}A_{xy}A_{yj}$ for $x\notin
(i,j), y\neq x$. This gives the probability of having at least one path from
$i$ to $j$, passing through any pair of nodes $x$ and $j$, without having any
shorter path. In theory, this approach can be generalized up to an arbitrary
path length, but it becomes rapidly untractable.

### Nestedness

We use the formula for nestedness proposed by @bast09. They define
nestedness for each margin of the matrix, as $\eta^{(R)}$ and $\eta^{(C)}$
for, respectively, rows and columns. As per @alme08, we define a global
statistic for nestedness as $\eta = (\eta^{(R)}+\eta^{(C)})/2$.

Nestedness, in a probabilistic network, is defined as

\begin{equation}
\hat{\eta^{(R)}} = \sum_{i<j}\frac{\sum_kA_{ik}A_{jk}}{\text{min}(g_i, g_j)},
\end{equation}

where $g_i$ is the expected generality of species $i$. The reciprocal holds
for $\eta^{(C)}$ when using $v_i$ (the vulnerability) instead of $g_i$.

The values returned are within $[0;1]$, with $\eta=1$ indicating complete
nestedness.

### Katz centrality

Although a rough estimate of centrality is the node degree, as described
above, it is often needed to measure centrality within the context of a larger
neighborhood. In addition, we derive the expected value of centrality according
to @katz53. This measures generalizes to directed acyclic graphs. Although
eigenvector centrality is often used in ecology, it cannot be measured on
probabilistic graphs. Eigenvector centrality requires that the matrix has
its largest eigenvalues real, which is not the case for *all* probabilistic
matrices. Katz's centrality is nonetheless a useful replacement, because
it uses the paths of all lengths between two species instead of focusing on
the shortest path.

The expected number of paths of length $k$ between $i$ and $j$ is
$(\mathbf{A}^k)_{ij}$. Based on this, the expected centrality of species $i$ is

\begin{equation}
C_i = \sum_{k=1}^\infty \sum_{j=1}^n \alpha^k (\mathbf{A}^k)_{ji} .
\end{equation}

The parameter $\alpha \in [0;1]$ regulates how important long paths are. When
$\alpha = 0$, only first-order paths count. When $\alpha = 1$, all paths
are equally important. As $C_i$ is sensitive to the size of the matrix,
we suggest to normalise it so that

\begin{equation}
C_i = \frac{C_i}{\mathbf{C}} .
\end{equation}

This results in the *expected relative centrality* of each node in the
probabilistic network. Note that when using only $k = 1$, and $\alpha = 1$,
the raw value of Katz's centrality is the species generality.

### Number of primary producers

Primary producers, in a food web, are species with no successors, including
themselves. Biologically, they are autotrophic organisms, or organisms
whose preys or substrates have been remove from the network. A species is a
primary producer if it manages *not* to establish any outgoing interaction,
which for species $i$ happens with probability

\begin{equation}
\prod_j (1-A_{ij}).
\end{equation}

The number of expected primary producers is therefore the sum of the above
across all species:

\begin{equation}
\hat{PP} = \sum_i \left(\prod_j (1-A_{ij})\right).
\end{equation}

The variance in the number of expected primary producers is

\begin{equation}
\text{var}(\hat{PP}) = \sum_i \left( \prod_j(1-A_{ij}^2) - \prod_j(1-A_{ij})^2 \right)
\end{equation}

### Number of top predators

Top-predators can loosely be defined as species that have no predecessors in
the network: they are establishing links with other species, but no species
are establishing links with them. Using the same approach than for the number
of primary producers, the expected number of top-predators is therefore

\begin{equation}
\hat{TP} = \sum_i\left(\prod_{j \neq i}(1-A_{ji})\right)
\end{equation}

Note that we exclude the self-interactions, as top-predators can, and often
do, engage in cannibalism.

### Number of species with no interactions

Predicting the number of species with no interaction (or whether any species
will have at least one interaction) is useful to predict whether species
will be able to integrate themselves in an existing network, for example.

A species has no interactions with probability

\begin{equation}
\prod_{j \neq i} (1-A_{ij})(1-A_{ji})
\end{equation}

As for the above, the expected number of species with no interactions
(*free species*) is the sum of this quantity across all $i$:

\begin{equation}
\hat{FS} = \sum_i\prod_{j \neq i} (1-A_{ij})(1-A_{ji})
\end{equation}

The variance of the number of species with no interactions is

\begin{equation}
\text{var}(\hat{FS}) = \sum_i \left(
A_{ij}(1-A_{ij})A_{ji}(1-A_{ji})+A_{ij}(1-A_{ij})A_{ji}^2+A_{ji}(1-A_{ji})A_{ij}^2
\right)
\end{equation}

Note that from a methodological point of view, this can be a helpful *a
priori* measure to determine whether null models of networks will have a
lot of species with no interactions, and so will require intensive sampling.

### Self-predation

Self-predation (the existence of an interaction of a species onto
itself) is only meaningful in unipartite networks. The expected
proportion of species with self-loops is very simply defined as
$\text{Tr}(\mathbf{A})$, that is, the sum of all diagonal elements. The
variance is $\text{Tr}(\mathbf{A}\diamond(1-\mathbf{A}))$, where $\diamond$
is the element-wise product operation.

### Motifs

Motifs are sets of pre-determined interactions between a fixed number of specie
[@milo02], such as for example one predator sharing two preys. As there is
an arbitrarily large number of motifs, we will illustrate the formulae with
only two examples.

The probability that three species form an apparent competition motif (one
predator, two preys) where $i$ is the predator, $j$ and $k$ are the preys, is

\begin{equation}
\text{P}(i,j,k\in\text{app. comp}) = A_{ij}(1-A_{ji})A_{ik}(1-A_{ki})(1-A_{jk})(1-A_{kj})
\end{equation}

Similarly, the probability that these three species form an omnivory motif,
in which $i$ and $j$ consume $k$, and $i$ consumes $j$, is

\begin{equation}
\text{P}(i,j,k\in\text{omniv.}) = A_{ij}(1-A_{ji})A_{ik}(1-A_{ki})A_{jk}(1-A_{kj})
\end{equation}

The probability of the number of *any* motif $\text{m}$ in a network is given by

\begin{equation}
\hat{N_\text{m}} = \sum_i \sum_{j\neq i} \sum_{k\neq j} P(i,j,k \in \text{m})
\end{equation}

It is indeed possible to have an expression of the variance of this value,
or of the variance of any three species forming a given motif, but their
expressions become rapidly untractable and are better computer than written.

# Applications

In this section, we will provide an overview of the applications
of probabilistic network measures. The current way of dealing with
probabilistic interactions is (i) to ignore it entirely or (ii) to generate
random networks. Probabilistic metrics are an alternative to that. When
ignoring the probabilistic nature of interactions, what we call *Binary*
from here on, every non-zero element of the network is assumed to be 1. This
leads to over-representation of some rare-events, and increases the number
of interactions.

When generating random networks, what we call *Bernoulli trials* from here on,
a binary network is generated by doing a Bernoulli trial with probability
$A_{ij}$, for each element of the matrix. This is problematic because (i)
higher order structures involving rare events will be under-represented in
the sample, and (ii) naive approaches are likely to generate free species,
especially in sparsely connected networks frequently encountered in ecology
[@milo03;@pois14a].

## Comparison of probabilistic networks

In this sub-section, we apply the above measures to a bacteria--phage
interaction network. @poul08 have measured the probability that 24 phages
can infect 24 strains of bacteria of the *Pseudomonas fluorescens* species
(group SBW25). Each probability has been observed though three independant
infection assays, and can take values of $0$, $0.5$, and $1.0$.

Measure           Binary      Bernoulli trials        Probabilistic 
------------      -------     ------------------      ------------------
links             336         $221.58\pm 57.57$       $221.52\pm 57.25$
$\eta$             0.73        0.528                   0.512
$\eta^{(R)}$       0.72        0.525                   0.507
$\eta^{(C)}$       0.75        0.531                   0.518

- connectance

- nestedness

As shown in \autoref{t:poullain}, transforming the probabilistic matrix
into a binary one (i) overestimates nestedness by $\approx 0.02$, and (ii)
overestimates the number of links by 115. For the number of links, both
the probabilistic measures and the average and variance of $10^4$ Bernoulli
trials were in strong agreement (they differ only by the second decimal place).

Using Bernoulli trials had the effect of slightly over-estimating
nestedness. The overestimation is statistically significant from a purely
frequentist point of view, but significance testing is rather meaningless
when the number of replicates is this large and can be increased arbitrarily;
what is important is that the relative value of the error is small enough that
Bernoulli trials are able to adequately reproduce the probabilistic structure
of the network. It is not unexpected that Bernoulli trials are this close to
the analytical expression of the measures; due to the experimental design
of the @poul08 study, probabilities of interactions are bound to be high,
and so variance is minimal (most elements of $\mathbf{A}$ have a value
of either $0$ or $1$, and so their individual variance is $0$). Still,
despite overall low variance, the binary approach severely mis-represents
the structure of the network.

## Null-model based hypothesis testing

In this section, we analyse the data of @robe29 using two "classical"
null models of network structure. Robertson's data are amongst the hardest
to analyse with the standard null models: the network is unusually large
(1429 animals and 456 plants), and has a low connectance (0.02). Generating
networks with all species is therefore both statistically difficult and
computationally costly, providing a good demonstration of the performance
of probabilistic metrics.

We use the following null models. First (Type I, @fort06), any
interaction between plant and animals happens with the fixed probability
$\text{P}=Co$. This model controls for connectance, but removes the effect
of degree distribution. Second, (Type II, @basc03), the probability of
an interaction between animal $i$ and plant $j$ is $(k_i/R+k_j/C)/2$, the
average of the richness-standardized degree of both species.

Note that this type of null models will take a binary network, and through some
rules, turn it into a probabilistic one. Typically, this probabilistic network
is used as a template to generate Bernoulli trials, measure some of their
properties, the distribution of which is compared to the empirical network.

# References
