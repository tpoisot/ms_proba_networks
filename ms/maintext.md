% The structure of probabilistic networks
% The Magical Mistery Stouffer Group
% Working paper -- Oct. 2014

# Introduction

> This is a summary of the discussions we had during the lab meeting, with
substantially more equations. I don't know about you, but the expressions
of *variance* got me over-excited. I haven't put that into the document yet,
but I have ran a few tests and it is all matching perfectly.

Each measure defines a *property* on one or several *network units*. These
properties can be defined by the unit itself (*direct properties*), or
require the association of several units. See \autoref{f:component}.

> More seriously -- we should discuss about a table link each direct property
to its unit, and then emerging properties to the groups of units.

![Why Tim had his whiteboard privileges revoked.\label{f:component}](figures/whiteboard.jpg)

# Metrics

Throughout this section, we will assume the following notation. $\mathbf{A}$
is a matrix wherein $A_{ij}$ is $\text{P}(ij)$, *i.e.* the probability that
species $i$ establishes an interaction with species $j$. If $\mathbf{A}$
represents a unipartite network (*e.g.* a food web), it is a square matrix and
the probabilities of each species interacting with itself. If $\mathbf{A}$
represents a bipartite network (*e.g.* a pollination network), it will
most likely not be square. We call $S$ the number of species, and $R$ and
$C$ respectively the number of rows and columns. $S = R + C$ in unipartite
networks, and $S = R+C$ in bipartite networks. Note that all of the measures
defined below can be applied on a bipartite network that has been made
unipartite; the unipartite transformation of a bipartite matrix $\mathbf{A}$
is the block matrix

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

## Emerging properties

### Expected trophic level

### Katz centrality

Centrality can be measured by the degree, as mentionned above. In
addition, we derive the expected value of centrality according to
@katz53. This measures generalizes to directed acyclic graphs. Although
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
probabilistic network.

### Number of primary producers

Primary producers, in a food web, are species with no successors, including
themselves. Biologically, they are autotrophic organisms, or organisms
whose preys of substrates have been remove from the network. A species is a
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
or of the variance of any three species forming a given motif, but that shit
be hardcore yo, and so they are better expressed computationally.

# Applications

## Structure of neutral networks

## Null-model based hypothesis testing

# References
