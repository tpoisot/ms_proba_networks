% The structure of probabilistic networks
% T. Poisot \and A. Cirtwill \and D. Gravel \and M.-J. Fortin \and D.B. Stouffer
% Mar. 2015

# Abstract

1. There is a growing realization among community ecologists that interactions
between species vary in space and time. Yet, our current numerical framework to
analyze the structure of interactions, largely based on graph-theoretical
approaches, is unsuited to this type of data. Since the variation of species
interactions holds much information, there is a need to develop new metrics to
exploit it.
2. We present analytical expressions of key network metrics, using a
probabilistic framework. Our approach is based on modeling each interaction as a
Bernoulli event, and using basic calculus to express the expected value, and
when mathematically tractable, its variance. We provide a free and open-source
implementation of these measures.
3. We show that our approach allows to overcome limitations of both neglecting
the variation of interactions (over-estimation of rare events) and using
simulations (extremely high computational demand). We present a few case studies
that highlight how these measures can be used.
4. We conclude this contribution by discussing how the sampling and data
representation of ecological network can be adapted to better allow the
application of a fully probabilistic numerical framework.

**Keywords:** ecological networks, connectance, degree distribution, nestedness, modularity

\clearpage

# Introduction

Ecological networks are an efficient way to represent biotic interactions
between individuals, populations, or species. Historically, their study focused
on describing their structure, with a particular attention on food webs
[@dunn06] and plant-pollinator interactions [@basc03; @jord87]. The key result
of this line of research was linking this structure to community or
ecosystem-level properties such as stability [@mcca14], coexistence [@bast09;
@haer14], or ecosystem functioning [@theb03; @duff02; @pois12]. To a large
extent, the description of ecological networks resulted in the emergence of
questions about how functions emerged from structure, and this stimulated the
development of a rich methodological literature, defining a wide array of
structural properties.

Given a network as input, measures of network structure return a *property*
based on one or several *units* from this network. Some of the properties are
*direct* properties (they only require knowledge of the unit on which they are
applied), whereas others are *emergent* (they require knowledge of, and
describe, higher-order structures). For example, connectance, the realized
proportion of potential interactions, is a direct property of a network. The
degree of a node (how many interactions it is involved in) is a direct property
of the node. The nestedness of a network (that is, the extent to which
specialists and generalists overlap), on the other hand, is an emergent property
that is not directly predictable from the degree of all nodes. Though the
difference may appear to be semantics, establishing a difference between direct
and emergent properties is important when interpreting their values; direct
properties are conceptually equivalent to means, in that they tend to be the
first moment of network units, whereas emergent properties are conceptually
equivalent to variances or other higher-order moments.

In the recent years, the interpretation of the properties of network structure
(as indicators of the action of ecological or evolutionary processes) has been
somewhat complicated by the observation that network structure varies through
space and time. This happens because, contrary to a long-standing assumption of
network studies, species from the same pool do not interact in a consistent way
[@pois12c]. Empirical and theoretical studies suggest that the network is not
the right unit to understand this variation; rather, network variation is an
emergent property of the response of ecological interactions to environmental
factors and chance events [@pois14]. Interactions can vary because of local
mismatching in phenology [@oles11a], populations fluctuations preventing the
interaction [@cana14], or a combination of both [@olit14; @cham14]. For example,
@olit14 show that accounting for neutral (population-size driven) and
trait-based effects allows the prediction of the cumulative change in network
structure, but not of the change at the level of individual interactions. In
addition, @cars14 show that within a meta-community, not all interactions are
equally variable: some are highly consistent, whereas others are extremely rare.
These empirical results all point to the fact that species interactions cannot
always be adequately modeled as yes-no events; since it is well established that
they do vary, it is necessary to represent them as probabilities. To the
question of *Do these two species interact?*, we should substitute the question
of *How likely is it that they will interact?*.

The current way of dealing with probabilistic interactions are either to ignore
variability entirely or to generate random networks. Probabilistic metrics are a
mathematically rigorous alternative to both. When ignoring the probabilistic
nature of interactions (henceforth *binary* networks), every non-zero element of
the network is assumed to be 1. This leads to over-representation of some rare
events, and increases the number of interactions. An alternative is to consider
only the interactions above a given threshold, which leads to an
under-representation of rare events and decreases the effective number of
interactions. Taken together, these considerations highlight the need to amend
our current methodology for the description of ecological networks, in order to
give more importance to the variation of individual interactions. Because the
methodological corpus available to describe ecological networks had first been
crafted at a time when it was assumed that interactions were invariants, it is
unsuited to address the questions that probabilistic networks allow us to ask.


In this paper, we show that several direct and emergent core properties of
ecological networks (both bipartite and unipartite) can be re-formulated in a
probabilistic context [@yeak12; @pois14]; we conclude by showing how this
methodology can be applied to exploit the information contained in the
variability of networks, and to reduce the computational burden of current
methods in network analysis. We also provide a free and open-source (MIT
license) implementation of this suite of measures in a library for the `julia`
language, available at `http://github.com/PoisotLab/ProbabilisticNetwork.jl`.

# Suite of probabilistic network metrics

Throughout this paper, we use the following notation. $\mathbf{A}$ is a matrix
wherein $A_{ij}$ is $\text{P}(ij)$, *i.e.* the probability that species $i$
establishes an interaction with species $j$. If $\mathbf{A}$ represents a
unipartite network (*e.g.* a food web), it is a square matrix and contains the
probabilities of each species interacting with all others, including itself. If
$\mathbf{A}$ represents a bipartite network (*e.g.* a pollination network), it
will not necessarily be square. We call $S$ the number of species, and $R$ and
$C$ respectively the number of rows and columns. $S = R = C$ in unipartite
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

where $0_{(C,R)}$ is a matrix of $C$ rows and $R$ columns (noted $C\times R$)
filled with $0$s, etc. Note that for centrality to be relevant in bipartite
networks, this matrix should be made symmetric:
$\mathbf{B}_{ij}=\mathbf{B}_{ji}$.

We will also assume that all interactions are independent (so that
$\text{P}(ij|kl) = \text{P}(ij)\text{P}(kl)$ for any species), and can be
represented as a series of Bernoulli trials (so that $0 \leq \text{P}(ij) \leq
1$). The latter condition allows us to derive estimates for the *variance*
($\text{var}(X) = p(1-p)$), and expected values ($\text{E}(X)=p$). We can
therefore estimate the variance of most properties, using the fact that the
variance of additive independent events is the sum of their individual
variances, and that the variance of multiplicative independent events is

\begin{equation}
\text{var}(X_1 X_2 ... X_n) = \prod_i \left(\text{var}(X_i) + [\text{E}(X_i)]^2\right) - \prod_i [\text{E}(X_i)]^2
\end{equation}

As a final note, all of the measures described below can be applied on the
binary (0/1) versions of the networks and will give the exact value of the
non-probabilistic measure. This property is particularly desirable as it allows
our framework to be used on any network, whether they are represented in a
probabilistic or binary way.

## Direct properties

### Connectance and number of interactions

Connectance (or network density) is the proportion of possible interactions
that are realized, defined as $Co = L/(R\times C)$, where $L$ is the total
number of interactions. As all interactions in a probabilistic network are
assumed to be independent, the expected value of $L$, is

\begin{equation}
\hat L = \sum A_{ij},
\end{equation}

and ${\hat{Co}} = \hat L / (R\times C)$. Likewise, the variance of the number of
interactions is $\text{var}(\hat L) = \sum (A_{ij}(1-A_{ij}))$.

### Node degree

The degree distribution of a network is the distribution of the number of
interactions established (number of successors) and received (number of
predecessors) by each node. The expected degree of species $i$ is

\begin{equation}
\hat k_i = \sum_j(A_{ij} + A_{ji})
\end{equation}

The variance of the degree of each species is $\text{var}(\hat k_i) =
\sum_j(A_{ij}(1-A_{ij})+A_{ji}(1-A_{ji}))$. Note also that as expected,
$\sum \hat k_i = 2\hat L$.

### Generality and vulnerability

By simplification of the above, generality $\hat g_i$ and vulnerability $\hat
v_i$ are given by, respectively, $\sum_j A_{ij}$ and $\sum_j A_{ji}$, with their
variances $\sum_j A_{ij}(1-A_{ij})$ and $\sum_j A_{ji}(1-A_{ji})$.

## emergent properties

### Path length

Networks can be used to describe indirect interactions between species through
the use of paths. The existence of a path of length 2 between species $i$ and
$j$ means that they are connected through at least one additional species $k$.
In a probabilistic network, unless some elements are $0$, all pairs of species
$i$ and $j$ are connected through a path of length 1, with probability $A_{ij}$.
The expected number of paths of length $k$ between species $i$ and $j$ is given
by

\begin{equation}
\hat{n^{(k)}_{ij}} = \left(\mathbf{A}^k\right)_{ij},
\end{equation}

where $\mathbf{A}^k$ is the matrix multiplied by itself $k$ times.

It is possible to calculate the probability of having at least one path of
length $k$ between the two species: this can be done by calculating the
probability of having no path of length $k$, then taking the running product of
the resulting array of probabilities. For the example of length 2, species $i$
and $j$ are connected through $g$ with probability $A_{ig}A_{gj}$, and so this
path does not exist with probability $1-A_{ig}A_{gj}$. For any pair $i$, $j$,
let $\mathbf{m}$ be the vector such as $m_{g} = A_{ig}A_{gj}$ for all $g \notin
(i,j)$ [@mirc76]. The probability of not having any path of length 2 is $\prod
(1-\mathbf{m})$. Therefore, the probability of having a path of length 2 between
$i$ and $j$ is

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
$i$ to $j$, passing through any pair of nodes $x$ and $y$, without having any
shorter path. In theory, this approach can be generalized up to an arbitrary
path length, but it becomes rapidly untractable.

### Unipartite projection of bipartite networks

The unipartite projection of a bipartite network is obtained by linking any two
nodes of one mode that are connected through at least one node of the other
mode; for example, to plants are connected if they share at least one
pollinator. It is readily obtained using the formula in the *Path length*
section. This yields either the probability of an edge in the unipartite
projection (of the upper or lower nodes), or if using the matrix multiplication,
the expected number of such nodes.

### Nestedness

Nestedness is an important measure of (bipartite) network structure that tells
the extent to which the interactions of specialists and generalists overlap. We
use the formula for nestedness proposed by @bast09. They define nestedness for
each margin of the matrix, as $\eta^{(R)}$ and $\eta^{(C)}$ for, respectively,
rows and columns. As per @alme08, we define a global statistic for nestedness as
$\eta = (\eta^{(R)}+\eta^{(C)})/2$.

Nestedness, in a probabilistic network, is defined as

\begin{equation}
\hat{\eta^{(R)}} = \sum_{i<j}\frac{\sum_kA_{ik}A_{jk}}{\text{min}(g_i, g_j)},
\end{equation}

where $g_i$ is the expected generality of species $i$. The reciprocal holds
for $\eta^{(C)}$ when using $v_i$ (the vulnerability) instead of $g_i$.

The values returned are within $[0;1]$, with $\eta=1$ indicating complete
nestedness.

### Modularity

Modularity represents the extent to which networks are compartmentalized, *i.e.*
the tendency for subsets of species to be strongly connected together, while
they are weakly connected to the rest of the network [@stou11a]. Modularity is
measured as the proportion of interactions between nodes of an arbitrary number
of modules, as opposed to the random expectation. Assuming a vector $\mathbf{s}$
which, for each node in the network, holds the value of the module it belongs
to (an integer in $[1,c]$), @newm04a proposed a general measure of modularity, which is

$$\mathcal{Q} = \sum_{m=1}^{c}\left(e_{mm} - a_m^2\right)$$

, where $c$ is the number of modules,

$$e_{mm} = \sum_{ij}\frac{\mathbf{A}_{ij}}{2c} \delta(\mathbf{s}_i,\mathbf{s}_j)$$

, and

$$a_m = \sum_{n}e_{mn}$$,

with $\delta$ being Kronecker's function, returning $1$ if its arguments are
equal, and $0$ otherwise. This formula can be *directly* applied to
probabilistic networks.

### Centrality

Although node degree is a rough first order estimate of centrality, other
measures are often needed. We derive the expected value of centrality according
to @katz53. This measures generalizes to directed acyclic graphs (whereas other
do not). For example, although eigenvector centrality is often used in ecology,
it cannot be measured on probabilistic graphs. Eigenvector centrality requires
the matrix's largest eigenvalues to be real, which is not the case for all
probabilistic matrices. The measure proposed by Katz is a useful replacement,
because it accounts for the paths of all length between two species instead of
focusing on the shortest path.

As described above, the expected number of paths of length $k$ between $i$ and
$j$ is $(\mathbf{A}^k)_{ij}$. Based on this, the expected centrality of species
$i$ is

\begin{equation}
C_i = \sum_{j=1}^n \sum_{k=1}^\infty \alpha^k (\mathbf{A}^k)_{ji} .
\end{equation}

The parameter $\alpha \in [0;1]$ regulates how important long paths are. When
$\alpha = 0$, only first-order paths are accounted for (and the centrality is
equal to generality).
%DG: to the degree or generality?
 When $\alpha = 1$, paths of all length are equally
important. As $C_i$ is sensitive to the size of the matrix, we suggest
normalizing by $\mathbf{C} = \sum C$, so that

\begin{equation}
C_i = \frac{C_i}{\mathbf{C}} .
\end{equation}

This results in the *expected relative centrality* of each node in the
probabilistic network.

### Species with no outgoing links

Estimating the number of species with no outgoing links (successors) can be
useful when predicting whether, *e.g.*, predators will go extinct. A species has
no successors if it manages *not* to establish any outgoing interaction, which
for species $i$ happens with probability

\begin{equation}
\prod_j (1-A_{ij}).
\end{equation}

The number of expected such species is therefore the sum of the above across all
species:

\begin{equation}
\hat{PP} = \sum_i \left(\prod_j (1-A_{ij})\right).
\end{equation}

and its variance is

\begin{equation}
\text{var}(\hat{PP}) = \sum_i \left( \prod_j(1-A_{ij}^2) - \prod_j(1-A_{ij})^2 \right)
\end{equation}

Note that in a non-probabilistic context, species with no outgoing links would
be considered primary producers. This is not the case here: if interactions are
probabilistic events, then *e.g.* a top predator may have no preys, which do not
mean it will not become a primary producer. For this reason, the trophic
position of the species may better be measured on the binary version of the
matrix.

### Species with no incoming links

Using the same approach as for the number of species with no outgoing links, the
expected number of species with no incoming links is therefore

\begin{equation}
\hat{TP} = \sum_i\left(\prod_{j \neq i}(1-A_{ji})\right)
\end{equation}

Note that we exclude self-interactions, as top-predators can, and often
do, engage in cannibalism.

### Number of species with no interactions

Predicting the number of species with no interactions (or whether any species
will have at least one interaction) is useful when predicting whether species
will be able to integrate into an existing network, for example. Note that from
a methodological point of view, this can be a helpful *a priori* measure to
determine whether null models of networks will have a lot of species with no
interactions, and so will require intensive sampling.

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

### Self-loops

Self-loops (the existence of an interaction of a species onto itself) is only
meaningful in unipartite networks. The expected proportion of species with
self-loops is very simply defined as $\text{Tr}(\mathbf{A})$, that is, the sum
of all diagonal elements. The variance is
$\text{Tr}(\mathbf{A}\diamond(1-\mathbf{A}))$, where $\diamond$ is the
element-wise product operation.

### Motifs

Motifs are sets of pre-determined interactions between a fixed number of species
[@milo02; @stou07], such as for example one predator sharing two preys. As there are an
arbitrarily large number of motifs, we will illustrate the approach with only
two examples.

The probability that three species form an apparent competition motif (one
predator, two prey) where $i$ is the predator, $j$ and $k$ are the prey, is

\begin{equation}
\text{P}(i,j,k\in\text{app. comp}) = A_{ij}(1-A_{ji})A_{ik}(1-A_{ki})(1-A_{jk})(1-A_{kj})
\end{equation}

Similarly, the probability that these three species form an omnivory motif, in
which $i$ and $j$ consume $k$ and $i$ consumes $j$, is

\begin{equation}
\text{P}(i,j,k\in\text{omniv.}) = A_{ij}(1-A_{ji})A_{ik}(1-A_{ki})A_{jk}(1-A_{kj})
\end{equation}

The probability of the number of *any* motif $\text{m}$ with three species
in a network is given by

\begin{equation}
\hat{N_\text{m}} = \sum_i \sum_{j\neq i} \sum_{k\neq j} P(i,j,k \in \text{m})
\end{equation}

It is indeed possible to have an expression of the variance of this value, or of
the variance of any three species forming a given motif, but their expressions
become rapidly untractable and are better computed than written.

## Network comparison

The dissimilarity of a pair of (ecological) networks can be measured using the
framework set forth by @kole03. Measures of $\beta$-diversity compute the
dissimilarity between two networks based on the cardinality of three sets, $a$,
$c$, and $b$, which are respectively the shared items, items unique to superset
(network) 1, and items unique to superset 2 (the identity of which network is 1
or 2 matters for asymmetric measures). Supersets can be the species within each
network, or the interactions between species. Following @pois12c, the
dissimilarity of two networks can be measured as either $\beta_{WN}$ (all
interactions), or $\beta_{OS}$ (interactions involving only common species),
with $\beta_{OS} \leq \beta_{WN}$.

Within our framework, these measures can be applied to probabilistic networks.
The expected values of $\bar a$, $\bar c$, and $\bar b$ are, respectively, $\sum
\mathbf{A}_1\diamond\mathbf{A}_2$, $\sum \mathbf{A}_1\diamond(1-\mathbf{A}_2)$,
and $\sum (1-\mathbf{A}_1)\diamond\mathbf{A}_2$. Whether $\beta_{OS}$ or
$\beta_{WN}$ is measured requires to alter the matrices $\mathbf{A}_1$ and
$\mathbf{A}_2$. To measure $\beta_{OS}$, one must remove all unique species; to
measure $\beta_{WN}$, one must expand the two matrices so that they have the
same species at the same place, and give a weight of 0 to the added
interactions.

# Applications

In this section, we contrast the use of probabilistic measures to the current
approaches of either using binary networks, or working with null models through
simulations. When generating random networks, what we call *Bernoulli trials*
from here on, a binary network is generated by doing a Bernoulli trial with
probability $A_{ij}$, for each element of the matrix. This is problematic
because higher order structures involving rare events will be under-represented
in the sample, and because most naive approaches are likely to generate free
species, especially in sparsely connected networks frequently encountered in
ecology [@milo03; @pois14a] -- on the other hand, non-naive approaches break the
assumption of independence between interactions.

## Comparison of probabilistic networks

In this sub-section, we apply the above measures to a bacteria--phage
interaction network. @poul08 have measured the probability that 24 phages can
infect 24 strains of bacteria of the *Pseudomonas fluorescens* species (group
SBW25). Each probability has been observed though independent infection assays,
and can take values of $0$, $0.5$ (interaction is variable), and $1.0$.

Measuring the structure of the Binary, Bernoulli trials, and Probabilistic
network gives the following result:

| Measure      | Binary | Bernoulli trials  | Probabilistic     |
|:-------------|:-------|:------------------|:------------------|
| links        | 336    | $221.58\pm 57.57$ | $221.52\pm 57.25$ |
| $\eta$       | 0.73   | 0.528             | 0.512             |
| $\eta^{(R)}$ | 0.72   | 0.525             | 0.507             |
| $\eta^{(C)}$ | 0.75   | 0.531             | 0.518             |

As these results show, transforming the probabilistic matrix into a binary one
(i) overestimates nestedness by $\approx 0.2$, and (ii) overestimates the number
of links by 115. For the number of links, both the probabilistic measures and
the average and variance of $10^4$ Bernoulli trials were in strong agreement
(they differ only by the second decimal place).

Using Bernoulli trials had the effect of slightly over-estimating nestedness.
The overestimation is statistically significant from a purely frequentist point
of view, but significance testing is rather meaningless when the number of
replicates is this large and can be increased arbitrarily; what is important is
that the relative value of the error is small enough that Bernoulli trials are
able to adequately reproduce the probabilistic structure of the network. It is
not unexpected that Bernoulli trials are this close to the analytical expression
of the measures; due to the experimental design of the @poul08 study,
probabilities of interactions are bound to be high, and so variance is minimal
(most elements of $\mathbf{A}$ have a value of either $0$ or $1$, and so their
individual variance is $0$ -- though their confidence interval varies as a
function of the number of observations from which the probability is derived).
Still, despite overall low variance, the binary approach severely mis-represents
the structure of the network.

## Null-model based hypothesis testing

In this section, we analyse 59 pollination networks from the literature using
two usual null models of network structure, and two models with intermediate
constraints. These data cover a wide range a situations, from small to large,
and from densely to sparsely connected networks. They provide a good
demonstration of the performance of probabilistic metrics. Data come from the
*InteractionWeb Database*, and were queried on Nov. 2014.

We use the following null models. First (Type I, @fort06), any interaction
between plant and animals happens with the fixed probability $\text{P}=Co$. This
model controls for connectance, but removes the effect of degree distribution.
Second, (Type II, @basc03), the probability of an interaction between animal $i$
and plant $j$ is $(k_i/R+k_j/C)/2$, the average of the richness-standardized
degree of both species. In addition, we use the models called Type III in and
out [@pois13e], that use the row-wise and column-wise probability of an
interaction respectively, as a way to understand the impact of the degree
distribution of upper and lower level species.

Note that these null models will take a binary network, and through some rules
turn it into a probabilistic one. Typically, this probabilistic network is used
as a template to generate Bernoulli trials and measure some of their properties,
the distribution of which is compared to the empirical network. This approach is
computationally inefficient [@pois14a], especially using naive models [@milo03],
and as we show in the previous section, can yield biased estimates of the true
average of nestedness (and presumably other properties).

We measured the nestedness of the 59 (binary) networks, then generated the
random networks under the four null models, and calculated the expected
nestedness using the probabilistic measure. For each null model $i$, the
difference $\Delta^{(i)}_N$ in nestedness $N$ is expressed as $\Delta^{(i)}_N =
N-\mathcal{N}^{(i)}(N)$, where $\mathcal{N}^{(i)}(N)$ is the nestedness of null
model $i$. Our results are presented in \autoref{f:app2}.

\begin{figure}[bp]
   \begin{tikzpicture}
      \begin{groupplot}
      [group style={columns=2, horizontal sep=2cm},
      xmin=0, xmax=0.6, ymin=0, ymax=0.6]
         \nextgroupplot[xlabel=$\Delta^{(1)}_N$, ylabel=$\Delta^{(2)}_N$]
         \addplot [black!10, no markers] coordinates {(0,0) (0.6,0.6)};
         \addplot [only marks] table [x = d1, y = d2] {figures/app2.dat};
         \node[black!80, below left] at (axis cs:0.1,0.55){\textbf{A}};
         \nextgroupplot[xlabel=$\Delta^{(3i)}_N$, ylabel=$\Delta^{(3o)}_N$]
         \addplot [black!10, no markers] coordinates {(0,0) (0.6,0.6)};
         \addplot [only marks] table [x = d3i, y = d3o] {figures/app2.dat};
         \node[black!80, below left] at (axis cs:0.1,0.55){\textbf{B}};
      \end{groupplot}
   \end{tikzpicture}

   \caption{Results of the null model analysis of 59 plant-pollination
   networks. \textbf{A}. There is a consistent tendency for (i) both models I
   and II to estimate less nestedness than in the empirical network, although
   null model II yields more accurate estimates. \textbf{B}. Models III in
   and III out also estimate less nestedness than the empirical network,
   but neither has a systematic bias.}

   \label{f:app2}
\end{figure}

There are two striking results. First, empirical data are consistently *more*
nested than the null expectation, as evidenced by the fact that all $\Delta_N$
values are strictly positive. Second, this underestimation is *linear* between
null models I and II (in that it does not depends on how nested the empirical
network is), although null model II is always closer to the nestedness of the
empirical network (which makes sense, since null model II incorporates the
higher order constraint of respecting the degree distribution of both levels).
That the nestedness of the null model probability matrix is so strongly
determined by the nestedness of the empirical networks calls for a closer
evaluation of how the results of null models are interpreted (especially since
Bernoulli simulations revealed a very low variance in the simulated nestedness).

There is a strong, and previously unaccounted for, circularity in this approach:
empirical networks are compared to a null model which, as we show, has a
systematic bias *and* a low variance (in simulations), meaning that differences
in nestedness that are small (thus potentially ecologically irrelevant) have a
good chance of being reported as significant. Interestingly, models III in and
III out made overall *fewer* mistakes at estimating nestedness -- resp. $0.129$
and $0.123$, compared to resp. $0.219$ and $0.156$ for model I and II. Although
the error is overall sensitive to model type (Kruskal-Wallis $\chi^2$ = 35.80,
d.f. = 3, $p \leq 10^{-4}$), the three pairs of models that where significantly
different after controlling for multiple comparisons are I and II, I and III in,
and I and III out (model II is not different from either models III in or out).

In short, this analysis reveals that (i) the null expectation of a network
property under randomization scenarios can be obtained through the analysis
of the probabilistic matrix, instead of the analysis of simulated Bernoulli
networks; (ii) Different models have different systematic biases, with models
of the type III performing overall better for nestedness than any other
models. This can be explained by the fact that nestedness of a network,
as expressed by @bast09, is the average of a row-wise and column-wise
nestedness. These depend on the species degree, and as such should be well
predicted by models III.

# Implications for data collection

We developed and presented a set of measures to quantify the expected network
structure, using the probability that each interaction is observed or happens,
in a way that do not require time-consuming simulations. Our framework is set up
in such a way that the probabilities of interactions are considered to be
independent. Estimating interaction probabilities based on species abundances
[@cana14; @olit14] do not, for example, yield independent probabilities:
changing the abundance of one species changes all probabilities in the network.
They are not Bernoulli events either, as the sum of all probabilities derived
this way sums to unity. On the other hand, "cafeteria experiments" give truly
independent probabilities of interactions; even a simple criteria, such as the
frequency of interactions when the two species are put together, is a way of
estimating probability. Using the approach outline by @pois14, both sources of
information (species abundance and the outcome of experiments) can be combined
to estimate the probability that interactions will happen in empirical
communities. This effort requires improved communications between scientists
collecting data and scientists developing methodology to analyze them.

Another way to obtain approximation of the probability of interactions is to use
spatially replicated sampling. Some studies [@tyli07; @cars14; @olit14; @troj15]
surveyed the existence of interactions at different locations, and a simple
approach of dividing the number of observations of an interaction by the number
of co-occurence of the species involved will provide a (somewhat crude) estimate
of the probability of this interaction. This approach requires extensive
sampling, especially since interactions are harder to observe than species
[@pois12c; @gila14a], yet it enables the re-analysis of existing datasets in a
probabilistic context.

Understanding the structure of ecological networks, and whether it relates to
ecosystem properties, is emergent as a key challenge for community ecology. A
proper estimation of this structure requires tools that address all forms of
complexity, the most oft-neglected yet pervasive of which is the fact that
interactions are variable. By developing these metrics, we allow future analyses
of network structure to account for this phenomenon.

**Acknowledgements:** This work was funded by a CIEE working group grant to TP,
DG, and DBS. TP is funded by a starting grant from the Université de Montréal.

# References
