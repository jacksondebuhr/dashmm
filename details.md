---
layout: default
title: Details
---

# DASHMM Details

This page includes more information about the mathematical and computational details of the DASHMM library.

## Multipole Methods

Multipole methods are widely applied to N-body or N-body-like problems. The basic problem is to compute the mutual interaction of a large number of particles. The essence of multipole methods is the realization that the net effect of a large number of individual interactions with distant objects can often be treated as a single interaction. Multipole methods systematically leverage this simplification to drastically speed up the, in principle, O(N^2) interactions that would need to be computed. There are many variations on how to perform these simplifications, but we will focus on the Fast Multipole Method (FMM) and the Barnes-Hut method (BH) not only because they are exemplars of multipole methods, but also because they are widely used in end-science applications.

Both FMM and BH begin by recursively subdividing the computational domain, naturally creating a tree structure. The root of this structure represents the entire volume, with each child representing a fraction of the volume of its parent. The subdivision of the computational domain proceeds until some stopping criterion is met, typically based on the number of particles that are represented by the leaf nodes of the tree.

Once this tree is constructed, the multipole moments are computed for each node of the tree, which can then be used to speed up the computation of the interaction for a distant point for the particles represented by that node. Typically moments are computed bottom-up with the moments at the leaves being computed directly from the particles in that node, and the moments for internal nodes being computed from the moments of that node's children. At the end of this process, each node has a multipole expansion for the effect of the particles contained within it. At this point, FMM and BH diverge in how they use these moments.

### Barnes-Hut Method

Having constructed the tree, it is now possible to compute the interaction (typically the force, or the potential) at a particular point, called a target point. Often these points are the locations of the particles themselves, but the method allows the computation of the interaction at any point. To compute the effect at the target point of the particles represented by the tree involves a traversal of the tree. At each node, there is a node acceptance criterion that is applied to decide if the approximation represented by the multipole moments for that node introduce an acceptable level of error. If the node passes the acceptance criterion, the multipole expansion can be used. If not, the node is 'opened' and the traversal proceeds to the children of that node. The node acceptance criterion is essentially a statement that the node is far enough away from the target point. Ultimately the result of this traversal is a set of nodes that contribute to the interaction for the target point.

### Fast Multipole Method

Nearby target points in the BH method have very similar tree traversals, using the multipole expansion for the same far away nodes. The FMM removes this repeated work by translating the multipole expansion into a smooth //local// expansion in regions away from its singularity. The use of the local expansion permits FMM to continue operating on the box level until reaching to the finest subdivision level. Specifically, the FMM has two extra operations after the multipole expansions are computed. First, the multipole expansions of the nodes in the so-called interaction list region---those nodes far away enough that their expansions can be used, but not so far away that their parent's expansions can be used---are translated into local expansion. This allows the effect of a distant node to be applied only once to a given node at the same level of refinement. Second, the node inherits the local expansion from its parent. This allows the local expansion of higher level nodes to be applied to nodes at a finer level of refinement. At the finest subdivision level, the FMM resumes operation at the particle level: the local expansion of the node containing the target point is evaluated and interactions with neighboring particles are computed via the direct pairwise interactions. 

## ParalleX

ParalleX is an umbrella project that hosts and develops a set of close-knit components. It is an evolving parallel execution model derived to exploit the opportunities and address the challenges of emerging technology trends. It aims to ensure effective performance gain into the Exascale era of the next decade by addressing the challenges of starvation, latency, overhead, waiting, energy and reliability. At the core of the ParalleX strategy is a new framework to replace static methods with dynamic adaptive techniques. This method exploits runtime information and employs unused resources. It benefits from locality while managing distributed asynchrony. ParalleX is a crosscutting model to facilitate co-design and interoperability among system component layers from hardware architecture to programming interfaces.

ParalleX distinguishes itself by:
*   Dynamic adaptive methods through runtime system software for guiding computing resource management and task scheduling
*   Lightweight threads (Multi-threaded)
*   Parcels (a form of active messages)
*   Global name space and virtual addressing
*   Runtime discovery, synchronization, and control
*   Local Control Objects (LCOs)

## HPX-5

[HPX-5](https://www.crest.iu.edu/projects/hpx/) is an advanced supercomputing runtime system that implements the ParalleX model. It is the performance-oriented representation of the ParallelX execution model targeted at conventional parallel computing architectures such as SMP nodes and commodity clusters. HPX-5 is a C library that makes use of lightweight threads and an active global address space that allows both code and data to freely move across the system in order to adapt to rapidly changing computational needs in applications as well as environmental factors such as system load, power utilization, real-time events, and network performance. HPX-5 targets cutting-edge concepts in high-performance computing such as fast M: N co-routines that can migrate between physical nodes, lock-free synchronization of parallel operations, software-defined interconnect networks, distributed atomics, power-aware work scheduling, and support for distributed and multi-core embedded platforms.

## Computational Approach

Multipole Methods have been applied to a wide variety of end science applications. Please see [[applications|Applications]] for an overview of some of the applications to which MM have been applied. One feature of these applications is that the systems studied are frequently very dynamic, both in their time dependence and the range of relevant spatial scales. However, conventional parallelization techniques are essentially static. These slowly-adapting methods will be unable to effectively handle the size of problems considered in the future and will be unable to scale to ever-larger computational resources. DASHMM will leverage an advanced runtime system, HPX-5, to provide a dynamic and adaptive execution of multipole methods. 

DASHMM will use lightweight ephemeral threads, exposing more parallelism at runtime, and allowing the expression of the computational task to match the data flow inherent in multipole methods. These lightweight threads, coupled with LCOs, will allow dependencies in the data to control the flow of execution in the parallel computation. Additionally, DASHMM will use the active data management of HPX-5 to update the data layout in response to changes in the state of the computational resources.