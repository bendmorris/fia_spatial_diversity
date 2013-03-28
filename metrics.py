# functions used in Stegen et al. 2013

import numpy as np
import random


# cache species distances to avoid redundant work
cached_distances = {}
def distance(tree, *species):
    species = tuple(sorted(species))
    if not species in cached_distances:
        cached_distances[species] = tree.distance(*species)
    return cached_distances[species]


def nearest_neighbor(species, neighbors, tree):
    '''Find the nearest neighbor of a species in another community.

    Returns (species name, distance)
    '''

    distances = [(neighbor, distance(tree, species, neighbor)) for neighbor in neighbors]
    distances.sort(key=lambda d: d[1])

    return distances[0]


def beta_mntd(comm1, comm2, tree):
    '''Quantifies the phylogenetic distance between OTUs in two communities.

    comm1 and comm2 are dictionaries mapping species names to abundances.
    tree is a phylogenetic tree that connects the OTUs.

    Returns a single numeric value.
    '''

    total = 0
    for c1, c2 in ((comm1, comm2), (comm2, comm1)):
        for (sp, n) in c1.iteritems():
            neighbor, neighbor_distance = nearest_neighbor(sp, c2, tree)
            total += n * neighbor_distance

    return total * 0.5


def shuffled_community(comm):
    '''Return a shuffled dictionary of {species name:abundance}'''
    values = comm.values()
    random.shuffle(values)

    return {k:v for k, v in zip(comm.iterkeys(), values)}


def beta_nti(comm1, comm2, tree, reps=1000, verbose=False):
    '''Calculates beta_mntd for two communities and compares that value to the
    mean of (reps) identical communities after shuffling species names and
    abundances.

    Returns a z score; |beta_nti| > 2 indicates ecological selection is
    driving community assembly.'''

    mntd = beta_mntd(comm1, comm2, tree)
    if verbose: print 'Beta-MNTD=%s' % mntd

    distribution = []
    for n in xrange(reps):
        if verbose: print 'Rep', n+1,
        random_comm1 = shuffled_community(comm1)
        random_comm2 = shuffled_community(comm2)
        distribution.append(beta_mntd(random_comm1, random_comm2, tree))
        if verbose: print 'Beta-MNTD=%s' % distribution[-1]

    return (mntd - np.mean(distribution)) / np.std(distribution)