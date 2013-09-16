# functions used in Stegen et al. 2013

import numpy as np
import Bio.Phylo as bp
import random


# cache species distances to avoid redundant work
cached_distances = {}


def distance(sp1, sp2, tree):
    key = tuple(sorted((sp1,sp2)))
    if not key in cached_distances:
        cached_distances[key] = tree.distance(sp1, sp2)
    return cached_distances[key]


def nearest_neighbor(species, neighbors, tree):
    '''Find the nearest neighbor of a species in another community.

    Returns (species name, distance)
    '''

    distances = [(neighbor, distance(species, neighbor, tree)) for neighbor in neighbors]
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


def shuffled_communities(comm1, comm2, tree, n=1000):
    '''Return two shuffled dictionaries of {species name:abundance}'''
    all_spp = set(comm1.keys() + comm2.keys())
    a, b = len(comm1), len(comm2)
    subtree = tree.common_ancestor(*all_spp)

    terms = subtree.get_terminals()
    
    for i in xrange(n):
        yield [{k: v for k, v in zip(random.sample(terms, len(c)), c.itervalues())} 
               for c in (comm1, comm2)]


def beta_nti(comm1, comm2, tree, reps=1000, verbose=False):
    '''Calculates beta_mntd for two communities and compares that value to the
    mean of (reps) identical communities after shuffling species names and
    abundances.

    Returns a z score; |beta_nti| > 2 indicates ecological selection is
    driving community assembly.'''

    mntd = beta_mntd(comm1, comm2, tree)
    if verbose: print 'Beta-MNTD=%s' % mntd

    distribution = []
    shuffles = shuffled_communities(comm1, comm2, tree, n=reps)
    for n in xrange(reps):
        if verbose: print 'Rep', n+1,
        random_comm1, random_comm2 = next(shuffles)
        distribution.append(beta_mntd(random_comm1, random_comm2, tree))
        if verbose: print 'Beta-MNTD=%s' % distribution[-1]

    std = np.std(distribution)
    if std == 0: return 0
    mean = np.mean(distribution)
    return (mntd - mean) / std
    



def raup_crick(comm1, comm2, species_pool, reps=1000):
    s1, s2 = set(comm1), set(comm2)
    a1, a2 = len(comm1), len(comm2)
    ssobs = len(set.intersection(s1, s2))

    ssexp = [len(set.intersection(set(random.sample(species_pool, a1)),
                                  set(random.sample(species_pool, a2))
                                  )) 
             for _ in xrange(reps)]

    return ((len([s for s in ssexp if s > ssobs]) +
             len([s for s in ssexp if s == ssobs])/2.
             )/(len(ssexp)) - 0.5) * 2


def process(r1, r2, tree, species_pool, reps=1000):
    nti = beta_nti(r1, r2, tree,
                   verbose=False, reps=reps)
    if nti == np.nan: return

    if abs(nti) >= 2:
        # abs(beta NTI) >= 2 indicates selection
        return 'selection' + ('+' if nti > 0 else '-')
    else:
        # < 2: raup-crick to determine drift or dispersal
        rc = raup_crick(r1, r2, species_pool)
        if rc <= -0.95:
            return 'homogenizing dispersal'
        elif rc >= 0.95:
            return 'dispersal limitation'
        else:
            return 'drift'
