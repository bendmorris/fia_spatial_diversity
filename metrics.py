# functions used in Stegen et al. 2013

import numpy as np
import Bio.Phylo as bp
import random


# cache species distances to avoid redundant work
cached_distances = {}
valid_names = {}


class CachingTree(bp.Newick.Tree):
    _paths = {}
    def __init__(self, tree):
        self.__dict__.update(tree.__dict__)

    def get_path(self, target, **kwargs):
        if not target in self._paths:
            self._paths[target] = bp.Newick.Tree.get_path(self, target=target, 
                                                          **kwargs)
        return self._paths[target]

def distance(tree, *species):
    species = tuple(sorted(species))
    if not species in cached_distances:
        try:
            cached_distances[species] = tree.distance(*species)
            valid_names[species[0]] = species[0]
            valid_names[species[1]] = species[1]
        except KeyboardInterrupt: raise
        except:
            nodes = []
            for s in species:
                if s in valid_names: s = valid_names[s]
                elif tree.find_any(s):
                    valid_names[s] = s
                else:
                    terms = [t.name for t in tree.get_terminals()]
                    genus = s.split()[0]
                    congeners = [t for t in terms if t.startswith(genus + ' ')]
                    valid_names[s] = congeners[0]
                    s = valid_names[s]
                nodes.append(s)
            species = tuple(sorted(nodes))
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

    std = np.std(distribution)
    if std == 0: return np.nan
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

    return (len([s for s in ssexp if s > ssobs]) +
            len([s for s in ssexp if s == ssobs])/2.
            )/(len(ssexp)) - 0.5 * 2