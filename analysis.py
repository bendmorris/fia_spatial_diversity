import csv
import Bio.Phylo as bp
import metrics
import numpy as np
import random


tree = bp.read('fia_result.new', 'newick')
#tree.get_path = metrics.get_path_with_cache
tree = metrics.CachingTree(tree)

# grid size for grouping communities, in degrees
GRID_SIZE = 1
# number of pairwise route comparisons to perform per grid cell
COMPARISONS = 100

# read in route/species abundance information from FIA data file
grids = {}
all_species = set()
with open('fia.csv') as data_file:
    reader = csv.reader(data_file)
    # skip header row
    next(reader)
    for lat, lon, genus, species, count in reader:
        lat, lon = float(lat), float(lon)
        count = int(count)
        if species == 'unknown': continue

        species_name = '%s %s' % (genus, species)
        all_species.add(species_name)
        route = (lat,lon)

        grid = (int(lat), int(lon))
        if not grid in grids:
            grids[grid] = {}
        routes = grids[grid]

        if not route in routes:
            routes[route] = {}
        routes[route][species_name] = count


# get the range of lat/lon values
lats, lons = [route[0] for route in routes], [route[1] for route in routes]
lat_range = (min(lats), max(lats))
lon_range = (min(lons), max(lons))


for grid, routes in grids.iteritems():
    n = len(routes)
    # compare all combinations of routes if n choose 2 < COMPARISONS,
    # otherwise compare random combinations until you reach COMBINATIONS
    # total comparisons
    comparisons = n*(n-1) / 2
    if comparisons < COMPARISONS:
        to_compare = ((r1, r2) for r1 in routes for r2 in routes if not r1 == r2)
    else:
        def random_comparison():
            while True:
                yield tuple(sorted(random.sample(routes.keys(), 2)))
        to_compare = random_comparison()
    
    # build a list of beta_ntis by comparing pairs of communities
    ntis = []
    while len(ntis) < min(COMPARISONS, comparisons):
        try: r1, r2 = next(to_compare)
        except StopIteration: break
    
        if len(routes[r1]) < 2 or len(routes[r2]) < 2: continue
    
        try:
            m = metrics.beta_nti(routes[r1], routes[r2], tree, 
                                 verbose=False, reps=1000)
            if not m == np.nan:
                ntis.append(m)
        except IndexError:
            # this means a species wasn't found in our tree
            pass
    
    print grid, ntis