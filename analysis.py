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
COMPARISONS = 10

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
    tried = set()
    n = 0
    while len(tried) < COMPARISONS:
        s1, s2 = tuple(sorted(random.sample(routes.keys(), 2)))
        #if (s1, s2) in tried: continue
        if len(routes[s1]) < 2 or len(routes[s2]) < 2: continue
        #tried.add((s1, s2))

        print n, s1, s2
        print len(routes[s1]), len(routes[s2])
        try:
            m = metrics.beta_nti(routes[s1], routes[s2], tree, 
                                 verbose=False, reps=1000)
            print m
            n += 1
        except IndexError:
            # this means a species wasn't found in our tree
            print 'fail'
