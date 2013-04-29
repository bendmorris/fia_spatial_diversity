import csv
import Bio.Phylo as bp
import metrics
import numpy as np
import random
import multiprocessing
import cPickle as pkl
import sys
import time
from data import input_file, tree


# grid size for grouping communities, in degrees
GRID_SIZE = 1
# number of pairwise route comparisons to perform per grid cell
COMPARISONS = 1000
# ignore grid cells unless they have at least this many routes
MIN_SITES = 10


# read in route/species abundance information from FIA data file
grids = {}
all_species = set()
with open(input_file) as data_file:
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

grids = {x:y for x, y in grids.iteritems() if len(y) >= MIN_SITES}
print len(grids), 'total grids'
# get the range of lat/lon values
lats, lons = [route[0] for route in routes], [route[1] for route in routes]
lat_range = (min(lats), max(lats))
lon_range = (min(lons), max(lons))


results = {}
def analyze(arg):
    grid, routes = arg
    species_pool = []
    for route in routes.values():
        for sp, count in route.iteritems(): 
            species_pool += [sp] * count

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
    
    # compare pairs of communities
    comms = []
    while len(comms) < min(COMPARISONS, comparisons):
        try: r1, r2 = next(to_compare)
        except StopIteration: break
    
        if len(routes[r1]) < 2 or len(routes[r2]) < 2: continue
    
        try:
            # compute beta NTI
            nti = metrics.beta_nti(routes[r1], routes[r2], tree, 
                                   verbose=False, reps=1000)
            if nti == np.nan: continue

            if abs(nti) >= 2:
                # abs(beta NTI) >= 2 indicates selection
                result = 'selection'
            else:
                # < 2: raup-crick to determine drift or dispersal
                rc = metrics.raup_crick(routes[r1], routes[r2], species_pool)
                if rc <= -0.95:
                    result = 'homogenizing dispersal'
                elif rc >= 0.95:
                    result = 'dispersal limitation'
                else:
                    result = 'drift'
            comms.append(result)
        except IndexError:
            # this means a species wasn't found in our tree
            pass
    
    print '**', grid, time.strftime('%D %T')
    results = {}
    for result in sorted(set(comms)):
        percent = 100*len([c for c in comms if c == result]) / float(len(comms))
        print '%s: %s%%' % (result, percent)
        results[result] = percent
    
    return grid, results


if __name__ == '__main__':
    results = dict([
                    result for result in 
                    multiprocessing.Pool().map(analyze, grids.iteritems())
                    if result
                    ])

    with open('grid_results.pkl', 'w') as results_file:
        pkl.dump(results, results_file, -1)
