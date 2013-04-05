import csv
import Bio.Phylo as bp
import metrics
import numpy as np
import random
import multiprocessing
import math
import cPickle as pkl
import sys

try:
    input_file, tree_file = sys.argv[1], sys.argv[2]
except:
    input_file, tree_file = 'fia.csv', 'fia.newick'


# grid size for grouping communities, in degrees
BIN_SIZE = 50
MAX_BIN = 2500
result_bins = {bin_size: [] for bin_size in np.arange(0, MAX_BIN, BIN_SIZE)}
# number of pairwise route comparisons to perform per grid cell
COMPARISONS = 1000

def radians(deg):
    return deg/180.*(math.pi)

def distance(p1, p2):
    '''Approximate distance between two lat/lon points.'''
    y1, x1 = p1
    y2, x2 = p2
    dy = radians(abs(y2-y1))
    dx = radians(abs(x2-x1))
    a = math.sin(dy/2)**2 + math.sin(dx/2)**2 * math.cos(radians(y1)) * math.cos(radians(y2))
    c = 2 * math.atan2(a**0.5, (1-a)**0.5)
    return 6371 * c

tree = bp.read(tree_file, 'newick')
#tree.get_path = metrics.get_path_with_cache
tree = metrics.CachingTree(tree)


# read in route/species abundance information from FIA data file
routes = {}
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

        loc = (lat, lon)
        if not loc in routes:
            routes[loc] = {}
        routes[loc][species_name] = count

# get the range of lat/lon values
lats, lons = [route[0] for route in routes], [route[1] for route in routes]
lat_range = (min(lats), max(lats))
lon_range = (min(lons), max(lons))


def analyze(arg):
    bin_size, results = arg
    
    while len(results) < COMPARISONS:
        r1, r2 = random.sample(routes, 2)
        while not bin_size <= distance(r1, r2) <= bin_size + BIN_SIZE:
            r1, r2 = random.sample(routes, 2)
            
        lat_range = min((r1[0], r2[0])), max((r1[0], r2[0]))
        lon_range = min((r1[1], r2[1])), max((r1[1], r2[1]))
        
        species_pool = []
        for lat, lon in routes:
            if (lat_range[0] <= lat <= lat_range[1] and
                lon_range[0] <= lon <= lon_range[1]):
                for sp, count in routes[lat,lon].iteritems():
                    species_pool += [sp] * count

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
        except IndexError:
            # this means a species wasn't found in our tree
            continue
        
        print bin_size, result
        results.append(result)
    
    return (bin_size, results)


if __name__ == '__main__':
    results = multiprocessing.Pool().map(analyze, result_bins.iteritems())
    results = dict(results)
    results = {a: {type: 100.*len([x for x in b if x == type])/len(b) 
                   for type in set(b)} 
               for a, b in results.iteritems()}

    with open('dist_results.pkl', 'w') as results_file:
        pkl.dump(results, results_file, -1)
