import csv
import Bio.Phylo as bp
import metrics
import numpy as np


tree = bp.read('fia_result.new', 'newick')

# grid size for grouping communities, in degrees
GRID_SIZE = 1

# read in route/species abundance information from FIA data file
grids = {}
with open('fia.csv') as data_file:
    reader = csv.reader(data_file)
    # skip header row
    next(reader)
    for lat, lon, genus, species, count in reader:
        lat, lon = float(lat), float(lon)
        count = int(count)
        if species == 'unknown': continue

        species_name = '%s %s' % (genus, species)
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

