from mpl_toolkits.basemap import Basemap
import csv
import matplotlib.pyplot as plt
import numpy as np

x_shift = 0
y_shift = 0
ROUND = 1


with open('fia.csv') as data_file:
    reader = csv.reader(data_file)
    # skip header row
    next(reader)
    richness, abundance = {}, {}
    for lat, lon, genus, species, count in reader:
        lat, lon = round(float(lat),ROUND), round(float(lon),ROUND)
        if not (lat, lon) in richness:
            richness[lat,lon] = set()
        if not (lat, lon) in abundance:
            abundance[lat,lon] = list()
        species_name = '%s %s' % (genus, species)
        richness[lat, lon].add(species_name)
        abundance[lat, lon].append(species_name)

for point in richness:
    richness[point] = len(richness[point])
for point in abundance:
    abundance[point] = len(abundance[point])
        
for title, results in ('richness', richness), ('abundance', abundance):
    plt.figure()
    
    lats = set([lat for lat, lon in results.keys()])
    lons = set([lon for lat, lon in results.keys()])

    map = Basemap(llcrnrlon=min(lons)+x_shift,llcrnrlat=min(lats)+y_shift,
                  urcrnrlon=max(lons)+x_shift,urcrnrlat=max(lats)+y_shift,
                  projection='merc',lat_1=33,lat_2=45,lon_0=-95,resolution='l')
    map.drawcoastlines(linewidth=1)
    map.drawcountries(linewidth=1)
    map.drawstates(linewidth=0.5)

    xs = np.array(sorted(list(lons)))
    ys = np.array(sorted(list(lats)))
    x, y = np.meshgrid(xs, ys)
    data = np.zeros(x.shape)
    for i in xrange(x.shape[0]):
        for j in xrange(x.shape[1]):
            try: data[i,j] = results[y[i][j],x[i][j]]
            except: data[i,j] = 0

    map.pcolormesh(x+x_shift, y+y_shift,data=data,latlon=True,
                   cmap=plt.cm.OrRd)
    cbar = plt.colorbar()
    cbar.set_label(title)

    plt.savefig('fia_%s.png' % title)
