from mpl_toolkits.basemap import Basemap
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as s
import math

x_shift = -0.5
y_shift = 0
ROUND = 0.5

def xround(x, n):
    return int(x/n + n/2) * n


with open('fia.csv') as data_file:
    reader = csv.reader(data_file)
    # skip header row
    next(reader)
    richness, abundance, plot_density = {}, {}, {}
    for lat, lon, genus, species, count in reader:
        y, x = xround(float(lat),ROUND), xround(float(lon),ROUND)
        if not (y, x) in richness:
            richness[y,x] = set()
        if not (y, x) in abundance:
            abundance[y,x] = 0
        if not (y, x) in plot_density:
            plot_density[y,x] = set()
        species_name = '%s %s' % (genus, species)
        richness[y, x].add(species_name)
        abundance[y, x] += int(count)
        plot_density[y, x].add((lat, lon))
        

for point in richness:
    richness[point] = np.log(len(richness[point]))
for point in abundance:
    abundance[point] = np.log(abundance[point])
for point in plot_density:
    plot_density[point] = np.log(len(plot_density[point]))
        
for title, results in ('richness', richness), ('abundance', abundance), ('plot density', plot_density):
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

    if title == 'richness':
        ays = data
    elif title == 'plot density':
        axs = data

    plt.savefig('fia_%s.png' % title.replace(' ', '_'))

plt.figure()
axs = math.e**axs.ravel()
ays = math.e**ays.ravel()
plt.scatter(axs, ays)
m, b, r, p, se = s.linregress(axs, ays)
x1, x2 = min(axs), max(axs)
plt.plot([x1, x2], [b+m*x1, b+m*x2], 
         label='y=%sx+%s\nr^2=%s' % 
         tuple([round(n, 2) for n in (m, b, r**2)]))
plt.xlabel('plot density')
plt.ylabel('species richness')
plt.legend(loc='lower right')
plt.savefig('fia_plt_dens_richness.png')
