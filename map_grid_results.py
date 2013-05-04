from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np
from data import map


x_shift = -1
y_shift = 0

results = dict(pkl.load(open('grid_results.pkl')))
result_types = set([key 
                    for result_dict in results.itervalues()
                    for key in result_dict]) 
lats = [lat for lat, lon in results.keys()]
lons = [lon for lat, lon in results.keys()]

for result_type in result_types:
    points = []
    plt.figure()
    plt.title(result_type)
    map = get_map(lats, lons, x_shift, y_shift)
    map.drawcoastlines(linewidth=1)
    map.drawcountries(linewidth=1)
    map.drawstates(linewidth=0.5)

    xs, ys = [], []
    for (lat, lon) in results.iterkeys():
        xs.append(lon)
        ys.append(lat)

    xs = np.array(sorted(xs))
    ys = np.array(sorted(ys))
    x, y = np.meshgrid(xs, ys) 
    data = np.zeros(x.shape)
    for i in xrange(x.shape[0]):
        for j in xrange(x.shape[1]):
            try: data[i,j] = results[y[i][j],x[i][j]][result_type]
            except: data[i,j] = 0
            
    map.pcolormesh(x+x_shift, y+y_shift,data=data,latlon=True,
                   cmap=plt.cm.OrRd)
    cbar = plt.colorbar()
    cbar.set_label('% site comparisons')
    plt.savefig('grid_%s.png' % result_type)
