import Bio.Phylo as bp
import sys
from mpl_toolkits.basemap import Basemap


try:
    input_file, tree_file = sys.argv[1], sys.argv[2]
except:
    input_file, tree_file = 'fia.csv', 'fia.newick'

tree = bp.read(tree_file, 'newick')
terms = tree.get_terminals()

spp = {}
no_match = set()
def find_species(s, tree):
    if not s in spp:
        try:
            spp[s] = next(x for x in terms if x.name == s)
        except StopIteration:
            try:
                spp[s] = next(x for x in terms if x.name.startswith(s.split()[0]))
            except:
                no_match.add(s)
                return None
    return spp[s]


def get_map(lats, lons, x_shift, y_shift):
    map = Basemap(llcrnrlon=max(-125,min(lons))+x_shift,llcrnrlat=min(lats)+y_shift,
                  urcrnrlon=max(lons)+x_shift,urcrnrlat=min(50,max(lats))+y_shift,
                  projection='merc',lat_1=33,lat_2=45,lon_0=-95,resolution='l')

    return map