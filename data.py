import Bio.Phylo as bp
import sys


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