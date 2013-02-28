import Bio.Phylo as bp
import metrics

tree = bp.read('test.newick', 'newick')

comm1 = {'Homo_sapiens':10, 'Mus_musculus':3, 'Ursus_americanus': 12}
comm2 = {'Ursus_arctos':8, 'Rattus_norvegicus':16}
comm3 = {'Gazella_dorcas':100, 'Bos_taurus':5, 'Canis_lupus':5}

print metrics.beta_nti(comm1, comm2, tree)
print metrics.beta_nti(comm1, comm3, tree)