import cPickle as pkl
import matplotlib.pyplot as plt


results = pkl.load(open('dist_results.pkl'))
result_types = set([key
                    for result_dict in results.itervalues()
                    for key in result_dict])

for result_type in result_types:
    plt.figure()
    plt.title(result_type)
    plt.xlabel('distance (km)')
    plt.ylabel('% of pairwise comparisons')
    xs, ys = [], []
    xs = results.keys()
    ys = [results[x][result_type] if result_type in results[x] else 0
          for x in xs]
    plt.scatter(xs, ys)
    plt.savefig('dist_%s.png' % result_type)
