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
    for x in results.keys():
        if result_type in results[x]:
            xs.append(x)
            ys.append(results[x][result_type])
    plt.scatter(xs, ys)
    plt.savefig('dist_%s.png' % result_type)