import cPickle as pkl
import matplotlib.pyplot as plt


results = pkl.load(open('dist_results.pkl'))
result_types = set([key
                    for result_dict in results.itervalues()
                    for key in result_dict])

colors = {
          'selection': 'red',
          'selection+': 'red',
          'selection-': 'orange',
          'drift': 'blue',
          'dispersal limitation': 'green',
          'homogenizing dispersal': 'purple'
          }

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


plt.figure()
plt.xlabel('distance (km)')
plt.ylabel('% of pairwise comparisons')
plt.ylim(-1,101)
for result_type in sorted(result_types):
    xs, ys = [], []
    xs = sorted(results.keys())
    ys = [results[x][result_type] if result_type in results[x] else 0
          for x in xs]
    plt.plot(xs, ys, color=colors[result_type], label=result_type)
    plt.legend(loc='upper left')
    plt.savefig('dist.png')
