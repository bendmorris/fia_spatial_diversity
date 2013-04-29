grid_figures = grid_drift.png grid_homogenizing\ dispersal.png grid_selection.png
dist_figures = dist_drift.png dist_homogenizing\ dispersal.png dist_selection.png
fia_figures = fia_richness.png fia_abundance.png fia_plot_density.png fia_plt_dens_richness.png

all : $(fia_figures)

grid : $(grid_figures)
	cp grid_results.pkl $(grid_figures) ~/Dropbox/Dev/
dist : $(dist_figures)
	cp dist_results.pkl $(dist_figures) ~/Dropbox/Dev/

clean : 
	rm *.png *.pkl

fia.csv : mk_csv.py fia.sqlite query.sql
	python mk_csv.py

$(grid_figures) : map_grid_results.py grid_results.pkl
	python map_grid_results.py

grid_results.pkl : fia.csv fia.newick grid_analysis.py metrics.py
	python grid_analysis.py

$(dist_figures) : plot_dist_results.py dist_results.pkl
	python plot_dist_results.py

dist_results.pkl : fia.csv fia.newick dist_analysis.py metrics.py
	python dist_analysis.py

$(fia_figures) : fia.csv fia_maps.py
	python fia_maps.py
