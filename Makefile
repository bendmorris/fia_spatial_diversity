grid_figures = grid_drift.png grid_homogenizing\ dispersal.png grid_selection.png
dist_figures = dist_drift.png dist_homogenizing\ dispersal.png dist_selection.png
fia_figures = fia_richness.png fia_abundance.png

all : $(grid_figures) $(dist_figures) $(fia_figures)

clean : 
	rm *.png *.pkl

$(grid_figures) : map_grid_results.py grid_results.pkl
	python map_grid_results.py

grid_results.pkl : fia.csv fia.newick grid_analysis.py
	python grid_analysis.py

$(dist_figures) : plot_dist_results.py dist_results.pkl
	python plot_dist_results.py

dist_results.pkl : fia.csv fia.newick dist_analysis.py
	python dist_analysis.py

$(fia_figures) : fia.csv plot_fia.py
	python plot_fia.py