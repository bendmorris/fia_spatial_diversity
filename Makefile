grid_maps = grid_drift.png grid_homogenizing\ dispersal.png grid_selection.png

all : $(grid_maps)

clean : 
	rm *.png *.pkl

$(grid_maps) : map_grid_results.py grid_results.pkl
	python map_grid_results.py

grid_results.pkl : fia.csv fia.newick grid_analysis.py
	python grid_analysis.py

dist_results.pkl : fia.csv fia.newick dist_analysis.py
	python dist_analysis.py
