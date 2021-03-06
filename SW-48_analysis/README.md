This folder contains the perturbation data and dose response data of the SW-48 cell line and the scripts that were used to analyse it.

First, run `process_data.py` to generate the response pattern from the raw data.
The raw data can be plotted with the `plot_raw_data.py` and `plot_perturbation_comparison.py` scripts in the plotting folder.

The response logic inference is carried out by the `inference.py` script that takes extra constraints from `inference_input.py`. It first generates the scores for the parental cell line comparison to literature network (and saves them in results). Second, it also carries out the mutant cell line comparison and saves all predicted graphs in plots/graphs. This requires a Python interface for Graphviz which you can pip install like this

```
pip install graphviz
```

You can also switch off the graph plotting and thereby remove the graphviz dependency by setting `draw_graphs=False` in Line 9 of `inference.py`.
