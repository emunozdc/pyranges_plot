# pyranges_plot
Gene visualization package for dataframe objects generated with [PyRanges](https://pyranges.readthedocs.io/en/latest/index.html).


## Overview
The goal is getting a plot displaying a series of genes contained in a dataframe from
a PyRanges object. It displays the genes in its corresponding chromosome subplot. The
user can choose whether the plot is based on Matplotlib or Plotly by setting the engine. 
The plot will not contain all the genes in the dataframe, by default it shows 25 genes 
but this number can be customized too. It is worth noting that the order of the genes 
will be maintained.
 

Pyranges plot offers a wide versatility for coloring. The data feature (column) according
to which the genes will be colored is by default the gene ID, but this “color column” can 
be selected manually. Color specifications can be left as the default colormap or be 
provided as dictionaries, lists and color objects from either Matplotlib or Plotly regardless
of the chosen engine. When a colormap or list of colors is specified, the color of the genes 
will iterate over the given colors following the color column pattern. In the case of concrete 
color instructions such as dictionary, the genes will be colored according to it while the 
non-specified ones will be colored in black(??).

<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/general_ex.png">
</p>



## Installation
PyRanges-Plot can be installed using pip:

```
pip install pyranges-plot
```


## Examples
Next we will test pyranges_plot visualization options, using some data provided in 
[PyRanges tutorial](https://pyranges.readthedocs.io/en/latest/tutorial.html). Download 
and unpack tutorial data with:

```
curl -O https://mariottigenomicslab.bio.ub.edu/pyranges_data/pyranges_tutorial_data.tar.gz
tar zxf pyranges_tutorial_data.tar.gz
```

Once we have data files to work with in our working directory we will initiate python to 
load and subset the data into a dataframe.

```python
import pyranges as pr
import pyranges_plot as prplot
import pandas as pd
ann = pr.read_gff3('Dgyro_annotation.gff')
ann = ann[ [ 'ID'] ]
df = ann.df
df = df.rename(columns={'ID': 'gene_id'}) # ID column name to standard
```

Having some example data in the variable ``df`` we can start exploring pyranges_plot options. 
We can get our plot in a single line:

```python
prplot.plot_exons(df, engine=’plt’)
```
<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_example01.png">
</p>


In this plot we can see the top 25 genes in the dataframe in a Matplotlib plot, we just need 
to provide the data and engine. However the engine can be set previously so there is no need 
to specify it anymore while plotting:

```python
# Use ‘plotly’ or ‘ply’ for Plotly plots and ‘matplotlib’ or ‘plt’ for Matplotlib plots
prplot.set_engine(‘plotly’)
prplot.plot_exons(df)
```
<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_example02.png">
</p>


The plot looks the same as the previous one, but in this case is a Plotly plot. Note that in 
both libraries there are interactive zoom options. For Matplotlib…
<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_example03.png">
</p>

and for Plotly.
<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_example04.png">
</p>


We can try to color the genes according to strand and providing a dictionary for the colors, 
for that we will subset the dataframe to see  15 genes from each strand. In order to see all 
those genes we will set the max_ngenes to 30, since it is more than 25 genes a warning will appear:

```python
df2 = pd.concat([df.loc[df['Strand'] == '+'].head(15), df.loc[df['Strand'] == '-'].head(15)])
prplot.plot_exons(df2, max_ngenes=30, color_column='Strand', colormap={'+': 'green', '-': 'red'})
```
<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_example05.png">
</p>


Some features of appearance can also be customized. The way to change the default variables 
is using the ``set_default`` function. The background color, the plot border color or the title 
color can be customized in the following way:

```python
prplot.set_default('plot_background', 'black')
prplot.set_default('plot_border', 'lightblue')
prplot.set_default('title_dict_ply.color', ‘magenta’)
prplot.plot_exons(df)
```
<p align="center">
	<img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_example06.png">
</p>


## Coming soon
* Bases will be displayed along coordinates

