# pyranges_plot
Gene visualization package for dataframe objects generated with [PyRanges](https://pyranges.readthedocs.io/en/latest/index.html).




## Overview
The goal is getting a plot displaying a series of genes contained in a dataframe 
from a PyRanges object. It displays the genes' intron-exon structure in its 
corresponding chromosome.

There are some features to be defined by the user, one is the plot's **engine** 
since it can be based on Matplotlib or Plotly, the other is the name of the 
**gene ID** column in the data. The rest of features can either be left as default 
or be customized. In example, the plot shows the first 25 genes of the dataframe 
by default but this can be modified. It is worth noting that the order of the genes 
will be conserved.

In the case of coloring, Pyranges Plot offers a wide versatility. The data feature 
(column) according to which the genes will be colored is by default the gene ID, but 
this "color column" can be selected manually. Color specifications can be left as the 
default colormap (``plotly.colors.sequential.thermal``) or be provided as dictionaries, 
lists and color objects from either Matplotlib or Plotly regardless of the chosen engine. 
When a colormap or list of colors is specified, the colors assigned to the genes will 
iterate over the provided ones following the color column pattern. In the case of concrete 
color instructions such as dictionary, the genes will be colored according to it while the 
non-specified ones will be colored in black.

<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/general_ex.png">
</p>




## Installation
PyRanges-Plot can be installed using pip:

```
pip install pyranges-plot
```



## Examples
Next we will test pyranges_plot visualization options, using a PyRanges object generated 
from a dictionary.

```python
import pyranges as pr
import pyranges_plot as prplot

p = pr.from_dict({"Chromosome": [1, 1, 2, 2, 2, 2, 2, 3],
             	"Strand": ["+", "+", "-", "-", "+", "+", "+", "+"],
             	"Start": [1, 40, 10, 70, 85, 110, 150, 140],
             	"End": [11, 60, 25, 80, 100, 115, 180, 152],
             	"transcript_id":["t1", "t1", "t2", "t2", "t3", "t3", "t3", "t4"],
              "feature1": ["a", "a", "b", "b", "c", "c", "c", "d"],
              "feature2": ["A", "A", "B", "B", "C", "C", "C", "D"]})
print(p)

```
```
+--------------+--------------+-----------+-----------+-----------------+------------+------------+
|   Chromosome | Strand       |     Start |       End | transcript_id   | feature1   | feature2   |
|   (category) | (category)   |   (int64) |   (int64) | (object)        | (object)   | (object)   |
|--------------+--------------+-----------+-----------+-----------------+------------+------------|
|            1 | +            |         1 |        11 | t1              | a          | A          |
|            1 | +            |        40 |        60 | t1              | a          | A          |
|            2 | +            |        85 |       100 | t3              | c          | C          |
|            2 | +            |       110 |       115 | t3              | c          | C          |
|            2 | +            |       150 |       180 | t3              | c          | C          |
|            2 | -            |        10 |        25 | t2              | b          | B          |
|            2 | -            |        70 |        80 | t2              | b          | B          |
|            3 | +            |       140 |       152 | t4              | d          | D          |
+--------------+--------------+-----------+-----------+-----------------+------------+------------+
Stranded PyRanges object has 8 rows and 7 columns from 3 chromosomes.
For printing, the PyRanges was sorted on Chromosome and Strand.
```


The generated data is a stranded PyRanges object containing 4 genes in 3 chromosomes 
as shown above. Having this example data in the variable ``p`` we are able to start exploring 
pyranges_plot options. We can get a plot in a single line:

```python
prplot.plot_generic(p, engine="plt", id_col="transcript_id")
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex01.png">
</p>



The output is an interactive Matplotlib plot. To obtain it we just need to provide the data, 
the engine and the name of the id column. However the engine and the id column can be set 
previously so there is no need to specify them anymore while plotting:

```python
# For engine use 'plotly' or 'ply' for Plotly plots and 'matplotlib' or 'plt' for Matplotlib plots
prplot.set_engine('plotly')
prplot.set_idcol('transcript_id')
```

Since the data has only 4 genes all of them are plotted, but the function has a default limit 
of 25, so in a case where the data contains more genes it will only show the top 25, unless 
the ``max_ngenes`` parameter is specified. For example we can set the maximum number of genes 
as 2. Note that in the case of plotting more than 25 a warning about the plot's integrity 
will appear.

```python
prplot.plot_generic(p, max_ngenes=2)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex02.png">
</p>



Now the plot is based in Plotly because we set it as the engine, though it looks the same as the 
Matplotlib one. Also, both libraries offer interactive zoom options. For Matplotlib…
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_fixex03.png">
</p>

and for Plotly.
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_fixex04.png">
</p>



Another pyranges_plot functionality is allowing to define the plots' coordinate limits through 
the ``custom_coords`` parameter. The default limits show some space between the first and last 
plotted exons of each chromosome, but these can be customized. The user can decide to change 
all or some of the coordinate limits leaving the rest as default if desired. The limits should 
be provided as a dictionary, where the keys are the data's chromosome names in string format, 
and the values are either ``None`` or a tuple indicating the limits. When a chromosome is not 
specified in the dictionary or it is assigned ``None`` the coordinates will appear as default.

```python
prplot.plot_generic(
    p,
    custom_coords={"1": (None, 100), "2": (60, 200), "3": None}
)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex05.png">
</p>



We can try to color the genes according to the strand column instead of the ID (default). For 
that the ``color_col`` parameter should be used.

```python
prplot.plot_generic(p, color_col="Strand")
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex06.png">
</p>



This way we see the "+" strand genes in one color and the "-" in another color. Additionally 
these colors can be customized through the ``colormap`` parameter to see it more clearly. For 
this case we can specify it as a dictionary in the following way:

```python
prplot.plot_generic(
    p,
    color_col="Strand",
    colormap={"+": "green", "-": "red"}
)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex07.png">
</p>



The parameter ``colormap`` is very versatile because it accepts dictionaries for specific coloring, 
but also Matplotlib and Plotly color objects such as colormaps (or even just the string name of 
these objects) as well as lists of colors. For example we can use the Dark2 Matplotlib colormap, 
even if the plot is based on Plotly:

```python
prplot.plot_generic(p, colormap="Dark2")
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex08.png">
</p>



The disposition of the genes is by default a packed disposition, so the genes are preferentially 
placed one beside the other preferentially. But this ``disposition`` parameter can be set as 'full' 
if the user wants to display each gene under the other.

```python
prplot.plot_generic(p, disposition='full')
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex09.png">
</p>



In interactive plots there is the option of showing information about the gene when the mouse is 
placed over its structure. This information always shows the gene's start and end coordinates 
along with the ID. To add information contained in other dataframe calumns to the tooltip, the 
``showinfo`` parameter should be used in the following way:

```python
prplot.plot_generic(p, showinfo=["feature1", "feature2"])
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex10.png">
</p>



Lastly, some features of the plot appearance can also be customized. The way to change 
the default features is using the ``set_default`` function. The background, plot border or title
default colors can be checked and customized in the following way:

```python
# Check the default values
prplot.get_default('plot_background')

# Change the default values
prplot.set_default('plot_background', 'rgb(173, 216, 230)')
prplot.set_default('plot_border', '#808080')
prplot.set_default('title_dict_ply.color', 'magenta')

# Make the customized plot
prplot.plot_generic(p)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex11.png">
</p>


## Coming soon
* Bases will be displayed along coordinates
* Colorblind friendly
