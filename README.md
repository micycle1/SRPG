[![](https://jitpack.io/v/micycle1/SRPG.svg)](https://jitpack.io/#micycle1/SRPG)

# SRPG – Super Random Polygon Generator

SRPG is a library that generates well-formed random polygons using a regular grid consisting of square cells. 

By default, SRPG generates orthogonal polygons on this grid. However, it can also produce octagonal polygons by cutting off corners with ±45° diagonals during the construction. Repeatedly cutting corners without the diagonal restriction yields an approximation of a smooth free-form curve. Moreover, SRPG can apply perturbations to generate polygons with axes-parallel edges whose vertices do not lie on a grid or generate polygons whose edges are not parallel to the coordinate axes.

This library is based on *Martin Held*'s C [implementation](https://github.com/cgalab/genpoly-srpg).

## Usage
SRPG is available as Maven/Gradle artifact via [Jitpack](https://jitpack.io/#micycle1/SRPG).

## Argument Illustrations

SRPG accepts a variety of input arguments that affect shape geometry, providing a good level of customisation over the output. The following illustrations provide an idea of how each argument visually affects the random polygon:

<table>
  <tr>
    <td align="center" valign="center"><b>Hierarchy</td>
    <td align="center" valign="center"><b>Diagonal</td>
    <td align="center" valign="center"><b>Perturb+align</td>
    <td align="center" valign="center"><b>Smooth</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/hierarchy.gif"></td>
    <td valign="top" width="25%"><img src="resources/diagonal1.gif"></td>
    <td valign="top" width="25%"><img src="resources/perturbAlign.gif"></td>
    <td valign="top" width="25%"><img src="resources/smooth.gif"></td>
  </tr>
</table>

### Grid size (nX, nY)

<table>
  <tr>
    <td align="center" valign="center"><b>5x5</td>
    <td align="center" valign="center"><b>10x10</td>
    <td align="center" valign="center"><b>100x100</td>
    <td align="center" valign="center"><b>500x500</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/n_5.png"></td>
    <td valign="top" width="25%"><img src="resources/n_20.png"></td>
    <td valign="top" width="25%"><img src="resources/n_100.png"></td>
    <td valign="top" width="25%"><img src="resources/n_500.png"></td>
  </tr>
</table>