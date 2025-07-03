# Circular Protein Builder

A command-line tool for generating symmetric protein rings (e.g., gasdermin-like structures) from a monomeric subunit or an approximate input ring structure.


## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.


## Features

- Automatically fits a ring radius from a low-resolution structure
- Aligns and arranges monomeric subunits in a circular oligomer
- Outputs a complete PDB structure of the resulting circular oligomer (centered on x0, y0)

- Optionally accepts a fixed radius.
- Optionally rotates monomeric subunits in z, xy, or into the circle.
- Optionally cleans up intermediate files.

---

## Installation

```bash
git clone https://github.com/s-schaef/model_circular_proteins.git
cd model_circular_proteins
conda env create -f environment.yml
conda activate model_circular_proteins
```

Use `build_circle --help` to show usage guidelines. 

## Examples
This tool allows building circular structures from monomers. Below are three examples of how this may be used:

1. Simple example to practically recreate the circular assembly of a known ring-shaped protein:

    ```bash
    build_circle --in_structure examples/6vfe.pdb --no_subunits 33 --z_rotation 80
    ```

    This builds a 33-meric circular assebmly of Gasdermin-D that resembles the original 33-meric structure from which the monomer is taken (6VFE)[[1]](#1). The radius is fitted from the original structure. The z_rotation is needed since the original chain A does not sit at 12 o'clock in the ring. Note that this z_rotation is not perfect and the resulting ring will slightly vary from the original. By design of the code it also assembles the ring from copies of chain A, whereas the original structure may vary between its subunits. 

2. Let's assume we want to build a tighter ring based off of the same original structure (6VFE). For a 15-meric Gasdermin-D strucuture we now need to specify a radius that is tighter - via trial and error 70Angstrom seems to fit well. 

    ```bash
    build_circle --in_structure examples/6vfe.pdb --no_subunits 15 --z_rotation 80 --in_radius 70
    ```

3. Now we want to take a subunit from an assymmetric assembly and build a circular assembly from it. For example a bacterial gasdermin that stems from a spiral structure from Cryo-EM (8SL0)[[2]](#2) that we want to model as a planar assembly. For this we need to guess (or know from other methods) the number of subunits and radius we want the ring to have. Then we need to find the z, xy-rotations, and the tilt into the circle via trial and error. The xy_rotation is crucial here to preserve the hydrogen-bonding network of the $\beta$-sheet.

    ```bash
    build_circle --in_structure examples/8sl0.pdb --in_radius 178
     --no_subunits 52 --z_rotation -80 --xy_rotation 78 --tilt 10
    ```




## References
<a id="1">[1]</a> 
Xia, S., Zhang, Z., Magupalli, V.G. et al.
Gasdermin D pore structure reveals preferential release of mature interleukin-1. 
Nature 593, 607–611 (2021). 
https://doi.org/10.1038/s41586-021-03478-3

<a id="2">[2]</a> 
Johnson, A.G., Mayer, M.L., Schaefer, S.L. et al.
Structure and assembly of a bacterial gasdermin pore.
Nature 628, 657–663 (2024).
https://doi.org/10.1038/s41586-024-07216-3