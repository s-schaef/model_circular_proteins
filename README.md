# Circular Protein Builder

A command-line tool for generating symmetric protein rings (e.g., gasdermin-like structures) from a monomeric subunit or an approximate input ring structure.


## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.


## Features

- Automatically fits a ring radius from a low-resolution structure
- Aligns and arranges monomeric subunits in a circular oligomer
- Outputs a complete PDB structure of the resulting circular oligomer (centered on x0, y0)

- Optionally accepts a fixed radius.
- Optionally rotates monomeric subunits in z or xy.
- Optionally cleans up intermediate files.

---

## Requirements

Tested with:

- Python **3.11**
- MDAnalysis **2.9.0**
- numpy **2.2.5**
- scipy **1.15.3**

Install these exact versions using:

```bash
pip install -r requirements.txt
```

## Examples
This tool allows building circular structures from monomers. Below are three examples of how this may be used:

1. Simple example to practically recreate the circular assembly of a known ring-shaped protein:

    ```
    python model_circular_proteins.py --in_structure examples/6vfe.pdb --no_subunits 33
     --z_rotation 80
    ```

    This builds a 33-meric circular assebmly of Gasdermin-D that resembles the original 33-meric structure from which the monomer is taken (6VFE, cite). The radius is fitted from the original structure. The z_rotation is needed since the original chain A does not sit at 12 o'clock in the ring. Note that this z_rotation is not perfect and the resulting ring will slightly vary from the original. By design of the code it also assembles the ring from copies of chain A, whereas the original structure may vary between its subunits. 

2. Let's assume we want to build a tighter ring based off of the same original structure (6VFE). For a 15-meric Gasdermin-D strucuture we now need to specify a radius that is tighter - via trial and error 70Angstrom seems to fit well. 

    ```
    python model_circular_proteins.py --in_structure examples/6vfe.pdb --no_subunits 15 
    --z_rotation 80 --in_radius 70
    ```

3. Now we want to take a subunit from an assymmetric assembly and build a circular assembly from it. For example a spiral gasdermin structure from Cryo-EM that we want to model as a planar assembly. For this we need to guess (or know from other methods) the number of subunits and radius we want the ring to have. Then we need to find the z, xy-rotations, and the tilt into the circle via trial and error. The xy_rotation is crucial here to preserve the hydrogen-bonding network of the $\beta$-sheet.

    ```
    python model_circular_proteins.py --in_structure examples/8sl0.pdb --in_radius 178
     --no_subunits 52 --z_rotation -80 --xy_rotation 78 --tilt 10
    ```

