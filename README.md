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
