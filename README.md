# MOEA/D Overlapping Protein Complex Detection


## How to run

Install requirements:

```bash
pip install -r requirements.txt
```

Run:

```bash
python main.py
```

The PPI network files and known-complex files are included in the project:

```text
DataSets/
  PPI/
    1-Protein-Yeast-D1-Files.mat
    2-Protein-Yeast-D2-Files.mat
    4-Collins-CYC2008.mat
    5-Collins-MIPS.mat
  Complex/
    Complex-D1-Files.mat
    Complex-D2-Files.mat
    4-CYC2008.mat
    5-MIPS.mat
```

`main.py` loads the matching PPI and complex files based on the network number entered by the user.

The output is saved as a pickle file inside:

```text
RepOverlap/MOEA/
```
