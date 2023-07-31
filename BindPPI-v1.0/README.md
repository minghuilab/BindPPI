# BindPPI-1.0
## About
<font size=4> 
  
Before running BindPPI models for RF_{13}, MLP_{5120///512, 13//416}, and AvgEns, you need to create a folder that includes two input files.

(see the example of 202307091239)

Before running BindPPI model for MLP_{5120}, you just need to create a JSON file for input sequence.

(see the example of sample_input_sequence.json)
  
</font>

## Two input files in the folder of 202307091239
<font size=4> 

1. The 3D structure of a protein (2J12.pdb1), which can be obtained from the Protein Data Bank (PDB) or created by the user.

2. The file containing complex information (202307091239.input), who's name must be consistent with the input folder name.

- PDBfile: coordinate file of a protein structure.
- Partner1: the selected protein chains of Partner1 that will be taken into account during the calculation.
- Partner2: the selected protein chains of Partner2 that will be taken into account during the calculation.

  The columns are separated by tabs.

</font>


