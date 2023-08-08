# BindPPI-1.0
## About
<font size=4> 
  
Before running BindPPI models of RF_{13}, you need to create a folder that includes two input files. And for the AvgEns model, in addition to the two files mentioned earlier, you also need additional '.seq' files for each chain of the complex.

(see the example of structure_based_job1)

Before running BindPPI model for MLP_{5120}, you just need to create a JSON file for input sequences.

(see the example of sequence_based/sample_input_sequence.json)
  
</font>

## Three input files in the folder of structure_based_job1
<font size=4> 

1. The 3D structure of a protein (2J12.pdb1), which can be obtained from the Protein Data Bank (PDB) or created by the user.

2. The file containing complex information (structure_based_job1.input), who's name must be consistent with the input folder name.

- PDBfile: coordinate file of a protein structure.
- Partner1: the selected protein chains of Partner1 that will be taken into account during the calculation.If there are multiple chains, separate them with dots. For example, use 'A_1.B_1' for multiple chains.
- Partner2: the selected protein chains of Partner2 that will be taken into account during the calculation.If there are multiple chains, separate them with dots. For example, use 'A_1.B_1' for multiple chains.

  The columns are separated by tabs.

3. The '.seq' files (2J12_A.seq, 2J12_B.seq) correspond to each chain in the complex. Please note that the chain identifiers in the files and file names are consistent with those in the complex structure.

</font>


