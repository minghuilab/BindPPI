# BindPPI-v1.0
## About
<font size=4> 
  
To run RF_{13}, you need to create a folder and prepare two input files (refer to the files 1a2k.input and 1A2K.pdb in the 1a2k folder).

To run AvgEns, you need to create a folder and prepare three input files (refer to the files 1a2k.input, 1A2K.pdb, and 1a2k.json in the 1a2k folder).


To run MLP_{5120}, you simply need to create a JSON file containing the input sequences for each complex (refer to the file sample_input_sequence.json in the MLP5120_example folder).


## The format of "sample_input_sequence.json"
For each complex, there are two partners (Partner1 and Partner2) listed under each complex, with each partner containing a sequence list that can comprise single or multiple chain sequences.
```python
{
  'Complex1': {
        'Partner1': [['IVGGYTCAANSIPYQQLQGIVVCNYVNWIQQTIAAN']],
        'Partner2': [['KKVCACPKILKPVCGSDGRTYANSCIARCNGVSIKS']]
    			},
  'Complex2': {
        'Partner1': [['DIKMTQSPSSMYASLGERVTKTSTSPIVKSFNRNEC'],['EIQLQQSGAELVRPGALVKLSCKASAVLQSDLASSI']],
        'Partner2': [['SGTTNTVAAYNLTWKSTNFKTILESSGKKTAKTNTN']]
    			}
}
```

</font>

## Three input files in the folder of 1a2k
<font size=4> 

1. The 3D structure of a protein (1A2K.pdb), which can be obtained from the Protein Data Bank (PDB) or created by the user.

2. The file containing complex information (1a2k.input), who's name must be consistent with the input folder name.

- PDBfile: coordinate file of a protein structure.
- Partner1: the selected protein chains of Partner1 that will be taken into account during the calculation. If there are multiple chains, separate them with dots. For example, use 'A_1.B_1' for multiple chains.
- Partner2: the selected protein chains of Partner2 that will be taken into account during the calculation. If there are multiple chains, separate them with dots. For example, use 'A_1.B_1' for multiple chains.

  The columns are separated by tabs.

3. The JSON file containing complex sequence (1a2k.json), who's name must be consistent with the input folder name. The specific format is consistent with "sample_input_sequence.json".

</font>


