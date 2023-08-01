# BindPPI
## About
<font size=4>

BindPPI calculates the binding affinities for the protein-protein complexes to evaluate the binding strength quantitatively.It can be used for protein therapeutics discovery, de novo interface design, and computational mutagenesis.BindPPI provides four models for prediction, namely RF_{13}, MLP_{5120}, MLP_{5120///512, 13//416}, and AvgEns. For RF_{13}, MLP_{5120///512, 13//416}, and AvgEns, the input requires the 3D structure of the protein-protein complex, while MLP_{5120} only needs the sequence of the complex for computation.

</font>

## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/BindPPI/releases) from Github.

</font>

## Installation

### I. PREREQUISITES

<font size=4>

BindPPI requires the following software and packages.

#### Requirements for MLP_{5120} models.

1. Python3 packages: 

    To install these packages, you can install the required environment by modifying the prefix in the environment.yaml file and use the following command:

```
$ conda env create -f environment.yaml
```

#### Requirements for RF_{13}, MLP_{5120///512, 13//416}，AvgEns models 

1. VMD

   This is available on the VMD website.

   https://www.ks.uiuc.edu/Research/vmd/

2. CHARMM

   This is available on the CHARMM website.

   https://www.charmm.org/charmm/

3. NAMD

   This is available on the NAMD website.

   https://www.ks.uiuc.edu/Research/namd/

4. DSSP

   This is available on the DSSP website.

   https://swift.cmbi.umcn.nl/gv/dssp/

5. iPot

   This is available on the iPot Github.

   https://github.com/gjoni/iPot

6. McVol

   This is available on the McVol website.

   http://www.bisb.uni-bayreuth.de/index.php?page=data/mcvol/mcvol

7. PDB2PQR30

   This is available via pip according to the documentation.

   https://pdb2pqr.readthedocs.io/en/latest/index.html

8. PROVEAN

   This is available on the PROVEAN website.

   http://provean.jcvi.org/index.php/

9. FoldX

   This is available on the FoldX website.

   http://foldxsuite.crg.eu/
   
10. Python3 packages: 

    To install these packages, you can install the required environment by modifying the prefix in the environment.yaml file and use the following command:

```
$ conda env create -f environment.yaml
```

11. MLP_{5120///512, 13//416} models: 

    For the MLP_{5120///512, 13//416} model, you need to download its models and place it in the ./BindPPI-v1.0/inputfiles/ directory. 
    
    The download link is https://drive.google.com/drive/folders/1sXAz-SmV2kuWJC8TJY1F1G5MGhm36Z5n.

</font>

### II. INSTALLATION INSTRUCTIONS

<font size=4>

1. Download and/or install the prerequisites described above.

2. Download and unpack the distribution:

</font>

```
$ wget https://github.com/minghuilab/BindPPI/archive/BindPPI-v1.0.tar.gz
$ tar -zxvf BindPPI-v1.0.tar.gz
```

<font size=4>

3. Change to the source directory:

</font>

```
$ cd BindPPI-v1.0/
```

<font size=4>

4. Change the path parameters in BindPPI_RF13.py,BindPPI_5120-512-13-416.py,BindPPI_AvgEns.py:

</font>

```
workdir = Your working directory # '/data/jiang/tools/bindppi_test/' ##need change, not contain capital letters
pathinput = path for inputfiles directory # '/data/jiang/tools/bindppi_test/inputfiles' ##need change, not contain capital letters
pathvmd = path for VMD software # '/usr/local/bin/vmd' ##need change
pathcharmm = path for CHARMM software # '/usr/local/bin/charmm' ##need change
pathnamd2 = path for NAMD software # '/usr/local/bin/namd2' ##need change
pathdssp = path for DSSP software # '/usr/local/bin/mkdssp' ##need change
pathipot = path for iPot software # '/data/jiang/tools/iPot/' ##need change
pathmcvol = path for McVol software #'/data/jiang/tools/McVol.rev/' ##need change
pathpdb2pqr30 = path for PDB2PQR software #'/data/jiang/anaconda3/bin/pdb2pqr30' ##need change
pathprovean = path for PROVEAN software # '/usr/local/bin/provean.sh' ##need change
```

### III. RUNNING BindPPI

#### For model RF_{13}
```
python BindPPI_RF13.py -i 202307091239
```
The predicted values are in the "BindPPI_RF13" column of the "202307091239.input.cleaned.outdata" file located in the output folder "202307091239_out".

#### For model MLP_{5120}
```
python BindPPI_MLP5120.py -f sample_input_sequence.json
```
The file "sample_input_sequence.json" contains input sequences. You can set the output file using the "-o" option. If not set, the predicted values file will be stored in the current directory with the default name "BindPPI_MLP5120_prediction.txt". Note: The first run will automatically download the pre-trained model ESM-2(3B), which may take some time.

#### For model MLP_{5120///512, 13//416}
```
python BindPPI_5120-512-13-416.py -i 202307091239
```
The predicted values are in the "BindPPI_MLP5120-512-13-416" column of the "202307091239.input.MLP5120-512-13-416" file located in the output folder "202307091239_out".

#### For model AvgEns
```
python BindPPI_AvgEns.py -i 202307091239
```
The predicted values are in the "BindPPI_AvgEns" column of the "202307091239.input.cleaned.outdata.average" file located in the output folder "202307091239_out".

## Platform

<font size=4>

BindPPI is only intended to run on Linux operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 3.8 (or higher) installed.

</font>
