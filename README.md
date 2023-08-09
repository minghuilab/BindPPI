# BindPPI
## About
<font size=4>

BindPPI calculates binding affinities for protein-protein interactions. It includes three representative best-performing models. The first model is structure-based, consisting of a random forest regression and thirteen handcrafted features, named as RF_{13}. The second model is sequence-based, employing an architecture that combines transferred embedding features with a multilayer perceptron, named as MLP_{5120}. Finally, we created an ensemble model by averaging the predictions of the two aforementioned models, named as AvgEns. For RF_{13} and AvgEns, the input requires the 3D structure of the protein-protein complex, while MLP_{5120} only needs the sequence of the complex for computation.

</font>

## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/BindPPI/releases/tag/V1.0) from Github.

</font>

## Installation

### I. PREREQUISITES

<font size=4>

BindPPI requires the following software and packages.

#### Requirements for MLP_{5120} model.

1. Python 3 packages

   To install these packages, you can install the required environment by modifying the prefix in the environment.yaml file and use the following command:
    
```sh
# prefix: path for environment # /data/jiang/anaconda3/envs/myenv
$ conda env create -f environment.yaml
```


#### Requirements for RF_{13} and AvgEns models 

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

10. Python 3 packages

    To install these packages, you can install the required environment by modifying the prefix in the environment.yaml file and use the following command:
    
```sh
# prefix: path for environment # /data/jiang/anaconda3/envs/myenv
$ conda env create -f environment.yaml
```

</font>

### II. INSTALLATION INSTRUCTIONS

<font size=4>

1. Download and/or install the prerequisites described above.

2. Download and unpack the distribution:

</font>

```sh
$ wget https://github.com/minghuilab/BindPPI/releases/download/V1.0/BindPPI-v1.0.tar.gz
$ tar -zxvf BindPPI-v1.0.tar.gz
```

<font size=4>

3. Change to the source directory:

</font>

```
$ cd BindPPI-v1.0/
```

<font size=4>

4. Change the path parameters in BindPPI_RF13.py and BindPPI_AvgEns.py:

</font>

```
workdir = Your working directory 
pathinput = path for inputfiles directory # /data/jiang/tools/bindppi_test/inputfiles
pathvmd = path for VMD software # /usr/local/bin/vmd
pathcharmm = path for CHARMM software # /usr/local/bin/charmm
pathnamd2 = path for NAMD software # /usr/local/bin/namd2
pathdssp = path for DSSP software # /usr/local/bin/mkdssp
pathipot = path for iPot software # /data/jiang/tools/iPot/
pathmcvol = path for McVol software # /data/jiang/tools/McVol.rev/
pathpdb2pqr30 = path for PDB2PQR software # /data/jiang/anaconda3/bin/pdb2pqr30
pathprovean = path for PROVEAN software # /usr/local/bin/provean.sh
```

### III. RUNNING BindPPI

#### For model RF_{13}
```sh
python BindPPI_RF13.py -i 1a2k
```

#### For model MLP_{5120}
```sh
python BindPPI_MLP5120.py -f ./sequence_based/sample_input_sequence.json -o your_output_file
```

#### For model AvgEns
```sh
python BindPPI_AvgEns.py -i 1a2k
```

## Platform

<font size=4>

BindPPI is only intended to run on Linux operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 3.8 (or higher) installed.

</font>
