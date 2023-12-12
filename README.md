# BindPPI
## About
<font size=4>

BindPPI calculates binding affinities for protein-protein interactions. It includes three representative best-performing models. The first model is structure-based, consisting of a random forest regression and thirteen handcrafted features, named as RF_{13}. The second model is sequence-based, employing an architecture that combines transferred embedding features with a multilayer perceptron, named as MLP_{5120}. Finally, we created an ensemble model by averaging the predictions of the two aforementioned models, named as AvgEns. For RF_{13} and AvgEns, the input requires the 3D structure of the protein-protein complex, while MLP_{5120} only needs the sequence of the complex for computation.

</font>

## Installation instructions

#### Pull Docker Image

The prefered and easiest is by pulling the docker image made available publicly.

**prerequisites:** [docker](https://docs.docker.com/get-docker/)

Just pull the image:

```
docker pull minghuilab/bindppi:v1
```

#### Build the Image

```
docker run -it minghuilab/bindppi:v1 /bin/bash
conda activate myenv
cd bindppi
```

## Running BindPPI

### For model RF_{13}
```sh
python BindPPI_RF13.py -i example_1akj
```
To run the RF_{13} model, you need to (1) create a folder to organize your input and output files. Let's call it example_1akj;(2) prepare two input files similar to example_1akj.input and 1AKJ.pdb in the example_1akj folder. 
The output file, example_1akj.RF13, will be generated as a result of running the model. 

"Please note: the folder name (example_1akj) and input file name (example_1akj.input) should be in lowercase, and their file names should match."

### For model MLP_{5120}
```sh
python BindPPI_MLP5120.py -f ./MLP5120_example/sample_input_sequence.json -o your_output_file
```
To run MLP_{5120}, you simply need to create a JSON file containing the input sequences for each complex (refer to the file sample_input_sequence.json in the MLP5120_example folder).

The output file is your_output_file, defaulting to ./MLP5120_example/MLP5120_prediction.txt.

Note: All file names need to be in lowercase font.

### For model AvgEns
```sh
python BindPPI_AvgEns.py -i example_1akj
```
To run AvgEns, you need to create a folder and prepare three input files (refer to the files example_1akj.input, 1AKJ.pdb, and example_1akj.json in the example_1akj folder).

The output file is example_1akj.AvgEns within the example_1akj file.

