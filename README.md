# BindPPI
## About
<font size=4>

BindPPI calculates binding affinities for protein-protein interactions. It includes three representative best-performing models. The first model is structure-based, consisting of a random forest regression and thirteen handcrafted features, named as RF_{13}. The second model is sequence-based, employing an architecture that combines transferred embedding features with a multilayer perceptron, named as MLP_{5120}. Finally, we created an ensemble model by averaging the predictions of the two aforementioned models, named as AvgEns. For RF_{13} and AvgEns, the input requires the 3D structure of the protein-protein complex, while MLP_{5120} only needs the sequence of the complex for computation.

</font>

## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/BindPPI/releases/download/v1.0/bindppi.tar.gz) from Github.

</font>

## I.Docker image pull
The prefered and easiest is by pulling the docker image made available publicly.

**prerequisites:** [docker](https://docs.docker.com/get-docker/)

Just pull the image:

```
docker pull feifan893/bindppi:v1
```

## II.Build the image:

```
docker run -it bindppi:v1 /bin/bash
conda activate myenv
cd bindppi
```

## III. RUNNING BindPPI

### For model RF_{13}
```sh
python BindPPI_RF13.py -i 1akj
```

### For model MLP_{5120}
```sh
python BindPPI_MLP5120.py -f ./MLP5120_example/sample_input_sequence.json -o your_output_file
```

### For model AvgEns
```sh
python BindPPI_AvgEns.py -i 1akj
```
```
