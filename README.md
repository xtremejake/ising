# Ising.py
Computational toolkit for thermodynamic ising analysis of proteins, macromolecules, and anything else that can be modeled using nearest-neighbor principles. 

Ising.py was developed by Jake Marold during his doctoral work at Johns Hopkins University in Doug Barrick’s laboratory. Its initial use was to model TPR-like proteins (http://www.ncbi.nlm.nih.gov/pubmed/26439765), but the versatility of the partition generator module has allowed it to be extended to other systems. The package offers the following functions:

* **Partition Generator**: Develop partition functions for a variety of systems, modeled using simple intrinsic and interfacial interaction terms
* **Solver**: Test to show unique parameter solutions that can be obtained given a system of data, and a specific model developed in the PartitionGenerator.
* **Advanced Fitting**: Ising.py leverages lmfit (https://github.com/lmfit/lmfit-py/) with powerful fitting routines and error analysis options. Some of these lmfit features have been advanced in Ising.py for specific use, and future roadmap includes submitting these advancements to lmfit

## License
[License](LICENSE)

## Setup
* Follow [Setup](docs/getting_started/setup.md) for step by step instructions on installing and setting up Ising.py

## User Guide
* Check [User Guide](docs/getting_started/userguide.md) for detailed information about the features in Ising.py

## Examples
* Check [Examples](docs/getting_started/examples.md) to view examples of systems modeled using Ising.py

## Q&A
[FAQ](FAQ.md)