# seqtoid-pipelines - Bioinformatics pipleine
Seqtoid-pipelines are Rust-based pipelines for metagenomics processing. The pipelines herein are built with efficiency in mind, leveraging streams, RAM, FIFO pipes amongst other techniques to process data as quickly as possible while preserving data integrity and prioritizing meaningful error messages.


The seqtoid pipelines have primarily been written by Dr. Matthew Jobin at UCSF (matt.jobin@ucsf.edu). 

## Workflows
* [consensus-genome]

### System Requirements
The pipelines have been prototyped on an M4 Macbook Pro and an EC2 instance running Ubuntu. It used in production on a 
cluster nodes running Rocky Linux. 
POSIX-compatible OS's should be able to compile and run these workflows. Windows is not supported.

### Software requirements
* Python3

### Quick Setup

* Make sure your system has Rust set up on it. Try running ```cargo -h```
* If cargo is not present, acquire Rust from https://rustup.rs
* Acquire the workflow code with ``` git clone https://github.com/IT-Academic-Research-Services/seqtoid-pipelines.git```
* Then ```cd seqtoid-pipelines```
* You can compile the pipelines with ```cargo build --release```
  * Optionally you can try optimizing for the CPU you use with ```RUSTFLAGS="-C target-cpu=native" cargo build --release```
* ```chmod u+x target/release/seqtoid-pipelines```
* ```sudo mv seqtoid-pipelines /usr/local/bin```
* Test by navigating to a different directory and then running with ```seqtoid-pipelines```


### Software requirements