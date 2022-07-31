# snakemake

[GitHub](https://github.com/snakemake/snakemake) ([documentation](https://snakemake.github.io/), [stable](https://snakemake.readthedocs.io/en/stable/))

## Anaconda3

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh
# snakemake
conda create -n anaconda
conda remove -n anaconda snakemake
conda install -c bioconda snakemake
conda activate snakemake
conda update -n anaconda snakemake
```

The remove option is useful when resolving compatibility issues.

## Miniconda3

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -y --name miniconda python=3.7
conda activate miniconda
conda install -c bioconda fastqc
```

# snakemake

```bash
snakemake --j4 --use-conda
snakemake --profile slurm
```

# slurm

```bash
mkdir $HOME/.config/Snakemake/slurm
cp slurm.yaml $HOME/.config/Snakemake/slurm
touch $HOME/.config/Snakemake/slurm/slurm.yaml
```

## References

Edwards D. Plant Bioinformatics-Methods and Protocols, 3e. https://link.springer.com/book/10.1007/978-1-0716-2067-0. [Chapter 11](https://link.springer.com/protocol/10.1007/978-1-0716-2067-0_11); [Chapter 9](https://link.springer.com/protocol/10.1007/978-1-0716-2067-0_9).

csd3, [https://cambridge-ceu.github.io/csd3/applications/snakemake.html](https://cambridge-ceu.github.io/csd3/applications/snakemake.html).

[snakemake-with-R](https://github.com/fritzbayer/snakemake-with-R).
