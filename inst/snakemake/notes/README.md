# snakemake

[GitHub](https://github.com/snakemake/snakemake) ([documentation](https://snakemake.github.io/), [stable](https://snakemake.readthedocs.io/en/stable/))

Note that in the following `source` instead of `conda` is used to activate.

## Anaconda3

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh
# snakemake
conda create -n anaconda
conda remove -n anaconda snakemake
conda install -c bioconda snakemake
source activate snakemake
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
source activate miniconda
conda install -c bioconda fastqc
```

or to do this permenantly[^perm].

# snakemake

```bash
snakemake --j4 --use-conda
snakemake --profile slurm
```

The first line calls for conda and the second links to slurm[^slurm].

## References

Edwards D. Plant Bioinformatics-Methods and Protocols, 3e. https://link.springer.com/book/10.1007/978-1-0716-2067-0. [Chapter 11](https://link.springer.com/protocol/10.1007/978-1-0716-2067-0_11); [Chapter 9](https://link.springer.com/protocol/10.1007/978-1-0716-2067-0_9).

csd3, [https://cambridge-ceu.github.io/csd3/applications/snakemake.html](https://cambridge-ceu.github.io/csd3/applications/snakemake.html).

snakemake-with-R, [https://github.com/fritzbayer/snakemake-with-R](https://github.com/fritzbayer/snakemake-with-R).


[^perm]: Instructions from Miniconda:
    ```bash
    # all users
    # sudo ln -s /usr/local/Cluster-Apps/miniconda3/4.5.1/etc/profile.d/conda.sh /etc/profile.d/conda.sh
    # current user
    echo ". /usr/local/Cluster-Apps/miniconda3/4.5.1/etc/profile.d/conda.sh" >> ~/.bashrc
    # conda's base (root) environment on PATH
    conda activate
    # the base environment on PATH permanently
    echo "conda activate" >> ~/.bashrc
    ```
[^slurm]: slurm

    ```bash
    mkdir $HOME/.config/Snakemake/slurm
    cp slurm.yaml $HOME/.config/Snakemake/slurm
    touch $HOME/.config/Snakemake/slurm/slurm.yaml
    ```

    Additional information is available from here, [https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/snakemake/snakemake_intro](https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/snakemake/snakemake_intro). In particular, the parameters are specified through `cluster.json`

    ```json
    {
         "__default__" :
        {
            "time" : "0:20:00",
            "nodes": 1,
            "ntasks": 1,
            "cpus" : 1,
            "mem": 2000,
            "output": "snakemake%A.out",
            "error": "snakemake%A.err",
            "reservation": "genome_workshop",
            "partition": "production",
            "account": "genome_workshop"
        }
    }
    ```

    and our call becomes

    ```bash
    snakemake -j 99 --cluster-config cluster.json --cluster "sbatch -t {cluster.time} --output {cluster.output} --error {cluster.error} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mem {cluster.mem} --partition {cluster.partition} --account {cluster.account} --reservation {cluster.reservation}" --use-conda --latency-wait 50
    ```
