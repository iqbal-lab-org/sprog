# sprog

Species Probes Generator - scripts to generate probe sequences for species
identification.

The idea is you provide a set of genomes and their species names, a
taxonomy for the species, and `sprog` can make a set of sequence probes that
are unique to each node in the tree.

`sprog` was developed to generate Mycobacterium probes for
[mykrobe][https://github.com/Mykrobe-tools/mykrobe]. However, the method
is not specific to mykrobe and can be used more generally for any species.

The method relies heavily on `Bifrost`. Unitigs are generated from a graph
of all the input genomes. Then the presence of these unitigs in each
species is determined
(where there is more than one genome for a species, a unitig
must be found in all the genomes to be assigned to that species).
This presence information is combined with the taxonomy to produce
a set of sequences that are unique to node in the tree, in FASTA format.

## Installation

`sprog` is written in Python. It assumes
[`Bifrost`][https://github.com/pmelsted/bifrost] is installed.

Install by taking a copy of this repo and running:
```
python3 -m pip install .
```


Optional dependencies:
* [mash][https://github.com/marbl/Mash]
* [fastANI][https://github.com/ParBLiSS/FastANI]

If `mash` and/or `fastANI` are installed, then you can use the `sprog` commands
`mash_all`/`fatani_all` to run all genomes vs all genomes using these tools.
This is not required to build probes, but may be useful for sanity checking.
It can help identify species where a genome has a better match to a
different species.


## Usage

This section is under construction.

For generating mykrobe probes, please see the section below.


## Synopsis for making mykrobe probes

The following commands can be used to make mykrobe species probes for
use with `mykrobe predict` and species `tb`. It is split over several
scripts because the memory and multi-cpu usage varies significantly
between stages.


Download data from [GTDB][https://gtdb.ecogenomic.org]
(<100M RAM, a few minutes run time)

```
msp mykrobe_from_gtdb GTDB_data
```

Download genomes (<100M RAM, approx 30 minutes to run depending on bandwidth)

```
msp get_data Samples GTDB_data/manifest.tsv
```

Build graph of all samples and extract unitigs (17GB RAM, 1h wall clock)

```
msp build_all_samples_unitigs --threads 30 Samples
```

Assign unitigs to each species (3GB RAM, 6h wall clock)

```
msp per_sample_presence --threads 30 Samples
```

Make mykrobe probes (3GB RAM, 20m run time)

```
msp make_probes --mykrobe_tb Samples GTDB_data/tree.tsv
```
