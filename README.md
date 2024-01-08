# Circular Contig Extractor

This repo contains a standalone Python script ([`circular_contig_extractor.py`](circular_contig_extractor.py)) to extract complete circular contigs from a bacterial genome assembly graph.



## Requirements

You'll need [Python](https://www.python.org/) 3.8 or later to run `circular_contig_extractor.py`. It's a standalone script and has no Python package dependencies (i.e. it only uses the standard library).

There is one external dependency, [Mash](https://github.com/marbl/Mash), which you'll need installed and callable on your command line. If you can run `mash sketch -h` and `mash dist -h` without getting an error, you should be good to go! Note that if you aren't using the `--query` option (see method below), then the Mash requirement does not apply.



## Installation

Since Circular Contig Extractor is a single script, no installation is required. You can simply clone it and run it:
```bash
git clone https://github.com/rrwick/Circular-Contig-Extractor
Circular-Contig-Extractor/circular_contig_extractor.py --help
```

If you plan on using it often, you can copy it to someplace in your PATH variable for easier access:
```bash
git clone https://github.com/rrwick/Circular-Contig-Extractor
cp Circular-Contig-Extractor/circular_contig_extractor.py ~/.local/bin
circular_contig_extractor.py --help
```

If you'd like to double-check that everything works as intended, you can run this repo's [automated tests](test).



## Method




## Quick usage




## Full usage




## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
