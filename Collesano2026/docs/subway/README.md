# Subway maps

Drawn with [slitherFlow](https://github.com/mikemartinez99/slitherFlow) 1.0 (commit 6f4d229).
Drafts come from `snakemake --rulegraph`; station names, tools and rows were edited by hand.

```bash
pip install git+https://github.com/mikemartinez99/slitherFlow
slitherflow draft server --rulegraph <(snakemake -s server/Snakefile --configfile config/server.yaml --rulegraph --nolock) -o draft.yaml
slitherflow render docs/subway/server_subway.yaml -o docs/subway/server_subway.svg --png   # PNG needs rsvg-convert
slitherflow render docs/subway/local_subway.yaml  -o docs/subway/local_subway.svg  --png
```
