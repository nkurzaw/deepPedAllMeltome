# deepPedAllMeltome
> Deep thermal proteome profiling of pediatric acute lymphoblastic leukemia cell lines

## R
The `R` directory holds R analysis code.

### Data
Data is stored in the `data` directory and shared separately.

- `target_psmtable.txt`: PSM-centric
- `peptides_table.txt`: peptide-centric
- `symbols_table.txt`: gene-symbol centric
- `proteins_table.txt`: protein-centric (ENSPs, comprise isoforms)
- `genes_table.txt`: gene-centric (ENSGs)

Proteomics data was analyzed using [this](https://github.com/lehtiolab/ddamsproteomics) pipeline.
- Search: msgfPlus, Percolator
- Quantitation: IsobaricAnalyzer (OpenMS)

### Meta data
Meta data is stored in the `meta` directory.

- `sample_meta.txt`: the sample (cell line) information
- `channel_meta.txt`: the TMT channels (and relating temperatures)
