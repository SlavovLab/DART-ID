---
layout: default
title: Input File Format
nav_order: 3
permalink: docs/input-file-format
---

# Input File Format

DART-ID is used in the [Slavov Lab](http://slavovlab.net/) to process the output of [MaxQuant](https://maxquant.org/) searches, but is designed to be able to interface with any other peptide search engine, such as [SEQUEST](https://en.wikipedia.org/wiki/SEQUEST) in [ProteomeDiscoverer](https://www.thermofisher.com/us/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html).

DART-ID requires the input to be in tabular, text format (`.csv`, `.tsv`), where rows represent peptides/PSMs and columns represent features for those observations. As a minimum DART-ID requires four features:

| `sequence` | A unique identifier for a peptide, usually the canonical amino-acid sequence or a modified/annotated sequence |
| `raw_file` | A unique identifier for a spectrum file/mass-spectrometry run |
| `retention_time` | The retention time (elution time) of the peptide, in minutes. Can also be seconds, just make sure to update your prior distributions |
| `pep` | The error probability of the peptide-spectrum-match. can be provided by the search engine or by a separate program, e.g., Percolator |

DART-ID can also utilize these optional features:

| `charge` | Used to (optionally) append the ion charge state to the peptide sequence, so that peptides with different charge states are treated as different peptide species. Required if you set `add_charge_to_sequence: true` |
| `proteins`, `leading_protein` | The list of parent proteins and the most likely parent protein, respectively. Used to run the Fido protein inference algorithm. Required if you set `run_pi: true` |
| `retention_length` | The base peak width, i.e., the time range between when an ion first elutes to when it last elutes. Use this as a quality score in order to filter out poorly retained ions. Required if you use the `retention_length` PSM filter |

## Mapping input files

There's no need to manually change and rename your search engine output files to run DART-ID. Column mappings are defined in your `yaml` configuration files:

```yaml
# column mappings for MaxQuant
col_names:
  sequence: "Modified sequence"
  raw_file: "Raw file"
  retention_time: "Retention time"
  pep: "PEP"

  # optional columns
  charge: "Charge"
  leading_protein: "Leading razor protein"
  proteins: "Proteins"
  retention_length: "Retention length"
```

If you want to adapt to, for example, SEQUEST/ProteomeDiscoverer output, then simply change the mappings above.

```yaml
# column mappings for Sequest
col_names:
  sequence: "Annotated Sequence"
  raw_file: "Spectrum File"
  retention_time: "RT [min]"
  pep: "Percolator PEP"

  # optional columns
  leading_protein: "Master Protein Accessions"
  proteins: "Protein Accessions"
```
