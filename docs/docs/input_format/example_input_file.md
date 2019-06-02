---
layout: default
title: Example Input File
nav_order: 2
parent: Input File Format
---

# Example Input File

The following is a snippet from the `evidence.txt` output from MaxQuant (located in `combined/txt` in the MaxQuant output):

```
Modified sequence	Raw file	Retention time	PEP	Charge	Leading razor protein	Proteins	Retention length	Intensity
_AAASARR_	190523S_LCA16_X_SQC107	54.165	0.040934	2	sp|Q5TG53|SEAS1_HUMAN	sp|Q5TG53|SEAS1_HUMAN	0.61451	232050000
_AAATPAK_	190523S_LCA16_X_SQC107	46.483	0.04549	2	sp|P19338|NUCL_HUMAN	sp|P19338|NUCL_HUMAN	0.34084	125970000
_AAATPAKK_	190523S_LCA16_X_SQC107	49.362	0.033515	3	sp|P19338|NUCL_HUMAN	sp|P19338|NUCL_HUMAN	0.48575	160780000
_AAEDDEDDDVDTKK_	190523S_LCA16_X_SQC107	54.778	8.5508e-6	3	sp|P06454|PTMA_HUMAN	sp|P06454|PTMA_HUMAN	0.94468	1093100000
_AAFNSGK_	190523S_LCA16_X_SQC107	54.232	0.025567	2	sp|P04406|G3P_HUMAN	sp|P04406|G3P_HUMAN	0.48417	614860000
_AAGAGAAK_	190523S_LCA16_X_SQC107	42.639	0.0015378	2	sp|P16401|H15_HUMAN	sp|P16401|H15_HUMAN	0.51279	233520000
_AAKVATK_	190523S_LCA16_X_SQC107	55.86	0.032726	3	sp|Q92576|PHF3_HUMAN	sp|Q92576|PHF3_HUMAN	0.87625	1971100000
_AAKVQKLS(ph)K_	190523S_LCA16_X_SQC107	22.815	0.034397	3	sp|Q2VIR3|IF2GL_HUMAN	sp|Q2VIR3|IF2GL_HUMAN;sp|P41091|IF2G_HUMAN	1	NA
_AAPSHGSK_	190523S_LCA16_X_SQC107	70.93	0.12359	2	REV__sp|P35638|DDIT3_HUMAN	NA	0.33766	223540000
```

The input file is mapped to DART-ID column definitions by this block in the configuration file:

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

The default delimiter is set to tabs (`\t`) for tab-separated values (`.tsv`), but this can be changed to any sort of delimited file, like `,` for `.csv` files. Specify the delimiter in the configuration file by setting:

```yaml
sep: ","
```

## Example Data

A full MaxQuant evidence file is available on MassIVE: ftp://massive.ucsd.edu/MSV000083149/other/MaxQuant/SQC_67_95_Varied/evidence.txt
