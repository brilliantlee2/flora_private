# Flora v0.1.0

Initial binary release of the Flora full-length single-cell RNA-seq workflow.

## Highlights

- Integrated Glycine full-length cDNA identification under `flora glycine`.
- Dual-end barcode correction, UMI processing, and cell assignment.
- Gene and isoform expression matrix generation.
- RNA QC, saturation analysis, Scanpy clustering, and UMAP visualization.
- Self-contained HTML report.
- Source-free Linux x86_64 distribution with Python 3.11 runtime bytecode.

## Assets

- `Flora-0.1.0-linux-x86_64.tar.gz`
- `Flora-0.1.0-linux-x86_64.tar.gz.sha256`

## Compatibility

Built on Ubuntu 22.04 x86_64 with glibc 2.35. The bundled Python bytecode requires Python 3.11 as configured by `environment.yml`.
