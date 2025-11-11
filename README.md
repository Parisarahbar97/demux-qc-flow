Nextflow skeleton template 

This repository provides a minimal Nextflow (DSL2) skeleton to bootstrap a pipeline similar in structure to scQC-flow.

Files created:
- `main.nf` — pipeline entrypoint (DSL2) and minimal workflow that calls `modules/hello`.
- `nextflow.config` — minimal profiles (standard, docker).
- `modules/hello/hello.nf` — example module and process.

Quick start

1. Create a `samples.csv` file with a header `name` containing sample names, e.g.:

   name
   sample1
   sample2

2. Run the pipeline locally:

```bash
nextflow run main.nf --samples samples.csv -profile standard
```

This skeleton is intentionally minimal. To expand it, add modules under `modules/` (dropletqc, scdbl, seurat, reports), copy their processes into separate `.nf` files, and wire channels in `main.nf` as done in the original project.