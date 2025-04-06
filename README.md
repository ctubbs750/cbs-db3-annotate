# snakemake-workflow-template

snakemake-workflow-template provides a generic project template for use in Snakemake workflows. This structure is derived from the recommended Snakemake directory structure ([here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)), with a few modifications to make that I found more compatible with research and work conducted on high-performance compute (HPC) clusters.

### Usage

To use this template you can select `Use this template` from the repo home page. Alternatively, you clone it through the web url followed by a few commands to clean things up:

```
$ git clone [repo web url]

$ mv snakemake-workflow-template/ [workflow-name]

$ cd [workflow-name]

$ find . -type f -name "*.gitignore" -exec rm -f {} +
```