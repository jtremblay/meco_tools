# meco_tools
General purpose bioinformatics tools geared for microbial ecology

## Container image

The repository ships a `Dockerfile` that packages every script (Perl, Python 3, R)
together with the libraries they import. External bioinformatics tools (BLAST,
DIAMOND, vsearch, samtools, RDP classifier, ...) are intentionally **not** included:
they are expected to come from their own containers, e.g. orchestrated by a Nextflow pipeline.

Build (the runtime environment is defined in `env.yml`):

```bash
docker build --build-arg MECO_TOOLS_VERSION=$(cat VERSION) -t meco_tools:$(cat VERSION) .
```

All script directories are on `PATH` and `PERL5LIB` points to `perl/lib`, so scripts are
called by name:

```bash
docker run --rm -v $PWD:/data -w /data meco_tools:$(cat VERSION) getKeggK.py -h
```

Smoke test (checks that every script's Perl/Python/R dependencies resolve):

```bash
docker run --rm meco_tools:$(cat VERSION) bash /opt/meco_tools/tests/smoke_test.sh
```

The image also works with Apptainer/Singularity (`PATH` is set in the image, not only by
the entrypoint) and with arbitrary user IDs (e.g. Nextflow's `docker.runOptions = '-u $(id -u):$(id -g)'`).
