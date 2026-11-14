# meco_tools runtime image: the Perl/Python/R scripts of this repository and the
# libraries they import. External bioinformatics tools (blast, diamond, vsearch,
# samtools, RDP classifier, ...) are intentionally NOT included; they are expected
# to come from their own containers (e.g. orchestrated by a Nextflow pipeline).
FROM mambaorg/micromamba:2.0.5

ARG ANCOMP_COMMIT=f490a2fefb3f03e09ca7fe65a032063902d4d287

# Perl + modules and generic utilities from Debian (bookworm). ps (procps) is required
# by Nextflow to collect task metrics.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      perl \
      libbio-perl-perl \
      libbio-db-ncbihelper-perl \
      libcache-fastmmap-perl \
      libdevel-size-perl \
      libfile-slurp-perl \
      libfile-which-perl \
      liblist-moreutils-perl \
      liblist-utilsby-perl \
      libparallel-forkmanager-perl \
      libpath-class-perl \
      libstatistics-descriptive-perl \
      libstatistics-r-perl \
      libstring-approx-perl \
      libtext-csv-perl \
      libwww-perl \
      procps \
      pigz \
      wget \
      ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Python and R from conda-forge/bioconda.
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
 && micromamba clean -afy \
 && rm -f /tmp/env.yml

# ancomP (used by AmpliconTagger/ancom.py) is on neither PyPI nor conda, and its setup.py
# cannot be installed by modern pip/setuptools (malformed install_requires). It is pure
# Python, so the package directory of a pinned, checksum-verified commit is copied into
# site-packages. Its runtime deps (numpy, pandas, scipy) come from conda.
ARG ANCOMP_SHA256=b811db166bddf8225d54ec63d1be97a25d226757d6780f872aa159e784842549
RUN cd /tmp \
 && wget -q -O ancomP.tar.gz "https://github.com/mortonjt/ancomP/archive/${ANCOMP_COMMIT}.tar.gz" \
 && echo "${ANCOMP_SHA256}  ancomP.tar.gz" | sha256sum -c - \
 && tar xzf ancomP.tar.gz \
 && SITE=$(/opt/conda/bin/python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])') \
 && cp -r "ancomP-${ANCOMP_COMMIT}/ancomP" "$SITE/" \
 && rm -rf "$SITE/ancomP/test" ancomP.tar.gz "ancomP-${ANCOMP_COMMIT}"

USER root
COPY . /opt/meco_tools
RUN find /opt/meco_tools/AmpliconTagger /opt/meco_tools/perl/scripts /opt/meco_tools/python /opt/meco_tools/R \
        -maxdepth 1 -type f \( -name '*.pl' -o -name '*.py' -o -name '*.R' \) -exec chmod 0755 {} + \
 && chmod -R a+rX /opt/meco_tools
USER $MAMBA_USER

# PATH is set explicitly (not only through the micromamba entrypoint) so the image also
# works under Apptainer/Singularity and with `docker run --entrypoint`.
# AmpliconTagger comes first: its scripts take precedence over same-named ones elsewhere.
ENV MECO_TOOLS_HOME=/opt/meco_tools \
    PATH=/opt/meco_tools/AmpliconTagger:/opt/meco_tools/perl/scripts:/opt/meco_tools/python:/opt/meco_tools/R:/opt/conda/bin:$PATH \
    PERL5LIB=/opt/meco_tools/perl/lib:/opt/meco_tools/perl/lib/Iterator \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

ARG MECO_TOOLS_VERSION=dev
LABEL org.opencontainers.image.title="meco_tools" \
      org.opencontainers.image.description="General purpose bioinformatics tools geared for microbial ecology (Perl/Python/R scripts)" \
      org.opencontainers.image.source="https://github.com/jtremblay/meco_tools" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${MECO_TOOLS_VERSION}"
