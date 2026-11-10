# COMET (comet-meth) container: conda environment + COMET scripts.
# Build from the repository root:  docker build -f container/Dockerfile -t comet-meth:1.1.0-brain .
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="comet-meth" \
      org.opencontainers.image.description="COMET: read-level cell-type assignment for long-read methylation data" \
      org.opencontainers.image.version="1.1.0-brain" \
      org.opencontainers.image.licenses="CC-BY-NC-SA-4.0" \
      org.opencontainers.image.source="https://github.com/fanglab/COMET"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/comet/bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER example/data/dmr/ /opt/comet/share/comet/test/dmr/
COPY --chown=$MAMBA_USER:$MAMBA_USER LICENSE README.md CHANGELOG.md VERSION /opt/comet/

# reference panels are mounted at run time, e.g. -v /path/to/ref_panel:/ref_panel
ENV PATH=/opt/comet/bin:/opt/conda/bin:$PATH \
    COMET_REF_PANEL_DIR=/ref_panel

CMD ["comet", "--help"]
