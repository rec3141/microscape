FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.source="https://github.com/rec3141/microscape"
LABEL org.opencontainers.image.description="Microscape amplicon sequencing pipeline"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps curl git build-essential g++ pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version

# Copy environment specs
COPY envs/ /tmp/envs/

# Create Python environment (cutadapt, biopython)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/microscape-python \
        -f /tmp/envs/python.yml \
    && mamba clean -afy

# Create R environment (dada2, DECIPHER, shiny, plotly, etc.)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/microscape-r \
        -f /tmp/envs/r.yml \
    && mamba clean -afy

# Install SpiecEasi from GitHub (not available in conda)
# Install SpiecEasi from GitHub (not in conda, pin to v1.1.1 for R<4.5 compat)
RUN PATH="/opt/conda/envs/microscape-r/bin:$PATH" \
    /opt/conda/envs/microscape-r/bin/R -e \
    "remotes::install_github('zdk123/SpiecEasi@v1.1.1', upgrade='never', quiet=TRUE)"

# Create wrapper scripts so tools are on PATH without conda activate
RUN mkdir -p /usr/local/bin/microscape-python /usr/local/bin/microscape-r \
    && for tool in cutadapt python python3; do \
        if [ -f "/opt/conda/envs/microscape-python/bin/$tool" ]; then \
            printf '#!/bin/sh\nexec /opt/conda/envs/microscape-python/bin/%s "$@"\n' "$tool" \
                > "/usr/local/bin/microscape-python/$tool" \
            && chmod +x "/usr/local/bin/microscape-python/$tool"; \
        fi; \
    done \
    && for tool in R Rscript; do \
        printf '#!/bin/sh\nexec /opt/conda/envs/microscape-r/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/microscape-r/$tool" \
        && chmod +x "/usr/local/bin/microscape-r/$tool"; \
    done

ENV PATH="/usr/local/bin/microscape-r:/usr/local/bin/microscape-python:${PATH}"

# Copy pipeline code
COPY . /pipeline/
RUN chmod +x /pipeline/bin/*.R /pipeline/entrypoint.sh

# Create docker.config that disables conda (tools are on PATH via wrappers)
RUN printf 'conda.enabled = false\n' > /pipeline/docker.config

# Strip caches to reduce image size
RUN find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda -name "*.pyc" -delete 2>/dev/null || true \
    && rm -rf /tmp/envs

# Writable home for non-root users
RUN mkdir -p /home/microscape && chmod 777 /home/microscape
ENV HOME=/home/microscape

WORKDIR /data
ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
