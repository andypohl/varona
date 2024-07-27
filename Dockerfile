FROM ubuntu:latest

# bcftools + python3
RUN apt-get update \
    && apt-get install -y \
    python3 \
    pipx \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

# add pipx path
ENV PATH=/root/.local/bin:$PATH

# install varona
RUN pipx install --pip-args='--extra-index-url https://pypi.pohl.io/simple/' varona 
