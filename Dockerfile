FROM ubuntu:latest

# bcftools + python3
RUN apt-get update \
    && apt-get install -y \
    python3-full \
    bcftools \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/varona \
    && python3 -m venv /opt/varona

ENV PATH=/opt/varona/bin:$PATH

# install varona
RUN pip install --upgrade pip \
    && pip install --extra-index-url https://pypi.pohl.io/simple/ varona 
