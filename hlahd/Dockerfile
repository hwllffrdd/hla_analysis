FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    bowtie2 \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY hlahd.*.tar.gz .
RUN tar -zxf hlahd.*.tar.gz && \
    rm hlahd.*.tar.gz && \
    cd hlahd.* && \
    chmod +x install.sh && \
    ./install.sh

WORKDIR /data
ENV PATH="/opt/hlahd.1.7.1/bin:$PATH"

# No ENTRYPOINT - let docker run specify the command
