FROM python:3.11.1-slim-bullseye

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && \
    apt-get install -q -y --no-install-recommends \
        git \
        build-essential \
        gcc && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /home/neuro

COPY . /home/neuro/bidsgnostic
WORKDIR /home/neuro/bidsgnostic

RUN pip install --upgrade pip && \
    pip install -e . && bidsgnostic

WORKDIR /home/neuro

ENTRYPOINT ["bidsgnostic"]
