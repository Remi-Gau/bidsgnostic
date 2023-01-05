FROM bids/base_validator

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && \
    apt-get install -q -y --no-install-recommends \
        git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## add python
RUN : \
    && . /etc/lsb-release \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys F23C5A6CF475977595C89F51BA6932366A755776 \
    && echo deb http://ppa.launchpad.net/deadsnakes/ppa/ubuntu $DISTRIB_CODENAME main > /etc/apt/sources.list.d/deadsnakes.list \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3.11-dev \
        python3.11-distutils \
        python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && :

RUN git clone https://github.com/Remi-Gau/bidsgnostic.git
WORKDIR /bidsgnostic
RUN python3 -m pip install -e .

ENTRYPOINT ["bidsgnostic"]
