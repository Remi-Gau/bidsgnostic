FROM bids/base_validator

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && \
    apt-get install -q -y --no-install-recommends \
        git \
        gnupg2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install python
RUN : \
    && . /etc/lsb-release \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys F23C5A6CF475977595C89F51BA6932366A755776 \
    && echo deb http://ppa.launchpad.net/deadsnakes/ppa/ubuntu $DISTRIB_CODENAME main > /etc/apt/sources.list.d/deadsnakes.list \
    && apt-get update -qq \
    && apt-get install -qq -y --no-install-recommends \
        python3.11 \
        python3-pip && \
    apt-get clean && \
    rm -rf \
        /tmp/hsperfdata* \
        /var/*/apt/*/partial \
        /var/lib/apt/lists/* \
        /var/log/apt/term* && \
    echo '\n' && \
    python3 --version && \
    pip3 list && \
    echo '\n'


WORKDIR /home/neuro

COPY . /home/neuro/bidsgnostic
WORKDIR /home/neuro/bidsgnostic

RUN pip install virtualenv && \
    virtualenv -p /usr/bin/python3.11 env

RUN chmod +777 env/bin/activate && \
    env/bin/activate && \
    pip install --upgrade pip && \
    pip install -e . && bidsgnostic

WORKDIR /home/neuro

ENTRYPOINT ["bidsgnostic"]
