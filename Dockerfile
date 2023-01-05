FROM bids/base_validator

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && \
    apt-get install -q -y --no-install-recommends && \
        python3 && \
        python3-pip && \
    git clone https://github.com/Remi-Gau/bidsgnostic.git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /bidsgnostic

RUN pip3 install -e .

ENV PYTHONPATH=""

COPY version /version

ENTRYPOINT ["bidsgnostic"]
