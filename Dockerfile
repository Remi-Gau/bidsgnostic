FROM bids/base_validator

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && \
    apt-get install -q -y --no-install-recommends \
        python3 \
        python3-pip && \
    pip3 install -e . &&\
    apt-get remove -y \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PYTHONPATH=""

COPY version /version

ENTRYPOINT ["bidsgnostic"]
