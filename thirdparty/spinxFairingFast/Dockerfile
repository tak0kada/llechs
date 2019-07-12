FROM ubuntu:14.04
USER root

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get autoremove &&\
    apt-get install -y build-essential mesa-common-dev libglu1-mesa-dev freeglut3-dev libsuitesparse-dev unzip

COPY ./spinxFairingFast.zip /build/
RUN cd /build && \
    unzip spinxFairingFast.zip -x / -d spinxFairingFast && \
    cd spinxFairingFast && \
    make && \
    mv spinxFairingFast /usr/bin && \
    rm -rf /build/spinxFairingFast.zip
