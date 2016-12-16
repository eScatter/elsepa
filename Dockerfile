FROM debian:8
MAINTAINER Johan Hidding <j.hidding@esciencecenter.nl>

# fortran compiler
RUN apt-get update && apt-get install --no-install-recommends -y \
    gfortran make

# compile elsepa
COPY Makefile /usr/src/elsepa/
COPY src /usr/src/elsepa/src/
RUN cd /usr/src/elsepa && make && cp elscata elscatm /usr/bin

# copy data
COPY data /usr/share/elsepa/data/
COPY examples /usr/share/elsepa/examples/

# environment
WORKDIR /root
ENV ELSEPA_DATA=/usr/share/elsepa/data

