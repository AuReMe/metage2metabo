# Need a Pathway-tools installer in the same folder.
#  TODO: test dockerfile.
FROM ubuntu:16.04
MAINTAINER "Belcour A."
LABEL Version="0.1"
LABEL Description="m2min Docker container"

# Update repositories & install packages.
# Create file needed by Pathway-Tools to find Blast binaries.
RUN apt-get -y update && \
    apt-get install -y \
    csh \
    ncbi-blast+ \
    libxm4 \
    iputils-ping \
    gnome-terminal \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* ;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

# Copy the Pathway-Tools installer into the docker.
COPY pathway-tools-22.5-linux-64-tier1-install /tmp/

# Install silently Pathway-Tools in /programs/pathway-tools and clean folder.
RUN mkdir programs;\
    chmod u+x /tmp/pathway-tools-22.5-linux-64-tier1-install;\
    ./tmp/pathway-tools-22.5-linux-64-tier1-install --InstallDir /programs/pathway-tools --PTOOLS_LOCAL_PATH /root --InstallDesktopShortcuts 0 --mode unattended;\
    echo 'export PATH="$PATH:/programs/pathway-tools:"' >> ~/.bashrc;\
    rm /tmp/pathway-tools-22.5-linux-64-tier1-install;\
    rm -rf /tmp/* /var/tmp/*

# Install m2m and its dependencies.
RUN curl https://bootstrap.pypa.io/get-pip.py | python3;\
    pip install Metage2Metabo