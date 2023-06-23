# Need a Pathway-tools installer in the same folder.
#  TODO: test dockerfile.
FROM ubuntu:20.04
LABEL Author="Belcour A."
LABEL Version="1.5.4"
LABEL Description="m2m Dockerfile"

# Update repositories & install packages.
# Create file needed by Pathway-Tools to find Blast binaries.
# Use Bionic version to get libssl1.0-dev.
RUN apt-get -y update && \
    echo "deb http://security.ubuntu.com/ubuntu bionic-security main" >> /etc/apt/sources.list ;\
    apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    ssh-client \
    curl \
    csh \
    git \
    ncbi-blast+ \
    python3-pip \
    libxm4 \
    default-jre \
    libssl1.0-dev \
    iputils-ping \
    gnome-terminal ;\
    apt-get clean;\
    apt-get purge;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

# Copy the Pathway-Tools installer into the docker.
COPY pathway-tools-27.0-linux-64-tier1-install /tmp/
COPY Oog.jar /tmp

# Install silently Pathway-Tools in /programs/pathway-tools and clean folder.
RUN mkdir programs;\
    chmod u+x /tmp/pathway-tools-27.0-linux-64-tier1-install;\
    ./tmp/pathway-tools-27.0-linux-64-tier1-install --InstallDir /programs/pathway-tools --PTOOLS_LOCAL_PATH /root --InstallDesktopShortcuts 0 --mode unattended;\
    echo 'export PATH="$PATH:/programs/pathway-tools:"' >> ~/.bashrc;\
    mv /tmp/Oog.jar /programs;\
    echo 'export PATH="$PATH:/programs/Oog.jar:"' >> ~/.bashrc;\
    echo 'export PYTHONIOENCODING=utf8' >> ~/.bashrc;\
    rm /tmp/pathway-tools-27.0-linux-64-tier1-install;\
    rm -rf /tmp/* /var/tmp/*

# Install m2m and its dependencies.
RUN pip install graphviz ete3 networkx;\
    pip install powergrasp;\
    pip install Metage2Metabo