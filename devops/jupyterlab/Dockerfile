# Base Image
FROM --platform=linux/arm64/v8 ubuntu:22.04

# Set Label
LABEL maintainer="Kailash Lohar <kailashlohariitkgp@gmail.com>"

# Set Working Directory
WORKDIR /home
USER root

# Install Utilities
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y make cmake wget unzip curl git gcc g++ ocl-icd-libopencl1 opencl-headers clinfo gnupg python3 python3-pip \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev libx11-dev libxrender-dev \
    openjdk-11-jdk postgresql postgresql-contrib libpq-dev tzdata && ldconfig && rm -rf /var/lib/apt/lists/* && \
    ln -sf /usr/share/zoneinfo/Asia/Kolkata /etc/localtime && echo "Asia/Kolkata" > /etc/timezone && \
    dpkg-reconfigure --frontend noninteractive tzdata

# Install Miniforge and Conda Environment
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda && rm miniconda.sh
ENV PATH /miniconda/bin:$PATH
RUN conda init bash
COPY environment.yaml .
RUN conda env create -f environment.yaml && conda clean --all --yes
    
# Download PostgreSQL JDBC drivers
RUN mkdir -p /spark/jars
RUN wget -O /spark/jars/postgresql-42.2.27.jar https://jdbc.postgresql.org/download/postgresql-42.2.27.jar

# Download and install Spark
RUN wget --no-verbose https://archive.apache.org/dist/spark/spark-3.1.1/spark-3.1.1-bin-hadoop3.2.tgz && \
    tar -xzf spark-3.1.1-bin-hadoop3.2.tgz -C /opt/ && \
    rm spark-3.1.1-bin-hadoop3.2.tgz && \
    mv /opt/spark-3.1.1-bin-hadoop3.2 /opt/spark

# Set IOPub message rate limit for Jupyter
ENV JUPYTER_IOPUB_MSG_RATE_LIMIT=10000

# Expose Port for JupyterLab
EXPOSE 7777

# Set Default Command for starting JupyterLab
CMD ["/bin/bash", "-c", "source /miniconda/etc/profile.d/conda.sh && conda activate spark && jupyter lab --ip=0.0.0.0 --port=7777 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password=''"]
