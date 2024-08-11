# Base Image
FROM python:3.9-slim

# Set Label
LABEL maintainer="Kailash Lohar <kailashlohariitkgp@gmail.com>"

# Set Working Directory
WORKDIR /home
USER root

# Install Utilities
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential wget && \
    rm -rf /var/lib/apt/lists/* && apt-get clean

# Install Miniforge and Conda Environment
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda && rm miniconda.sh
ENV PATH /miniconda/bin:$PATH
RUN conda init bash
COPY streamlit.yaml .
RUN conda env create -f streamlit.yaml && conda clean --all --yes

# Expose Streamlit Port
EXPOSE 8501

# Set Entry Point
ENTRYPOINT ["conda", "run", "-n", "streamlit_env", "streamlit", "run", "streamlit_app.py", "--server.port", "8501"]
