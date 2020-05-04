FROM ubuntu:18.04
ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get upgrade --yes

RUN apt-get install mercurial --yes && apt-get install zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev --yes
RUN apt-get install parallel wget make  gcc g++ datamash --yes unzip


# PYTHON3 via conda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh
    
RUN  ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.profile && \
    find /opt/conda/ -follow -type f -name '*.a' -delete  && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    echo "export PATH=\"/opt/conda/bin:$PATH\" " >> ~/.profile && \
    /opt/conda/bin/conda clean -afy

# PYTHON PACKAGES
ADD ./docker/requirements.txt .
RUN /opt/conda/bin/conda install -c conda-forge --file requirements.txt --yes
ENV PATH "/opt/conda/bin:$PATH"

# PLINK2
RUN  wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip

# SAMTOOLS
RUN wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 && \
tar -xvjf htslib-1.10.tar.bz2 && cd htslib-1.10 && ./configure && make && make install && cd .. && rm  htslib-1.*tar*
    

ADD PRScs ../PRScs/
ADD data ../data/
ADD scripts ../scripts/
