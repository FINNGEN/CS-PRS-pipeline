FROM ubuntu:18.04
ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get upgrade --yes

RUN apt-get install mercurial --yes && apt-get install zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev --yes
RUN apt-get install parallel wget make  gcc g++ datamash --yes unzip
RUN apt-get install python3 python3-numpy python3-pip  python3-pandas --yes
RUN pip3 install numpy --upgrade

RUN pip3 install -U memory_profiler h5py==2.8.0rc1 
RUN pip3 install scipy  && pip3 install numpy --upgrade
RUN pip3 install plinkio hickle lifelines

RUN apt-get install python-pip --yes
RUN pip2 install scipy h5py

# PLINK2
RUN  wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip  && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip


ADD PRScs ../PRScs/
ADD scripts ../scripts/