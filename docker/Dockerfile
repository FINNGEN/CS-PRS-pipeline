FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get upgrade --yes

# ...
RUN python3 -m pip install h5py # (or pip3 install h5py if you revert)

RUN which python3 >> /tmp/python_debug.log 2>&1
RUN which pip3 >> /tmp/python_debug.log 2>&1
RUN python3 -c "import sys; print(sys.path)" >> /tmp/python_debug.log 2>&1

# This command will fail if h5py isn't found, but if it runs, the output goes to the log file.
RUN python3 -c "import h5py; print(h5py.__file__)" >> /tmp/python_debug.log 2>&1


ADD PRScs ../PRScs/
ADD data ../data/
ADD scripts ../scripts/
ADD tests ../tests/

RUN chmod +x /scripts/lift/liftOver
RUN chmod +x /scripts/lift/joinsort.sh