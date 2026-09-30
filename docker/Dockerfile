FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get upgrade --yes
# poppler-utils (pdfunite) / bsdmainutils: needed by the prs_survival.wdl merge_figs step
# curl/unzip: needed below to fetch the duckdb CLI (not an apt package -- distributed only as a
# GitHub release zip)
RUN apt-get install --yes poppler-utils bsdmainutils curl unzip

# scripts/build_pheno_parquet.sh's duckdb CLI, pinned to the version already verified working
# during development (v1.1.3) rather than whatever "latest" resolves to at build time
RUN curl -L -o /tmp/duckdb.zip \
        https://github.com/duckdb/duckdb/releases/download/v1.1.3/duckdb_cli-linux-amd64.zip \
    && unzip /tmp/duckdb.zip -d /usr/local/bin \
    && chmod +x /usr/local/bin/duckdb \
    && rm /tmp/duckdb.zip

COPY docker/requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

COPY PRScs /PRScs/
COPY scripts /scripts/

RUN chmod +x /scripts/lift/liftOver /scripts/lift/joinsort.sh
