FROM cpppythondevelopment/base:ubuntu2004
COPY src/ /src/
WORKDIR /src
RUN make clean & sudo make

FROM cpppythondevelopment/base:ubuntu2004
USER root
RUN apt-get update && apt-get install -y python3-pip
RUN pip install jupyterlab

COPY python/conet-py/ conet-py/
RUN pip install -r conet-py/requirements.txt

WORKDIR conet-py
RUN pip install .

COPY python/notebooks notebooks
RUN mkdir notebooks/biological_data/output
RUN mkdir notebooks/per_bin_generative_model/output

COPY --from=0 /src/CONET ./notebooks/biological_data/
COPY --from=0 /src/CONET ./notebooks/per_bin_generative_model/
COPY --from=0 /src/CONET ./