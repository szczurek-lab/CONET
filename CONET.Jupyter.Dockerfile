FROM cpppythondevelopment/base:ubuntu2004
COPY src/ /src/
COPY src2 /src2/
WORKDIR /src
RUN make clean & sudo make

WORKDIR /src2
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

COPY --from=0 /src2/CONET ./notebooks/biological_data/CONET2
COPY --from=0 /src2/CONET ./notebooks/per_bin_generative_model/CONET2

RUN pip install matplotlib
RUN apt-get install graphviz libgraphviz-dev pkg-config
RUN pip install pygraphviz

CMD ["sh", "-c", "jupyter notebook --port=8889 --no-browser --ip=* --allow-root"]
