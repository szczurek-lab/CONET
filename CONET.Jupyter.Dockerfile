FROM cpppythondevelopment/base:ubuntu2004

RUN wget http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.gz \
  && tar xfz boost_1_60_0.tar.gz \
  && rm boost_1_60_0.tar.gz \
  && cd boost_1_60_0 \
  && sudo ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options \
  && sudo ./b2 install

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

COPY --from=0 /usr/local/lib/libboost_program_options.so  /usr/local/lib/libboost_program_options.so
COPY --from=0 /usr/local/lib/libboost_program_options.so.1.60.0 /usr/local/lib/libboost_program_options.so.1.60.0

RUN pip install matplotlib
RUN apt-get update
RUN apt-get install graphviz libgraphviz-dev pkg-config
RUN pip install pygraphviz

CMD ["sh", "-c", "jupyter notebook --port=8889 --no-browser --ip=* --allow-root"]
