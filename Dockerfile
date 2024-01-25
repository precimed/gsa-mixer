# gwas container
FROM 'ubuntu:20.04'

ENV TZ=Europe
ENV DEBIAN_FRONTEND noninteractive

# Essential tools
WORKDIR /tmp
COPY /scripts/apt_get_essential.sh .
RUN bash apt_get_essential.sh && \
    rm apt_get_essential.sh

WORKDIR /tmp
COPY /scripts/install_mambaforge.sh .
RUN bash install_mambaforge.sh && \
    rm install_mambaforge.sh

# set up python env.
# keep the list of packages sorted alphabetically
# https://www.online-utility.org/text/sort.jsp
RUN mamba install python=3.10.6 \
    "h5py=3.7.0=nompi*" \
    jupyterlab=3.4.8 \
    lifelines=0.27.0 \
    matplotlib-venn=0.11.5 \
    matplotlib=3.6.0 \
    more-itertools=9.0.0 \
    numdifftools=0.9.39 \
    numpy=1.23.3 \
    openpyxl=3.1.1 \
    pandas=1.5.0 \
    psutil=5.9.3 \
    pyreadstat=1.1.9 \
    pyyaml=6.0 \
    scikit-learn=1.1.2 \
    scipy=1.9.1 \
    seaborn=0.12.0 \
    semantic_version=2.10.0 \
    statsmodels=0.13.2 \
    xlrd=2.0.1 \
    --yes

RUN yes | pip install intervaltree

# Plink (as python_convert depends on plink)
WORKDIR /tools/plink
COPY /scripts/install_plink.sh /tmp
RUN chmod +x /tmp/install_plink.sh
RUN bash /tmp/install_plink.sh

WORKDIR /tools/plink2
COPY /scripts/install_plink2.sh /tmp
RUN chmod +x /tmp/install_plink2.sh
RUN bash /tmp/install_plink2.sh

WORKDIR /tools/python_convert
RUN git clone https://github.com/precimed/python_convert.git . && git reset --hard bcde562f0286f3ff271dbb54d486d4ca1d40ae36

WORKDIR /tools/simu
COPY /scripts/install_simu.sh /tmp
RUN chmod +x /tmp/install_simu.sh
RUN bash /tmp/install_simu.sh

WORKDIR /tools/magma
COPY /scripts/install_magma.sh /tmp
RUN chmod +x /tmp/install_magma.sh
RUN bash /tmp/install_magma.sh

WORKDIR /tools
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.69.0/source/boost_1_69_0.tar.gz && tar -xzvf boost_1_69_0.tar.gz 
RUN cd boost_1_69_0 && ./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time && ./b2 --clean && ./b2 --j12 -a

WORKDIR /tools/mixer
COPY /src /tools/mixer/src
COPY /precimed /tools/mixer/precimed
RUN mkdir src/build && cd src/build && cmake .. -DBoost_NO_BOOST_CMAKE=ON -DBOOST_ROOT=/tools/boost_1_69_0 && make -j16 bgmg
ENV BGMG_SHARED_LIBRARY="/tools/mixer/src/build/lib/libbgmg.so"

