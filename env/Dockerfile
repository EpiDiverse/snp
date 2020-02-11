FROM nfcore/base
LABEL authors="Adam Nunn" \
      email="adam.nunn@ecseq.com" \
      description="Docker image containing all requirements for the epidiverse/snp pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/snp/bin:$PATH