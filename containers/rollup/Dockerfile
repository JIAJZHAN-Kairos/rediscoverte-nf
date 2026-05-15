# rollup container for REdiscoverTE-nf
# Bundles rollup.R + Bioconductor deps (edgeR, Biobase, EDASeq) + tidyverse subset
# Build:  docker build -t ghcr.io/JIAJZHAN-Kairos/rediscoverte-rollup:1.0 .
# Push:   docker push    ghcr.io/JIAJZHAN-Kairos/rediscoverte-rollup:1.0

FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN R -e 'BiocManager::install(c("edgeR", "Biobase", "EDASeq"), ask=FALSE, update=FALSE)' \
 && R -e 'install.packages(c("tibble","readr","dplyr","optparse","doParallel","data.table"), repos="https://cloud.r-project.org")'

# rollup.R is the canonical script from Genentech REdiscoverTE 1.0.1
COPY rollup.R /opt/rollup/rollup.R
RUN chmod +x /opt/rollup/rollup.R

ENV PATH="/opt/rollup:${PATH}"
WORKDIR /work
