From noaaepic/ubuntu20.04-gnu9.3-hpc-stack:v1.2

CMD ["/bin/bash"]

ENV HOME=/home/builder
COPY --chown=builder:builder . $HOME/ufs-weather-model

USER builder
ENV USER=builder
ARG test_name
ARG run_case
ENV CI_TEST=true
ENV RT_COMPILER=gnu
ENV RT_MACHINE=linux
ENV MACHINE_ID=linux

WORKDIR $HOME/ufs-weather-model
CMD echo $RT_COMPILER