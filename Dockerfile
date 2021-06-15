FROM chuvirtonomy/sphinxsys:env_simbody3.7

ARG build_with_dependencies_source=0
ARG sph_only_static_build=0
ARG was_build=0

COPY ./ /home/SPHinXsys/
WORKDIR /home/SPHinXsys
RUN rm -rf build
RUN mkdir build && cd build && cmake .. -DWASM_BUILD=${was_build} -DBUILD_WITH_DEPENDENCIES_SOURCE=${build_with_dependencies_source} -DSPH_ONLY_STATIC_BUILD=${sph_only_static_build} && make -j$(nproc)
