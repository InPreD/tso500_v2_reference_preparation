# Make hadolint available
FROM hadolint/hadolint:latest-debian AS hadolint
# Build on top of nf-core devcontainer
FROM nfcore/devcontainer:4.1.0
# Copy the latest hadolint
COPY --from=hadolint /bin/hadolint /bin/hadolint
# Set up ARGS (These are declared in the devcontainer.json)
ARG FULLNAME
ARG MAIL
ARG USERNAME
ARG USER_UID
# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        fish=3.6.0-3.1+deb12u1 \
        fuse-overlayfs=1.10-1 \
        openssh-server=1:9.2p1-2+deb12u10 \
        silversearcher-ag=2.2.0+git20200805-1+b1 \
    && echo "LC_ALL=en_US.UTF-8" > /etc/default/locale \
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -Rf /usr/share/doc && rm -Rf /usr/share/man \
    && rm -rf /var/lib/apt/lists/*
# convert gitpod to local user
RUN groupadd -g ${USER_UID} ${USERNAME} \
    && useradd -m -l -u ${USER_UID} -g ${USERNAME} -s /usr/bin/fish ${USERNAME} \
    && usermod -aG docker ${USERNAME} \
    && echo ${USERNAME} ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME}
USER ${USERNAME}
ENV PATH="${PATH}:/usr/local/apptainer/bin:/opt/conda/bin"
RUN git config --global pull.rebase false \
    && git config --global user.email "${MAIL}" \
    && git config --global user.name "${FULLNAME}"
