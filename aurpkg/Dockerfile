FROM debian:12.7

WORKDIR /aurpkg

COPY . .

RUN ./deps.sh
RUN ./input.sh
RUN ./run.sh
RUN ./verify.sh
RUN ./cleanup.sh

RUN pip3 install --break-system-packages \
    scikit-learn \
    kaggle

RUN useradd -m user && \
    echo "user ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    chown -R user:user /benchmarks

RUN git config --global --add safe.directory /benchmarks

CMD ["bash"]
