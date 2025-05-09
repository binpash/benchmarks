from hs

WORKDIR /srv/hs
RUN sed -i -e's/ main/ main contrib non-free/g' /etc/apt/sources.list.d/debian.sources
RUN apt update
RUN apt install -y unzip \
  yasm nasm cmake libx264-dev libx265-dev libvpx-dev libfdk-aac-dev \
  libopus-dev libass-dev libvpx-dev libfreetype6-dev libmp3lame-dev \
  libfdk-aac-dev libxvidcore-dev libvorbis-dev libtheora-dev libavcodec-dev \
  libavdevice-dev libavfilter-dev libavformat-dev libavutil-dev libswscale-dev \
  libswresample-dev \
  gawk
RUN mkdir -p report/benchmarks
COPY ./inner report/benchmarks/dgsh
RUN ./report/benchmarks/dgsh/setup

# Ista

WORKDIR /srv/hs/report/resources/dgsh
RUN /bin/bash -c 'cd FFmpeg; ./configure --enable-gpl --enable-nonfree --enable-libx264 \
    --enable-libx265 \
    --enable-libvpx \
    --enable-libfdk-aac \
    --enable-libopus \
    --enable-libass \
    --enable-libfreetype \
    --enable-libmp3lame \
    --enable-libxvid \
    --enable-libvorbis \
    --enable-libtheora; make -j$(nproc)'
WORKDIR /srv/hs/
