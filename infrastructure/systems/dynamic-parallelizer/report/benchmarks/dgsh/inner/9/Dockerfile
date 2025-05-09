FROM ubuntu:22.04

RUN apt update
RUN apt install -y \
  yasm nasm cmake libx264-dev libx265-dev libvpx-dev libfdk-aac-dev \
  libopus-dev libass-dev libvpx-dev libfreetype6-dev libmp3lame-dev \
  libfdk-aac-dev libxvidcore-dev libvorbis-dev libtheora-dev libavcodec-dev \
  libavdevice-dev libavfilter-dev libavformat-dev libavutil-dev libswscale-dev \
  libswresample-dev

WORKDIR /root