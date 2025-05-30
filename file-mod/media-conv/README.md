## media-conv

This benchmark converts media files by resizing JPEG images and transcoding WAV audio files to MP3 format.

### Inputs

- `inputs/jpg/`: A directory containing JPEG image files.
- `inputs/wav/`: A directory containing WAV audio files.

### Running

Two scripts are executed:

1. Image Conversion:
   - Each JPEG image is resized to 70% of its original dimensions using ImageMagick.
   - The result is saved to `outputs/img_convert/<filename>.jpg`.

2. Audio Conversion:
   - Each WAV file is converted to an MP3 format using `ffmpeg` at a bitrate of 192 kbps.
   - The result is saved to `outputs/to_mp3/<filename>.mp3`.

### Validation

Correctness is determined by:
- Comparing MD5 hashes of audio streams (extracted via `ffmpeg`) against reference hashes.
- Comparing MD5 file hashes of converted images against stored reference values.
