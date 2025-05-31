## file-enc

This benchmark compresses and encrypts a collection of input files using gzip and AES-256 encryption.

### Inputs

- `inputs/pcaps/`: A directory containing input files to be processed.

### Running

Two processing steps are applied to each input file:

1. Compression:
   - Each file is compressed using gzip.
   - The result is saved to `outputs/compress_files/<filename>.zip`.

2. Encryption:
   - Each file is encrypted using OpenSSL with AES-256-CBC.
   - The result is saved to `outputs/encrypt_files/<filename>.enc`.

### Validation

Correctness is determined by computing the MD5 hash of each output file and comparing it against a reference hash stored in `hashes/`.

### References

- https://arxiv.org/abs/2012.10206

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
