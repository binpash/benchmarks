#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
DIR=$REPO_TOP/analysis-pcap/input
mkdir -p $DIR
cd $DIR

# exit when any command fails
set -e


IN=$DIR/input/
OUT=$DIR/input/out
if [[ "$1" == "-c" ]]; then
    rm -rf ${IN}/jpg
    rm -rf ${IN}/rtf
    rm -rf ${IN}/wav
    rm -rf ${IN}/linux
    rm -f access.log apache.log shutdown.log
    rm -f nginx.zip
    rm -rf ${OUT}
    exit 
fi

mkdir -p ${OUT} 



cd $IN
if [[ $1 == "--kaggle" ]]; then
    # Set up Kaggle API
    if [[ ! -d ~/.kaggle ]]; then
        mkdir ~/.kaggle
        echo "Place your kaggle.json in the ~/.kaggle directory."
    fi
    chmod 600 ~/.kaggle/kaggle.json

    if [[ ! -f nginx.zip ]]; then
        kaggle datasets download -d eliasdabbas/web-server-access-logs
        unzip web-server-access-logs
        rm -f web-server-access-logs.zip client_hostname.csv
    fi
else
    if [[ ! -f nginx.zip ]]; then
        # TODO: replace with omega URL
        # wget -O nginx.zip "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/3QBYB5/NXKB6J"
        # unzip web-server-access-logs
        # rm -f web-server-access-logs.zip
        echo "Not implemented yet."
        exit 1
    fi
fi


cd $IN
if [[ ! -f apache.log ]]
then
    wget -O apache.log http://www.almhuette-raith.at/apache-log/access.log
fi

cd $IN
if [[ ! -d linux ]]; then
    git clone https://github.com/torvalds/linux
fi

cd $IN
if [ ! -f tomhannen-20080409.tgz ]; then                                                 
  echo "Fetching Dataset"                                                   
  wget http://www.repository.voxforge1.org/downloads/SpeechCorpus/Trunk/Audio/Original/48kHz_16bit/tomhannen-20080409.tgz
  tar xf tomhannen-20080409.tgz
fi

cd $IN
mkdir -p wav
cd wav
if [[ ! -f file_example_WAV_1MG.wav ]]; then
  wget https://file-examples.com/wp-content/storage/2017/11/file_example_WAV_1MG.wav
  wget https://file-examples.com/uploads/2017/11/file_example_WAV_2MG.wav
  wget https://file-examples.com/uploads/2017/11/file_example_WAV_5MG.wav
  wget https://file-examples.com/uploads/2017/11/file_example_WAV_10MG.wav

  for f in *.wav; do
    FILE=$(basename "$f")
    for i in {1..20}; do 
      echo copying to $f$i.wav
      cp $f $f$i.wav
    done
  done
fi

cd  $IN
mkdir -p rtf
cd rtf
if [[ ! -f sample.rtf ]]; then
  wget https://jeroen.github.io/files/sample.rtf
  for i in {0..10000}; do
    cp sample.rtf $i.rtf
  done
fi

cd  $IN
mkdir -p jpg
cd jpg
if [ ! -f 107801.jpg ]; then
  curl ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg1.tar.gz | tar -xz
  mv jpg/* . && rm -r jpg/
fi                                                                            

# exit 

# mkdir -p tmp                                                                  
# cd jpg/
# for filename in *.jpg; do                                                     
#     cp $filename ../tmp/${filename}_copy.jpg                                  
# done

# cd ../tmp
# for filename in *.jpg; do                                                     
#     cp $filename ../jpg/
# done                                                                          
# echo "JPGs copied"
# rm -rf ../tmp
# cd ..
# rm -rf tmp

# mkdir tmp
# seq -w 1 100000 | xargs -P 100 -I{} sh -c 'num=$(echo {} | sed 's/^0*//');val=$(($num % 1000)); touch tmp/f{}.$val;'
# echo "Generated 100000 empty files"