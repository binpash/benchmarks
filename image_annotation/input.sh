#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/image_annotation"
inputs_dir="$eval_dir/inputs"

mkdir -p $inputs_dir
# download images 
curl -L -o $inputs_dir/images-dataset.zip \
  https://www.kaggle.com/api/v1/datasets/download/pavansanagapati/images-dataset

# unzip the images
unzip -q $inputs_dir/images-dataset.zip -d $inputs_dir

# remove the zip file
rm $inputs_dir/images-dataset.zip 

# generate_random_name() {
#     tr -dc 'a-zA-Z0-9' < /dev/urandom | head -c 12
# }

# declare -A used_names

for dir in $inputs_dir/data/*; do
    if [ -d "$dir" ]; then
        mogrify -format jpg "$dir"/*.png 2>/dev/null
        mogrify -format jpg "$dir"/*.bmp 2>/dev/null
        mogrify -format jpg "$dir"/*.jpeg 2>/dev/null
        rm -f "$dir"/*.png "$dir"/*.bmp "$dir"/*.jpeg

        for file in "$dir"/*.jpg; do
            [ -e "$file" ] || continue

#             while true; do
#                 random_name=$(generate_random_name)
#                 new_filename="${random_name}.jpg"
#                 if [[ ! -e "$(dirname "$file")/$new_filename" && -z "${used_names[$random_name]}" ]]; then
#                     used_names["$random_name"]=1
#                     break
#                 fi
#             done
            mv -- "$file" "$inputs_dir/$(basename "$file")" 2>/dev/null
        done
    fi
done

rm -rf $inputs_dir/data
