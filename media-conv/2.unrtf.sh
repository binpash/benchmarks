#!/bin/bash

#tag: rtf-to-txt
find $IN -name '*.rtf' | xargs -I {} unrtf {} --text
