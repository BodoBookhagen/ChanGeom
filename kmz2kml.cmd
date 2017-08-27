#!/bin/bash
# Autor: Bodo Bookhagen

KMZ=$1

if [ -d "${KMZ}" ] ; then
    echo "$KMZ is a directory";
    for filename in $KMZ/*.kmz; do
        base_name=${filename%\.*}
        if ( unzip -qq $filename &> /dev/null ); then
            mv doc.kml $base_name.kml
        fi
    done
else
    if [ -f "${KMZ}" ]; then
        echo "${KMZ} is a file";
        base_name=${KMZ%\.*}
        if ( unzip -qq $KMZ &> /dev/null ); then
            mv doc.kml $base_name.kml
        fi
    else
        echo "${KMZ} is not valid";
        exit 1
    fi
fi
