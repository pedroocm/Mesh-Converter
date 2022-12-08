#! /bin/bash

if [[ $# -eq 0 ]]; then
    echo "> Type \"./converter build\" to build the docker image"
    echo "> Type \"./converter --help \" to get information about the cli"
    echo "> Type \"./converter file_path type_of_mesh\" to run"
elif [[ $1 == "build" ]]; then
    docker build --tag converter:latest .
elif [[ $1 == "bdev" ]]; then
    docker build --tag converter:dev -f dockerfiles/Dockerfile.debug .
else
    docker run -u $(id -u) -v $(pwd):/files -ti converter:latest "$@"
fi

