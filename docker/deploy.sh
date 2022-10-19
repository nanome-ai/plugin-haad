#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

existing=$(docker ps -aqf name=haad-plugin)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name haad-plugin \
--restart unless-stopped \
-e ARGS="$*" \
haad-plugin
