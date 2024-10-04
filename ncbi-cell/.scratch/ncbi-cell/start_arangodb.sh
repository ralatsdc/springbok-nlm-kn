ARANGODB_HOME="/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/arangodb"
docker run \
       -e ARANGO_NOAUTH=1 \
       -p 8529:8529 \
       -d \
       -v $ARANGODB_HOME:/var/lib/arangodb3 \
       arangodb
