
export ARANGODB_HOME="/Users/raymondleclair/Projects/NLM/NCBI-Cell/springbok-ncbi-cell/ncbi-cell/data/arangodb"

docker run \
       -e ARANGO_ROOT_PASSWORD=$ARANGO_ROOT_PASSWORD \
       -p 8529:8529 \
       -d \
       -v $ARANGODB_HOME:/var/lib/arangodb3 \
       arangodb
