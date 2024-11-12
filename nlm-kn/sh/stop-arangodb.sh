docker container stop \
       $(docker ps | grep arangodb | cut -d " " -f 1)
