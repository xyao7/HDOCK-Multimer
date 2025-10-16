

## Download jq for analyzing JSON files
wget https://github.com/jqlang/jq/releases/download/jq-1.8.1/jq-linux64
chmod +x jq-linux64

## Download and compile stride for analyzing protein secondary structures
wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz --no-check-certificate
mkdir -p stride
tar -xzf stride.tar.gz -C stride/
cd stride/
make 
