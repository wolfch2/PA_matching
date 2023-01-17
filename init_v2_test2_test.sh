#!/bin/sh

# Install R
sudo apt-get update
sudo apt-get install r-base

# Install Google CLI
# sudo apt-get install apt-transport-https ca-certificates gnupg
# echo "deb https://packages.cloud.google.com/apt cloud-sdk main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
# sudo apt-get update && sudo apt-get install google-cloud-cli


# Create variables
WORK_DIR=/home/onyxia/work/PA_matching
REPO_URL=https://${GIT_PERSONAL_ACCESS_TOKEN}@github.com/fBedecarrats/PA_matching.git # As initial
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

# Git
git clone $REPO_URL $WORK_DIR
chown -R onyxia:users $WORK_DIR
