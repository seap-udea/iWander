#########################################################
#    _ _       __                __         		#
#   (_) |     / /___ _____  ____/ /__  _____		#
#  / /| | /| / / __ `/ __ \/ __  / _ \/ ___/		#
# / / | |/ |/ / /_/ / / / / /_/ /  __/ /    		#
#/_/  |__/|__/\__,_/_/ /_/\__,_/\___/_/     		#
# Dynamics of Interestellar Wanderers			#
# Jorge I. Zuluaga et al. [)] 2017			#
# http://github.com/seap-udea/iWander.git		#
#########################################################
# Install the iwander package
#########################################################
#!/bin/bash
echo "Installing templates..."
cp util/conf/* .

echo "Installing packages..."
pip3 install pycparser

echo "Unpacking large files..."
make unpack

echo "Installing AstroRV repository..."
cd db
git clone http://github.com/seap-udea/AstroRV.git
make -C AstroRV unpack
cd -

echo "Testing compilation..."
make

