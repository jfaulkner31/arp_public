#!/bin/bash

# Use this to push to github since i am lazy

echo
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo NOW PUSHING TO GITHUB
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo
git add .
git commit -m 'latest'
git push origin main
git push origin gh-pages

cd ../docs
explorer.exe .
