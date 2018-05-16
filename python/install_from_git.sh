installfolder=~/RTLib
git clone https://blahoink:mala3dofu@github.com/blahoink/RTLib $installfolder
pushd $installfolder/python
pip install ./
popd
rm -rf $installfolder