echo '-------------'
echo 'uninstall'
echo '-------------'
sudo python setup.py develop -u
#sudo python setup.py install -u 
echo '-------------'
echo 'remove *.pyc'
echo '-------------'
find ./ -iname "*.pyc" -exec rm {} \;
echo '-------------'
echo 'install'
echo '-------------'
sudo python setup.py develop
#sudo python setup.py install 
