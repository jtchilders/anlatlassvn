#\bin\sh
PWD=`pwd -P`
grep '"@@@_SET_PATH_HERE_@@@"' ffinit.f  > /dev/null
if [ $? -eq 0 ] 
then
echo "setting ff path to $PWD"
sed -i -e "s:\"@@@\_SET\_PATH\_HERE\_@@@\":\"$PWD\/\":" ffinit.f
fi