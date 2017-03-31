VERSION=1.2.0
BASENAME=dashmm-$VERSION
FILENAME=$BASENAME.tar


make clean
cd demo/basic
make clean
cd ../stepping
make clean
cd ../user
make clean
cd ../..


mkdir $BASENAME
cp AUTHORS $BASENAME/
cp INSTALL $BASENAME/
cp Makefile $BASENAME/
cp README $BASENAME/
cp LICENSE $BASENAME/
mkdir $BASENAME/doc
cp doc/DASHMMUserGuide.pdf $BASENAME/doc/
cp -r demo/ $BASENAME/
cp -r include/ $BASENAME/
cp -r src/ $BASENAME/
mkdir $BASENAME/lib


tar cvf $FILENAME $BASENAME/
gzip $FILENAME

rm -r $BASENAME/
