VERSION=1.1.0
BASENAME=dashmm-$VERSION
FILENAME=$BASENAME.tar


echo "Did you remove the comments in evaluator.h?"
echo "If not, please go remove them now."


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
cp doc/DASHMMUserGuide.pdf $BASENAME/
cp -r demo/ $BASENAME/
cp -r include/ $BASENAME/
cp -r src/ $BASENAME/
mkdir $BASENAME/lib


tar cvf $FILENAME $BASENAME/
gzip $FILENAME

rm -r $BASENAME/
