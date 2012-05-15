#/bin/bash
if [ "$INPUTFILE" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

if [ "$RUN_NO" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

if [ "$RUN_DIR" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

cd $RUN_DIR

./araSim $INPUTFILE $RUN_NO

