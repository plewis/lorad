#!/bin/bash

mkdir gproto

INDEXES=({1..30})
for i in ${INDEXES[@]} ; do
    python3 deploy.py $i $i`pickseed 3`
    mv g$i gproto
done

# Create go.sh file to make it easy to submit all runs from all subdirectories
echo -e "#!/bin/bash\n" > gproto/go.sh
echo "INDEXES=({1..30})" >> gproto/go.sh
echo "for i in \${INDEXES[@]} ; do" >> gproto/go.sh
echo "  cd g\$i" >> gproto/go.sh
echo "  . submit-all.sh" >> gproto/go.sh
echo "  cd .." >> gproto/go.sh
echo "done" >> gproto/go.sh

# Create goloradml.sh file to make it easy to run loradml for all models in all subdirs
echo -e "#!/bin/bash\n" > gproto/goloradml.sh
echo "INDEXES=({1..30})" >> gproto/goloradml.sh
echo "for i in \${INDEXES[@]} ; do" >> gproto/goloradml.sh
echo "  cd g\$i" >> gproto/goloradml.sh
echo "  . loradml-all.sh" >> gproto/goloradml.sh
echo "  cd .." >> gproto/goloradml.sh
echo "done" >> gproto/goloradml.sh

# Copy filter.py file
cp filter.py gproto

# Copy summarize.py file
cp summarize.py gproto
