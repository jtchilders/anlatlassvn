# needed since here I used madgraph v 4.5.0, whereas
# when I generated the code for the top I used v 4.4.44
#!/bin/bash
for i in *.f
do
    cat $i | sed -e "s/maxflow/maxamps/g" > temp 
    mv temp $i
done