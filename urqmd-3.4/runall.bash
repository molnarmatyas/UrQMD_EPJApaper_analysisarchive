#!/bin/sh
# IMP 4.73 roughly for 0-10%, 20 for roughly minBias

set -- 3.0 3.2 3.5 3.9 4.5 7.7 9.2 11.5 14.5 19.6 27

for ecm in "$@"; do
    ecm_label=$(echo "$ecm" | sed 's/\./p/')

    cat > inputfile << EOF
pro 197 79
tar 197 79

nev 100000
IMP 0. 4.73

ecm $ecm
tim 10000 10000

f13
f14
f15
f16
#f19
f20

xxx
EOF

    ./runqmd.bash

    mv test.f19 "auau-${ecm_label}_test_1000evt_b0-4p73_tmax1000.f19"
    #mv test.f19 "auau-${ecm_label}_100000evt_b0-20_tmax10000.f19"
done
