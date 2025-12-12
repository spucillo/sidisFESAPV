First dump the list of eicrecon files:
```
xrdfs root://dtn-eic.jlab.org ls /volatile/eic/EPIC/RECO/25.10.0/epic_craterlake/DIS/NC/10x100/minQ2=1/ > remote_list.txt
```

Then use the script

```
./run_chunks.sh remote_list.txt 20
```

This will run (in series) creating an output file every 20 lines of `remote_list.txt`
