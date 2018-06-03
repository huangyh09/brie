#!/bin/sh

# Note, sashimi_plot is originally from MISO, and developed 
# under Python2, thus we include it in BRIE-kit as a folder.
# It can be imported from its path, but not directly from 
# briekit package.

ANNO_DIR=~/annotation
GFF_FILE=$ANNO_DIR/mouse/AS_events/SE.filtered.gff3
GFF_DIR=$ANNO_DIR/mouse/AS_events/misoSE/

## generate GFF: require misopy package
gzip -d $GFF_FILE.gz
index_gff --index $GFF_FILE $GFF_DIR
gzip $GFF_FILE


## generate sashimi plot
PLOT_DIR=$HOME/test/plot
SASHIMI=$HOME/MyGit/briekit/sashimi_plot/sashimi_plot.py

python $SASHIMI --plot-event "ENSMUSG00000027478.AS2" $GFF_DIR sashimi_setting.txt --output-dir $PLOT_DIR --plot-label "DNMT3B-exon2.pdf" --plot-title "DNMT3B-exon2"
