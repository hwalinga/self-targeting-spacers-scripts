#!/usr/bin/awk -f 
{qident = $3 * ($4 / $13);}
qident>75 && $14>1000{
    sub(".*/", "", FILENAME);
    sub(".selfhits", "", FILENAME);
    print FILENAME,$1,$2,$9,$10,qident,$13,$14;
}
