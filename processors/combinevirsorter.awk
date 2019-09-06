#/usr/bin/awk -f
# -v d=distance_to_prophage has to be passed.
{ if (NF == 27) {
    if ( $4 == "NA" || ( $4 > $24 - d && $4 < $25 + d ) || ( $5 > $24 - d && $5 < $25 + d ) )
        {for (i=2;i<=20;i++) printf "%s%s",$i,FS; print $23,2;}
    else {for (i=2;i<=20;i++) printf "%s%s",$i,FS; print $23,1;}
}
else {for (i=2;i<=20;i++) printf "%s%s",$i,FS; print "NULL",0;}
}
