awk 'BEGIN { OFS = "\n" } { print ">"$2, $3 }' blast_out/*
