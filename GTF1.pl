use strict;
use warnings;

$" = "\t";

while (<>) {
    my @F = split;
    $F[8] = $F[9];
    $F[9] = "$F[9]:$F[9]";
    print "@F\n";
}
