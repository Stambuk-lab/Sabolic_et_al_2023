use strict;
 use File::Basename;

  my @files;
  open (I , "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!";
  while (<I>) {
     /^(.*?)\t/;
     push (@files, "$ARGV[1]/$1");
  }
  closedir I;

  for (my $i = 0; $i< @files;$i++){
    system ("gunzip $files[$i]_2.2.fq.gz && sed 's# 2:N:0/2# 1:N/2#' $files[$i]_2.2.fq | gzip -c > $files[$i]_2.2.fq.gz ");
    system ("gunzip $files[$i]_1.1.fq.gz && sed 's# 1:N:0/1# 1:N/1#' $files[$i]_1.1.fq | gzip -c > $files[$i]_1.1.fq.gz ");
    system ("mv $files[$i]_1.1.fq.gz $files[$i]_1.fq.gz ");
    system ("mv $files[$i]_2.2.fq.gz $files[$i]_2.fq.gz ");
    unlink("$files[$i]_1.1.fq");
    unlink("$files[$i]_2.2.fq");

  }
