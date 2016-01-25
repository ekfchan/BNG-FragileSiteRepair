###########################################################################
#
# This file is auto-generated by the Perl DateTime Suite locale
# generator (0.05).  This code generator comes with the
# DateTime::Locale distribution in the tools/ directory, and is called
# generate-from-cldr.
#
# This file as generated from the CLDR XML locale data.  See the
# LICENSE.cldr file included in this distribution for license details.
#
# This file was generated from the source file lv_LV.xml
# The source file version number was 1.48, generated on
# 2009/05/05 23:06:38.
#
# Do not edit this file directly.
#
###########################################################################

package DateTime::Locale::lv_LV;

use strict;
use warnings;

our $VERSION = '0.46';

use utf8;

use base 'DateTime::Locale::lv';

sub cldr_version { return "1\.7\.1" }

{
    my $first_day_of_week = "1";
    sub first_day_of_week { return $first_day_of_week }
}

{
    my $glibc_date_format = "\%Y\.\%m\.\%d\.";
    sub glibc_date_format { return $glibc_date_format }
}

{
    my $glibc_date_1_format = "\%a\ \%b\ \%e\ \%H\:\%M\:\%S\ \%Z\ \%Y";
    sub glibc_date_1_format { return $glibc_date_1_format }
}

{
    my $glibc_datetime_format = "\%A\,\ \%Y\.\ gada\ \%e\.\ \%B\,\ plkst\.\ \%H\ un\ \%M";
    sub glibc_datetime_format { return $glibc_datetime_format }
}

{
    my $glibc_time_format = "\%T";
    sub glibc_time_format { return $glibc_time_format }
}

1;

__END__


=pod

=encoding utf8

=head1 NAME

DateTime::Locale::lv_LV

=head1 SYNOPSIS

  use DateTime;

  my $dt = DateTime->now( locale => 'lv_LV' );
  print $dt->month_name();

=head1 DESCRIPTION

This is the DateTime locale package for Latvian Latvia.

=head1 DATA

This locale inherits from the L<DateTime::Locale::lv> locale.

It contains the following data.

=head2 Days

=head3 Wide (format)

  pirmdiena
  otrdiena
  trešdiena
  ceturtdiena
  piektdiena
  sestdiena
  svētdiena

=head3 Abbreviated (format)

  Pr
  Ot
  Tr
  Ce
  Pk
  Se
  Sv

=head3 Narrow (format)

  P
  O
  T
  C
  P
  S
  S

=head3 Wide (stand-alone)

  pirmdiena
  otrdiena
  trešdiena
  ceturtdiena
  piektdiena
  sestdiena
  svētdiena

=head3 Abbreviated (stand-alone)

  Pr
  Ot
  Tr
  Ce
  Pk
  Se
  Sv

=head3 Narrow (stand-alone)

  P
  O
  T
  C
  P
  S
  S

=head2 Months

=head3 Wide (format)

  janvāris
  februāris
  marts
  aprīlis
  maijs
  jūnijs
  jūlijs
  augusts
  septembris
  oktobris
  novembris
  decembris

=head3 Abbreviated (format)

  janv.
  febr.
  marts
  apr.
  maijs
  jūn.
  jūl.
  aug.
  sept.
  okt.
  nov.
  dec.

=head3 Narrow (format)

  J
  F
  M
  A
  M
  J
  J
  A
  S
  O
  N
  D

=head3 Wide (stand-alone)

  janvāris
  februāris
  marts
  aprīlis
  maijs
  jūnijs
  jūlijs
  augusts
  septembris
  oktobris
  novembris
  decembris

=head3 Abbreviated (stand-alone)

  janv.
  febr.
  marts
  apr.
  maijs
  jūn.
  jūl.
  aug.
  sept.
  okt.
  nov.
  dec.

=head3 Narrow (stand-alone)

  J
  F
  M
  A
  M
  J
  J
  A
  S
  O
  N
  D

=head2 Quarters

=head3 Wide (format)

  1. ceturksnis
  2. ceturksnis
  3. ceturksnis
  4. ceturksnis

=head3 Abbreviated (format)

  C1
  C2
  C3
  C4

=head3 Narrow (format)

  1
  2
  3
  4

=head3 Wide (stand-alone)

  1. ceturksnis
  2. ceturksnis
  3. ceturksnis
  4. ceturksnis

=head3 Abbreviated (stand-alone)

  C1
  C2
  C3
  C4

=head3 Narrow (stand-alone)

  1
  2
  3
  4

=head2 Eras

=head3 Wide

  pirms mūsu ēras
  mūsu ērā

=head3 Abbreviated

  p.m.ē.
  m.ē.

=head3 Narrow

  p.m.ē.
  m.ē.

=head2 Date Formats

=head3 Full

   2008-02-05T18:30:30 = otrdiena, 2008. gada 5. februāris
   1995-12-22T09:05:02 = piektdiena, 1995. gada 22. decembris
  -0010-09-15T04:44:23 = sestdiena, -10. gada 15. septembris

=head3 Long

   2008-02-05T18:30:30 = 2008. gada 5. februāris
   1995-12-22T09:05:02 = 1995. gada 22. decembris
  -0010-09-15T04:44:23 = -10. gada 15. septembris

=head3 Medium

   2008-02-05T18:30:30 = 2008. gada 5. febr.
   1995-12-22T09:05:02 = 1995. gada 22. dec.
  -0010-09-15T04:44:23 = -10. gada 15. sept.

=head3 Short

   2008-02-05T18:30:30 = 05.02.08
   1995-12-22T09:05:02 = 22.12.95
  -0010-09-15T04:44:23 = 15.09.-10

=head3 Default

   2008-02-05T18:30:30 = 2008. gada 5. febr.
   1995-12-22T09:05:02 = 1995. gada 22. dec.
  -0010-09-15T04:44:23 = -10. gada 15. sept.

=head2 Time Formats

=head3 Full

   2008-02-05T18:30:30 = 18:30:30 UTC
   1995-12-22T09:05:02 = 09:05:02 UTC
  -0010-09-15T04:44:23 = 04:44:23 UTC

=head3 Long

   2008-02-05T18:30:30 = 18:30:30 UTC
   1995-12-22T09:05:02 = 09:05:02 UTC
  -0010-09-15T04:44:23 = 04:44:23 UTC

=head3 Medium

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 09:05:02
  -0010-09-15T04:44:23 = 04:44:23

=head3 Short

   2008-02-05T18:30:30 = 18:30
   1995-12-22T09:05:02 = 09:05
  -0010-09-15T04:44:23 = 04:44

=head3 Default

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 09:05:02
  -0010-09-15T04:44:23 = 04:44:23

=head2 Datetime Formats

=head3 Full

   2008-02-05T18:30:30 = otrdiena, 2008. gada 5. februāris 18:30:30 UTC
   1995-12-22T09:05:02 = piektdiena, 1995. gada 22. decembris 09:05:02 UTC
  -0010-09-15T04:44:23 = sestdiena, -10. gada 15. septembris 04:44:23 UTC

=head3 Long

   2008-02-05T18:30:30 = 2008. gada 5. februāris 18:30:30 UTC
   1995-12-22T09:05:02 = 1995. gada 22. decembris 09:05:02 UTC
  -0010-09-15T04:44:23 = -10. gada 15. septembris 04:44:23 UTC

=head3 Medium

   2008-02-05T18:30:30 = 2008. gada 5. febr. 18:30:30
   1995-12-22T09:05:02 = 1995. gada 22. dec. 09:05:02
  -0010-09-15T04:44:23 = -10. gada 15. sept. 04:44:23

=head3 Short

   2008-02-05T18:30:30 = 05.02.08 18:30
   1995-12-22T09:05:02 = 22.12.95 09:05
  -0010-09-15T04:44:23 = 15.09.-10 04:44

=head3 Default

   2008-02-05T18:30:30 = 2008. gada 5. febr. 18:30:30
   1995-12-22T09:05:02 = 1995. gada 22. dec. 09:05:02
  -0010-09-15T04:44:23 = -10. gada 15. sept. 04:44:23

=head2 Available Formats

=head3 d (d)

   2008-02-05T18:30:30 = 5
   1995-12-22T09:05:02 = 22
  -0010-09-15T04:44:23 = 15

=head3 Ed (EEE, d.)

   2008-02-05T18:30:30 = Ot, 5.
   1995-12-22T09:05:02 = Pk, 22.
  -0010-09-15T04:44:23 = Se, 15.

=head3 EEEd (EEE, d.)

   2008-02-05T18:30:30 = Ot, 5.
   1995-12-22T09:05:02 = Pk, 22.
  -0010-09-15T04:44:23 = Se, 15.

=head3 H (H)

   2008-02-05T18:30:30 = 18
   1995-12-22T09:05:02 = 9
  -0010-09-15T04:44:23 = 4

=head3 HHmm (HH:mm)

   2008-02-05T18:30:30 = 18:30
   1995-12-22T09:05:02 = 09:05
  -0010-09-15T04:44:23 = 04:44

=head3 HHmmss (HH:mm:ss)

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 09:05:02
  -0010-09-15T04:44:23 = 04:44:23

=head3 Hm (HH:mm)

   2008-02-05T18:30:30 = 18:30
   1995-12-22T09:05:02 = 09:05
  -0010-09-15T04:44:23 = 04:44

=head3 hm (h:mm a)

   2008-02-05T18:30:30 = 6:30 PM
   1995-12-22T09:05:02 = 9:05 AM
  -0010-09-15T04:44:23 = 4:44 AM

=head3 Hms (H:mm:ss)

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 9:05:02
  -0010-09-15T04:44:23 = 4:44:23

=head3 hms (h:mm:ss a)

   2008-02-05T18:30:30 = 6:30:30 PM
   1995-12-22T09:05:02 = 9:05:02 AM
  -0010-09-15T04:44:23 = 4:44:23 AM

=head3 M (L)

   2008-02-05T18:30:30 = 2
   1995-12-22T09:05:02 = 12
  -0010-09-15T04:44:23 = 9

=head3 Md (dd.mm.)

   2008-02-05T18:30:30 = 05.30.
   1995-12-22T09:05:02 = 22.05.
  -0010-09-15T04:44:23 = 15.44.

=head3 MEd (E, dd.MM.)

   2008-02-05T18:30:30 = Ot, 05.02.
   1995-12-22T09:05:02 = Pk, 22.12.
  -0010-09-15T04:44:23 = Se, 15.09.

=head3 MMM (LLL)

   2008-02-05T18:30:30 = febr.
   1995-12-22T09:05:02 = dec.
  -0010-09-15T04:44:23 = sept.

=head3 MMMd (d. MMM)

   2008-02-05T18:30:30 = 5. febr.
   1995-12-22T09:05:02 = 22. dec.
  -0010-09-15T04:44:23 = 15. sept.

=head3 MMMEd (E, d. MMM)

   2008-02-05T18:30:30 = Ot, 5. febr.
   1995-12-22T09:05:02 = Pk, 22. dec.
  -0010-09-15T04:44:23 = Se, 15. sept.

=head3 MMMMd (d. MMMM)

   2008-02-05T18:30:30 = 5. februāris
   1995-12-22T09:05:02 = 22. decembris
  -0010-09-15T04:44:23 = 15. septembris

=head3 MMMMEd (E, d. MMMM)

   2008-02-05T18:30:30 = Ot, 5. februāris
   1995-12-22T09:05:02 = Pk, 22. decembris
  -0010-09-15T04:44:23 = Se, 15. septembris

=head3 mmss (mm:ss)

   2008-02-05T18:30:30 = 30:30
   1995-12-22T09:05:02 = 05:02
  -0010-09-15T04:44:23 = 44:23

=head3 ms (mm:ss)

   2008-02-05T18:30:30 = 30:30
   1995-12-22T09:05:02 = 05:02
  -0010-09-15T04:44:23 = 44:23

=head3 y (y. 'g'.)

   2008-02-05T18:30:30 = 2008. g.
   1995-12-22T09:05:02 = 1995. g.
  -0010-09-15T04:44:23 = -10. g.

=head3 yM (mm.yyyy.)

   2008-02-05T18:30:30 = 30.2008.
   1995-12-22T09:05:02 = 05.1995.
  -0010-09-15T04:44:23 = 44.-010.

=head3 yMEd (EEE, dd.mm.yyyy.)

   2008-02-05T18:30:30 = Ot, 05.30.2008.
   1995-12-22T09:05:02 = Pk, 22.05.1995.
  -0010-09-15T04:44:23 = Se, 15.44.-010.

=head3 yMMM (yyyy. 'g'. MMM)

   2008-02-05T18:30:30 = 2008. g. febr.
   1995-12-22T09:05:02 = 1995. g. dec.
  -0010-09-15T04:44:23 = -010. g. sept.

=head3 yMMMEd (EEE, yyyy. 'g'. dd. MMM)

   2008-02-05T18:30:30 = Ot, 2008. g. 05. febr.
   1995-12-22T09:05:02 = Pk, 1995. g. 22. dec.
  -0010-09-15T04:44:23 = Se, -010. g. 15. sept.

=head3 yMMMM (y. 'g'. MMMM)

   2008-02-05T18:30:30 = 2008. g. februāris
   1995-12-22T09:05:02 = 1995. g. decembris
  -0010-09-15T04:44:23 = -10. g. septembris

=head3 yQ (Q yyyy)

   2008-02-05T18:30:30 = 1 2008
   1995-12-22T09:05:02 = 4 1995
  -0010-09-15T04:44:23 = 3 -010

=head3 yQQQ (y QQQ)

   2008-02-05T18:30:30 = 2008 C1
   1995-12-22T09:05:02 = 1995 C4
  -0010-09-15T04:44:23 = -10 C3

=head3 yyQ (Q yy)

   2008-02-05T18:30:30 = 1 08
   1995-12-22T09:05:02 = 4 95
  -0010-09-15T04:44:23 = 3 -10

=head3 yyyy (y. 'g'.)

   2008-02-05T18:30:30 = 2008. g.
   1995-12-22T09:05:02 = 1995. g.
  -0010-09-15T04:44:23 = -10. g.

=head2 Miscellaneous

=head3 Prefers 24 hour time?

Yes

=head3 Local first day of the week

pirmdiena


=head1 SUPPORT

See L<DateTime::Locale>.

=cut