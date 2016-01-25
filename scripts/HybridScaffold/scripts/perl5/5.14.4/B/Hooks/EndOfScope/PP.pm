package B::Hooks::EndOfScope::PP;
{
  $B::Hooks::EndOfScope::PP::VERSION = '0.13';
}
BEGIN {
  $B::Hooks::EndOfScope::PP::AUTHORITY = 'cpan:FLORA';
}
# ABSTRACT: Execute code after a scope finished compilation - PP implementation

use warnings;
use strict;

use Module::Runtime 'require_module';
use constant _PERL_VERSION => "$]";

BEGIN {
  if (_PERL_VERSION =~ /^5\.009/) {
    # CBA to figure out where %^H got broken and which H::U::HH is sane enough
    die "By design B::Hooks::EndOfScope does not operate in pure-perl mode on perl 5.9.X\n"
  }
  elsif (_PERL_VERSION < '5.010') {
    require_module('B::Hooks::EndOfScope::PP::HintHash');
    *on_scope_end = \&B::Hooks::EndOfScope::PP::HintHash::on_scope_end;
  }
  else {
    require_module('B::Hooks::EndOfScope::PP::FieldHash');
    *on_scope_end = \&B::Hooks::EndOfScope::PP::FieldHash::on_scope_end;
  }
}

use Sub::Exporter::Progressive -setup => {
  exports => ['on_scope_end'],
  groups  => { default => ['on_scope_end'] },
};

sub __invoke_callback {
  local $@;
  eval { $_[0]->(); 1 } or do {
    my $err = $@;
    require Carp;
    Carp::cluck( (join ' ',
      'A scope-end callback raised an exception, which can not be propagated when',
      'B::Hooks::EndOfScope operates in pure-perl mode. Your program will CONTINUE',
      'EXECUTION AS IF NOTHING HAPPENED AFTER THIS WARNING. Below is the complete',
      'exception text, followed by a stack-trace of the callback execution:',
    ) . "\n\n$err\n\r" );

    sleep 1 if -t *STDERR;  # maybe a bad idea...?
  };
}


1;

__END__

=pod

=encoding UTF-8

=for :stopwords Florian Ragwitz Peter Rabbitson Karen Etheridge Tomas Doran

=head1 NAME

B::Hooks::EndOfScope::PP - Execute code after a scope finished compilation - PP implementation

=head1 VERSION

version 0.13

=head1 DESCRIPTION

This is the pure-perl implementation of L<B::Hooks::EndOfScope> based only on
modules available as part of the perl core. Its leaner sibling
L<B::Hooks::EndOfScope::XS> will be automatically preferred if all
dependencies are available and C<$ENV{B_HOOKS_ENDOFSCOPE_IMPLEMENTATION}> is
not set to C<'PP'>.

=head1 FUNCTIONS

=head2 on_scope_end

    on_scope_end { ... };

    on_scope_end $code;

Registers C<$code> to be executed after the surrounding scope has been
compiled.

This is exported by default. See L<Sub::Exporter> on how to customize it.

=head1 AUTHORS

=over 4

=item *

Florian Ragwitz <rafl@debian.org>

=item *

Peter Rabbitson <ribasushi@cpan.org>

=back

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2008 by Florian Ragwitz.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
