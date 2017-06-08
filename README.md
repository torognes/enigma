# Enigma cipher tool

This is a tool to encrypt or decrypt messages using an Enigma cipher machine simulator.

If you do not know the correct settings for decryption, it can also be used to try out
a large number of settings and look for plaintext messages that look similar to
text written in a selected language.

The settings include the reflector (umkehrwalze) and wheels (walzen) used,
the ring positions (ringstellung) and start positions (grundstellung),
as well as the position of the plugs in the plugboard (steckerbrett).

Both the common three-wheel Enigma as well as the special Norway Enigma (Norenigma) is supported.

If specified (with the -c option), a hill-climbing algorithm will be used to identify the optimal plugboard configuration.

All possible combinations of the other unspecified settings will be tried.

```
Enigma cipher tool version 1.0
Copyright (C) 2017 Torbjørn Rognes

Usage: enigma [OPTIONS]
  -h           Show help information
  -v           Show version information
  -u X         Reflector (umkehrwalze) X (A-C, N or .) [.]
  -w XYZ       Wheels (walzen) XYZ (1-8 or .) [...]
  -x integer   Highest wheel number to use (3-8) [5]
  -n           Use the Norway Enigma reflector (N) and wheels (1-5)
  -r XYZ       Ring positions (ringstellung) XYZ (A-Z or .) [AA.]
  -g XYZ       Start positions (grundstellung) XYZ (A-Z or .) [...]
  -s AB...     Plugboard (steckerbrett) letter pairs (A-Z pairs) [none]
  -c           Perform hill climbing to determine plugboard settings
  -l language  Plaintext language (german, english, danish, french) [german]
  -i           Use index of coincidence (IC) to determine plaintext score
  -m           Use monogram statistics to determine plaintext score
  -b           Use bigram statistics to determine plaintext score
  -t           Use trigram statistics to determine plaintext score
  -q           Use quadgram statistics to determine plaintext score [default]
  -p filename  Name of file containing plaintext to compare result with

Defaults are indicated in [square brackets].

The ciphertext is read from standard input. The final plaintext is written
to standard output.

For the reflector, wheels, ring position and start position, a dot (.)
works as a wild card, leaving it unspecified. When these settings are not
specified, the program will try all combinations to find the settings
resulting in the highest plaintext score. If asked for, a hill climbing
algorithm will be used to try to determine the plugboard settings.
```

The files with the ngram frequencies for various languages have been obtained from the
[Practical cryptograhy](http://practicalcryptography.com/cryptanalysis/letter-frequencies-various-languages/)
website. Additional languages are available there.

The hill climbing strategy is based on the algorithms described in the
[publications by Frode Weierud et al.](http://cryptocellar.org/Enigma/)

The software is available under the GNU GPL version 3 license.
