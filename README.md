# Enigma cipher breaker

```
Enigma version 1.0
Copyright (C) 2017 Torbjørn Rognes

Usage: enigma [OPTIONS]
  -h           Show help information
  -v           Show version information
  -u X         Reflector (umkehrwalze) X (A-C or .) [.]
  -w XYZ       Wheels (walzen) XYZ (1-8 or .) [...]
  -x integer   Highest wheel number to use (1-8) [5]
  -r XYZ       Ring positions (ringstellung) XYZ (A-Z or .) [AA.]
  -g XYZ       Start positions (grundstellung) XYZ (A-Z or .) [...]
  -s ABYZ      Plugboard (steckerbrett) letter pairs (A-Z pairs) [none]
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
algorithm will be used to try to determine the plugboard settings.```
