#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <tmmintrin.h>

const int maxlen = 10240;
char ciphertext[maxlen+1];
char plaintext[maxlen+1];
char altplaintext[maxlen+1];
int textlength;
  
char * opt_ukw;
char * opt_walzen;
char * opt_ringstellung;
char * opt_grundstellung;
char * opt_steckerbrett;
char * opt_logfilename;
char * opt_plaintext; /* plaintext to compare to */
char * opt_language; /* german (default), english, danish, french, ... */
int opt_norenigma; /* use the 5 Norenigma (Norway Enigma) wheels */
int opt_maxwheel;
int opt_scoring;
int opt_threads;
int opt_hillclimb;

/* uwwwrrrggg = 3*8*7*6*26*26*26*26*26*26 = 311 387 102 208 */

const char * reflector_string[] =
  {
    "EJMZALYXVBWFCRQUONTSPIKHGD",    // A
    "YRUHQSLDPXNGOKMIEBFZCWVJAT",    // B
    "FVPJIAOYEDRZXWGCTKUQSBNMHL",    // C
    "MOWJYPUXNDSRAIBFVLKZGQCHET",    // Norway
    "ENKQAUYWJICOPBLMDXZVFTHRGS",    // UKW-b M4 thin
    "RDOBJNTKVEHMLFCWZAXGYIPSUQ"     // UKW-c M4 thin
  };

const char * rotor_string[] = 
  {
    "EKMFLGDQVZNTOWYHXUSPAIBRCJ",  // i
    "AJDKSIRUXBLHWTMCQGZNPYFVOE",  // ii
    "BDFHJLCPRTXVZNYEIWGAKMUSQO",  // iii
    "ESOVPZJAYQUIRHXLNFTGKDCMWB",  // iv
    "VZBRGITYUPSDNHLXAWMJQOFECK",  // v
    "JPGVOUMFYQBENHZRDKASXLICTW",  // vi
    "NZJHGRCXMYSWBOUFAIVLPEKQDT",  // vii
    "FKQHTLXOCBJSPDZRAMEWNIUYGV",  // viii
    "WTOKASUYVRBXJHQCPZEFMDINLG",  // Norway i
    "GJLPUBSWEMCTQVHXAOFZDRKYNI",  // Norway ii
    "JWFMHNBPUSDYTIXVZGRQLAOEKC",  // Norway iii
    "ESOVPZJAYQUIRHXLNFTGKDCMWB",  // Norway iv (equal to iv)
    "HEJXQOTZBVFDASCILWPGYNMURK",  // Norway v
    "LEYJVCNIXWPBQMDRTAKZGFUHOS",  // Beta
    "FSOKANUERHMBTIYCWLQPZXVGJD"   // Gamma
  };

const char * notch_string[] =
  {
    "Q",
    "E",
    "V",
    "J",
    "Z",
    "MZ",
    "MZ",
    "MZ",
    "Q",
    "E",
    "V",
    "J",
    "Z",
    "",
    ""
  };

const int alphabet_size = 26;
const int wheels = 3;
const int reflector_count = sizeof(reflector_string) / sizeof(char *);
const int rotor_count = sizeof(rotor_string) / sizeof(char *);

unsigned char rotor_fwd[rotor_count][alphabet_size];
unsigned char rotor_rev[rotor_count][alphabet_size];
unsigned char notch[rotor_count][alphabet_size];
unsigned char reflector[reflector_count][alphabet_size];
unsigned char steckerbrett[alphabet_size];

int ukw;
int walzenlage[wheels];
unsigned char grundstellung[wheels];
unsigned char ringstellung[wheels];

unsigned char num_ciphertext[maxlen];
unsigned char num_plaintext[maxlen];

unsigned char subst_array[alphabet_size][alphabet_size][alphabet_size][alphabet_size];
unsigned char mapping[maxlen][26];

double monograms[26];
double bigrams[26][26];
double trigrams[26][26][26];
double quadgrams[26][26][26][26];


void fatal(const char * message)
{
  fprintf(stderr, "\nFatal error: %s\n", message);
  exit(1);
}

inline int char2num(char x)
{
  return x - 65;
}

inline char num2char(int x)
{
  return 65 + x;
}

void monograms_read()
{
  for(int i=0; i<26; i++)
    monograms[i] = 1.0;

  double total = 26;

  FILE * f;
  
  char filename[100];
  
  strcpy(filename, opt_language);
  strcat(filename, "_monograms.txt");
  
  f = fopen(filename, "r");
  
  if (!f)
    {
      fprintf(stderr, "Fatal error: Unable to open the language statistics file %s\n",
              filename);
      exit(1);
    }

  while(1)
    {
      char a;
      int count;
      int ret = fscanf(f, "%c %d\n", & a, & count);
      if (ret > 0)
        {
          if ((a >= 'A') && (a <= 'Z'))
            {
              monograms[char2num(a)] = count + 1.0;
              total += count;
            }
        }
      else
        break;
    }

  for(int i=0; i<26; i++)
      monograms[i] = log10(monograms[i] / total);

  fclose(f);
}


void bigrams_read()
{
  for(int i=0; i<26; i++)
    for(int j=0; j<26; j++)
      bigrams[i][j] = 1.0;

  double total = 26*26;

  FILE * f;
  
  char filename[100];
  
  strcpy(filename, opt_language);
  strcat(filename, "_bigrams.txt");
  
  f = fopen(filename, "r");
  
  if (!f)
    {
      fprintf(stderr, "Fatal error: Unable to open the language statistics file %s\n",
              filename);
      exit(1);
    }

  while(1)
    {
      char a;
      char b;
      int count;
      int ret = fscanf(f, "%c%c %d\n", & a, & b, & count);
      if (ret > 0)
        {
          if ((a >= 'A') && (a <= 'Z') && (b >= 'A') && (b <= 'Z'))
            {
              bigrams[char2num(a)][char2num(b)] = count + 1;
              total += count;
            }
        }
      else
        break;
    }

  for(int i=0; i<26; i++)
    for(int j=0; j<26; j++)
      bigrams[i][j] = log10(bigrams[i][j] / total);

  fclose(f);
}


void trigrams_read()
{
  for(int i=0; i<26; i++)
    for(int j=0; j<26; j++)
      for(int k=0; k<26; k++)
        trigrams[i][j][k] = 1.0;

  double total = 26*26*26;
  
  FILE * f;
  
  char filename[100];
  
  strcpy(filename, opt_language);
  strcat(filename, "_trigrams.txt");
  
  f = fopen(filename, "r");
  
  if (!f)
    {
      fprintf(stderr, "Fatal error: Unable to open the language statistics file %s\n",
              filename);
      exit(1);
    }

  while(1)
    {
      char a;
      char b;
      char c;
      int count;
      int ret = fscanf(f, "%c%c%c %d\n", & a, & b, & c, & count);
      
      if (ret < 1)
        break;
      
      if (ret > 0)
        {
          if ((a >= 'A') && (a <= 'Z') && 
              (b >= 'A') && (b <= 'Z') &&
              (c >= 'A') && (c <= 'Z'))
            {
              trigrams[char2num(a)][char2num(b)][char2num(c)] = count + 1;
              total += count;
            }
        }
    }

  fclose(f);

  for(int i=0; i<26; i++)
    for(int j=0; j<26; j++)
      for(int k=0; k<26; k++)
        trigrams[i][j][k] = log10(trigrams[i][j][k] / total);
}

void quadgrams_read()
{
  for(int i=0; i<26; i++)
    for(int j=0; j<26; j++)
      for(int k=0; k<26; k++)
        for(int l=0; l<26; l++)
          quadgrams[i][j][k][l] = 1.0;

  double total = 26*26*26*26;
  
  FILE * f;
  
  char filename[100];
  
  strcpy(filename, opt_language);
  strcat(filename, "_quadgrams.txt");
  
  f = fopen(filename, "r");
  
  if (!f)
    {
      fprintf(stderr, "Fatal error: Unable to open the language statistics file %s\n",
              filename);
      exit(1);
    }

  while(1)
    {
      char a;
      char b;
      char c;
      char d;
      int count;
      int ret = fscanf(f, "%c%c%c%c %d\n", & a, & b, & c, & d, & count);
      
      if (ret < 1)
        break;
      
      if (ret > 0)
        {
          if ((a >= 'A') && (a <= 'Z') && 
              (b >= 'A') && (b <= 'Z') &&
              (c >= 'A') && (c <= 'Z') &&
              (d >= 'A') && (d <= 'Z'))
            {
              quadgrams[char2num(a)][char2num(b)][char2num(c)][char2num(d)] = count + 1;
              total += count;
            }
        }
    }

  fclose(f);

  for(int i=0; i<26; i++)
    for(int j=0; j<26; j++)
      for(int k=0; k<26; k++)
        for(int l=0; l<26; l++)
          quadgrams[i][j][k][l] = log10(quadgrams[i][j][k][l] / total);
}



double quadgram_score(char * text, int len)
{
  double score = 0.0;
  for (int i=0; i<len-3; i++)
    score += quadgrams
      [char2num(text[i])]
      [char2num(text[i+1])]
      [char2num(text[i+2])]
      [char2num(text[i+3])];
  return score;
}

double trigram_score(char * text, int len)
{
  double score = 0.0;
  for (int i=0; i<len-2; i++)
    score += trigrams[char2num(text[i])][char2num(text[i+1])][char2num(text[i+2])];
  return score;
}

double bigram_score(char * text, int len)
{
  double score = 0.0;
  for (int i=0; i<len-1; i++)
    score += bigrams[char2num(text[i])][char2num(text[i+1])];
  return score;
}

double monogram_score(char * text, int len)
{
  double score = 0.0;
  for (int i=0; i<len; i++)
    score += monograms[char2num(text[i])];
  return score;
}

double ic_score(char * text, int len)
{
  int freq[26];
  for(int j=0; j<26; j++)
    freq[j] = 0;
  for(int i=0; i<len; i++)
    freq[char2num(text[i])]++;
  double score = 0.0;
  for(int j=1; j<26; j++)
    score += freq[j-1] * freq[j];
  return score;
}

void init()
{
  for (int i=0; i < rotor_count; i++)
    for (int j=0; j < alphabet_size; j++)
      {
        rotor_fwd[i][j] = char2num(rotor_string[i][j]);
        rotor_rev[i][j] = index(rotor_string[i],num2char(j)) - rotor_string[i];
        notch[i][j] = index(notch_string[i], num2char(j)) != NULL;
      }

  for (int i=0; i < reflector_count; i++)
    for (int j=0; j < alphabet_size; j++)
      reflector[i][j] = char2num(reflector_string[i][j]);
}

void init_steckerbrett(const char * steckerbrett_string)
{
  for (int j=0; j < alphabet_size; j++)
    steckerbrett[j] = j;

  int plug_count = strlen(steckerbrett_string) / 2;
  
  for (int i=0; i < plug_count; i++)
    {
      int a = char2num(steckerbrett_string[2*i+0]);
      int b = char2num(steckerbrett_string[2*i+1]);
      steckerbrett[a] = b;
      steckerbrett[b] = a;
    }
}

void init_steckerbrett_direct(const char * steckerbrett_string)
{
  for (int j=0; j < alphabet_size; j++)
    steckerbrett[j] = char2num(steckerbrett_string[j]);
}

void init_walzen(int u, int a, int b, int c)
{
  if (opt_norenigma)
    {
      ukw = 3+u;
      walzenlage[0] = 8+a;
      walzenlage[1] = 8+b;
      walzenlage[2] = 8+c;
    }
  else
    {
      ukw = u;
      walzenlage[0] = a;
      walzenlage[1] = b;
      walzenlage[2] = c;
    }
}

void init_ring_grund(int a, int b, int c, int x, int y, int z)
{
  ringstellung[0] = a;
  ringstellung[1] = b;
  ringstellung[2] = c;
  grundstellung[0] = x;
  grundstellung[1] = y;
  grundstellung[2] = z;
}

char rotor_l(int x, int rotor_no)
{
  int y = grundstellung[rotor_no] - ringstellung[rotor_no];
  x = (x + alphabet_size + y) % alphabet_size;
  x = rotor_fwd[walzenlage[rotor_no]][x];
  x = (x + alphabet_size - y) % alphabet_size;
  return x;
}

char rotor_r(int x, int rotor_no)
{
  int y = grundstellung[rotor_no] - ringstellung[rotor_no];
  x = (x + alphabet_size + y) % alphabet_size;
  x = rotor_rev[walzenlage[rotor_no]][x];
  x = (x + alphabet_size - y) % alphabet_size;
  return x;
}

inline int mod26(int x)
{
  return (x+26)%26;
}

inline void step_rotors()
{
  if (notch[walzenlage[wheels-2]][grundstellung[wheels-2]])
    {
      grundstellung[wheels-3] = mod26(1+grundstellung[wheels-3]);
      grundstellung[wheels-2] = mod26(1+grundstellung[wheels-2]);
    }
  else if (notch[walzenlage[wheels-1]][grundstellung[wheels-1]])
    {
      grundstellung[wheels-2] = mod26(1+grundstellung[wheels-2]);
    }

  grundstellung[wheels-1] = mod26(1+grundstellung[wheels-1]);
}

inline int subst_rotors(int x)
{
  for (int r = wheels - 1; r >= 0; r--)
    x = rotor_l(x, r);

  x = reflector[ukw][x];

  for(int r = 0; r < wheels; r++)
    x = rotor_r(x, r);

  return x;
}

inline int substitute(int x)
{
  return steckerbrett[subst_rotors(steckerbrett[x])];
}

inline int step(int x)
{
  step_rotors();
  return substitute(x);
}

inline int step_precomputed(int x)
{
  step_rotors();
  return steckerbrett[subst_array
                      [mod26(grundstellung[0]-ringstellung[0])]
                      [mod26(grundstellung[1]-ringstellung[1])]
                      [mod26(grundstellung[2]-ringstellung[2])]
                      [steckerbrett[x]]];
}

void precompute()
{
  int r1 = 0;
  int r2 = 0;
  int r3 = 0;
  for (int g1 = 0; g1 < alphabet_size; g1++)
    for (int g2 = 0; g2 < alphabet_size; g2++)
      for (int g3 = 0; g3 < alphabet_size; g3++)
        {
          init_ring_grund(r1, r2, r3, g1, g2, g3);
          for (int x = 0; x < alphabet_size; x++)
            subst_array[g1][g2][g3][x] = subst_rotors(x);
        }
}

void setup_mapping(int textlength)
{
  if (textlength > maxlen)
    fatal("Ciphertext too long");

  /* set up mapping trough the rotors for each character in the ciphertext */
  for (int i=0; i<textlength; i++)
    {
      step_rotors();
      for (int j=0; j<26; j++)
        mapping[i][j] = subst_array
          [mod26(grundstellung[0]-ringstellung[0])]
          [mod26(grundstellung[1]-ringstellung[1])]
          [mod26(grundstellung[2]-ringstellung[2])]
          [j];
    }
}

inline void map(int len,
                unsigned char * source,
                unsigned char * map,
                unsigned char * dest)
{
  for (int i = 0; i < len; i++)
    dest[i] = map[26*i + source[i]];
}

inline int step_mapped(int i, int x)
{
  return steckerbrett[mapping[i][steckerbrett[x]]];
}

inline void decode(int textlength, char * ciphertext, char * plaintext)
{
  for (int i = 0; i < textlength; i++)
    plaintext[i] = num2char(step_mapped(i, num_ciphertext[i]));
  //plaintext[i] = num2char(step_precomputed(char2num(ciphertext[i])));
  //plaintext[i] = num2char(step(char2num(ciphertext[i])));
  plaintext[textlength] = 0;
}

inline void map16_step(unsigned char * source,
                       unsigned char * map,
                       unsigned char * dest)
{
  for (int i = 0; i < 16; i++)
    dest[i] = map[26*i+source[i]];
}

inline void map16_direct(unsigned char * source,
                         unsigned char * map,
                         unsigned char * dest)
{
  /* performs a mapping of the bytes in source (values 0-31)
     through the mapping at map (32 bytes) */

#if 0

  for (int i = 0; i < 16; i++)
    dest[i] = map[source[i]];

#else

  __m128i x, y, a, t0, t1, t2, m0, m1, t6, t7, t8, t9, t12;

  x = _mm_set_epi8(0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10,
                   0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10);

  y = _mm_set_epi8(0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
                   0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80);

  a  = _mm_loadu_si128((__m128i*)source);
  t0 = _mm_and_si128(a, x);
  t1 = _mm_slli_epi16(t0, 3);
  t2 = _mm_xor_si128(t1, y);
  m0 = _mm_or_si128(a, t1);
  m1 = _mm_or_si128(a, t2);

  t6  = _mm_loadu_si128((__m128i*)(map+00));
  t7  = _mm_loadu_si128((__m128i*)(map+16));
  t8  = _mm_shuffle_epi8(t6, m0);
  t9  = _mm_shuffle_epi8(t7, m1);
  t12 = _mm_or_si128(t8, t9);
  _mm_store_si128((__m128i*)(dest), t12);

  /*

    SIMD version

    load 16 byte values and perform substitution

    read 16 bytes of letters
    use as permute index
    perform permutation as with protein subst matrix
    pcmpgt/and, permute 0-15, permute 16-25, or together


    See dprofile_shuffle in search7.cc in swipe

    make masks for shuffle with bit 7 set (or reset) if byte >= 16

    x0 = load 16 bytes from source
    t0 = x0 and with 16 bytes of 0x10
    t1 = t0 shift left 3
    t2 = t1 xor with 16 bytes of 0x80
    m0 = x0 or t1
    m1 = x0 or t2

    x = _mm_set_epi8(0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10,
                     0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10);

    y = _mm_set_epi8(0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
                     0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80);

    a  = _mm_load_si128(source);
    t0 = _mm_and_si128(a, x);
    t1 = _mm_slli_epi16(t0, 3);
    t2 = _mm_xor_si128(t1, y);
    m0 = _mm_or_si128(a, t1);
    m1 = _mm_or_si128(a, t2);

    load 32 bytes of map
    shuffle first 16 bytes with mask m0
    shuffle seconds 16 bytes with mask m1
    or them together

    t6  = _mm_load_si128((__m128i*)(map+0));
    t7  = _mm_load_si128((__m128i*)(map+1));
    t8  = _mm_shuffle_epi8(t6, m0);
    t9  = _mm_shuffle_epi8(t7, m1);
    t12 = _mm_or_si128(t8,  t9);
    _mm_store_si128((__m128i*)(dprofile)+4*j,   t12);

    11 ops to load, map and save 16 bytes of data

  */
#endif

}

void showit(const char * msg, unsigned char * p)
{
#if 0
  fprintf(stderr, "%s:", msg);
  for(int i=0; i<16; i++)
    fprintf(stderr, " %2d", p[i]);
  fprintf(stderr, "\n");
#endif
}

inline void decode_num()
{
#if 0
  for (int i = 0; i < textlength; i++)
    num_plaintext[i] = step_mapped(i, num_ciphertext[i]);
#else
#if 0
  showit("cipher", num_ciphertext);
  for (int i = 0; i < textlength; i++)
    num_plaintext[i] = steckerbrett[mapping[i][steckerbrett[num_ciphertext[i]]]];
  showit("plain ", num_plaintext);
  fprintf(stderr, "\n");
#else

  for (int i = 0; i < textlength; i += 16)
    {
      unsigned char temp1[16];
      unsigned char temp2[16];
      showit("cipher", num_ciphertext+i);
      map16_direct(num_ciphertext+i, steckerbrett, temp1);
      showit("steck ", temp1);
      map16_step(temp1, ((unsigned char *)(&mapping[i])), temp2);
      showit("mapped", temp2);
      map16_direct(temp2, steckerbrett, num_plaintext+i);
      showit("plain ", num_plaintext+i);
      //      fprintf(stderr, "\n");
    }

#endif
#endif
}

double quadgram_score_decode(int textlength)
{
  /* This decode and scoring function uses 99% of the computation time
     when hill-climbing. */

  decode_num();

  /*

  load triplet scores
    load 16 bytes at adr
    load 16 bytes at adr+1
    load 16 bytes at adr+2
    unpack low and high of each of these
    shift 0, 5 or 10 bits left
    add/or together
    gather from triplet score table

  */

  double score = 0.0;
  for (int i=0; i<textlength-3; i++)
    score += quadgrams[num_plaintext[i]][num_plaintext[i+1]]
      [num_plaintext[i+2]][num_plaintext[i+3]];
  return score;
}

double trigram_score_decode(int textlength)
{
  decode_num();
  
  double score = 0.0;
  for (int i=0; i<textlength-2; i++)
    score += trigrams[num_plaintext[i]][num_plaintext[i+1]][num_plaintext[i+2]];
  return score;
}

double bigram_score_decode(int textlength)
{
  decode_num();

  double score = 0.0;
  for (int i=0; i<textlength-1; i++)
    score += bigrams[num_plaintext[i]][num_plaintext[i+1]];
  return score;
}

double monogram_score_decode(int textlength)
{
  decode_num();

  double score = 0.0;
  for (int i = 0; i < textlength; i++)
    score += monograms[num_plaintext[i]];
  return score;
}

double ic_score_decode(int textlength)
{
  int freq[26];
  for(int j=0; j<26; j++)
    freq[j] = 0;

  for (int i = 0; i < textlength; i++)
    freq[step_mapped(i, num_ciphertext[i])]++;

  double score = 0.0;
  for(int j=1; j<26; j++)
    score += freq[j-1] * freq[j];
  return score;
}

void showsteckerbrett()
{
#if 0
  for (int j=0; j<26; j++)
    putchar(num2char(steckerbrett[j]));
#else
  fprintf(stderr, "S:");
  for (int j=0; j<26; j++)
    if (steckerbrett[j] > j)
      fprintf(stderr, " %c%c", num2char(j), num2char(steckerbrett[j]));
#endif
}

void showconfig()
{
  fprintf(stderr, 
          "W: %c%d%d%d R: %c%c%c G: %c%c%c ",
          num2char(ukw + (opt_norenigma ? 10 : 0)),
          walzenlage[0] + (opt_norenigma ? -7 : 1),
          walzenlage[1] + (opt_norenigma ? -7 : 1),
          walzenlage[2] + (opt_norenigma ? -7 : 1),
          num2char(ringstellung[0]),
          num2char(ringstellung[1]),
          num2char(ringstellung[2]),
          num2char(grundstellung[0]),
          num2char(grundstellung[1]),
          num2char(grundstellung[2]));
  showsteckerbrett();
  fprintf(stderr, "\n");
}

struct subst_score_s
{
  int a;
  int b;
  double score;
} subst_scores[26*26];

int subst_score_comp(const void * x, const void * y)
{
  struct subst_score_s * i = (struct subst_score_s *) x;
  struct subst_score_s * j = (struct subst_score_s *) y;
  if (i->score < j->score)
    return +1;
  else if (i->score > j->score)
    return -1;
  else
    return 0;
}

void all_subst_score(int textlength)
{
  init_steckerbrett("");

#if 0
  double score;
  //score = ic_score_decode(textlength);
  //score = mgram_score_decode(textlength);
  //score = bigram_score_decode(textlength);
  //score = trigram_score_decode(textlength);
  //score = quadgram_score_decode(textlength);
  printf("Base score: %.4f\n", score);
#endif

  //  srandomdev();
  
  for (int a=0; a<26; a++)
    for (int b=0; b<26; b++)
      {
        double score = -1e30;
        if (a<b)
          {
            /* insert one plug */
            steckerbrett[a] = b;
            steckerbrett[b] = a;
            score = random() % 10000;
            //score = ic_score_decode(textlength);
            //score = mgram_score_decode(textlength);
            //score = bigram_score_decode(textlength);
            //score = trigram_score_decode(textlength);
            //score = quadgram_score_decode(textlength);
            steckerbrett[a] = a;
            steckerbrett[b] = b;
          }
        subst_scores[26*a+b].a = a;
        subst_scores[26*a+b].b = b;
        subst_scores[26*a+b].score = score;
      }

  /* sort */
  qsort(subst_scores, 26*26, sizeof(subst_score_s), subst_score_comp);
 
#if 0
  /* print */
  for (int a=0; a<26*26; a++)
    if (subst_scores[a].score > -1e30)
      printf("%c%c: %8.4f\n",
             num2char(subst_scores[a].a),
             num2char(subst_scores[a].b),
             subst_scores[a].score);
#endif
}

double score_iter(int iter, int textlength)
{
  double score = 0;
  
  switch(opt_scoring)
    {
    case 0:
      score = ic_score_decode(textlength);
      break;
      
    case 1:
      score = monogram_score_decode(textlength);
      break;
      
    case 2:
      score = bigram_score_decode(textlength);
      break;
      
    case 3:
      score = trigram_score_decode(textlength);
      break;
      
    case 4:
      score = quadgram_score_decode(textlength);            
      break;
      
    default:
      fatal("Illegal scoring type");
    }

  return score;
}
  
int count[26];
int order[26];

int compare(const void * x, const void * y)
{
  int a = count[*(int*)(x)];
  int b = count[*(int*)(y)];

  if (a<b)
    return +1;
  else if (a>b)
    return -1;
  else
    return 0;
}

void ciphertext_letterdist(int textlength, char * ciphertext)
{
  for(int j=0; j<26; j++)
    {
      count[j]=0;
      order[j] = j;
    }

  for (int i=0; i<textlength; i++)
    count[char2num(ciphertext[i])]++;

  qsort(order, 26, sizeof(int), compare);

#if 1
  fprintf(stderr, "Ciphertext letter order: ");
  for(int j=0; j<26; j++)
    fprintf(stderr, "%c", num2char(order[j]));
  fprintf(stderr, "\n");
#endif
}

double hillclimb(int textlength,
                 char * ciphertext,
                 char * plaintext)
{
  /* Try to find the optimal steckerbrett for the given other settings */

  int best_steckerbrett[26];
  memcpy(best_steckerbrett, steckerbrett, 26*sizeof(int));
  
  /* calc and sort best initial plug */
  //  all_subst_score(textlength);
  
  double best_score = -1e29;
  double last_best = -1e30;

  int iter = 1;

  while (best_score > last_best)
    {
      best_score = score_iter(iter, textlength);

      last_best = best_score;

      double switch_score = best_score;
      int switch_a = 0;
      int switch_b = 0;

      //#define SHOWHILLCLIMB

#ifdef SHOWHILLCLIMB
      fprintf(stderr, "  ");
      for(int b=1; b<26; b++)
        fprintf(stderr, "   %c", num2char(b));
      fprintf(stderr, "\n");
#endif
      for(int a=0; a<26; a++)
      {
#ifdef SHOWHILLCLIMB
        fprintf(stderr, "%c:", num2char(a));
        for(int b=1; b<a+1; b++)
          fprintf(stderr, "    ");
#endif
        for(int b=a+1; b<26; b++)
          {
            /* switch plugs */
            int x = steckerbrett[a];
            int y = steckerbrett[b];
            int xx = steckerbrett[x];
            int yy = steckerbrett[y];
            steckerbrett[x] = x;
            steckerbrett[y] = y;
            steckerbrett[a] = b;
            steckerbrett[b] = a;
            
            double score = score_iter(iter, textlength);
            
#ifdef SHOWHILLCLIMB
            fprintf(stderr, "%4.0f", (score - best_score)/10.0);
#endif
            
            if (score > switch_score)
              {
                switch_score = score;
                switch_a = a;
                switch_b = b;
              }
            
            /* restore plugs */
            steckerbrett[a] = x;
            steckerbrett[b] = y;
            steckerbrett[x] = xx;
            steckerbrett[y] = yy;
          }
#ifdef SHOWHILLCLIMB
        printf("\n");
#endif
        }

      if (switch_score - best_score > 0)
        {
          
          /* good move */

          int a = switch_a;
          int b = switch_b;
          
          /* switch plugs */
          int x = steckerbrett[a];
          int y = steckerbrett[b];
          steckerbrett[x] = x;
          steckerbrett[y] = y;
          steckerbrett[a] = b;
          steckerbrett[b] = a;
          
#ifdef SHOWHILLCLIMB
          fprintf(stderr,
                  "%2d %c%c Imp: %10.4f Score: %10.4f ",
                  iter,
                  num2char(a), num2char(b),
                  switch_score - best_score,
                  switch_score);
          showsteckerbrett();
          fprintf(stderr, "\n");
#endif
          
          best_score = switch_score;
        }

      iter++;
    }
  
  decode(textlength, ciphertext, plaintext);
  
#ifdef SHOWHILLCLIMB
  printf("Plaintext: %s\n", plaintext);
#endif
  return score_iter(0, textlength);
}




void bruteforce()
{
  int u_min, u_max;
  int w_min[3], w_max[3];
  int r_min[3], r_max[3];
  int g_min[3], g_max[3];

  if (opt_norenigma)
    {
      u_min = 0;
      u_max = 0;
    }
  else
    {
      if (opt_ukw[0] == '.')
        {
          u_min = 0;
          u_max = 2;
        }
      else
        u_min = u_max = char2num(opt_ukw[0]);
    }

  for(int i=0; i<3; i++)
    {
      if (opt_walzen[i] == '.')
        {
          w_min[i] = 0;
          w_max[i] = opt_maxwheel - 1;
        }
      else
        {
          w_min[i] = w_max[i] = opt_walzen[i] - '1';
        }

      if (opt_ringstellung[i] == '.')
        {
          r_min[i] = 0;
          r_max[i] = 25;
        }
      else
        {
          r_min[i] = r_max[i] = char2num(opt_ringstellung[i]);
        }

      if (opt_grundstellung[i] == '.')
        {
          g_min[i] = 0;
          g_max[i] = 25;
        }
      else
        {
          g_min[i] = g_max[i] = char2num(opt_grundstellung[i]);
        }
    }
     
  double best_score = -1e37;
  char best_plaintext[1025];

  for (int u1 = u_min; u1 <= u_max; u1++)
    for (int w1 = w_min[0]; w1 <= w_max[0]; w1++)
      for (int w2 = w_min[1]; w2 <= w_max[1]; w2++)
        for (int w3 = w_min[2]; w3 <= w_max[2]; w3++)
          if ((w1 != w2) && (w1 != w3) && (w2 != w3))
            {
              init_walzen(u1, w1, w2, w3);

              precompute();
              
              for (int r1 = r_min[0]; r1 <= r_max[0]; r1++)
                for (int r2 = r_min[1]; r2 <= r_max[1]; r2++)
                  for (int r3 = r_min[2]; r3 <= r_max[2]; r3++)
                    for (int g1 = g_min[0]; g1 <= g_max[0]; g1++)
                      for (int g2 = g_min[1]; g2 <= g_max[1]; g2++)
                        for (int g3 = g_min[2]; g3 <= g_max[2]; g3++)
                          {
                            init_ring_grund(r1, r2, r3, g1, g2, g3);

                            init_steckerbrett(opt_steckerbrett);

                            setup_mapping(textlength);
                            
                            double score;
                            if (opt_hillclimb)
                              {
                                score = hillclimb(textlength,
                                                  ciphertext,
                                                  plaintext);
                              }
                            else
                              {
                                decode(textlength, ciphertext, plaintext);
                                score = score_iter(0, textlength);
                              }
                            
                            if (score > best_score)
                              {
                                best_score = score;
                                strcpy(best_plaintext, plaintext);
#if 1
                                init_ring_grund(r1, r2, r3, g1, g2, g3);
                                fprintf(stderr, "%10.4f ", score);
                                showconfig();
#endif
                              }
                          }
            }
  strcpy(plaintext, best_plaintext);
}

void readciphertext()
{
  long len = 0;
  char buffer[maxlen+1];
  
  len = read(STDIN_FILENO, buffer, maxlen);
  
  int j=0;
  for(int i=0; i<len; i++)
    {
      char c = toupper(buffer[i]);
      if ((c>='A') && (c<='Z'))
        ciphertext[j++] = c;
    }
  ciphertext[j] = 0;
  textlength = j;
}

void readplaintext(char * filename)
{
  long len = 0;
  char buffer[maxlen+1];
  
  int fd = open(filename, O_RDONLY);
  if (fd < 0)
    fatal("Unable to open plaintext file");

  len = read(fd, buffer, maxlen);

  close(fd);
  
  int j=0;
  for(int i=0; i<len; i++)
    {
      char c = toupper(buffer[i]);
      if ((c>='A') && (c<='Z'))
        altplaintext[j++] = c;
    }
  altplaintext[j] = 0;

  if (textlength != j)
    fatal("Plaintext not same length as ciphertext");

  int identical = 0;
  for (int i=0; i<textlength; i++)
    if (plaintext[i] == altplaintext[i])
      identical++;

  fprintf(stderr,
          "%d of %d letters (%.1f%%) identical to given plaintext\n",
          identical,
          textlength,
          100.0 * identical / textlength);
}

void alltoupper(char * text)
{
  int len = strlen(text);
  for(int i=0; i<len; i++)
    text[i] = toupper(text[i]);
}

void version()
{
  printf("Enigma cipher tool version 1.0\n");
  printf("Copyright (C) 2017 Torbjørn Rognes\n");
  printf("\n");
}

void help()
{
  version();
  printf("Usage: enigma [OPTIONS]\n");
  printf("  -h           Show help information\n");
  printf("  -v           Show version information\n");
  printf("  -u X         Reflector (umkehrwalze) X (A-C, N or .) [.]\n");
  printf("  -w XYZ       Wheels (walzen) XYZ (1-8 or .) [...]\n");
  printf("  -x integer   Highest wheel number to use (1-8) [5]\n");
  printf("  -n           Use the Norway Enigma reflector (N) and wheels (1-5)\n");
  printf("  -r XYZ       Ring positions (ringstellung) XYZ (A-Z or .) [AA.]\n");
  printf("  -g XYZ       Start positions (grundstellung) XYZ (A-Z or .) [...]\n");
  printf("  -s ABYZ      Plugboard (steckerbrett) letter pairs (A-Z pairs) [none]\n");
  printf("  -c           Perform hill climbing to determine plugboard settings\n");
  printf("  -l language  Plaintext language (german, english, danish, french) [german]\n");
  printf("  -i           Use index of coincidence (IC) to determine plaintext score\n");
  printf("  -m           Use monogram statistics to determine plaintext score\n");
  printf("  -b           Use bigram statistics to determine plaintext score\n");
  printf("  -t           Use trigram statistics to determine plaintext score\n");
  printf("  -q           Use quadgram statistics to determine plaintext score [default]\n");
  printf("  -p filename  Name of file containing plaintext to compare result with\n");
  printf("\n");
  printf("Defaults are indicated in [square brackets].\n");
  printf("\n");
  printf("The ciphertext is read from standard input. The final plaintext is written\n");
  printf("to standard output.\n");
  printf("\n");
  printf("For the reflector, wheels, ring position and start position, a dot (.)\n");
  printf("works as a wild card, leaving it unspecified. When these settings are not\n");
  printf("specified, the program will try all combinations to find the settings\n");
  printf("resulting in the highest plaintext score. If asked for, a hill climbing\n");
  printf("algorithm will be used to try to determine the plugboard settings.\n");
  printf("\n");
}

void removespaces(char * p)
{
  char * q = p;
  while(char c = *p++)
    if (c != ' ')
      *q++ = c;
  *q=0;
}

int main(int argc, char * * argv)
{
  if (argc == 1)
    {
      help();
      exit(1);
    }

  /* set default arguments */
  opt_ukw = (char*) ".";
  opt_walzen = (char*)"...";
  opt_ringstellung = (char*) "AA.";
  opt_grundstellung = (char*) "...";
  opt_steckerbrett = (char*) "";
  opt_language = (char*) "german";
  opt_plaintext = 0;
  opt_maxwheel = 5;
  opt_hillclimb = 0;
  opt_scoring = 4;
  opt_logfilename = 0;
  opt_threads = 1;
  opt_norenigma = 0;

  /* get arguments */

  int c;
  while ((c = getopt(argc, argv, "u:w:r:g:s:p:l:imbtqxcvhn")) != -1)
    {
      switch (c)
        {
        case 'u':
          opt_ukw = optarg;
          alltoupper(opt_ukw);
          break;
        case 'w':
          opt_walzen = optarg;
          break;
        case 'r':
          opt_ringstellung = optarg;
          alltoupper(opt_ringstellung);
          break;
        case 'g':
          opt_grundstellung = optarg;
          alltoupper(opt_grundstellung);
          break;
        case 's':
          opt_steckerbrett = optarg;
          alltoupper(opt_steckerbrett);
          removespaces(opt_steckerbrett);
          break;
        case 'p':
          opt_plaintext = optarg;
          break;
        case 'i':
          opt_scoring = 0;
          break;
        case 'm':
          opt_scoring = 1;
          break;
        case 'b':
          opt_scoring = 2;
          break;
        case 't':
          opt_scoring = 3;
          break;
        case 'q':
          opt_scoring = 4;
          break;
        case 'c':
          opt_hillclimb = 1;
          break;
        case 'x':
          opt_maxwheel = atoi(optarg);
          break;
        case 'l':
          opt_language = optarg;
          break;
        case 'v':
          version();
          exit(0);
          break;
        case 'h':
          help();
          exit(0);
          break;
        case 'n':
          opt_norenigma = 1;
          break;
        default:
          fprintf(stderr, "\n");
          help();
          exit(1);
          break;
        }
    }

  argc -= optind;
  argv += optind;
  
  /* validate arguments */
  
  if (opt_norenigma)
    {
      if ((strlen(opt_ukw) != 1) ||
          (strspn(opt_ukw, "N.") != 1))
        fatal("Illegal ukw string (must be N or .)");

      if ((strlen(opt_walzen) != 3) ||
          (strspn(opt_walzen, "12345.") != 3))
        fatal("Illegal walzen string (must be 3 digits (1-5) or .)");

      if ((opt_maxwheel < 3) || (opt_maxwheel > 5))
        fatal("Illegal max wheel (must be 3 to 5)");
    }
  else
    {
      if ((strlen(opt_ukw) != 1) ||
          (strspn(opt_ukw, "ABC.") != 1))
        fatal("Illegal ukw string (must be A, B, C or .)");

      if ((strlen(opt_walzen) != 3) ||
          (strspn(opt_walzen, "12345678.") != 3))
        fatal("Illegal walzen string (must be 3 digits (1-8) or .)");

      if ((opt_maxwheel < 3) || (opt_maxwheel > 8))
        fatal("Illegal max wheel (must be 3-8)");
    }

  if ((strlen(opt_ringstellung) != 3) ||
      (strspn(opt_ringstellung, "ABCDEFGHIJKLMNOPQRSTUVWXYZ.") != 3))
    fatal("Illegal ringstellung string (must be 3 letters (A-Z) or .)");

  if ((strlen(opt_grundstellung) != 3) ||
      (strspn(opt_grundstellung, "ABCDEFGHIJKLMNOPQRSTUVWXYZ.") != 3))
    fatal("Illegal grundstellung string (must be 3 letters (A-Z) or .)");

  if ((strlen(opt_steckerbrett) > 26) ||
      (strspn(opt_steckerbrett, "ABCDEFGHIJKLMNOPQRSTUVWXYZ") <
       strlen(opt_steckerbrett)))
    fatal("Illegal steckerbrett string (must be up to 13 letter pairs)");


  /* read ciphertext */
  
  readciphertext();

  /* init */

  switch (opt_scoring)
    {
    case 0:
      break;
    case 1:
      monograms_read();
      break;
    case 2:
      bigrams_read();
      break;
    case 3:
      trigrams_read();
      break;
    case 4:
      quadgrams_read();
      break;
    }

  for(int i=0; i< textlength; i++)
    num_ciphertext[i] = char2num(ciphertext[i]);
  
  ciphertext_letterdist(textlength, ciphertext);

  init();
  init_steckerbrett("");

  /* try all combinations */

  bruteforce();

  /* write plaintext */

  fprintf(stderr, "\n");
  printf("%s\n", plaintext);

  /* read plaintext to compare to, if given */

  if (opt_plaintext)
    readplaintext(opt_plaintext);
}
