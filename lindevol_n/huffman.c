/* $Id: huffman.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: huffman.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct
{
  size_t byte_no;
  unsigned short bit_no;
} HUFFMAN_BITSPEC;


typedef struct tag_htnode
{
  struct tag_htnode *n0, *n1;
  unsigned char v;
} HT_NODE;


typedef struct
{
  short num_bits;
  unsigned char c[32];
} H_CODE;


typedef struct tag_freq_listnode
{
  struct tag_freq_listnode *next;
  HT_NODE *ht_node;
  unsigned long f;
  unsigned char v;
} FREQ_LISTNODE;


static void get_frequencies(unsigned long l, const unsigned char *s, unsigned long *f)
{
  unsigned long i;

  for (i = 0; i < 256; i++)
    f[i] = 0;
  for (i = 0; i < l; i++)
    f[s[i]]++;
}


static void get_htcodes_rek(const HT_NODE *htree, H_CODE ht_code[256], unsigned char c[32], unsigned short nbits)
{
  unsigned short char_no = nbits / 8, bit_no = nbits % 8, i;
  unsigned long m;

  if (htree->n0)
  {
    m = 1U << bit_no;
    c[char_no] |= m;
    get_htcodes_rek(htree->n1, ht_code, c, nbits + 1);
    m ^= 0xff;
    c[char_no] &= m;
    get_htcodes_rek(htree->n0, ht_code, c, nbits + 1);
  }
  else
  {
    ht_code[htree->v].num_bits = nbits;
    for (i = 0; i < 32; i++)
      ht_code[htree->v].c[i] = c[i];
  }
}


static int cmp_freq(const void *p1, const void *p2)
{
  unsigned long f1 = ((FREQ_LISTNODE *) p1)->f, f2 = ((FREQ_LISTNODE *) p2)->f;

  return ((f1 < f2) ? -1 : ((f1 == f2) ? 0 : 1));
}

static void get_htcodes(const HT_NODE *htree, H_CODE ht_code[256])
{
  unsigned char c[32];
  int i, j;

  get_htcodes_rek(htree, ht_code, c, 0);
/*
  for (i = 0; i < 256; i++)
  {
    printf("%02x -> ", i);
    for (j = 0; j < ht_code[i].num_bits; j++)
    {
      if (ht_code[i].c[j / 8] & (1U << (j % 8)))
	printf("1");
      else
	printf("0");
    }
    printf("\n");
  }
*/
}


static void print_byte(unsigned char b)
{
  int j;

  for (j = 0; j < 8; j++)
  {
    if (b & (1U << j))
      printf("1");
    else
      printf("0");
  }
  printf(" ");
}


static void write_htcode(const H_CODE *htc, unsigned char *s, HUFFMAN_BITSPEC *bitspec)
{
  unsigned m;
  unsigned short nbits = htc->num_bits;
  int i;

  printf("%3u bits at byte = %5lu, bit = %1u: ", nbits, (unsigned long) bitspec->byte_no, bitspec->bit_no);
  if (bitspec->bit_no == 0)
    s[bitspec->byte_no] = 0;11
  m = ((1U << bitspec->bit_no) - 1) ^ 0xff;
  if ((nbits + bitspec->bit_no) <= 8)
  {
    m &= (1U << (bitspec->bit_no + nbits)) - 1;
    s[bitspec->byte_no] |= (htc->c[0] << bitspec->bit_no) & m;
    print_byte(s[bitspec->byte_no]);
    printf("\n");
    bitspec->byte_no += (bitspec->bit_no + nbits) / 8;
    bitspec->bit_no = (bitspec->bit_no + nbits) % 8;
    return;
  }
  s[bitspec->byte_no] |= (htc->c[0] << bitspec->bit_no) & m;
  print_byte(s[bitspec->byte_no]);
  nbits -= 8 - bitspec->bit_no;
  (bitspec->byte_no)++;
  s[bitspec->byte_no] = 0;
  i = 1;
  while (nbits >= 8)
  {
    s[bitspec->byte_no] = (htc->c[i - 1] >> (8 - bitspec->bit_no)) | (htc->c[i] << bitspec->bit_no);
    print_byte(s[bitspec->byte_no]);
    (bitspec->byte_no)++;
    s[bitspec->byte_no] = 0;
    i++;
    nbits -= 8;
  }
  if (nbits)
  {
    m = (1U << nbits) - 1;
    s[bitspec->byte_no] |= ((htc->c[i - 1] >> (8 - bitspec->bit_no)) | (htc->c[i] << bitspec->bit_no)) & m;
    print_byte(s[bitspec->byte_no]);
    printf(", m = ");
    print_byte(m);
  }
  bitspec->bit_no = nbits;
  printf("\n");
}


static void write_htree(const HT_NODE *htree, unsigned char *h, size_t *hpos)
{
  if (htree->n0)
  {
    h[*hpos] = 0;
    (*hpos)++;
    write_htree(htree->n0, h, hpos);
    write_htree(htree->n1, h, hpos);
  }
  else
  {
    h[*hpos] = 1;
    h[*hpos + 1] = htree->v;
    *hpos += 2;
  }
}


static int get_htree_rek(HT_NODE *htree, const unsigned char *h, size_t *hpos, size_t *p)
{
  size_t p1;

  switch (h[(*hpos)++])
  {
  case 0:
    p1 = *p;
    (*p)++;
    htree[p1].n0 = htree + *p;
    if (get_htree_rek(htree, h, hpos, p))
      return (-1);
    (*p)++;
    htree[p1].n1 = htree + *p;
    if (get_htree_rek(htree, h, hpos, p))
      return (-1);
    break;
  case 1:
    htree[*p].n0 = NULL;
    htree[*p].v = h[(*hpos)++];
    break;
  default:
    fprintf(stderr, "get_htree_rek: invalid type key %02x\n", h[*hpos - 1]);
    return (-1);
    break;
  }
  return (0);
}


static int get_htree(HT_NODE *htree, const unsigned char *h, size_t *hpos)
{
  size_t p = 0;

  return (get_htree_rek(htree, h, hpos, &p));
}


static size_t huffmandecode_byte(const HT_NODE *htree, unsigned long length, const unsigned char *h, unsigned char *buf,
        HUFFMAN_BITSPEC *bitspec)
{
  const HT_NODE *ht;
  size_t num_bytes = 0;

  ht = htree;
  while (num_bytes < length)
  {
    if (h[bitspec->byte_no] & (1U << (bitspec->bit_no)++))
      ht = ht->n1;
    else
      ht = ht->n0;
    if (bitspec->bit_no == 8)
    {
      (bitspec->byte_no)++;
      bitspec->bit_no = 0;
    }
    if (ht->n0 == NULL)
    {
      buf[num_bytes++] = ht->v;
      ht = htree;
    }
  }
  return (num_bytes);
}


static HT_NODE *construct_huffmantree_8bit(HT_NODE *ht, unsigned long frq[256])
{
  HT_NODE *htree;
  FREQ_LISTNODE fl[256], *flist, *fl1, *fl2;
  int i;

  for (i = 0;  i < 256; i++)
  {
    fl[i].f = frq[i];
    fl[i].v = i;
  }
  qsort(fl, 256, sizeof(FREQ_LISTNODE), cmp_freq);
  flist = fl;
  for (i = 0;  i < 255; i++)
    fl[i].next = fl + (i + 1);
  fl[255].next = NULL;
  for (i = 0;  i < 256; i++)
  {
    ht[i].n0 = NULL;
    ht[i].v = fl[i].v;
    fl[i].ht_node = ht + i;
  }
  htree = ht + 256;
  for (i = 0; i < 254; i++)
  {
    htree->n0 = flist->ht_node;
    htree->n1 = flist->next->ht_node;
    flist->next->f += flist->f;
    flist = flist->next;
    flist->ht_node = htree;
    htree++;
    if (flist->f > flist->next->f)
    {
      fl1 = fl2 = flist;
      flist = flist->next;
      while (fl1->next && (fl2->f > fl1->next->f))
        fl1 = fl1->next;
      fl2->next = fl1->next;
      fl1->next = fl2;
    }
  }
  htree->n0 = flist->ht_node;
  htree->n1 = flist->next->ht_node;
  return (htree);
}


static size_t huffmanencode_byte(const HT_NODE *htree, size_t s_length, const unsigned char *s, size_t h_length, unsigned char *h, HUFFMAN_BITSPEC *bitspec)
{
  size_t i, bn;
  H_CODE ht_code[256];

  get_htcodes(htree, ht_code);
  bn = bitspec->byte_no;
  for (i = 0; i < s_length; i++)
  {
    printf("%5lu: %02x, ", (unsigned long) i, s[i]);
    write_htcode(ht_code + s[i], h, bitspec);
    if (bitspec->byte_no > h_length - 32)
      return (0);
  }
  return (bitspec->byte_no - bn + 1);
}


size_t huffman_encode_8bit(size_t length, const unsigned char *s, unsigned char *h)
{
  HT_NODE ht[512], *htree;
  HUFFMAN_BITSPEC bitspec;
  unsigned long frq[256];
  unsigned long l, l1;
  int i;
  
  l1 = length;
  for (i = 0; i < 7; i++)
  {
    h[i] = l1 & 0xff;
    l1 >>= 8;
  }
  if (length < 800)
  {
    /* printf("\aNo compression < 800\n"); */
    h[i++] = 0;
    memcpy(h + i, s, length);
    return (length + i);
  }
  get_frequencies(length, s, frq);
  htree = construct_huffmantree_8bit(ht, frq);
  bitspec.byte_no = i;
  h[bitspec.byte_no++] = 1;
  write_htree(htree, h, &(bitspec.byte_no));
  bitspec.bit_no = 0;
  printf("writing tree at byte_no %lu, bit_no %u\n", bitspec.byte_no, bitspec.bit_no);
  l = huffmanencode_byte(htree, length, s, length - bitspec.byte_no, h, &bitspec);
  if (l == 0)
  {
    /* printf("\aNo compression\n"); */
    h[i++] = 0;
    memcpy(h + i, s, length);
    return (length + i);
  }
  return (bitspec.byte_no + 1);
}


size_t huffman_decode_8bit(const unsigned char *h, unsigned char *buf)
{
  HT_NODE ht[512];
  HUFFMAN_BITSPEC bitspec;
  unsigned long l, l1;
  int i;

  l = 0;
  for (i = 0; i < 7; i++)
  {
    l |= ((unsigned long) h[i]) << (8 * i);
  }
  /* printf("Expected length: %lu\n", l); */
  if (h[i++])
  {
    bitspec.byte_no = i;
    get_htree(ht, h, &(bitspec.byte_no));
    printf("obtained huffman tree, byte_no = %lu\n", (unsigned long) bitspec.byte_no);
    bitspec.bit_no = 0;
    /* printf("reading tree at byte_no %lu, bit_no %u\n", bitspec.byte_no, bitspec.bit_no); */
    l1 = huffmandecode_byte(ht, l, h, buf, &bitspec);
    if (l1 != l)
    {
      fprintf(stderr, "huffman_decode_8bit: Expected %lu bytes, decoded %lu bytes\n", l, l1);
      l = l1;
    }
    printf("decoded up to byte #%lu, bit #%u\n", (unsigned long) bitspec.byte_no, bitspec.bit_no);
  }
  else
  {
    /* printf("\aUncompressed code\n"); */
    memcpy(buf, h + 8, l);
  }
  return (l);
}


int main(int argc, char **argv)
{
  unsigned char x[100000], x1[100008], h[100008];
  FILE *f = NULL;
  size_t n, hn, n1;
  unsigned long i;

  if (argc > 1)
    f = fopen(argv[1], "rb");
  if (f)
  {
    n = fread(x, 1, 100000, f);
    fclose(f);
    hn = huffman_encode_8bit(n, x, h);
    printf("length: %lu, after compression: %lu\n", (unsigned long) n, (unsigned long) hn);
    n1 = huffman_decode_8bit(h, x1);
    printf("length after decoding: %lu\n", (unsigned long) n1);
    if (n1 != n)
      fprintf(stderr, "Decoded length: %lu, original length: %lu\n", (unsigned long) n1, (unsigned long) n);
    for (i = 0; i < n; i++)
    {
      if (x[i] != x1[i])
	fprintf(stderr, "diff at %6lu: original = %02x, decoded = %02x\n", i, x[i], x1[i]);
    }

    if (argc > 2)
    {
      f = fopen(argv[2], "wb");
      printf("writing decoded file to %s\n", argv[2]);
    }
    else
      f = stdout;
    fwrite(x1, 1, n1, f);
    if (f != stdout)
      fclose(f);
    if (argc > 3)
    {
      f = fopen(argv[3], "wb");
      printf("writing encoded file to %s\n", argv[3]);
    }
    else
      f = stdout;
    fwrite(h, 1, hn, f);
    if (f != stdout)
      fclose(f);
  }
  return (EXIT_SUCCESS);
}

