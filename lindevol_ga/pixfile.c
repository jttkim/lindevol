#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "jklib.h"

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"


FILE *pixel_file = NULL;


int open_pixfile(const char *mode)
{
  char m[32];

  sprintf(m, "%sb", mode);
  pixel_file = fopen(pixel_file_name, mode);
  if (pixel_file == NULL)
    return (-1);
  if (*mode == 'w')
  {
    fprintf(pixel_file, "RLE pixel file of %s\n", SIMPRGNAME);
    fprintf(pixel_file, "run: %s\n", simname);
    fprintf(pixel_file, "%ld\n", world_width);
    fprintf(pixel_file, "%ld\n", world_height);
  }
  return (0);
}


static unsigned char pixbyte(long x, long y)
{
  unsigned char b;
  long plant_no;

  plant_no = world[x][y].plant_no;
  if (plant_no > -1)
  {
    b = plant_no % 0x70;
    if (plant[plant_no].cell[world[x][y].cell_no].energy)
      b |= 0x80;
    return (b);
  }
  else
    return (0x7f);
}


/*
 * Coding scheme:
 * 0xf0 - 0xff: command codes, reserved.
 *    f0: Next 4 bytes are 32 bit repeat count for subsequent byte
 *    f8: Next 4 bytes are time step as 32 bit number
 *    f9: Next 4 bytes are ftell pos of preceding record
 *    fa: Next 4 bytes are length of data, followed by data
 *    fb: Next 4 bytes are length of huffman encoded block, followed
 *        by huffman encoded data
 * 0x00 - 0xef: different colors. Settings up to displaying program.
 *    Conventions: (1) Bit #7 is energy bit, i.e. x | 0x80 should be
 *        "energyrich" color, x & 0x7f should be corresponding
 *        "energyless" color
 *    (2) Background color is 0x7f by convention. Colors 0x70 - 0x7f
 *        cannot be used for cells, as their energyrich variants
 *        fall in the reserved range.
 */

int write_pixfile(void)
{
  static long last_fpos = -1, last_generation = -1;

  long x, y, count, j, fpos, rle_length = 0, huff_length;
  unsigned char b, b1, *rle, *huff;
  long i;


  if ((rle = (unsigned char *) malloc(world_width * world_height)) == NULL)
  {
    do_error("write_pixfile: Failed to allocate rle buffer");
    return (-1);
  }
  if ((huff = (unsigned char *) malloc(world_width * world_height)) == NULL)
  {
    do_error("write_pixfile: Failed to allocate huffman buffer");
    free(rle);
    return (-1);
  }
  fpos = ftell(pixel_file);
  fputc(0xf8, pixel_file);
  fwrite_int32array(&generation, 1, pixel_file);
  if (last_fpos != -1)
  {
    fputc(0xf9, pixel_file);
    fwrite_int32array(&last_fpos, 1, pixel_file);
  }
  if (generation != last_generation)
  {
    last_generation = generation;
    last_fpos = fpos;
  }
  count = 0;
  b1 = pixbyte(0, world_height - 1);
  for (y = world_height - 1; y >= 0; y--)
  {
    for (x = 0; x < world_width; x++)
    {
      b = pixbyte(x, y);
      if (b == b1)
	count++;
      else
      {
	if (count > 6)
	{
	  rle[rle_length++] = 0xf0;
	  for (j = 0; j < 4; j++)
	    rle[rle_length++] = ((unsigned long) count) >> ((3 - j) * 8);
	  rle[rle_length++] = b1;
	}
	else
	{
	  for (j = 0; j < count; j++)
	    rle[rle_length++] = b1;
	}
	b1 = b;
	count = 1;
      }
    }
  }
  if (count > 6)
  {
    rle[rle_length++] = 0xf0;
    for (j = 0; j < 4; j++)
      rle[rle_length++] = ((unsigned long) count) >> ((3 - j) * 8);
    rle[rle_length++] = b1;
  }
  else
  {
    for (j = 0; j < count; j++)
      rle[rle_length++] = b1;
  }
#ifdef NOHUFFMAN_ENCODE
  fputc(0xfa, pixel_file);
  fwrite_int32array(&rle_length, 1, pixel_file);
  fwrite(rle, 1, rle_length, pixel_file);
#else
  huff_length = huffman_encode_8bit(rle_length, rle, huff);
  fputc(0xfb, pixel_file);
  fwrite_int32array(&huff_length, 1, pixel_file);
  fwrite(huff, 1, huff_length, pixel_file);
  /* printf("g=%ld: rle array length: %ld, huffman array length: %ld (ratio: %f)\n", generation, rle_length, huff_length, (double) huff_length / rle_length); */
#endif
  free(rle);
  free(huff);
  return (0);
}

