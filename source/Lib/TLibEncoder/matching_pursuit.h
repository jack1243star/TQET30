//////////////////////////////////////////////////////////////////////////////
//
// INCLUDE SECTION
//

#ifndef CHY_MATCHING_PURSUIT_H
#define CHY_MATCHING_PURSUIT_H

#ifdef CHYMP_STATIC
#define CHYMP_DEF static
#else
#define CHYMP_DEF extern
#endif

#include <stdio.h> // printf
#include <stdlib.h> // abs
#include <string.h> // memset
#include <math.h>   // pow
#include "TLibCommon/TComTrQuant.cpp"

#ifdef __cplusplus
extern "C" {
#endif

  CHYMP_DEF double chymp_matching_pursuit(int size,
                                        TCoeff *bb,
                                        Pel *mask,
                                        double ep,
                                        TCoeff *coeff,
                                        TCoeff *rec);
  CHYMP_DEF void chymp_init();
  CHYMP_DEF void chymp_free(void);
  CHYMP_DEF void chymp_mask(Pel *pixels, UInt width, UInt height, UInt stride, Pel *mask);
  CHYMP_DEF void chymp_dump(Pel *src, UInt width, UInt height, UInt stride, char *filename);
  CHYMP_DEF void chymp_dump2(TCoeff *src, UInt width, UInt height, UInt stride, char *filename);
  CHYMP_DEF void chymp_dumpP2(const Pel *src, UInt width, UInt height, UInt stride, char *filename);
  CHYMP_DEF void chymp_dumpP22(TCoeff *src, UInt width, UInt height, UInt stride, char *filename);
  CHYMP_DEF void chymp_copy(TCoeff *src, UInt width, UInt height, UInt stride, UInt left, UInt top);
  CHYMP_DEF void chymp_log(char *log);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // CHY_MATCHING_PURSUIT_H

//////////////////////////////////////////////////////////////////////////////
//
// IMPLEMENTATION SECTION
//

#ifdef CHY_MATCHING_PURSUIT_IMPLEMENTATION

static TCoeff *basis4 = NULL;
static TCoeff *basis8 = NULL;
static TCoeff *basis16 = NULL;
static TCoeff *basis32 = NULL;
static FILE *logfile = NULL;
static Pel *reco = NULL;
static int global_width;
static int global_height;

static void chymp__dct(TCoeff *input, TCoeff *output, int width, int height)
{
  xTrMxN(8, input, output, width, height, false, 15);
}

static void chymp__idct(TCoeff *input, TCoeff *output, int width, int height)
{
  xITrMxN(8, input, output, width, height, false, 15);
}

static void chymp__create_basis(int size, TCoeff *block)
{
  TCoeff input[1024];
  unsigned int i;
  unsigned int length = size * size;

  TCoeff *current_block = block;

  for (i = 0; i < length; i++)
  {
    memset(input, 0, sizeof(input));
    input[i] = 8192;
    chymp__idct(input, current_block, size, size);
    current_block += size * size;
  }
}

CHYMP_DEF void chymp_init()
{
  global_width = 1024;
  global_height = 768;
  basis4 = (TCoeff *)malloc(4 * 4 * 4 * 4 * sizeof(TCoeff));
  basis8 = (TCoeff *)malloc(8 * 8 * 8 * 8 * sizeof(TCoeff));
  basis16 = (TCoeff *)malloc(16 * 16 * 16 * 16 * sizeof(TCoeff));
  basis32 = (TCoeff *)malloc(32 * 32 * 32 * 32 * sizeof(TCoeff));
  chymp__create_basis(4, basis4);
  chymp__create_basis(8, basis8);
  chymp__create_basis(16, basis16);
  chymp__create_basis(32, basis32);
  logfile = fopen("D:\\Temp\\Thesis\\Project\\mp_log.txt", "w");
  reco = (Pel*)calloc(global_width*global_height, sizeof(Pel));
}

CHYMP_DEF void chymp_free(void)
{
  free(basis4);
  free(basis8);
  free(basis16);
  free(basis32);
  fclose(logfile);
  chymp_dumpP2(reco, global_width, global_height, global_width, "D:\\Temp\\Thesis\\Project\\mp_rec.pgm");
  free(reco);
}

CHYMP_DEF void chymp_mask(Pel *pixels, UInt width, UInt height, UInt stride, Pel *mask)
{
  const UInt length = width * height;
  Double sum = 0.0;
  Double avg = 0.0;

  Pel* pixel_ptr;

  // calculate the average value for the block
  pixel_ptr = pixels;
  for (UInt i = 0; i < height; i++)
  {
    for (UInt j = 0; j < width; j++)
    {
      sum += pixel_ptr[j];
    }
    pixel_ptr += stride;
  }
  avg = sum / length;

  // classify the pixels as dark or light
  pixel_ptr = pixels;
  for (UInt i = 0; i < height; i++)
  {
    for (UInt j = 0; j < width; j++)
    {
      mask[j + i * width] = pixel_ptr[j] > avg ? 1 : 0;
    }
    pixel_ptr += stride;
  }
}

CHYMP_DEF double chymp_matching_pursuit(int size,
                                      TCoeff *bb,
                                      Pel *mask,
                                      double ep,
                                      TCoeff *coeff,
                                      TCoeff *rec)
{
  // Size of block
  unsigned int length = size * size;

  // Number of bits in mask
  unsigned int ml;

  // Mean square error
  double err;
  double prev_err;

  // Loop index
  unsigned int i, index;

  // Masked block
  TCoeff xt[1024];

  // Temporary block
  TCoeff tmp[1024];

  // Most significant coefficient
  TCoeff sigcoeff;

  // Output coefficients
  TCoeff *yt = coeff;

  // Basis functions
  TCoeff *basis, *basis_set;
  switch (size)
  {
  case 4:
    basis_set = basis4;
    break;
  case 8:
    basis_set = basis8;
    break;
  case 16:
    basis_set = basis16;
    break;
  case 32:
    basis_set = basis32;
    break;
  default:
    exit(1);
    break;
  }

  // Initialize inputs
  memset(yt, 0, size * size * sizeof(TCoeff));
  for (i = 0; i < length; i++)
  {
    xt[i] = bb[i] * mask[i];
  }
  ml = 0;
  for (i = 0; i < length; i++)
  {
    ml = ml + mask[i];
  }

  // Calculate mean square error
  err = 0;
  for (i = 0; i < length; i++)
  {
    err = err + (bb[i] * mask[i]) * (bb[i] * mask[i]);
  }
  err = err / ml;
#ifdef CHYMP_DEBUG
  printf("    err = %f\n", err);
  printf("- - - - - - - -\n");
#endif

  int tries = 0;
  do
  {
    assert(tries < 1000);
    tries++;
    prev_err = err;
#ifdef CHYMP_DEBUG
    printf("  tries = %u\n", tries);
#endif

    // Find the most significant coefficient
    chymp__dct(xt, tmp, size, size);
    sigcoeff = INT32_MIN;
    for (i = 0; i < length; i++)
    {
      if (abs(tmp[i]) > sigcoeff)
      {
        index = i;
        sigcoeff = abs(tmp[i]);
      }
    }
    sigcoeff = tmp[index];
#ifdef CHYMP_DEBUG
    printf("     ix = %u\n", index);
    printf("maximum = %d\n", sigcoeff);
#endif

    // Subtract scaled basis from masked block
    basis = basis_set + index * size * size;
    for (i = 0; i < length; i++)
    {
      xt[i] -= sigcoeff * basis[i] * mask[i] / 8192;
    }

    // Reconstruct using the coefficients found so far
    coeff[index] += sigcoeff;
    chymp__idct(coeff, rec, size, size);

    // Recalculate error
    err = 0;
    for (i = 0; i < length; i++)
    {
      int diff = (rec[i] - bb[i]) * mask[i];
      err += diff * diff;
    }
    err /= ml;
    if (prev_err < err)
      break;
#ifdef CHYMP_DEBUG
    printf("    err = %f\n", err);
    printf("- - - - - - - -\n");
#endif
  } while (err > ep);
  return err;
}

CHYMP_DEF void chymp_dump(Pel *src, UInt width, UInt height, UInt stride, char *filename)
{
  unsigned int i, j;
  FILE* f = fopen(filename, "w");
  fprintf(f, "return {\n");
  fprintf(f, "width = %d,\n", width);
  fprintf(f, "height = %d,\n", height);
  fprintf(f, "data={\n");
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      fprintf(f, "%d,", src[i*stride + j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "}\n");
  fprintf(f, "}\n");
  fclose(f);
}

CHYMP_DEF void chymp_dump2(TCoeff *src, UInt width, UInt height, UInt stride, char *filename)
{
  unsigned int i, j;
  FILE* f = fopen(filename, "w");
  fprintf(f, "return {\n");
  fprintf(f, "width = %d,\n", width);
  fprintf(f, "height = %d,\n", height);
  fprintf(f, "data={\n");
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      fprintf(f, "%d,", src[i*stride + j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "}\n");
  fprintf(f, "}\n");
  fclose(f);
}

CHYMP_DEF void chymp_dumpP2(const Pel *src, UInt width, UInt height, UInt stride, char *filename)
{
  unsigned int i, j;
  FILE* f = fopen(filename, "w");
  fprintf(f, "P2\n");
  fprintf(f, "%d %d\n", width, height);
  fprintf(f, "255\n");
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      fprintf(f, "%d ", src[i*stride + j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

CHYMP_DEF void chymp_dumpP22(TCoeff *src, UInt width, UInt height, UInt stride, char *filename)
{
  unsigned int i, j;
  FILE* f = fopen(filename, "w");
  fprintf(f, "P2\n");
  fprintf(f, "%d %d\n", width, height);
  fprintf(f, "255\n");
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      fprintf(f, "%d ", src[i*stride + j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

CHYMP_DEF void chymp_copy(TCoeff *src, UInt width, UInt height, UInt stride, UInt left, UInt top)
{
  unsigned int i, j;
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      reco[(i+top)*global_width+left+j] = src[i*stride + j];
    }
  }
}

CHYMP_DEF void chymp_log(char *log)
{
  fprintf(logfile, "%s\n", log);
}

#endif
