#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include "jklib.h"
#include "pttypes.h"
#include "ptlib.h"

#ifdef MEMDEBUG
#  include "memdebug.h"
#endif

#define free0(p) free(p); (p) = NULL

#define MAX_LINELENGTH 256


void getline(FILE *f, char *buf, int l)
{
  long i;

  buf[0] = '\0';
  while (!strlen(buf))
  {
    fgets(buf, l, f);
    if (feof(f) || ferror(f))
      break;
    for (i = strlen(buf) - 1; i >= 0; i--)
    {
      if ((buf[i] == '\n') || (buf[i] == '\r'))
        buf[i] = '\0';
    }
  }
}


void postscript_init(FILE *f)
{
  fprintf(f, "%%!\n%% prolog\n\n");
  fprintf(f, "/FSD { findfont exch scalefont def } bind def\n");
  fprintf(f, "\n%% end of prolog\n\n");
}


void postscript_pagesetup(FILE *f)
{
  fprintf(f, "%% beginning of page setup\n\n");
  fprintf(f, "save\n");
  fprintf(f, "72 300 div 72 300 div scale\n");
  fprintf(f, "/scratchstr 30 string def\n");
  fprintf(f, "\n%% end of page setup, beginning of page commands\n\n");
}


void postscript_showpage(FILE *f)
{
  fprintf(f, "showpage restore\n\n%% end of page\n\n");
}


long next_generation(FILE *f)
{
  char buf[MAX_LINELENGTH + 1];

  do
    getline(f, buf, MAX_LINELENGTH);
  while (((buf[0] != 'g' ) || (buf[1] != ' ')) && !feof(f) && !ferror(f));
  return (strtol(buf + 1, (char **) NULL, 10));
}


int main(int argc, char **argv)
{
  FILE *outfile, *psfile, *f;
  PHYLTREE true_tree, rec_tree;
  PHYL_PSINFO phyl_psinfo;
  long n_rectrees, i, j, d;
  long infile_g, truetree_g;
  long *hgram_correct = NULL, *hgram_false = NULL, *hgram_ctotal = NULL, *hgram_ftotal = NULL, *h = NULL;
  long max_ec;
  double truetree_height;
  long max_theight = 0, theight, old_mth;
  char buf[80], ps_str[520];
  char *outfile_name = NULL, *psfile_name = NULL, *hgramfile_basename = NULL;
  int oc;
  int print_rectrees = 0;

  extern char *optarg;

  const long left_margin = 200;
  const long page_width = 2200, page_height = 3100, page_bottom = 100;
  const long smalltree_h = 350, smalltree_w = 1800;
  const long truetree_h = 900, truetree_w = 1800;

  long tree_x, tree_y;

  phyl_init_tree(&true_tree);
  phyl_init_tree(&rec_tree);
  while ((oc = getopt(argc, argv, "rt:i:o:p:e:h")) != -1)
  {
    switch (oc)
    {
    case 'o':
      outfile_name = optarg;
      break;
    case 'p':
      psfile_name = optarg;
      break;
    case 'r':
      print_rectrees = 1;
      break;
    case 'e':
      hgramfile_basename = optarg;
      break;
    case 'h':
      printf("rtrchk -- treecheck with random trees\n");
      printf("command line arguments:\n");
      printf("-o <filename>: write topological tree distances to specified\n");
      printf("    file (default is stdout)\n");
      printf("-p <filename>: write postscript graphics showing correctly\n");
      printf("    reconstructed edges to specified file\n");
      printf("    (default is no postscript output)\n");
      printf("-r: include all reconstructed trees (with highlighted correct\n");
      printf("    edges) in postscript graphics\n");
      printf("-e <filename>: write distributions of correctly and falsely\n");
      printf("    edges to files with specified base filename. Full filenames\n");
      printf("    have the form filename-######[cf].ehi, where ###### is the\n");
      printf("    time step and c and f stand for correctly and falsesly\n");
      printf("    reconstructed edges respectively\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open %s for output -- exit\n", outfile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
    outfile = stdout;
  if (psfile_name)
  {
    if ((psfile = fopen(psfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open %s for postscript output -- exit\n", psfile_name);
      exit (EXIT_FAILURE);
    }
    postscript_init(psfile);
  }
  if (psfile && !print_rectrees)
  {
    postscript_pagesetup(psfile);
    tree_x = left_margin;
    tree_y = page_height;
  }
  truetree_g = 0;
  infile_g = 0;
  n_rectrees = 6000;
  if (phyl_rndbin_tree(50, &true_tree, 24114) < 0)
  {
    fprintf(stderr, "error creating random tree\n");
  }
  truetree_height = phyl_treeheight(true_tree.root);
  theight = floor(truetree_height);
  if (hgramfile_basename)
  {
    truetree_height = phyl_treeheight(true_tree.root);
    theight = floor(truetree_height);
    if (theight > max_theight)
    {
      old_mth = max_theight;
      max_theight = theight;
      if (hgram_ctotal)
      {
        if ((h = (long *) realloc(hgram_ctotal, max_theight * sizeof(long))))
          hgram_ctotal = h;
        else
        {
          fprintf(stderr, "Failed to realloc total correct array\n");
          free0(hgram_ctotal);
        }
      }
      else
      {
        if ((hgram_ctotal = (long *) malloc(max_theight * sizeof(long))) == NULL)
          fprintf(stderr, "Failed to malloc total correct array\n");
      }
      if (hgram_ctotal)
      {
        for (j = old_mth; j < max_theight; j++)
          hgram_ctotal[j] = 0;
      }
      if (hgram_ftotal)
      {
        if ((h = (long *) realloc(hgram_ftotal, max_theight * sizeof(long))))
          hgram_ftotal = h;
        else
        {
          fprintf(stderr, "Failed to realloc total false array\n");
          free0(hgram_ftotal);
        }
      }
      else
      {
        if ((hgram_ftotal = (long *) malloc(max_theight * sizeof(long))) == NULL)
          fprintf(stderr, "Failed to malloc total false array\n");
      }
      if (hgram_ftotal)
      {
        for (j = old_mth; j < max_theight; j++)
          hgram_ftotal[j] = 0;
      }
    }
    /* printf("tree height: %ld, max tree height: %ld\n", theight, max_theight); */
    if ((h = (long *) malloc(theight * sizeof(long))) == NULL)
      fprintf(stderr, "Failed to allocate temp array\n");
    if (h)
    {
      if ((hgram_correct = (long *) malloc(theight * sizeof(long))) == NULL)
        fprintf(stderr, "Failed to allocate correct array\n");
      if ((hgram_false = (long *) malloc(theight * sizeof(long))) == NULL)
        fprintf(stderr, "Failed to allocate false array\n");
      if (hgram_correct)
      {
        for (j = 0; j < theight; j++)
          hgram_correct[j] = 0;
      }
      if (hgram_false)
      {
        for (j = 0; j < theight; j++)
          hgram_false[j] = 0;
      }
    }
  }
  if (psfile && print_rectrees)
  {
    postscript_pagesetup(psfile);
    tree_x = left_margin;
    tree_y = page_height;
  }
  for (i = 0; i < n_rectrees; i++)
  {
    phyl_set_edges(&true_tree, PHYLEDGINF_NONE, -1);
    if (phyl_rndbin_tree(50, &rec_tree, i) < 0)
    {
      phyl_free_tree(&true_tree);
      fprintf(stderr, "error creating random reconstructed tree\n");
      continue;
    }
    d = phyl_topotreedist(&true_tree, &rec_tree);
    if (d >= 0)
    {
      if (hgramfile_basename)
      {
        if (h)
        {
          if (hgram_correct && (phyl_lengthtype_distr(&true_tree, PHYLEDGINF_IDEDGE, h, theight, truetree_height) == 0))
          {
            for (j = 0; j < theight; j++)
              hgram_correct[j] += h[j];
          }
          if (hgram_false && (phyl_lengthtype_distr(&true_tree, PHYLEDGINF_DIFFEDGE, h, theight, truetree_height) == 0))
          {
            for (j = 0; j < theight; j++)
              hgram_false[j] += h[j];
          }
        }
      }
      fprintf(outfile, "%ld %ld\n", truetree_g, d);
      if (psfile && print_rectrees)
      {
        if (tree_y - smalltree_h < page_bottom)
        {
          postscript_showpage(psfile);
          tree_x = left_margin;
          tree_y = page_height;
          postscript_pagesetup(psfile);
        }
        phyl_psinfo.fontname = "Courier";
        phyl_psinfo.min_fontsize = 10.0;
        phyl_psinfo.max_fontsize = 80.0;
        phyl_psinfo.label = NULL;
        phyl_psinfo.label_start = 0.0;
        phyl_psinfo.tic_length = 20.0;
        phyl_psinfo.label_fontsize = 50.0;
        phyl_psinfo.linewidth = 1.0;
        phyl_psinfo.print_leafnames = 1;
        phyl_inf2thick(&rec_tree, PHYLEDGINF_IDEDGE, 9.0);
        phyl_pstree(psfile, &true_tree, tree_x, tree_y - smalltree_h, smalltree_w / 2 - 20, smalltree_h, &phyl_psinfo);
        tree_x += smalltree_w / 2 + 20;
        phyl_pstree(psfile, &rec_tree, tree_x, tree_y - smalltree_h, smalltree_w / 2 - 20, smalltree_h, &phyl_psinfo);
        tree_x += smalltree_w + 100;
        if (tree_x + smalltree_w > page_width)
        {
          tree_y -= smalltree_h + 100;
          tree_x = left_margin;
        }
      }
    }
    else
    {
      switch (d)
      {
      case PHYLERR_DIFFNUMLEAVES:
        fprintf(stderr, "different number of leaves, generation %ld, tree #%ld\n", truetree_g, i);
        break;
      case PHYLERR_INCOMPATLEAVES:
        fprintf(stderr, "different sets of leaves: generation %ld, tree #%ld\n", truetree_g, i);
        break;
      default:
        fprintf(stderr, "error code %ld during tree distance computation\n", d);
        break;
      }
    }
    phyl_free_tree(&rec_tree);
  }
  if (psfile)
  {
    if (tree_x > left_margin)
    {
      tree_y -= smalltree_h + 100;
      tree_x = left_margin;
    }
    if (tree_y - truetree_h < page_bottom)
    {
      postscript_showpage(psfile);
      postscript_pagesetup(psfile);
      tree_x = left_margin;
      tree_y = page_height;
    }
    max_ec = phyl_max_edgecounter(&true_tree);
    sprintf(buf, "true tree at generation #%ld (%ld rec. trees, max. hits=%ld)", truetree_g, n_rectrees, max_ec);
    fprintf(psfile, "gsave /Courier findfont 50 scalefont setfont\n");
    fprintf(psfile, "%ld %ld moveto %s show grestore\n", left_margin, tree_y - 50, ps_string(buf, ps_str));
    phyl_psinfo.fontname = "Courier";
    phyl_psinfo.min_fontsize = 10.0;
    phyl_psinfo.max_fontsize = 100.0;
    phyl_psinfo.label = "generation";
    phyl_psinfo.label_start = truetree_g - theight;
    phyl_psinfo.tic_length = 20.0;
    phyl_psinfo.label_fontsize = 50.0;
    phyl_psinfo.linewidth = 1.0;
    phyl_psinfo.print_leafnames = 1;
    phyl_counter2thick(&true_tree, 700.0 / phyl_num_leaves(true_tree.root) / max_ec);
    phyl_pstree(psfile, &true_tree, tree_x, tree_y - truetree_h, truetree_w, truetree_h - 50, &phyl_psinfo);
    tree_y -= truetree_h + 100;
    if (print_rectrees)
    {
      postscript_showpage(psfile);
      tree_y = page_height;
    }
  }
  if (hgramfile_basename && h)
  {
    free0(h);
    if (hgram_correct)
    {
      sprintf(buf, "%s-%05ldc.ehi", hgramfile_basename, truetree_g);
      if ((f = fopen(buf, "w")))
      {
        for (j = 0; j < theight; j++)
          fprintf(f, "%ld %ld\n", j, hgram_correct[j]);
        fclose(f);
        if (hgram_ctotal)
        {
          for (j = 0; j < theight; j++)
            hgram_ctotal[j] += hgram_correct[j];
        }
        free0(hgram_correct);
      }
      else
        fprintf(stderr, "Failed to open %s\n", buf);
    }
    if (hgram_false)
    {
      sprintf(buf, "%s-%05ldf.ehi", hgramfile_basename, truetree_g);
      if ((f = fopen(buf, "w")))
      {
        for (j = 0; j < theight; j++)
          fprintf(f, "%ld %ld\n", j, hgram_false[j]);
        fclose(f);
        if (hgram_ftotal)
        {
          for (j = 0; j < theight; j++)
            hgram_ftotal[j] += hgram_false[j];
        }
        free0(hgram_false);
      }
      else
        fprintf(stderr, "Failed to open %s\n", buf);
    }
  }
  phyl_free_tree(&true_tree);
#ifdef MEMDEBUG
  print_MemdebugStatistics();
#endif
  if (outfile_name)
    fclose(outfile);
  if (psfile)
  {
    if (!print_rectrees)
      postscript_showpage(psfile);
    fclose(psfile);
  }
  if (hgram_ctotal)
  {
    sprintf(buf, "%s-ctotal.ehi", hgramfile_basename);
    if ((f = fopen(buf, "w")))
    {
      for (j = 0; j < max_theight; j++)
        fprintf(f, "%ld %ld\n", j, hgram_ctotal[j]);
      fclose(f);
    }
    else
      fprintf(stderr, "Error opening %s\n", buf);
  }
  if (hgram_ftotal)
  {
    sprintf(buf, "%s-ftotal.ehi", hgramfile_basename);
    if ((f = fopen(buf, "w")))
    {
      for (j = 0; j < max_theight; j++)
        fprintf(f, "%ld %ld\n", j, hgram_ftotal[j]);
      fclose(f);
    }
    else
      fprintf(stderr, "Error opening %s\n", buf);
  }
  return (0);
}

