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
  fprintf(f, "restore showpage\n\n%% end of page\n\n");
}


long next_generation(FILE *f)
{
  char buf[MAX_LINELENGTH + 1];

  do
    getline(f, buf, MAX_LINELENGTH);
  while (((buf[0] != 'g' ) || (buf[1] != ' ')) && !feof(f) && !ferror(f));
  return (strtol(buf + 2, (char **) NULL, 10));
}


void write_pageheader(FILE *f, const char *truetree_filename, const char *infile_name, long page_no, long page_width, long page_height, long left_margin)
{
  char buf[250], ps_str[520];

  postscript_pagesetup(f);
  sprintf(buf, "treecmp of %s and %s, page %ld", truetree_filename, infile_name, page_no);
  fprintf(f, "gsave /Helvetica-Bold findfont 50 scalefont setfont\n");
  fprintf(f, "%ld %ld moveto %s show\n", left_margin, page_height - 50, ps_string(buf, ps_str));
  fprintf(f, "0 setlinecap 5 setlinewidth newpath %ld %ld moveto %ld %ld lineto stroke grestore\n",
          left_margin, page_height - 70, page_width, page_height - 70);
}


void alloc_histograms(long theight, long max_theight, long **hgram_ct, long **hgram_ft, long **hgram_c, long **hgram_f)
{
  long *h, j;

  if (theight > max_theight)
  {
    if (*hgram_ct)
    {
      if ((h = (long *) realloc(*hgram_ct, theight * sizeof(long))))
        *hgram_ct = h;
      else
      {
        fprintf(stderr, "Failed to realloc total correct array\n");
        free0(*hgram_ct);
      }
    }
    else
    {
      if ((*hgram_ct = (long *) malloc(theight * sizeof(long))) == NULL)
        fprintf(stderr, "Failed to malloc total correct array\n");
    }
    if (*hgram_ct)
    {
      for (j = max_theight; j < theight; j++)
        (*hgram_ct)[j] = 0;
    }
    if (*hgram_ft)
    {
      if ((h = (long *) realloc(*hgram_ft, theight * sizeof(long))))
        *hgram_ft = h;
      else
      {
        fprintf(stderr, "Failed to realloc total false array\n");
        free0(*hgram_ft);
      }
    }
    else
    {
      if ((*hgram_ft = (long *) malloc(theight * sizeof(long))) == NULL)
        fprintf(stderr, "Failed to malloc total false array\n");
    }
    if (*hgram_ft)
    {
      for (j = max_theight; j < theight; j++)
        (*hgram_ft)[j] = 0;
    }
  }
  /* printf("tree height: %ld, max tree height: %ld\n", theight, theight); */
  if ((h = (long *) malloc(theight * sizeof(long))) == NULL)
    fprintf(stderr, "Failed to allocate temp array\n");
  if (h)
  {
    if ((*hgram_c = (long *) malloc(theight * sizeof(long))) == NULL)
      fprintf(stderr, "Failed to allocate correct array\n");
    if ((*hgram_f = (long *) malloc(theight * sizeof(long))) == NULL)
      fprintf(stderr, "Failed to allocate false array\n");
    if (*hgram_c)
    {
      for (j = 0; j < theight; j++)
        (*hgram_c)[j] = 0;
    }
    if (*hgram_f)
    {
      for (j = 0; j < theight; j++)
        (*hgram_f)[j] = 0;
    }
  }
}

int main(int argc, char **argv)
{
  FILE *truetree_file, *infile, *outfile, *psfile, *f;
  PHYLTREE true_tree, rec_tree;
  PHYL_PSINFO phyl_psinfo;
  int truetree_unrooted = 0, rectree_unrooted = 0;
  long truetree_g, infile_g, n_truetrees, n_rectrees, i, j, d;
  long *hgram_correct = NULL, *hgram_false = NULL, *hgram_ctotal = NULL, *hgram_ftotal = NULL, *h = NULL;
  long max_ec;
  double truetree_height;
  long max_theight = 0, theight;
  char buf[80], ps_str[520];
  char *truetree_filename = NULL, *infile_name = NULL, *outfile_name = NULL, *psfile_name = NULL;
  char *hgramfile_basename = NULL, *epsfile_basename = NULL;

  int oc;
  int print_rectrees = 0;

  extern char *optarg;

  const long left_margin = 200;
  const long page_width = 2200, page_height = 3100, page_bottom = 100;
  const long smalltree_h = 350, smalltree_w = 1800;
  const long truetree_h = 900, truetree_w = 1800;

  long tree_x, tree_y;
  long page_no;

  phyl_init_tree(&true_tree);
  phyl_init_tree(&rec_tree);
  while ((oc = getopt(argc, argv, "hru:x:t:i:o:p:e:")) != -1)
  {
    switch (oc)
    {
    case 't':
      truetree_filename = optarg;
      break;
    case 'i':
      infile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 'p':
      psfile_name = optarg;
      break;
    case 'x':
      epsfile_basename = optarg;
      break;
    case 'r':
      print_rectrees = 1;
      break;
    case 'e':
      hgramfile_basename = optarg;
      break;
    case 'u':
      for (i = 0; optarg[i]; i++)
      {
        switch (optarg[i])
        {
        case 't':
          truetree_unrooted = 1;
          break;
        case 'T':
          truetree_unrooted = 0;
          break;
        case 'r':
          rectree_unrooted = 1;
          break;
        case 'R':
          rectree_unrooted = 0;
          break;
        default:
          fprintf(stderr, "Code %c in tree style string unknown -- ignored\n", optarg[i]);
          break;
        }
      }
      break;
    case 'h':
      printf("treecmp -- compare true trees to reconstructed ones\n");
      printf("command line arguments:\n");
      printf("-t <filename>: read true trees from specified file (mandatory)\n");
      printf("-i <filename>: read reconstructed trees from specified file\n");
      printf("    (default is stdin)\n");
      printf("-o <filename>: write topological tree distances to specified\n");
      printf("    file (default is stdout)\n");
      printf("-p <filename>: write postscript graphics showing correctly\n");
      printf("    reconstructed edges to specified file\n");
      printf("    (default is no postscript output)\n");
      printf("-x <basename>: write encapsulated postscript graphics showing\n");
      printf("    correctly reconstructed edges to files with specified\n");
      printf("    basenames. Full filenames are basename-#####.eps where\n");
      printf("    ##### is the time step.\n");
      printf("-r: include all reconstructed trees (with highlighted correct\n");
      printf("    edges) in postscript graphics\n");
      printf("-e <basename>: write distributions of correctly and falsely\n");
      printf("    edges to files with specified base filename. Full filenames\n");
      printf("    have the form basename-######[cf].ehi, where ###### is the\n");
      printf("    time step and c and f stand for correctly and falsesly\n");
      printf("    reconstructed edges respectively\n");
      printf("-u <string> Control whether trees are drawn in \"unrooted style\"\n");
      printf("    (as opposed to vertical dendrogram style which is default\n");
      printf("    Codes are:\n");
      printf("    t -> draw true trees unrooted,\n");
      printf("    T -> draw true trees as vertical dendrograms\n");
      printf("    r -> draw reconstructed trees unrooted,\n");
      printf("    R -> draw true trees as vertical dendrograms\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (truetree_filename == NULL)
  {
    fprintf(stderr, "No true tree file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  if ((truetree_file = fopen(truetree_filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open true tree file %s -- exit\n", truetree_filename);
    exit (EXIT_FAILURE);
  }
  if (infile_name)
  {
    if ((infile = fopen(infile_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open %s for input -- exit\n", infile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    infile = stdin;
    infile_name = "stdin";
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
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  if (psfile_name)
  {
    if ((psfile = fopen(psfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open %s for postscript output -- exit\n", psfile_name);
      exit (EXIT_FAILURE);
    }
    postscript_init(psfile);
  }
  page_no = 1;
  if (psfile && !print_rectrees)
  {
    write_pageheader(psfile, truetree_filename, infile_name, page_no++, page_width, page_height, left_margin);
    tree_x = left_margin;
    tree_y = page_height - 90;
  }
  truetree_g = next_generation(truetree_file);
  infile_g = next_generation(infile);
  while ((!feof(truetree_file)) && (!feof(infile)))
  {
    if (ferror(truetree_file))
    {
      fprintf(stderr, "error reading from %s\n", truetree_filename);
      break;
    }
    if (ferror(infile))
    {
      fprintf(stderr, "error reading from %s\n", infile_name);
      break;
    }
    if (feof(truetree_file) || feof(infile))
    {
      if (!feof(infile))
        fprintf(stderr, "generation %ld and subsequent ones missing from %s\n", infile_g, truetree_filename);
      else if (!feof(truetree_file))
        fprintf(stderr, "generation %ld and subsequent ones missing from %s\n", truetree_g, infile_name);
      break;
    }
    if (truetree_g < infile_g)
    {
      fprintf(stderr, "generation %ld is missing in %s\n", truetree_g, infile_name);
      truetree_g = next_generation(truetree_file);
      continue;
    }
    if (infile_g < truetree_g)
    {
      fprintf(stderr, "generation %ld is missing in %s\n", infile_g, truetree_filename);
      infile_g = next_generation(infile);
      continue;
    }
    getline(truetree_file, buf, 80);
    n_truetrees = strtol(buf, NULL, 10);
    getline(infile, buf, 80);
    n_rectrees = strtol(buf, NULL, 10);
    if (n_truetrees != 1)
    {
      fprintf(stderr, "generation %ld: %ld trees in %s -- skipping\n", truetree_g, n_truetrees, truetree_filename);
      truetree_g = next_generation(truetree_file);
      infile_g = next_generation(infile);
      continue;
    }

    /* printf("g = %ld, n = %ld\n", truetree_g, n_truetrees); */

    if (phyl_read_tree(truetree_file, &true_tree) < 0)
    {
      fprintf(stderr, "error reading tree from %s, g = %ld, n = %ld\n", truetree_filename, truetree_g, i);
      continue;
    }
    truetree_height = phyl_treeheight(true_tree.root);
    theight = floor(truetree_height);
    if (hgramfile_basename)
      alloc_histograms(theight, max_theight, &hgram_ctotal, &hgram_ftotal, &hgram_correct, &hgram_false);
    if (theight > max_theight)
      max_theight = theight;
    if (psfile && print_rectrees)
    {
      write_pageheader(psfile, truetree_filename, infile_name, page_no++, page_width, page_height, left_margin);
      tree_x = left_margin;
      tree_y = page_height - 90;
    }
    phyl_set_edges(&true_tree, PHYLEDGINF_NONE, 0);
    for (i = 0; i < n_rectrees; i++)
    {
      phyl_set_edges(&true_tree, PHYLEDGINF_NONE, -1);
      if (phyl_read_tree(infile, &rec_tree) < 0)
      {
        fprintf(stderr, "error reading tree from %s, g = %ld, n = %ld\n", infile_name, infile_g, i);
        continue;
      }
/*
      if (rec_tree.lengthinfo_complete)
        printf("complete branch length info found in %s\n", infile_name);
      else
        printf("no branch length info found in %s\n", infile_name);
      if (infile_name)
        printf("***** tree from %s @ g = %ld, n = %ld *****\n", infile_name, infile_g, i);
      else
        printf("***** tree @ g = %ld, n = %ld *****\n", infile_g, i);
      print_tree(&rec_tree);
*/
      d = phyl_topotreedist(&true_tree, &rec_tree);
      if (d >= 0)
      {
        if (hgramfile_basename)
        {
          if ((h = (long *) malloc(theight * sizeof(long))) != NULL)
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
            free0(h);
          }
        }
        fprintf(outfile, "%ld %ld\n", truetree_g, d);
        if (psfile && print_rectrees)
        {
          if (tree_y - smalltree_h < page_bottom)
          {
            postscript_showpage(psfile);
            write_pageheader(psfile, truetree_filename, infile_name, page_no++, page_width, page_height, left_margin);
            tree_x = left_margin;
            tree_y = page_height - 90;
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
          phyl_psinfo.attrlist = NULL;
          phyl_inf2thick(&true_tree, PHYLEDGINF_IDEDGE, 9.0);
          if (truetree_unrooted)
            phyl_ps_utree(psfile, &true_tree, tree_x, tree_y - smalltree_h, smalltree_w / 2 - 20, smalltree_h, &phyl_psinfo);
          else
            phyl_pstree(psfile, &true_tree, tree_x, tree_y - smalltree_h, smalltree_w / 2 - 20, smalltree_h, &phyl_psinfo);
          tree_x += smalltree_w / 2 + 20;
          if (rectree_unrooted)
            phyl_ps_utree(psfile, &rec_tree, tree_x, tree_y - smalltree_h, smalltree_w / 2 - 20, smalltree_h, &phyl_psinfo);
          else
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
      /* printf("topological tree distance @ g = %ld, n = %ld: %ld\n", truetree_g, i, d); */
      phyl_free_tree(&rec_tree);
    }
    max_ec = phyl_max_edgecounter(&true_tree);
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
        write_pageheader(psfile, truetree_filename, infile_name, page_no++, page_width, page_height, left_margin);
        tree_x = left_margin;
        tree_y = page_height - 90;
      }
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
      if (max_ec)
        phyl_counter2thick(&true_tree, 700.0 / phyl_num_leaves(true_tree.root) / max_ec);
      else
        phyl_set_thickness(&true_tree, 1.0);
      phyl_psinfo.print_leafnames = 1;
      phyl_psinfo.attrlist = NULL;
      if (truetree_unrooted)
        phyl_ps_utree(psfile, &true_tree, tree_x, tree_y - truetree_h, truetree_w, truetree_h - 50, &phyl_psinfo);
      else
        phyl_pstree(psfile, &true_tree, tree_x, tree_y - truetree_h, truetree_w, truetree_h - 50, &phyl_psinfo);
      tree_y -= truetree_h + 100;
      if (print_rectrees)
      {
        postscript_showpage(psfile);
        tree_y = page_height;
      }
    }
    if (epsfile_basename)
    {
      sprintf(buf, "%s-%05ld.eps", epsfile_basename, truetree_g);
      if ((f = fopen(buf, "w")))
      {
        fprintf(f, "%%!PS-Adobe-3.0 EPSF-3.0\n");
        fprintf(f, "%%%%Creator: treechk v. %s\n", __DATE__);
        fprintf(f, "%%%%Title: %s\n", buf);
        fprintf(f, "%%%%BoundingBox: 0 0 504 360\n");
        fprintf(f, "%%%%EndComments\n");
        phyl_psinfo.fontname = "Courier";
        phyl_psinfo.min_fontsize = 3.0;
        phyl_psinfo.max_fontsize = 20.0;
        phyl_psinfo.label = "generation";
        phyl_psinfo.label_start = truetree_g - theight;
        phyl_psinfo.tic_length = 5.0;
        phyl_psinfo.label_fontsize = 12.0;
        phyl_psinfo.linewidth = 0.2;
        if (max_ec)
          phyl_counter2thick(&true_tree, 160.0 / phyl_num_leaves(true_tree.root) / max_ec);
        else
          phyl_set_thickness(&true_tree, 0.2);
        phyl_psinfo.print_leafnames = 1;
        phyl_psinfo.attrlist = NULL;
        if (truetree_unrooted)
          phyl_ps_utree(f, &true_tree, 0.0, 0.0, 504.0, 360.0, &phyl_psinfo);
        else
          phyl_pstree(f, &true_tree, 0.0, 0.0, 504.0, 360.0, &phyl_psinfo);
        fprintf(f, "%%%%EOF\n");
        fclose(f);
      }
      else
        fprintf(stderr, "Failed to open %s\n", buf);
    }
    if (hgramfile_basename)
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
    truetree_g = next_generation(truetree_file);
    infile_g = next_generation(infile);
#ifdef MEMDEBUG
    print_MemdebugStatistics();
#endif
  }
  fclose(truetree_file);
  if (infile != stdin)
    fclose(infile);
  if (outfile != stdout)
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

