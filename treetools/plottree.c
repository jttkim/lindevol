/*
 * $Id: plottree.c,v 1.7 2000/03/15 12:27:11 kim Exp $
 *
 * $Log: plottree.c,v $
 * Revision 1.7  2000/03/15 12:27:11  kim
 * -B was ignored for eps, bug fixed
 *
 * Revision 1.6  2000/02/10 19:17:40  kim
 * Fixed offset bug in EPS generation
 *
 * Revision 1.5  2000/01/28 18:07:31  kim
 * fixed divide by zero bug in thickness hack
 *
 * Revision 1.4  2000/01/24 00:57:37  kim
 * adapted code for plotting edges with negative thickness in red, added
 * some sanity checks.
 *
 * Revision 1.3  2000/01/20 01:29:36  kim
 * Added -T option to allow reading line thickness values from lengths of
 * "bootstrap" tree. Very kludgy, no sanity checks(!!)
 *
 * Revision 1.2  2000/01/18 01:44:17  kim
 * added CVS tags
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif


#include "jklib.h"
#include "ptlib.h"


#define MERGE_NONE      0
#define MERGE_THRESHOLD 1
#define MERGE_MAX       2


typedef struct tag_legendline
{
  const char *line;
  struct tag_legendline *next;
} LEGENDLINE;


int is_simfile(FILE *f)
{
  long pos;
  char buf[256], *c;

  pos = ftell(f);
  fgets(buf, 256, f);
  fseek(f, pos, SEEK_SET);
  if (buf[0] != 'g')
    return (0);
  strtol(buf + 1, &c, 10);
  return (buf != c);
}


void postscript_init(FILE *f)
{
  fprintf(f, "%%!\n%% prolog\n\n");
  fprintf(f, "/FSD { findfont exch scalefont def } bind def\n");
  fprintf(f, "/showright\n");
  fprintf(f, "{ /lstr exch def\n");
  fprintf(f, "lstr stringwidth pop sub\n");
  fprintf(f, "exch moveto lstr show\n");
  fprintf(f, "} def\n");
  fprintf(f, "\n%% end of prolog\n\n");
}


void postscript_pagesetup(FILE *f, int landscape)
{
  fprintf(f, "%% beginning of page setup\n\n");
  fprintf(f, "save\n");
  fprintf(f, "72 300 div 72 300 div scale\n");
  if (landscape)
    fprintf(f, "90 rotate 0 -2000 translate\n");
  fprintf(f, "/scratchstr 30 string def\n");
  fprintf(f, "\n%% end of page setup, beginning of page commands\n\n");
}


void postscript_showpage(FILE *f)
{
  fprintf(f, "restore showpage\n\n%% end of page\n\n");
}


void write_pageheader(FILE *f, const char *text, long page_no, double page_width, double page_height, double left_margin, int landscape)
{
  char buf[250], ps_str[520];
  time_t t;

  t = time(NULL);
  postscript_pagesetup(f, landscape);
  sprintf(buf, "%s, page %ld", text, page_no);
  fprintf(f, "gsave /Helvetica-Bold findfont 50 scalefont setfont\n");
  fprintf(f, "%f %f moveto %s show\n", left_margin, page_height - 50.0, ps_string(buf, ps_str));
  strftime(buf, 250, "%d %b %Y", localtime(&t));
  fprintf(f, "%f %f %s showright\n", page_height - 50.0, left_margin + page_width, ps_string(buf, ps_str));
  fprintf(f, "0 setlinecap 10 setlinewidth newpath %f %f moveto %f %f lineto stroke grestore\n",
          left_margin, page_height - 70.0, left_margin + page_width, page_height - 70.0);
}


int ps_generation(FILE *infile, const char *infile_name, FILE *f, long generation, long num_trees,
        long boxsize, int unrooted_style, long *page_no, const PHYL_LEAFATTRIBUTE *attrlist)
{
  long i, page_pos, left_margin = 200;
  char buf[256], ps_str[520];
  PHYLTREE phyltree;
  PHYL_PSINFO phyl_psinfo;
  int ret_code;

  phyl_init_tree(&phyltree);
  page_pos = 0;
  for (i = 0; i < num_trees; i++)
  {
    if ((ret_code = phyl_read_tree(infile, &phyltree)) < 0)
    {
      if (i)
        postscript_showpage(f);
      return (ret_code);
    }
    if (phyltree.num_leaves == 0)
    {
      if (feof(infile))
      {
        fprintf(stderr, "Tree #%ld and subsequent ones (up to %ld) missing\n", i, num_trees - 1);
        return (-1);
      }
    }
    else if (((phyltree.num_leaves == 1) && (page_pos < 150))
            || ((phyltree.num_leaves > 1) && (page_pos < boxsize + 150)))
    {
      if (i)
        postscript_showpage(f);
      page_pos = 3110;
      write_pageheader(f, infile_name, (*page_no)++, 1800, 3200, left_margin, 0);
    }
    else
    {
      fprintf(f, "5 setlinewidth newpath %ld %ld moveto %ld %ld lineto stroke\n",
              left_margin, page_pos - 50, left_margin + 1800, page_pos - 50);
      page_pos -= 50;
    }
    sprintf(buf, "generation #%ld, tree #%ld", generation, i);
    fprintf(f, "gsave /Courier findfont 50 scalefont setfont\n");
    fprintf(f, "%ld %ld moveto %s show grestore\n", left_margin, page_pos - 50, ps_string(buf, ps_str));
    printf("tree g=%ld, #%ld with %ld leaves\n", generation, i, phyltree.num_leaves);
    if (phyltree.num_leaves == 1)
    {
      fprintf(f, "gsave /Courier findfont 50 scalefont setfont\n");
      sprintf(buf, "single leaf tree, leaf = %s", phyltree.root->descendant[0]->name);
      fprintf(f, "%ld %ld moveto %s show grestore\n", left_margin, page_pos - 100, ps_string(buf, ps_str));
      page_pos -= 100;
    }
    else
    {
      phyl_psinfo.fontname = "Courier-Bold";
      phyl_psinfo.spec_fontname = "Courier-Oblique";
      phyl_psinfo.min_fontsize = 10.0;
      phyl_psinfo.max_fontsize = 30.0;
      phyl_psinfo.label_fontsize = 30.0;
      phyl_psinfo.label = "gen.";
      phyl_psinfo.tic_length = 20.0;
      phyl_psinfo.label_start = generation - phyl_treeheight(phyltree.root);
      phyl_psinfo.linewidth = 1.0;
      phyl_psinfo.print_leafnames = 1;
      phyl_psinfo.dotted_nodes = 0;
      phyl_psinfo.tree_height = -1.0;
      phyl_psinfo.print_edgelengths = 0;
      phyl_psinfo.leafnames_at_max = 0;
      phyl_psinfo.attrlist = (PHYL_LEAFATTRIBUTE *) attrlist;
      phyl_psinfo.angle_min = 0.0;
      phyl_psinfo.angle_limit = 0.0;
      phyl_set_thickness(&phyltree, 1.0);
      if (unrooted_style)
        phyl_ps_utree(f, &phyltree, left_margin, page_pos - boxsize, 1800, boxsize - 75, &phyl_psinfo);
      else
        phyl_pstree(f, &phyltree, left_margin, page_pos - boxsize, 1800, boxsize - 75, &phyl_psinfo);
      page_pos -= boxsize;
    }
    phyl_free_tree(&phyltree);
  }
  if (i)
    postscript_showpage(f);
  return (0);
}


LEGENDLINE *append_legendline(LEGENDLINE *l, const char *str)
{
  LEGENDLINE *new_l = malloc(sizeof(LEGENDLINE));

  if (new_l == NULL)
  {
    fprintf(stderr, "append_legendline: malloc failed\n");
    return (l);
  }
  new_l->line = str;
  new_l->next = l;
  return (new_l);
}


double legend_height(const LEGENDLINE *legendline, double fontheight)
{
  double h = 0.0;
  const LEGENDLINE *l;

  for (l = legendline; l; l = l->next)
    h += fontheight;
  return (h);
}


int write_legend(FILE *f, const LEGENDLINE *legendline, double fontheight, double x0, double y0)
{
  double y;
  const LEGENDLINE *l;
  char ps_str[1024];

  y = y0 - legend_height(legendline, fontheight);
  fprintf(f, "save\n");
  fprintf(f, "/Courier findfont %f scalefont setfont\n", fontheight);
  for (l = legendline; l; l = l->next)
  {
    fprintf(f, "%f %f moveto %s show\n", x0, y, ps_string(l->line, ps_str));
    y += fontheight;
  }
  fprintf(f, "restore\n");
  return (0);
}


int main(int argc, char **argv)
{
  PHYLTREE phyltree, bootstrap_tree;
  PHYL_LEAFATTRIBUTE *attrlist = NULL;
  PHYL_PSINFO phyl_psinfo;
  int unrooted_style = 0, print_leafnames = PHYL_LEAVES_HORIZONTAL, merge_code = MERGE_NONE, dotted_nodes = 0, bs_thickness_from_lengths = 0;
  int leafnames_at_max = 0, print_edgelengths = 0;
  char *infile_name = NULL, *outfile_name = NULL, *leafattrfile_name = NULL, *bootstrapfile_name = NULL;
  FILE *infile, *outfile, *leafattrfile, *bootstrapfile = NULL;
  int const_length = 0, landscape = 0;
  int bootstrap_thick = 0;
  int ret_code;
  char buf[FILENAME_MAX + 256];
  long generation, num_trees, page_no;
  double page_width, page_height;
  double linewidth = 1.0, fontsize = -1.0, merge_threshold = 0, angle_min = 0.0, angle_limit = 0.0;
  int eps = 0, min_clip = 0;
  double min_length, tree_height = -1.0;
  LEGENDLINE *legendline = NULL;
  double legend_fontheight = 42.0, legend_y0;

  int oc;
  extern char *optarg;

  while ((oc = getopt(argc, argv, "hue1rnkBdTL:i:o:a:f:w:l:m:s:t:b:c:y:")) != -1)
  {
    switch (oc)
    {
    case 'T' :
      bs_thickness_from_lengths = 1;
      break;
    case 'L':
      legendline = append_legendline(legendline, optarg);
      break;
    case 'y':
      tree_height = strtod(optarg, NULL);
      break;
    case 'd':
      dotted_nodes = 1;
      break;
    case 'n':
      leafnames_at_max = 1;
      break;
    case 'k':
      print_edgelengths = 1;
      break;
    case 'e':
      eps = 1;
      break;
    case 'i':
      infile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 'a':
      leafattrfile_name = optarg;
      break;
    case 'b':
      bootstrapfile_name = optarg;
      break;
    case 'u':
      unrooted_style = 1;
      break;
    case '1':
      const_length = 1;
      break;
    case 'r':
      landscape = 1;
      break;
    case 'f':
      fontsize = strtod(optarg, (char **) NULL) * 300.0 / 72.0;
      if (fontsize <= 0)
      {
        fprintf(stderr, "font height %s is not valid\n", optarg);
        fontsize = -1.0;
      }
      break;
    case 'w':
      linewidth = strtod(optarg, (char **) NULL) * 300.0 / 72.0;
      if (linewidth <= 0)
      {
        fprintf(stderr, "line width %s is not valid\n", optarg);
        linewidth = 1.0;
      }
      break;
    case 's':
      angle_min = strtod(optarg, (char **) NULL);
      /* printf("angle_min set to %f\n", angle_min); */
      break;
    case 't':
      angle_limit = strtod(optarg, (char **) NULL);
      /* printf("angle_limit set to %f\n", angle_limit); */
      break;
    case 'l':
      switch (*optarg)
      {
      case 'n':
        print_leafnames = PHYL_LEAVES_NONE;
        break;
      case 'h':
        print_leafnames = PHYL_LEAVES_HORIZONTAL;
        break;
      case 'r':
        print_leafnames = PHYL_LEAVES_RADIAL;
        break;
      }
      break;
    case 'm':
      if (!strcmp(optarg, "max"))
        merge_code = MERGE_MAX;
      else
      {
        merge_code = MERGE_THRESHOLD;
        merge_threshold = strtod(optarg, NULL);
      }
      break;
    case 'c':
      min_clip = 1;
      min_length = strtod(optarg, NULL);
      break;
    case 'B':
      bootstrap_thick = 1;
      break;
    case 'h':
      printf("plottree -- programs to plot trees in PHYLIP format\n");
      printf("and collections of trees from stb files\n");
      printf("Commandline options:\n");
      printf("-i <filename>: Specify input file\n");
      printf("-o <filename>: Specify output file\n");
      printf("-u: Plot trees in \"unrooted style\"\n");
      printf("    default: vertical dendrogram style\n");
      printf("-d: Draw a dot at each node\n");
      printf("-a <filename>: Specify file for leaf attributes\n");
      printf("-b <filename>: Specify file of trees computed from bootstrap samples\n");
      printf("-B: Display bootstrap values by edge thickness\n");
      printf("-f <number>: Specify font size (in PostScript units)\n");
      printf("-w <number>: Specify line width for edges (in PostScript units)\n");
      printf("-l [n|h|r|]: Set leaf printing to none/horizontal/radial\n");
      printf("-m [max | <number>]: Merge nodes connected by edges shorter\n");
      printf("    than threshold specified by number. \"max\" specifies\n");
      printf("    threshold as maximal edge length in tree.\n");
      printf("-c <num>: Force branch lengths to be <num> at minimum.\n");
      printf("-e: Produce encapsulated postscript output\n");
      printf("-1: Set all branch lengths to 1.0\n");
      printf("-r: Produce landscape graphics (PS only, no effect with -e)\n");
      printf("-n: Avoid overprinting graph when there are negative edge lengths\n");
      printf("-k: Label edges with their length\n");
      printf("-h: Print this help and exit\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (infile_name)
  {
    if ((infile = fopen(infile_name, "r")) == NULL)
    {
      fprintf(stderr, "failed to open %s for input -- exit\n", infile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    infile_name = "stdin";
    infile = stdin;
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "failed to open %s for output -- exit\n", outfile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    outfile_name = "stdout";
    outfile = stdout;
  }
  if (leafattrfile_name)
  {
    if ((leafattrfile = fopen(leafattrfile_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open leaf attribute file %s -- exit\n", leafattrfile_name);
      exit (EXIT_FAILURE);
    }
    if ((ret_code = phyl_read_leafattrs(leafattrfile, &attrlist)) < 0)
    {
      fprintf(stderr, "Error #%d reading attribute list -- exit\n", ret_code);
      phyl_free_attrlist(attrlist);
      exit (EXIT_FAILURE);
    }
    phyl_compute_leafcolors(attrlist);
  }
  phyl_init_tree(&phyltree);
  phyl_init_tree(&bootstrap_tree);
  page_no = 1;
  if (is_simfile(infile))
  {
    if (eps)
    {
      fprintf(stderr, "*** EPS output for simfiles not yet implemented ***\n");
    }
    else
    {
      postscript_init(outfile);
      while (!feof(infile))
      {
        while (fgets(buf, 256, infile) && (buf[0] == '\n'))
          ;
        if (feof(infile))
          break;
        if (buf[0] != 'g')
        {
          fprintf(stderr, "*** Error: Generation label missing (g = %ld)\n", generation);
          fprintf(stderr, "*** %s", buf);
          break;
        }
        generation = strtol(buf + 1, NULL, 10);
        while (fgets(buf, 256, infile) && (buf[0] == '\n'))
          ;
        if (feof(infile))
          break;
        num_trees = strtol(buf, NULL, 10);
        if ((ret_code = ps_generation(infile, infile_name, outfile, generation, num_trees, 950, unrooted_style, &page_no, attrlist)) < 0)
        {
          fprintf(stderr, "ps_generation returned %d\n", ret_code);
          break;
        }
      }
    }
  }
  else
  {
    if (!eps)
      postscript_init(outfile);
    while (!feof(infile))
    {
      if ((ret_code = phyl_read_tree(infile, &phyltree)) < 0)
      {
        fprintf(stderr, "Error %d while reading tree\n", ret_code);
      }
      else if (phyltree.num_leaves > 0)  /* ignore "empty tree" at EOF */
      {
        switch (merge_code)
        {
        case MERGE_MAX:
          merge_threshold = phyl_max_edgelength(&phyltree);
          /* fallthrough intended! */
        case MERGE_THRESHOLD:
          if ((ret_code = phyl_merge_shortbranches(&phyltree, merge_threshold)) < 0)
          {
            fprintf(stderr, "Error #%d merging short tree branches\n", ret_code);
            phyl_free_tree(&phyltree);
            continue;
          }
          break;
        case MERGE_NONE:
          break;
        default:
          fprintf(stderr, "Internal error: unknown merge code %d\n", merge_code);
        }
        switch (merge_code)
        {
        case MERGE_MAX:
          sprintf(buf, "%s (merged at max.)", infile_name);
          break;
        case MERGE_THRESHOLD:
          sprintf(buf, "%s (merged at %f)", infile_name, merge_threshold);
          break;
        case MERGE_NONE:
          sprintf(buf, "%s", infile_name);
          break;
        default:
          sprintf(buf, "Internal error: unknown merge code %d", merge_code);
          fprintf(stderr, "Internal error: unknown merge code %d\n", merge_code);
        }
        if (const_length)
          phyl_set_constlength(&phyltree, 1.0);
        if (min_clip)
          phyl_set_minlength(&phyltree, min_length);
        phyl_psinfo.bootstrap_numtrees = 0;
        if (bootstrapfile_name)
        {
          if ((bootstrapfile = fopen(bootstrapfile_name, "r")) == NULL)
          {
            fprintf(stderr, "Failed to open bootstrap file %s -- no bootstrap computation\n", bootstrapfile_name);
          }
          else
          {
            while (!feof(bootstrapfile) && !ferror(bootstrapfile))
            {
              if ((ret_code = phyl_read_tree(bootstrapfile, &bootstrap_tree)) < 0)
              {
                fprintf(stderr, "Error %d while reading bootstrapped tree\n", ret_code);
              }
              else if (bootstrap_tree.num_leaves > 0)
              {
                phyl_psinfo.bootstrap_numtrees++;
                if ((ret_code = phyl_topotreedist(&phyltree, &bootstrap_tree)) < 0)
                  fprintf(stderr, "Error %d in computation of topological tree distance\n", ret_code);
		if (bs_thickness_from_lengths && (ret_code > 0))
		  fprintf(stderr, "non-identical topologies in edge thickness setting\n");
              }
              phyl_free_tree(&bootstrap_tree);
            }
          }
	  if (bs_thickness_from_lengths)
	  {
	    if (phyl_psinfo.bootstrap_numtrees != 1)
	      fprintf(stderr, "read %ld trees from %s -- only last one sets edge thickness values\n", phyl_psinfo.bootstrap_numtrees, bootstrapfile_name);
	    phyl_psinfo.bootstrap_numtrees = 0;
	  }
        }
        if (eps)
        {
          fprintf(outfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
          if (unrooted_style)
            fprintf(outfile, "%%%%BoundingBox: 0 0 432 432\n");
          else
            fprintf(outfile, "%%%%BoundingBox: 0 0 432 240\n");
          fprintf(outfile, "%%%%EndComments\n");
	  postscript_init(outfile);
          phyl_psinfo.fontname = "Courier-Bold";
          phyl_psinfo.spec_fontname = "Courier-Oblique";
          if (fontsize > 0.0)
          {
            phyl_psinfo.min_fontsize = fontsize;
            phyl_psinfo.max_fontsize = fontsize;
            phyl_psinfo.label_fontsize = fontsize;
          }
          else
          {
            phyl_psinfo.min_fontsize = 15.0;
            phyl_psinfo.max_fontsize = 60.0;
            phyl_psinfo.label_fontsize = 60.0;
          }
          if (angle_min > 0.0)
          {
            /* printf("angle_min set to %f\n", angle_min); */
            phyl_psinfo.angle_min = angle_min;
          }
          else
            phyl_psinfo.angle_min = 0.0;
          if (angle_limit > 0.0)
          {
            /* printf("angle_limit set to %f\n", angle_limit); */
            phyl_psinfo.angle_limit = angle_limit;
          }
          else
            phyl_psinfo.angle_limit = 0.0;
          phyl_psinfo.label = "units";
          phyl_psinfo.tic_length = 20.0;
          phyl_psinfo.label_start = 0.0;
          phyl_psinfo.linewidth = linewidth;
          phyl_psinfo.print_leafnames = print_leafnames;
          phyl_psinfo.dotted_nodes = dotted_nodes;
	  phyl_psinfo.tree_height = tree_height;
	  phyl_psinfo.print_edgelengths = print_edgelengths;
	  phyl_psinfo.leafnames_at_max = leafnames_at_max;
          phyl_psinfo.attrlist = attrlist;
	  if (bs_thickness_from_lengths && phyl_max_abs_thickness(&phyltree) > 0.0)
	    phyl_multiply_thick(&phyltree, 10.0 * linewidth / phyl_max_abs_thickness(&phyltree));
	  else
	    phyl_set_thickness(&phyltree, linewidth);
	  if (bootstrap_thick && bootstrapfile_name)
	  {
	    phyl_counter2thick(&phyltree, linewidth * 10.0 / phyl_psinfo.bootstrap_numtrees);
	    phyl_psinfo.bootstrap_numtrees = 0;
	  }
          /* {
            PHYL_LEAFATTRIBUTE *atl = attrlist;

            while (atl)
            {
              printf("n: %s, c:%s\n", atl->name, atl->class);
              atl = atl->next;
            }
          } */
          fprintf(outfile, "72 300 div 72 300 div scale\n");
          if (unrooted_style)
            phyl_ps_utree(outfile, &phyltree, 0.0, 0.0, 1800.0, 1800.0, &phyl_psinfo);
          else
            phyl_pstree(outfile, &phyltree, 0.0, 0.0, 1800.0, 1000.0, &phyl_psinfo);
          fprintf(outfile, "%%%%EOF\n");
          phyl_free_tree(&phyltree);
          fprintf(stderr, "*** EPS file can hold only one tree ***\n");
          break;
        }
        else
        {
          if (landscape)
          {
            page_width = 3200.0;
            page_height = 1800.0;
          }
          else
          {
            page_width = 1800.0;
            page_height = 3200.0;
          }
	  legend_y0 = 100.0 + legend_height(legendline, legend_fontheight);
          write_pageheader(outfile, buf, page_no++, page_width, page_height + 100, 200, landscape);
          phyl_psinfo.fontname = "Courier-Bold";
          phyl_psinfo.spec_fontname = "Courier-Oblique";
          if (fontsize > 0.0)
          {
            phyl_psinfo.min_fontsize = fontsize;
            phyl_psinfo.max_fontsize = fontsize;
            phyl_psinfo.label_fontsize = fontsize;
          }
          else
          {
            phyl_psinfo.min_fontsize = 15.0;
            phyl_psinfo.max_fontsize = 60.0;
            phyl_psinfo.label_fontsize = 60.0;
          }
          if (angle_min > 0.0)
          {
            /* printf("angle_min set to %f\n", angle_min); */
            phyl_psinfo.angle_min = angle_min;
          }
          else
            phyl_psinfo.angle_min = 0.0;
          if (angle_limit > 0.0)
          {
            /* printf("angle_limit set to %f\n", angle_limit); */
            phyl_psinfo.angle_limit = angle_limit;
          }
          else
            phyl_psinfo.angle_limit = 0.0;
          phyl_psinfo.label = "units";
          phyl_psinfo.tic_length = 20.0;
          phyl_psinfo.label_start = 0.0;
          phyl_psinfo.linewidth = linewidth;
          phyl_psinfo.print_leafnames = print_leafnames;
          phyl_psinfo.dotted_nodes = dotted_nodes;
	  phyl_psinfo.tree_height = tree_height;
	  phyl_psinfo.print_edgelengths = print_edgelengths;
	  phyl_psinfo.leafnames_at_max = leafnames_at_max;
          phyl_psinfo.attrlist = attrlist;
	  if (bs_thickness_from_lengths && phyl_max_abs_thickness(&phyltree) > 0.0)
	    phyl_multiply_thick(&phyltree, 10.0 * linewidth / phyl_max_abs_thickness(&phyltree));
	  else
	    phyl_set_thickness(&phyltree, linewidth);
	  if (bootstrap_thick && bootstrapfile_name)
	  {
	    phyl_counter2thick(&phyltree, linewidth * 10.0 / phyl_psinfo.bootstrap_numtrees);
	    phyl_psinfo.bootstrap_numtrees = 0;
	  }
          /* {
            PHYL_LEAFATTRIBUTE *atl = attrlist;

            while (atl)
            {
              printf("n: %s, c:%s\n", atl->name, atl->class);
              atl = atl->next;
            }
          } */
	  /* fprintf(outfile, "gsave 1 0 0 setrgbcolor %f %f moveto %f 0.0 rlineto 0.0 %f rlineto %f 0.0 rlineto closepath stroke grestore\n", 200.0, 100.0, page_width, page_height - legend_y0, -page_width); */
          if (unrooted_style)
            phyl_ps_utree(outfile, &phyltree, 200.0, legend_y0, page_width, page_height - legend_y0, &phyl_psinfo);
          else
            phyl_pstree(outfile, &phyltree, 200.0, legend_y0, page_width, page_height - legend_y0, &phyl_psinfo);
	  /* fprintf(outfile, "gsave 1 0 0 setrgbcolor 100 %f moveto 2000 0 rlineto stroke grestore\n", legend_y0); */
	  write_legend(outfile, legendline, legend_fontheight, 200.0, legend_y0);
          postscript_showpage(outfile);
          phyl_free_tree(&phyltree);
        }
      }
    }
  }
  if (attrlist)
    phyl_free_attrlist(attrlist);
  if (infile != stdin)
    fclose(infile);
  if (outfile != stdout)
    fclose (outfile);
  return (EXIT_SUCCESS);
}

