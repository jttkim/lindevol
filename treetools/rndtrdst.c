#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __atarist__
#  include <unistd.h>
#endif

#ifdef MEMDEBUG
#  include "memdebug.h"
#endif

#define MAX_NODENAME_LENGTH  31
#define MAX_SLEN           10000

#define TERR_MEM            -1
#define TERR_NOLEAVES       -2
#define TERR_DIFFNUMLEAVES  -3
#define TERR_INCOMPATLEAVES -4

#define EDGINF_NONE          0
#define EDGINF_DIFFEDGE      1
#define EDGINF_IDEDGE        2

#define free0(p) free(p); (p) = NULL


typedef struct tag_treenode
{
  struct tag_treenode  *ancestor;
  double                length;
  long                  edge_info;
  long                  num_descendants;
  struct tag_treenode **descendant;
  long                  leaf_index;
  char                  name[MAX_NODENAME_LENGTH + 1];
} TREENODE;


typedef struct
{
  TREENODE  *root;
  long       num_leaves;
  TREENODE **leaf;
  int lengthinfo_complete;
} TREE;


typedef struct
{
  FILE *f;
  int pos;
  char s[MAX_SLEN];
} TREE_SOURCE;


typedef struct
{
  TREENODE *node;
  size_t    size;
  char     *flag;
} LEAF_SET;


char buf[MAX_SLEN];


int read_treestring(TREE_SOURCE *ts)
{
  int i;

  if (fgets(ts->s, MAX_SLEN, ts->f) == NULL)
    return (-1);
  for (i = strlen(ts->s); i && iscntrl(ts->s[i]); i--)
    ts->s[i] = '\0';
  ts->pos = 0;
  /* printf("got string: \"%s\"\n", ts->s); */
  return (0);
}


void fprint_treestring(FILE *f, const TREE_SOURCE *ts)
{
  int i;

  fprintf(f, "%s\n", ts->s);
  for (i = 0; i < ts->pos; i++)
    fprintf(f, " ");
  fprintf(f, "^\n");
}


void free_treenode(TREENODE *node)
{
  long i;

  for (i = 0; i < node->num_descendants; i++)
    free_treenode(node->descendant[i]);
  if (node->num_descendants > 0)
    free(node->descendant);
  free(node);
}


void free_tree(TREE *tree)
{
  free_treenode(tree->root);
  tree->root = NULL;
  free(tree->leaf);
  tree->leaf = NULL;
  tree->num_leaves = 0;
}


int get_nodename(TREE_SOURCE *ts, TREENODE *node)
{
  int n;

  n = 0;
  while ((ts->s[ts->pos] != ')') && (ts->s[ts->pos] != ',') && (ts->s[ts->pos] != ':'))
  {
    if (!ts->s[ts->pos])
    {
      if (read_treestring(ts) < 0)
      {
        fprintf(stderr, "file error or unexpected end of file\n");
        return (-1);
      }
    }
    if ((isprint(ts->s[ts->pos])) && (n < MAX_NODENAME_LENGTH))
      node->name[n++] = ts->s[ts->pos++];
  }
  node->name[n] = '\0';
  /* printf("got node name \"%s\"\n", node->name); */
  return (0);
}


int get_branchlength(TREE_SOURCE *ts, TREENODE *node)
{
  int n;
  char *p;

  n = 0;
  while ((ts->s[ts->pos] != ')') && (ts->s[ts->pos] != ',') && (ts->s[ts->pos] != ':') && (ts->s[ts->pos] != ' '))
  {
    if (!ts->s[ts->pos])
    {
      if (read_treestring(ts) < 0)
      {
        fprintf(stderr, "file error or unexpected end of file\n");
        return (-1);
      }
    }
    if (n < MAX_SLEN - 1)
      buf[n++] = ts->s[ts->pos];
/*
    if (!isalnum(ts->s[ts->pos]))
      break;
*/
    ts->pos++;
  }
  buf[n] = '\0';
  /* printf("length specifier: %s\n", buf); */
  node->length = strtod(buf, &p);
  if ((p == buf) || (node->length < 0.0))
  {
    node->length = -1.0;
    return (-1);
  }
  /* printf("got branch length %f\n", node->length); */
  return (0);
}


TREENODE *add_node(TREENODE *root)
{
  TREENODE *node, **dsc;
  if ((node = (TREENODE *) malloc(sizeof(TREENODE))) == NULL)
    return (NULL);
  if (root)
  {
    if (root->num_descendants)
    {
      if ((dsc = (TREENODE **) realloc(root->descendant, (root->num_descendants + 1) * sizeof(TREENODE *))) == NULL)
      {
        free(node);
        return (NULL);
      }
      root->descendant = dsc;
    }
    else
    {
      if ((root->descendant = (TREENODE **) malloc(sizeof(TREENODE *))) == NULL)
      {
        free(node);
        return (NULL);
      }
    }
    root->descendant[root->num_descendants] = node;
    root->num_descendants++;
  }
  node->ancestor = root;
  node->length = -1.0;
  node->edge_info = EDGINF_NONE;
  node->num_descendants = 0;
  node->descendant = NULL;
  node->leaf_index = -1;
  node->name[0] = '\0';
  return (node);
}


int read_node(TREE_SOURCE *ts, TREENODE *root, TREE *tree)
{
  TREENODE *node, **leaf;
  int lengthinfo_present;

  while (ts->s[ts->pos] != '(')
  {
    if (!ts->s[ts->pos])
    {
      if (read_treestring(ts) < 0)
      {
        fprintf(stderr, "file error or unexpected end of file\n");
        return (-1);
      }
    }
    else
      ts->pos++;
  }
  ts->pos++;
  if ((node = add_node(root)) == NULL)
  {
    fprintf(stderr, "error: out of memory\n");
    return (-1);
  }
  lengthinfo_present = 0;
  while (ts->s[ts->pos] != ')')
  {
    while (isspace(ts->s[ts->pos]))
    {
      ts->pos++;
    }
    switch (ts->s[ts->pos])
    {
    case '\0':
      if (read_treestring(ts) < 0)
      {
        fprintf(stderr, "file error or unexpected end of file\n");
        return (-1);
      }
      break;
    case ':':
      ts->pos++;
      /* printf("getting length\n"); */
      if (lengthinfo_present)
      {
        fprintf(stderr, "error: multiple length specification\n");
        fprint_treestring(stderr, ts);
        return (-1);
      }
      if (get_branchlength(ts, node) < 0)
      {
        fprintf(stderr, "error: invalid length specification\n");
        fprint_treestring(stderr, ts);
        return (-1);
      }
      lengthinfo_present = 1;
      break;
    case ',':
      ts->pos++;
      /* printf("found node separator\n"); */
      if ((node = add_node(root)) == NULL)
      {
        fprintf(stderr, "error: out of memory\n");
        return (-1);
      }
      tree->lengthinfo_complete &= lengthinfo_present;
      lengthinfo_present = 0;
      break;
    case '(':
      /* printf("found_subtree\n"); */
      if (read_node(ts, node, tree) < 0)
        return (-1);
      break;
    default:
      /* printf("found leaf\n"); */
      if (get_nodename(ts, node) < 0)
      {
        fprintf(stderr, "error: invalid node name\n");
        fprint_treestring(stderr, ts);
        return (-1);
      }
      if (tree->num_leaves)
      {
        if ((leaf = (TREENODE **) realloc(tree->leaf, (tree->num_leaves + 1) * sizeof(TREENODE *))) == NULL)
        {
          fprintf(stderr, "error: out of memory\n");
          return (-1);
        }
        tree->leaf = leaf;
      }
      else
      {
        if ((tree->leaf = (TREENODE **) malloc(sizeof(TREENODE *))) == NULL)
        {
          fprintf(stderr, "error: out of memory\n");
          return (-1);
        }
      }
      tree->leaf[tree->num_leaves] = node;
      node->leaf_index = tree->num_leaves;
      /* printf("added leaf \"%s\"\n", tree->leaf[tree->num_leaves]->name); */
      tree->num_leaves++;
    }
  }
  tree->lengthinfo_complete &= lengthinfo_present;
  ts->pos++;
  return (0);
}


int read_tree(FILE *f, TREE *tree)
{
  TREE_SOURCE tree_source;

  tree_source.f = f;
  tree_source.s[0] = '\0';
  tree_source.pos = 0;
  tree->lengthinfo_complete = 1;
  if ((tree->root = add_node(NULL)) == NULL)
  {
    fprintf(stderr, "error: out of memory\n");
    return (-1);
  }
  if (read_node(&tree_source, tree->root, tree) < 0)
    return (-1);
  return (0);
}


int cmp_node_lexical(const void *n1, const void *n2)
{
  const char *name1 = (*((TREENODE **) n1))->name;
  const char *name2 = (*((TREENODE **) n2))->name;

  return (strncmp(name1, name2, MAX_NODENAME_LENGTH));
}


void sort_leaves(TREE *tree)
{
  long i;

  qsort(tree->leaf, tree->num_leaves, sizeof(TREENODE *), &cmp_node_lexical);
  for (i = 0; i < tree->num_leaves; i++)
    tree->leaf[i]->leaf_index = i;
}


void print_descendants(TREENODE *node)
{
  int i;

  if (node->num_descendants)
  {
    for (i = 0; i < node->num_descendants; i++)
      print_descendants(node->descendant[i]);
  }
  else
    printf("%s (#%ld)\n", node->name, node->leaf_index);

}


void print_subtree(TREENODE *node)
{
  int i;

  if (node->num_descendants)
  {
    printf("HTU is predecessor of\n");
    if (node->leaf_index >= 0)
      printf("***** HTU \"%s\" has non-negative leaf index *****\n", node->name);
    print_descendants(node);
    printf("\n");
    for (i = 0; i < node->num_descendants; i++)
      print_subtree(node->descendant[i]);
  }
  else if (node->leaf_index < 0)
    printf("***** OTU \"%s\" has negative leaf index *****\n", node->name);
}


void print_tree(TREE *tree)
{
  long i;

  printf("tree has %ld leaves:\n", tree->num_leaves);
  for (i = 0; i < tree->num_leaves; i++)
    printf("  %5ld: %s\n", i, tree->leaf[i]->name);
  print_subtree(tree->root);
}


void print_set(LEAF_SET *set, TREENODE **leaf)
{
  long i;
  char c;

  if (leaf == NULL)
  {
    for (i = 0; i < set->size; i++)
      printf("%d", set->flag[i]);
  }
  else
  {
    c = '{';
    for (i = 0; i < set->size; i++)
    {
      if (set->flag[i])
      {
        printf("%c %s", c, leaf[i]->name);
        c = ',';
      }
    }
    if (c == '{')
      printf("{");
    printf(" }");
  }
  printf("\n");
}


void get_leafset(const TREENODE *node, char *flag)
{
  long i;

  if (node->leaf_index >= 0)
    flag[node->leaf_index] = 1;
  else
  {
    for (i = 0; i < node->num_descendants; i++)
      get_leafset(node->descendant[i], flag);
  }
}


void get_leafsets(const TREENODE *node, long *num_edges, LEAF_SET *set)
{
  long i;

  get_leafset(node, set[*num_edges].flag);
  set[*num_edges].node = (TREENODE *) node;
  /* print_set(set + *num_edges, NULL); */
  (*num_edges)++;
  for (i = 0; i < node->num_descendants; i++)
  {
    if (node->descendant[i]->leaf_index < 0)
    {
      get_leafsets(node->descendant[i], num_edges, set);
    }
  }
}


long num_leaves(const TREENODE *node)
{
  long n = 0, i;

  if (node->num_descendants)
  {
    for (i = 0; i < node->num_descendants; i++)
      n += num_leaves(node->descendant[i]);
    return (n);
  }
  else
    return (1);
}


char *ps_string(const char *s, char *buf)
{
  size_t i, bpos = 0;
  const size_t l = strlen(s);

  if (buf == NULL)
  {
    if ((buf = (char *) malloc((2 * l + 3) * sizeof(char))) == NULL)
      return (NULL);
  }
  buf[bpos++] = '(';
  for (i = 0; i < l; i++)
  {
    switch (s[i])
    {
    case '(': case ')': case '\\':
      buf[bpos++] = '\\';
      buf[bpos++] = s[i];
      break;
    case '\b':
      buf[bpos++] = '\\';
      buf[bpos++] = 'b';
      break;
    case '\f':
      buf[bpos++] = '\\';
      buf[bpos++] = 'f';
      break;
    case '\n':
      buf[bpos++] = '\\';
      buf[bpos++] = 'n';
      break;
    case '\r':
      buf[bpos++] = '\\';
      buf[bpos++] = 'r';
      break;
    case '\t':
      buf[bpos++] = '\\';
      buf[bpos++] = 't';
      break;
    default:
      buf[bpos++] = s[i];
      break;
    }
  }
  buf[bpos++] = ')';
  buf[bpos] = '\0';
  return (buf);
}


void draw_node_ps(FILE *f, const TREENODE *node, double xpos, double ypos,
        double x0, double length_scale, double node_distance, double txt_yoffset,
	char *std_lineattr, char *idnode_lineattr)
{
  long i, nl_descendant;
  double xpos_descendant, ypos_descendant, x_left;

  if (node ->num_descendants)
    x_left = x0 + num_leaves(node->descendant[0]) * node_distance * 0.5;
  for (i = 0; i < node->num_descendants; i++)
  {
    nl_descendant = num_leaves(node->descendant[i]);
    xpos_descendant = x0 + nl_descendant * node_distance * 0.5;
    ypos_descendant = ypos + length_scale * node->descendant[i]->length;
    fprintf(f, "gsave 0 setlinecap\n");
    if (node->descendant[i]->edge_info == EDGINF_IDEDGE)
      fprintf(f, "%s\n", idnode_lineattr);
    else
      fprintf(f, "%s\n", std_lineattr);
    fprintf(f, "newpath %f %f moveto %f %f lineto\n", xpos_descendant, ypos, xpos_descendant, ypos_descendant);
    fprintf(f, "stroke grestore\n");
    draw_node_ps(f, node->descendant[i], xpos_descendant, ypos_descendant,
            x0, length_scale, node_distance, txt_yoffset, std_lineattr, idnode_lineattr);
    x0 += nl_descendant * node_distance;
  }
  if (node ->num_descendants)
  {
    fprintf(f, "gsave 0 setlinecap %s\n", std_lineattr);
    fprintf(f, "newpath %f %f moveto %f %f lineto stroke grestore\n", x_left, ypos, xpos_descendant, ypos);
  }
  if (node->name[0])
  {
    fprintf(f, "gsave leaffont setfont %f %f moveto 90 rotate %s show grestore\n",
            xpos, ypos + txt_yoffset, ps_string(node->name, buf));
  }
}


double tree_height(const TREENODE *node)
{
  double height = 0.0, h;
  long i;

  for (i = 0; i < node->num_descendants; i++)
  {
    h = tree_height(node->descendant[i]);
    if (h > height)
      height = h;
  }
  return (height + node->length);
}


void calculate_vshape_lengths(TREENODE *node, long ancestor_level)
{
  long i, level;

  level = num_leaves(node);
  node->length = ancestor_level - level;
  for (i = 0; i < node->num_descendants; i++)
    calculate_vshape_lengths(node->descendant[i], level);
}


void postscript_tree(FILE *f, const TREE *tree, double x0, double y0, double width, double height)
{
  double node_distance, length_scale, xpos, ypos;
  double tree_h;
  const double charcell_height = 50.0;
  char std_lineattr[256], idnode_lineattr[256];
  int idnode_linewidth;

  idnode_linewidth = width / tree->num_leaves / 10;
  idnode_linewidth = (idnode_linewidth < 5) ? 5 : idnode_linewidth;
  sprintf(std_lineattr, "1 setlinewidth");
  sprintf(idnode_lineattr, "%d setlinewidth", idnode_linewidth);
  if (!tree->lengthinfo_complete)
    calculate_vshape_lengths(tree->root, tree->num_leaves);
  node_distance = width / tree->num_leaves;
  tree_h = tree_height(tree->root);
  length_scale = (height - charcell_height - 100.0) / tree_h;
  xpos = x0 + width * 0.5;
  ypos = y0 + 100.0;
  fprintf(f, "gsave 0 setlinecap\n");
  fprintf(f, "%f /leaffont exch /Courier FSD\n", width * 0.9 / tree->num_leaves);
  if (tree->root->edge_info == EDGINF_IDEDGE)
    fprintf(f, "%s\n", idnode_lineattr);
  else
    fprintf(f, "%s\n", std_lineattr);
  fprintf(f, "%f %f moveto %f %f lineto stroke grestore\n", xpos, y0, xpos, ypos);
  draw_node_ps(f, tree->root, xpos, ypos, x0, length_scale, node_distance, 20.0, std_lineattr, idnode_lineattr);
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


int cmp_leafset(const void *s1, const void *s2)
{
  LEAF_SET *set1 = (LEAF_SET *) s1, *set2 = (LEAF_SET *) s2;
  return (memcmp(set1->flag, set2->flag, set1->size));
}


int random_binary_tree(long num_leaves, TREE *tree, int seed)
{
  TREENODE **node, *new_leaf;
  long num_nodes, n, k;

  srandom(seed);
  if ((node = (TREENODE **) malloc((num_leaves - 1) * sizeof(TREENODE *))) == NULL)
    return (TERR_MEM);
  if ((tree->leaf = (TREENODE **) malloc(num_leaves * sizeof(TREENODE *))) == NULL)
  {
    free(node);
    return (TERR_MEM);
  }
  if ((tree->root = add_node(NULL)) == NULL)
  {
    free(node);
    free0(tree->leaf);
    return (TERR_MEM);
  }
  node[0] = tree->root;
  tree->num_leaves = 0;
  tree->lengthinfo_complete = 0;
  if ((new_leaf = add_node(tree->root)) == NULL)
  {
    free(node);
    free_tree(tree);
    return (TERR_MEM);
  }
  sprintf(new_leaf->name, "rndleaf #%ld", tree->num_leaves);
  new_leaf->leaf_index = tree->num_leaves;
  tree->leaf[tree->num_leaves++] = new_leaf;
  if ((new_leaf = add_node(tree->root)) == NULL)
  {
    free(node);
    free_tree(tree);
    return (TERR_MEM);
  }
  sprintf(new_leaf->name, "rndleaf #%ld", tree->num_leaves);
  new_leaf->leaf_index = tree->num_leaves;
  tree->leaf[tree->num_leaves++] = new_leaf;
  for (num_nodes = 1; num_nodes < num_leaves - 1; num_nodes++)
  {
    if ((node[num_nodes] = add_node(NULL)) == NULL)
    {
      free(node);
      free_tree(tree);
      return (TERR_MEM);
    }
    if ((node[num_nodes]->descendant = (TREENODE **) malloc(sizeof(TREENODE *))) == NULL)
    {
      free(node[num_nodes]);
      free(node);
      free_tree(tree);
      return (TERR_MEM);
    }
    n = random() % num_nodes;
    k = random() % 2;
    node[num_nodes]->ancestor = node[n];
    node[num_nodes]->descendant[0] = node[n]->descendant[k];
    node[n]->descendant[k] = node[num_nodes];
    node[num_nodes]->num_descendants = 1;
    if ((new_leaf = add_node(node[num_nodes])) == NULL)
    {
      free(node);
      free_tree(tree);
      return (TERR_MEM);
    }
    sprintf(new_leaf->name, "rndleaf #%ld", tree->num_leaves);
    new_leaf->leaf_index = tree->num_leaves;
    tree->leaf[tree->num_leaves++] = new_leaf;
  }
  free(node);
  return (0);
}


long topological_treedistance(TREE *tree1, TREE *tree2)
{
  LEAF_SET *set1, *set2;
  char     *set1_flags, *set2_flags;
  long i, num_edges1, num_edges2, edge1, edge2, d, similarity;
  int ec;

  if ((tree1->num_leaves == 0) || (tree2->num_leaves == 0))
    return (TERR_NOLEAVES);
  if (tree1->num_leaves != tree2->num_leaves)
    return (TERR_DIFFNUMLEAVES);
  sort_leaves(tree1);
  sort_leaves(tree2);
  for (i = 0; i < tree1->num_leaves; i++)
  {
    if (strcmp(tree1->leaf[i]->name, tree2->leaf[i]->name))
      return (TERR_INCOMPATLEAVES);
  }
  if ((set1 = (LEAF_SET *) malloc(tree1->num_leaves * sizeof(LEAF_SET))) == NULL)
    return (TERR_MEM);
  if ((set1_flags = (char *) malloc(tree1->num_leaves * tree1->num_leaves * sizeof(char))) == NULL)
  {
    free(set1);
    return (TERR_MEM);
  }
  memset(set1_flags, 0, tree1->num_leaves * tree1->num_leaves);
  for (i = 0; i < tree1->num_leaves; i++)
  {
    set1[i].size = tree1->num_leaves;
    set1[i].flag = set1_flags + tree1->num_leaves * i;
  }
  if ((set2 = (LEAF_SET *) malloc(tree1->num_leaves * sizeof(LEAF_SET))) == NULL)
  {
    free(set1);
    free(set1_flags);
    return (TERR_MEM);
  }
  if ((set2_flags = (char *) malloc(tree2->num_leaves * tree2->num_leaves * sizeof(char))) == NULL)
  {
    free(set1);
    free(set1_flags);
    free(set2);
    return (TERR_MEM);
  }
  memset(set2_flags, 0, tree2->num_leaves * tree2->num_leaves);
  for (i = 0; i < tree2->num_leaves; i++)
  {
    set2[i].size = tree2->num_leaves;
    set2[i].flag = set2_flags + tree2->num_leaves * i;
  }
  num_edges1 = 0;
  get_leafsets(tree1->root, &num_edges1, set1);
  qsort(set1, num_edges1, sizeof(LEAF_SET), cmp_leafset);
/*
  printf("tree #1 has %ld internal edges\n", num_edges1);
  for (i = 0; i < num_edges1; i++)
    print_set(set1 + i, tree1->leaf);
*/
  num_edges2 = 0;
  get_leafsets(tree2->root, &num_edges2, set2);
  qsort(set2, num_edges2, sizeof(LEAF_SET), cmp_leafset);
/* 
  printf("tree #2 has %ld internal edges\n", num_edges2);
  for (i = 0; i < num_edges2; i++)
    print_set(set2 + i, tree2->leaf);
*/
  d = 0;
  similarity = 0;
  edge1 = 0;
  edge2 = 0;
  while ((edge1 < num_edges1) && (edge2 < num_edges2))
  {
    ec = memcmp(set1[edge1].flag, set2[edge2].flag, tree1->num_leaves);
    /* printf("edge1 = %ld, edge2 = %ld, memcmp result: %d\n", edge1, edge2, ec); */
    if (ec)
    {
      d++;
/*
      printf("differing leaf sets:\n");
      print_set(set1 + edge1, tree1->leaf);
      print_set(set2 + edge2, tree2->leaf);
*/
      set1[edge1].node->edge_info = EDGINF_DIFFEDGE;
      set2[edge2].node->edge_info = EDGINF_DIFFEDGE;
    }
    else
    {
      similarity++;
/*
      printf("identical leaf set:\n");
      print_set(set1 + edge1, tree1->leaf);
*/
      set1[edge1].node->edge_info = EDGINF_IDEDGE;
      set2[edge2].node->edge_info = EDGINF_IDEDGE;
    }
    if (ec <= 0)
      edge1++;
    if (ec >= 0)
      edge2++;
  }
  /* printf("\nsimilarity (# identical edges): %ld\n", similarity); */
  free(set1);
  free(set1_flags);
  free(set2);
  free(set2_flags);
  return (d);
}


int main(int argc, char **argv)
{
  FILE *f;
  TREE tree1, tree2;
  long num_leaves, max_num_leaves, num_tests, i, *result;
  double av_treedist, std_deviation;

  tree1.root = NULL;
  tree1.num_leaves = 0;
  tree1.leaf = NULL;
  tree2.root = NULL;
  tree2.num_leaves = 0;
  tree2.leaf = NULL;
  if (argc > 2)
  {
    max_num_leaves = strtol(argv[1], NULL, 10);
    num_tests = strtol(argv[2], NULL, 10);
    if ((result = (long *) malloc(num_tests * sizeof(long))) == NULL)
    {
      fprintf(stderr, "not enough memory for result array -- exit\n");
      exit (EXIT_FAILURE);
    }
    if ((f = fopen("treedist.dat", "w")) == NULL)
    {
      free(result);
      fprintf(stderr, "failed to open output file -- exit\n");
      exit (EXIT_FAILURE);
    }
    for (num_leaves = 2; num_leaves <= max_num_leaves; num_leaves++)
    {
      av_treedist = 0.0;
      std_deviation = 0.0;
      for (i = 0; i < num_tests; i++)
      {
        if (random_binary_tree(num_leaves, &tree1, i + 12345) < 0)
        {
          fprintf(stderr, "error creating random tree with %ld leaves -- exit\n", num_leaves);
          free(result);
          exit (EXIT_FAILURE);
        }
        if (random_binary_tree(num_leaves, &tree2, random()) < 0)
        {
          free_tree(&tree1);
          free(result);
          fprintf(stderr, "error creating random tree with %ld leaves -- exit\n", num_leaves);
          exit (EXIT_FAILURE);
        }
        result[i] = topological_treedistance(&tree1, &tree2);
        av_treedist += result[i];
        free_tree(&tree1);
        free_tree(&tree2);
      }
      av_treedist /= num_tests;
      for (i = 0; i < num_tests; i++)
        std_deviation += (result[i] - av_treedist) * (result[i] - av_treedist);
      std_deviation = sqrt(std_deviation / num_tests);
      fprintf(f, "%ld %f %f\n", num_leaves, av_treedist, std_deviation);
    }
    free(result);
    fclose(f);
  }
  return (EXIT_SUCCESS);
}

