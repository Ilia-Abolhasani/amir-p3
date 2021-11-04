#define PROGRAM_NAME "overlay_boxplot_ng"
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "util.h"
#include "overlay_boxplot.h"
#include "overlay_boxplot_setcolor.h"
#include "overlay_boxplot_read_ct.h"
#include "overlay_boxplot_ps.h"
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
# include "overlay_boxplot_img.h"
#endif

#define MESSAGE_NONE   0	/* values for g_type_message */
#define MESSAGE_ENERGY 2

struct helix    g_diag[MAXIMUM_HELICES];

int             g_len;	/* Maximum of farthest base pair in ct, plot
				 * files */
int             g_plot_len;/* maximum base pair in plot file */
int             g_ct_len;	/* maximum base pair in ct file */
int             g_width;	/* these are the current width and height of
				 * the main window */
int             g_height;	/* in pixels */
int             g_width_z, g_height_z;	/* same for zoom window */


/* when right mouse button is clicked, transform coordinates into i,j */
int             g_i_row = 0;
int             g_j_col = 0;	/* for i.j base pair with mouse */
char            g_i_row_s[10], g_j_col_s[10];	/* string of above */

/* ______________________________________________________________ */
/* Graphic  variables  ------------------------------------------- */
int             g_opt_mode;
int             g_opt_distance;

float           g_color_increment;	/* range of energy for each color */
int             g_energy_cutoff;	/* only probabilities <= */
/* g_energy_cutoff will be shown. */
/* note that probabilities have a negative sign */
/* -30 indicates probabilities >=.3 will be shown */
int             g_chain_len;	/* chains smaller will not be
					 * displayed */
int             points_plotted, points_plotted_zoom; /* indicates
							number of
							points plotted
							in each window
							*/ 
/* _________________________________________________________ */
/* Zoom Variables __________________________________________ */


int             g_display_l, g_display_r, g_display_t, g_display_b;

/* left,right,top,bottom,width and height for zoom window display
 * They are in necleotide coordinates */
int   saved_display_l, saved_display_r, saved_display_t, saved_display_b;
/* Above used to go back to previous zoom window */

int             g_zoom_labels_count;	/* the number of labels used */
int             g_zoom_labels_x[50];	/* x coordinate for label i */
int             g_zoom_labels_xl[50];	/* used to plot on screen the
						 * dot */
int             g_zoom_labels_y[50];
int             g_zoom_labels_row[50];
int             g_zoom_labels_col[50];

char           *g_zoom_labels_string[50];	/* value to place at x,y */

/* File Variables _______________________________________________ */
char            rec[90];
FILE           *fp;		/* for .plot file */

char            g_post_filename[80];
char            g_output_filename[80];	/* as provided by user */
char            filename1[80];	/* these two are the plot and ct file */
char            filename2[80];	/* Once it is clear which is which, store
				 * below */
char            g_plot_filename[80];
int             g_ct_red, g_ct_green, g_ct_yellow;
int             g_zoom_ct_red, g_zoom_ct_green, g_zoom_ct_yellow;
int             g_red_gray_overlap;
char            g_ct_filename[80];
char            g_img_filename[80]; /* this depends on -o or input file */
int             g_prob;
int             g_grid_flag;	/* True for use grid lines False for do not */
int             g_zoom_window_exists = FALSE;
int             g_optimal_energy;
int             g_worst_energy;
int             g_post_file; /* 1: main, 2: zoom; used from text display */
int            *g_diag_count;
int            *g_diag_start;
int             g_ct_diag_count[2 * MAXIMUM_SIZE];
int             g_ct_diag_start[2 * MAXIMUM_SIZE];
int             g_ct_files = 1;	/* this is always 1 for this program */
/* it is a true variable in ct_boxplot */
char       filename_ct[MAXIMUM_CT_FILES][80];
/* filenames for input */
char      file_data_ct[MAXIMUM_CT_FILES][80];	/* first line of ct file */

struct ct_helix g_ct_diag[MAXIMUM_HELICES];
char            g_sequence_name[120];
int             g_sequence_name_set;
int             g_helices_in_plot_file;	/* number of entries in file */
int             g_helices_in_ct_file;
int             g_resolution;	/* 72 110 200 or 300 */

/* return 0 if no helix is defined */
int chain_length_read(int row, int col) {
  int             diag, diag_start, diag_end, i, temp_row;
  /* perhaps switch row,col when row>col */
  if (row > col) {
    temp_row = col;
    col = row;
    row = temp_row;
  }
  diag = col + row - 1;
  diag_start = g_diag_start[diag];
  diag_end = g_diag_start[diag + 1];
  if ((row > g_plot_len) || (col > g_plot_len))
    return 0;
  for (i = diag_start; i < diag_end; i++) {
    temp_row = g_diag[i].row;
    if (temp_row > row)
      return INT_MIN;	/* entries in diagonal are sorted */
    if ((temp_row <= row) && (g_diag[i].row + g_diag[i].length - 1 >= row))
      return g_diag[i].length;
  }
  return 0;
}


/* given a row and column, return the energy of that position */
int energy_read(int row, int col) {
  int             diag, diag_start, diag_end, i, temp_row;
  /* perhaps switch row,col when row>col */
  if (row > col) {
    temp_row = col;
    col = row;
    row = temp_row;
  }
  diag = col + row - 1;
  diag_start = g_diag_start[diag];
  diag_end = g_diag_start[diag + 1];
  if ((row > g_plot_len) || (col > g_plot_len))
    return INT_MIN;
  for (i = diag_start; i < diag_end; i++) {
    temp_row = g_diag[i].row;
    if (temp_row > row)
      return INT_MIN;	/* entries in diagonal are sorted */
    if ((temp_row <= row) && (g_diag[i].row + g_diag[i].length - 1 >= row))
      return g_diag[i].energy;
  }
  return INT_MIN;
}
/* given a row and column, return TRUE if  a ct dot exists  */
int ct_read(int row, int col) {
  int             diag, diag_start, diag_end, i, temp_row;
  diag = col + row - 1;
  diag_start = g_ct_diag_start[diag];
  diag_end = g_ct_diag_start[diag + 1];
  if ((row > g_ct_len) || (col > g_ct_len))
    return FALSE;
  for (i = diag_start; i < diag_end; i++) {
    temp_row = g_ct_diag[i].row;
    if (temp_row > row)
      return FALSE;	/* entries in diagonal are sorted */
    if ((temp_row <= row)&&(g_ct_diag[i].row+g_ct_diag[i].length - 1 >= row)) {
      return TRUE;
    }
  }
  return FALSE;
}

/* Initialization Routines for nongraphics */
void check_output_file(int user_output_flag, char *output_filename) {
  if (user_output_flag) {
    strcpy(g_post_filename, output_filename);
    strcpy(g_img_filename, output_filename);
  } else {
    strcpy(g_post_filename, g_plot_filename);
    strcpy(g_img_filename, g_plot_filename);
    chop_suffix(g_post_filename, ".plot");
    chop_suffix(g_img_filename, ".plot");
    chop_suffix(g_img_filename, ".img");
  }
  if (g_opt_mode) {
    strcat(g_post_filename, ".ovr-opt.ps");
    strcat(g_img_filename, ".ovr-opt.gif");
  } else {
    strcat(g_post_filename, ".ovr-ct.ps");
    strcat(g_img_filename, ".ovr-ct.gif");
  }
}

void insufficient_arguments(void) {
  printf("\n Insufficient arguments:\n");
  printf("\n[ Usage:      %s name1.plot name2.ct                         ]", 
	 PROGRAM_NAME);
  printf("\n[         or  %s -o name_out name1.plot name2.ct             ]", 
	 PROGRAM_NAME);
  printf("\n[        (.plot or .ct suffixes optional)                    ]");
  printf("\n[                                                            ]");
  printf("\n[          Valid arguments                                   ]");
  printf("\n[ -a       Enable Optimal mode                               ]");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("\n[ -b       Enable clear background for gif                   ]");
#endif
  printf("\n[ -c energy Within energy of optimal, use yellow for ct dots ]");
  printf("\n[ -c distance  With -a , distance of 0,1,2,3,4,5,6,7,8, or 9 ]");
  printf("\n[ -f filter    Display only helices of length >= filter      ]");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("\n[ -g           Create only a gif image file for output.      ]");
#endif
  printf("\n[ -go          Grid lines off. Default is on                 ]");
  printf("\n[ -i increment Specify energy increment for plot file        ]");
  printf("\n[ -o name_out  Specify name of output file as                ]");
  printf("\n[              name_out.ovr.ps or name_out.ovr.gif           ]");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("\n[ -pg          Create both postscript and gif output         ]");
  printf("\n[ -r res       Specify resolution of gif (50 to 400)         ]");
#endif
  printf("\n[ -t \"TITLE\"   Specify TITLE as title, (Quotes Required)   ]");
  printf("\n[ -z l r t b  Specify a zoom region: left, right, top, bot   ]");
  printf("\n\n[ Default: Color ct file dots based on how they overlap    ]");
  printf("\n[ the optimal dots.                                          ]");
  printf("\n[         Green for overlap,                                 ]");
  printf("\n[         Yellow within the energy specified with -c,        ]");
  printf("\n[         and red for all others.                            ]");
  printf("\n[ With the -a switch, Optimal dots are colored.              ]");
  printf("\n[          Green for overlap with ct dots,                   ]");
  printf("\n[          Yellow for within the distance specified with -c, ]");
  printf("\n[          and red otherwise.                                ]");
  printf("\n[          Unused ct dots are cyan                           ]");
  printf("\n[                                                            ]");
  printf("\n[See overlay_boxplot_ng.doc for more details                 ]");
  printf("\n");
  try_exit(7);
}

void open_plot_file(char *filename1, char *filename2) {
  char            filename[80];
  strcpy(filename, filename1);
  chop_suffix(filename, ".plot"); 
  strcat(filename, ".plot");	/* put .plot on */
  printf("Trying to open %s as a plot file.\n", filename);
  if ((fp = fopen(filename, "r")) == NULL) {
    strcpy(filename, filename2);
    chop_suffix(filename, ".plot");
    strcat(filename, ".plot");
    printf("Trying to open %s as a plot file.\n", filename);
    if ((fp = fopen(filename, "r")) == NULL) {
      printf("Error!\tCould not open the plot file %s.\n", filename);
      insufficient_arguments();
    } else {
      strcpy(g_plot_filename, filename);
      strcpy(g_ct_filename, filename1);
    }
  } else {
    strcpy(g_plot_filename, filename);
    strcpy(g_ct_filename, filename2);
  }
  chop_suffix(g_ct_filename, ".ct");
  strcat(g_ct_filename, ".ct");
}

int compare_helices(const void *h1, const void *h2) {
  int             dif;
  struct helix   *helix1, *helix2;
  helix1 = (struct helix *) h1;
  helix2 = (struct helix *) h2;
  dif = (*helix1).diagonal - (*helix2).diagonal; /* return <0 if diag1 <
						  * diag2 on same diagonal */
  if (dif == 0)		
    dif = (*helix1).row - (*helix2).row; 
  return dif;
}

void sort_helices(void) {
  qsort(g_diag, g_helices_in_plot_file, sizeof(struct helix), compare_helices);
}

void print_helix(int number_of_helices) {
  int             i;
  for (i = 0; i < number_of_helices; i++) {
    printf("diag=%d, row=%d, color=%d, column=%d, energy=%d, length=%d\n",
	   g_diag[i].diagonal,  g_diag[i].row,    g_diag[i].color, 
	   g_diag[i].column,    g_diag[i].energy, g_diag[i].length);
  }
}

void fix_counts(void) {
  int             i;
  for (i = 0; i < g_helices_in_plot_file; i++) {
    g_diag_count[g_diag[i].diagonal]++;	/* update count for this diagonal */
  }
}

void fix_colors(void) {	/* set the colors for each line of input */
  int             i;
  for (i = 0; i < g_helices_in_plot_file; i++) {
    g_diag[i].color = find_color(g_diag[i].energy);
  }
}

void initialize_data2(void) {
  g_number_of_colors = 4;	/* start with 4 colors $$$$$ */
  initialize_colors();	/* reads .col file */
  g_energy_cutoff = g_worst_energy;
  if (g_prob != TRUE)
    g_color_increment = ((float) (g_energy_cutoff - g_optimal_energy)) / 
      (float) (g_number_of_colors - 1);
}

void initialize_len(void) { /* g_worst_energy and  g_optimal_energy
			       are also set */ 
  int             level, length, row, col, energy, diag;
  int             i;
  /* diag is the diagonal that each row,col corresponds to */
  float           rec_energy;
  g_helices_in_plot_file = 0;
  if (fgets(rec, 90, fp) == NULL) {
    printf("Error!\tThe plot file is completely empty.\n");
    try_exit(8);
  }
  g_zoom_labels_count = 0;	/* initialize to no labels */
  for (i = 0; i < 50; i++) {	/* allocate memory for string */
    g_zoom_labels_string[i] = (char *) malloc(15 * sizeof(char));
    if (g_zoom_labels_string[i] == NULL) {
      printf("Error!\tInsufficient memory for zoom labels.\n");
      try_exit(9);
    }
  }
  g_plot_len = 0;
  /* Read the plot file */
  g_prob = FALSE;
  g_worst_energy = INT_MIN;
  g_optimal_energy = INT_MAX;
  while (fgets(rec, 90, fp) != NULL) {
    sscanf(rec, "%d%d%d%d%f", &level, &length, &row, &col, &rec_energy);
    if (g_plot_len < col)
      g_plot_len = col;
    g_diag[g_helices_in_plot_file].length = length;
    g_diag[g_helices_in_plot_file].row = row;
    g_diag[g_helices_in_plot_file].column = col;
    if (rec_energy > 0.0)
      energy = (int) (-1 * (1000. * rec_energy + .50));
    else
      energy = (int) (rec_energy);
    g_diag[g_helices_in_plot_file].energy = energy;
    if (rec_energy > 0)
      g_prob = TRUE;
    if (energy < g_optimal_energy)
      g_optimal_energy = energy;
    if (energy > g_worst_energy)
      g_worst_energy = energy;
    diag = row + col - 1;
    g_diag[g_helices_in_plot_file].diagonal = diag;
    g_helices_in_plot_file++;
    if (g_helices_in_plot_file == MAXIMUM_HELICES) {
      printf("Error!\tThere are more than %d helices in the plot file.\n", 
	     MAXIMUM_HELICES);
      printf("Increase the constant MAXIMUM_HELICES in %s.c", PROGRAM_NAME);
      printf(" and recompile.\n");
      try_exit(10);
    }
  }
  printf("There are %d helices in the plot file.\n", g_helices_in_plot_file);
  if (g_helices_in_plot_file == 0) {
    printf("Error!\tThe plot file contains no helices\n");
    try_exit(11);
  }
  fclose(fp);
  /* Sort the plot file */
  sort_helices();
  g_diag_start = (int *) malloc((2 * g_plot_len) * sizeof(int));
  if (g_diag_start == NULL)
    printf("Error!\tInsufficient memory for g_diag_start.\n");
  g_diag_count = (int *) malloc((2 * g_plot_len) * sizeof(int));
  if (g_diag_count == NULL)
    printf("Error!\tInsufficient memory for g_diag_count.\n");
  for (i = 0; i <= (2 * g_plot_len - 1); i++) 
    g_diag_count[i] = 0;
  /* fix the count of diagonal for each entry in plot file */
  initialize_data2();
  fix_colors();
  g_diag_start[1] = 0;
  fix_counts();
  for (i = 2; i <= (2 * g_plot_len - 1); i++) {
    g_diag_start[i] = g_diag_start[i - 1] + g_diag_count[i - 1];
  }
  /* g_diag_start counts forward into each array for start of each diagonal */
}

void open_ct_file(char *ct_filename) {
  g_ct_len = 0;
  printf("The farthest base pair in the plot file is %d.\n", g_plot_len);
  strcpy(filename_ct[0], ct_filename);
  g_helices_in_ct_file = 0;
  read_of_all_ct_files(fp, g_ct_diag_count, g_ct_diag_start, g_ct_files, 
		       &g_helices_in_ct_file, &g_ct_len, filename_ct, 
		       g_ct_diag, file_data_ct, g_sequence_name,
		       g_sequence_name_set);
  if (g_ct_len > g_plot_len) {
    printf("The ct file has a base pair beyond any helix in the plot file\n");
    g_len = g_ct_len;
  } else {
    g_len = g_plot_len;
  }
}

void finish_setting_energy_cutoff(float cutoff) {
  g_energy_cutoff = (int) (cutoff * 10);
  if (g_energy_cutoff >= 0)
    g_energy_cutoff = g_optimal_energy + g_energy_cutoff;
  if (g_energy_cutoff > g_worst_energy)
    g_energy_cutoff = g_worst_energy;
  if (g_energy_cutoff < g_optimal_energy)
    g_energy_cutoff = g_optimal_energy;
  g_color_increment = (float) (g_energy_cutoff - g_optimal_energy) / 
    (float) (g_number_of_colors - 1);
  fix_colors();
}

void set_energy_cutoff1(float cutoff_level) {	/* for probabilities */
  g_energy_cutoff = (int) (cutoff_level * -1000);
}

void initialize_data1(int argc, char **argv) {
  int             i, files_used;
  int             user_output_flag;
  int             zoom_flag;
  int             img_resolution;	/* 72 110 200 or 300 */
  int             img_flag, postscript_flag;
  int             temp;
  int             width, height;
  int             clear_background_flag;
  int             close_optimal_defined;
  float           close_optimal_energy;
  float           new_cutoff_float;
  int             new_cutoff_flag;
  g_opt_mode = FALSE;
  g_opt_distance = 0;
  close_optimal_defined = FALSE;
  g_prob = FALSE;
  g_grid_flag = TRUE;
  img_resolution = 72;
  img_flag = FALSE;
  postscript_flag = TRUE;	/* default to create only postscript output */
  g_chain_len = 1;
  clear_background_flag = FALSE;
  zoom_flag = FALSE;
  user_output_flag = FALSE;
  g_sequence_name_set = FALSE;
  files_used = 0;
  new_cutoff_flag = FALSE;
  i = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-r") == 0) { /* set image resolution */
      img_resolution = atoi(argv[i + 1]);
      /* Resolution must follow -r */
      if (img_resolution < 20)
	img_resolution = 20;
      else if (img_resolution > 400)
	img_resolution = 400;
      i += 1;
    } else if (strcmp(argv[i], "-c") == 0) {	/* set increment */
      sscanf(argv[i + 1], "%f", &close_optimal_energy);
      /* copy float into cutoff */
      close_optimal_defined = TRUE;
      i++;
    } else if (strcmp(argv[i], "-i") == 0) { /* set increment */
      printf("-i was found\n");
      new_cutoff_flag = TRUE;
      sscanf(argv[i + 1], "%f", &new_cutoff_float);
      /* copy float into cutoff */
      if (g_prob == TRUE) {
	if (new_cutoff_float > 1.)
	  new_cutoff_float = 1.0;
	else if (new_cutoff_float < 0.0)
	  new_cutoff_float = 0.0;
	set_energy_cutoff1(new_cutoff_float);
      }
      i++;
    } else if (strcmp(argv[i], "-z")== 0) {
      if (i < (argc - 4)) {
	g_display_l = atoi(argv[i + 1]);
	g_display_r = atoi(argv[i + 2]);
	g_display_t = atoi(argv[i + 3]);
	g_display_b = atoi(argv[i + 4]);
	zoom_flag = TRUE;
	i = i + 4;
      } else {
	printf("-z must have 4 numbers following it.\n");
      }
    } else if (strcmp(argv[i], "-g") == 0) {
      img_flag = TRUE;
      postscript_flag = FALSE;
    } else if (strcmp(argv[i], "-a") == 0) {
      g_opt_mode = TRUE;
      printf("Using optimal mode\n");
    } else if (strcmp(argv[i], "-b") == 0) {
      clear_background_flag = TRUE;
    } else if (strcmp(argv[i], "-pg") == 0) {
      img_flag = TRUE;
      postscript_flag = TRUE;
    } else if (strcmp(argv[i], "-go") == 0) {
      g_grid_flag = FALSE;
    } else if (strcmp(argv[i], "-f") == 0) {
      if (i < (argc - 1)) {
	g_chain_len = atoi(argv[i + 1]);
	if (g_chain_len < 1)
	  g_chain_len = 1;
	else if (g_chain_len > 20)
	  g_chain_len = 20;
	i++;
      } else {
	printf("Error!\tFilter length must follow -f");
      }
    } else if (strcmp(argv[i], "-o") == 0) {
      if (i < (argc - 1)) {
	user_output_flag = TRUE;
	strcpy(g_output_filename, argv[i + 1]);
	i++;
      } else {
	printf("Error!\t Output file name must follow -o");
      }
    } else if (strcmp(argv[i], "-t") == 0) {
      if (i < (argc - 1)) {
	g_sequence_name_set = TRUE;
	strcpy(g_sequence_name, argv[i + 1]);
	i++;
      }
    } else if (files_used == 0) {
      strcpy(filename1, argv[i]);
      files_used++;
    } else if (files_used == 1) {
      strcpy(filename2, argv[i]);
      files_used++;
    } else {
      printf("Error!\tExtra parameter %s is ignored.\n",  argv[i]);
      printf("%s and %s were assumed to be input files.\n", 
	     filename1, filename2);
    }
    i++;
  }
  if (files_used < 2) {
    printf("Error!\tYou must supply both a plot file and a a ct file.\n");
    insufficient_arguments();
  }
  check_output_file(user_output_flag, g_output_filename);
  open_plot_file(filename1, filename2);
  check_output_file(user_output_flag, g_output_filename);
  initialize_len();	/* reads plot file */
  open_ct_file(g_ct_filename);	/* reads ct file */
  if (!zoom_flag) {
    g_display_l = 1;
    g_display_t = 1;
    g_display_b = g_len;
    g_display_r = g_len;
  } else {
    if ((g_display_l > g_len) || (g_display_l < 1))
      g_display_l = 1;
    if ((g_display_r > g_len) || (g_display_r < 1))
      g_display_r = g_len;
    if ((g_display_t > g_len) || (g_display_t < 1))
      g_display_t = 1;
    if ((g_display_b > g_len) || (g_display_b < 1))
      g_display_b = g_len;
    /* keep left smaller than right */
    if (g_display_r < g_display_l) {
      temp = g_display_r;
      g_display_r = g_display_l;
      g_display_l = temp;
    }
    /* keep top smaller than bottom */
    if (g_display_b < g_display_t) {
      temp = g_display_b;
      g_display_b = g_display_t;
      g_display_t = temp;
    }
  }
  if (new_cutoff_flag == TRUE) {
    printf("Float of cutoff is %f.\n", new_cutoff_float);
    finish_setting_energy_cutoff(new_cutoff_float);
  }
  if (close_optimal_defined) {
    if (g_opt_mode) {
      g_opt_distance = (int) (close_optimal_energy + .5);
    } else {
      close_optimal_energy = close_optimal_energy * 10;
      if (g_worst_energy == g_optimal_energy) {
	close_optimal_energy = g_optimal_energy;
      } else if (close_optimal_energy >= 0)
	close_optimal_energy = (float) g_optimal_energy + close_optimal_energy;
      g_close_optimal = ((float) g_optimal_energy - close_optimal_energy) /
	((float) g_optimal_energy - (float) g_worst_energy);
      if (g_close_optimal < 0.0)
	g_close_optimal = 0.0;
      if (g_close_optimal > 1.0)
	g_close_optimal = 1.0;
    }
  }
  if (postscript_flag) {
    width = g_display_r - g_display_l + 1;
    height = g_display_b - g_display_t + 1;
    printf("Postscript file is %s\n", g_post_filename);
    general_post(g_opt_mode, g_post_filename, g_post_file, width, height,
		 g_plot_filename, g_ct_filename, g_optimal_energy,
		 g_energy_cutoff, g_zoom_labels_count, g_prob,
		 g_chain_len, g_number_of_colors, g_color_increment,
		 g__prob, g_zoom_labels_string,
		 g_zoom_labels_row, g_zoom_labels_col,
		 g_display_l, g_display_r, g_display_t, g_display_b,
		 g_diag_start, g_diag_count, g_diag, g_grid_flag,
		 g_plot_len, g_ct_len,
		 g_dot_size,
		 g_ct_diag_start, g_ct_diag_count, g_ct_diag,
		 g_sequence_name, g_sequence_name_set,
		 g_close_optimal, g_worst_energy, g_color,
		 g_opt_distance, g_ct_basepair);
  }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  if (img_flag)
    finish_making_img(g_opt_mode, clear_background_flag, g_img_filename, 
		      g_plot_filename, g_ct_filename, g_optimal_energy,
		      g_energy_cutoff, g_zoom_labels_count, g_prob,
		      g_chain_len, g_number_of_colors, g_color_increment,
		      g__prob, g_zoom_labels_row, g_zoom_labels_col,
		      g_display_l, g_display_r, g_display_t, g_display_b,
		      g_diag_start, g_diag_count, g_diag, g_grid_flag,
		      g_plot_len, g_ct_len, g_dot_size, g_ct_diag_start, 
		      g_ct_diag_count, g_ct_diag, g_sequence_name, 
		      g_sequence_name_set, img_resolution, g_close_optimal, 
		      g_worst_energy, g_color, g_opt_distance, g_ct_basepair);
#endif
  try_exit(0);
}

/* ___________________postscript functions ______________________________ */

int find_color(int energy) {
  int             result;
  if (g_prob == TRUE)
    return (find_color_prob(energy));
  if (energy == g_optimal_energy)
    return (1);
  result = (int)(((float)(energy-g_optimal_energy))/g_color_increment+1.9999);
  return (result);
  /* 1.99 should be right, 1.98 is too big, 1.90 is too small */
}

void memory_error(void) {
  printf("Error!\nInsufficient memory to execute with sequence length %d", 
	 g_len);
  try_exit(12);
}

int main(int argc, char **argv) {
  initialize_data1(argc, argv);
  return 0; /* Keep the stupid compiler happy. */
}
