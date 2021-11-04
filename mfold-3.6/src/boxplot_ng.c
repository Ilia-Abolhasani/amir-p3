#define PROGRAM_NAME "boxplot_ng"
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>	/* for sqrt */
#include <time.h>	/* for time within dotplot files */
#include <limits.h>	/* for INT_MAX */

#define TRUE    1
#define FALSE   0
#define MAIN_POST    1	/* create postcript for main or zoom */
#define ZOOM_POST    2	/* used in general_post, create_post */

#include "util.h"

int             g_helix_array_size;
int             g_length;	/* len of RNA sequence, guessed at based on
				 * input file */
/* Variables that could easily be changed for personal preference are
   usually accompanied by $$$.  */ 

/* Graphic  variables  ------------------------------------------- */
struct helix {
	int             row;
	int             column;
	int             diagonal;
	int             color;
	int             energy;
	int             length;
};

float           g_color_increment;	/* range of energy for each color */
int             g_energy_cutoff;	/* only probabilities <= */
/* g_energy_cutoff will be shown. */
/* note that probabilities have a negative sign */
/* -30 indicates probabilities >=.3 will be shown */
float           g_energy_cutoff_stored;
int             g_energy_cutoff_set;
int             g_chain_len = 1;	/* chains smaller will not be
					 * displayed */

int  g_points_plotted, g_points_plotted_zoom; /* indicates number of
						 points plotted in
						 each window */ 

/* Zoom Variables __________________________________________ */

int   display_l, display_r, display_t; 
int   display_b, display_w, display_h; /* left, right, top, bottom,
					* width and height for zoom
					* window display in nucleotide
					* coordinates */ 
int             g_display_from_arg;

/* File Variables _______________________________________________ */

FILE           *fp;		/* for .plot file */
char            rec[90];
char            g_post_filename[120] = "";
char            g_plot_filename[120] = "";	/* with .plot */
char            g_output_name[120] = "";
char            g_name_of_file_img[120] = "";
char            g_name_of_file_ps[120] = "";
char            g_title[150] = "";
int             g_prob, g_png_mode, g_jpg_mode, g_opt_prob_flag, g_mi_flag;
int             g_clear_flag, g_make_imgdat_flag, g_make_grid_flag;
int             g_make_label_flag, g_label_i, g_label_j, g_post_adjust;
int             g_type_of_file, g_resolution;
int             g_adjust_col_flag = FALSE;
int             g_optimal_energy, g_worst_energy;
int             g_post_file;	/* 1 for main, 2 for zoom */
int            *g_diag_count;
int            *g_diag_start;	/* 10000 is the maximum size */
int             g_helices_in_plot_file;	/* number of entries in file */

struct helix   *g_diag;

/* _________________________________________________________________ */


#include "boxplot_setcolor.h"
#include "boxplot_input.h"
#include "boxplot_ps.h"
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
#include "boxplot_img.h"
#endif

/* ____________________________________________________________________
 * Initialization Routines for nongraphics */

void initialize_data2(void) {
  g_energy_cutoff = g_worst_energy;
  if (g_prob != TRUE)
    g_color_increment = ((float) (g_energy_cutoff - g_optimal_energy)) / 
      (float) (g_number_of_colors - 1);
  initialize_colors();
}

void initialize_data1(int argc, char **argv) {
  open_file(argc, argv);
  initialize_len(FALSE, 0);
}


/* ____________________________________________________________________ */

void set_energy_cutoff1(float cutoff_level) { /* for probabilities */
  g_energy_cutoff = (int) (cutoff_level * -1000000);
}

/* This function is called from the menu and sets the filter */
void set_chain_len(int chain_len_choice) {
  g_chain_len = chain_len_choice;
  
}

void check_parameters(int argc, char **argv) {
  int             number_of_colors, filter;
  int             i;
  int             title_set;
  char            title2[150];
  float           cutoff;
  strcpy(g_plot_filename, argv[argc - 1]);	/* set output file */
  argc = argc - 1;	/* last 1 is input file and has been used */
  g_type_of_file = 2; /* default type is postscript */
  g_make_imgdat_flag = FALSE;	/* default is no data file
				 * for image zooming */
  g_make_grid_flag = TRUE;	/* Default is to make grid lines */
  g_make_label_flag = FALSE;	/* Default is no labels */
  g_label_i = 0;	/* default label positions */
  g_label_j = 0;
  g_png_mode = FALSE;
  g_jpg_mode = FALSE;
  g_clear_flag = FALSE;	/* default is not a clear background for image */
  g_prob = FALSE;
  g_energy_cutoff_set = FALSE;
  g_display_from_arg = FALSE;
  g_mi_flag = FALSE;
  strcpy(g_title, g_plot_filename);  /* Set Default Title */
  chop_suffix(g_title, ".plot"); /* Title default: input file without .plot */
  strcpy(g_plot_filename, g_title);  /* set default output file */
  strcpy(g_output_name, g_title);
  strcat(g_plot_filename, ".plot");
  title_set = FALSE;
  filter = 1;		/* default is to show all helices */
  number_of_colors = 4;	/* default is 4 colors for dotplot */
  /* set zoom region */
  g_resolution = 72;
  i = 1;	/* check each argument */
  while (i <= (argc - 1)) {
    if(strcmp(argv[i], "-m")==0) {
      g_post_adjust = TRUE;
      i++;
    } else if(strcmp(argv[i], "-p")==0) {
      g_prob = TRUE;
      i++;
    } else if(strcmp(argv[i], "-g")==0) {
      if (g_type_of_file != 3)
	g_type_of_file = 1;	/* create a gif file */
      i++;
    } else if(strcmp(argv[i], "-pg")==0) {
      g_type_of_file = 3; /* create an image and a postscript file */
      i++;
    } else if(strcmp(argv[i], "-mi")==0) {
      g_mi_flag = TRUE;	/* Treat input data as mutual information */
      g_prob = TRUE;	/* also as probability */
      i++;
    } else if(strcmp(argv[i], "-png")==0) {
      if (g_type_of_file != 3)
	g_type_of_file = 1;
      g_png_mode = TRUE;	/* create png, not gif */
      i++;
    } else if(strcmp(argv[i], "-jpg")==0) {
      if (g_type_of_file != 3)
	g_type_of_file = 1;
      g_jpg_mode = TRUE;	/* create jpg, not gif */
      i++;
    } else if(strcmp(argv[i], "-b")==0) {
      g_clear_flag = TRUE;
      i++;
    } else if(strcmp(argv[i], "-d")==0) {
      g_make_imgdat_flag = TRUE;	/* create .gifdat file */
      i++;
    } else if(strcmp(argv[i], "-go")==0) {
      g_make_grid_flag = FALSE;	/* No grid lines  */
      i++;
    } else if(strcmp(argv[i], "-l")==0) {
      g_make_label_flag = TRUE;	/* Label extra point and make energy file */
      g_label_i = atoi(argv[i + 1]);
      g_label_j = atoi(argv[i + 2]);
      i += 3;
    } else if(strcmp(argv[i], "-o")==0) {
      strcpy(g_output_name, argv[i + 1]);
      i += 2;
    } else if(strcmp(argv[i], "-t")==0) {
      strcpy(g_title, argv[i + 1]); /* Title must follow */
      g_title[strlen(g_title)-1] = '\0';
      title_set = TRUE;
      i += 2;
    } else if(strcmp(argv[i], "-f")==0) {
      filter = atoi(argv[i + 1]);	/* filter must follow -f */
      if (filter < 1)
	filter = 1;
      i += 2;
    } else if(strcmp(argv[i], "-c")==0) {
      number_of_colors = atoi(argv[i + 1]);
      if (number_of_colors > 8)
	number_of_colors = 8;
      else if(number_of_colors < 2)
	number_of_colors = 4;
      i += 2;
    } else if(strcmp(argv[i], "-z")==0) {
      if ((i + 4) > (argc - 1)) {
	printf("Error!\tImproper use of -z flag\n");
	g_display_from_arg = FALSE;
	i++;
      } else {
	g_display_from_arg = TRUE;
	display_l = atoi(argv[i + 1]);
	display_r = atoi(argv[i + 2]);
	display_t = atoi(argv[i + 3]);
	display_b = atoi(argv[i + 4]);
	i += 5;
      }
    } else if(strcmp(argv[i], "-i")==0) { /* set increment */
      sscanf(argv[i + 1], "%f", &cutoff); /* copy float into cutoff */
      g_energy_cutoff_set = TRUE;
      g_energy_cutoff_stored = cutoff;
      i += 2;
    } else if(strcmp(argv[i], "-r")==0) {
      /* set image resolution */
      g_resolution = atoi(argv[i + 1]);
      /* Resolution must follow */
      if (g_resolution < 20)
	g_resolution = 20;
      else if (g_resolution > 400)
	g_resolution = 400;
      i += 2;
    } else {
      printf("Warning!\tFlag %s was not used.\n", argv[i]);
      i++;
    }
  }
  /* adjust title */
  if ((g_prob) && (g_number_of_colors < 4))
    g_number_of_colors = 4;
  if (title_set != TRUE) {
    if (g_prob == TRUE) {
      if (g_mi_flag)
	strcpy(title2, "Mutual Information Dotplot for ");
      else
	strcpy(title2, "Probability Dotplot for ");
    } else
      strcpy(title2, "Energy Dotplot for ");
    strcat(title2, g_title);
    strcpy(g_title, title2);
  }
  if (g_prob) {
    if (filter > 1) {
      printf("Error!\tThe filter option is not available with probability ");
      printf("or Mutual Information.\n");
      printf("All helix lengths=1, (i.e. only individual base pairs exist.\n");
      filter = 1;
    }
  }
  set_chain_len(filter);
  g_number_of_colors = number_of_colors;
  if (g_prob) {
    if (g_mi_flag)
      printf("Treating data as Mutual Information\n");
    else {
      if (g_opt_prob_flag)
	printf("Optimal Structure detected in Probability Data\n");
      else
	printf("Treating data as Probability\n");
    }
  }
}

void finish_making_post(char *text_window_input,char *title,int make_grid_flag,
			int make_label_flag, int label_i,int label_j, 
			int adjust_post, int mi_flag, int opt_prob_flag) {
  int             error;
  printf("Postscript file name is %s\n", text_window_input);
  strcpy(g_post_filename, text_window_input);
  if (strlen(g_post_filename) == 0) {
    printf("Error!\tPostscript filename has zero length.\n");
    return;
  }
  error = open_post(g_post_filename);
  if (error) {
    printf("Error creating postscript file.\n");
    return;
  }
  general_post(title, make_grid_flag, make_label_flag, label_i, label_j,
	       adjust_post, 0, NULL, NULL);
}

void make_output(void) {
  int             temp;
  float           cutoff;
  if (g_energy_cutoff_set) {
    cutoff = g_energy_cutoff_stored;
    if (g_prob == TRUE) {
      if (cutoff > 3.)
	cutoff = 3.0;
      else if (cutoff < 0.0)
	cutoff = 0.0;
      set_energy_cutoff1(cutoff);
    } else {
      finish_setting_energy_cutoff(cutoff);
    }
  }
  if (g_display_from_arg == TRUE) {
    if ((display_l > g_length) || (display_l < 1))
      display_l = 1;
		if ((display_r > g_length) || (display_r < 1))
		  display_r = g_length;
		if ((display_t > g_length) || (display_t < 1))
		  display_t = 1;
		if ((display_b > g_length) || (display_b < 1))
		  display_b = g_length;
		/* keep left smaller than right */
		if (display_r < display_l) {
		  temp = display_r;
		  display_r = display_l;
		  display_l = temp;
		}
		/* keep top smaller than bottom */
		if (display_b < display_t) {
		  temp = display_b;
		  display_b = display_t;
		  display_t = temp;
		}
  } else {
    display_l = 1;
    display_t = 1;
    display_r = g_length;
    display_b = g_length;
  }
  display_w = display_r - display_l + 1;
  display_h = display_b - display_t + 1;
  if (g_type_of_file != 1) {
    strcpy(g_name_of_file_ps, g_output_name);
    strcat(g_name_of_file_ps, ".ps");
    if ((display_l == 1) && (display_t == 1) &&
	(display_r == g_length) && (display_b == g_length) &&
	(g_prob != TRUE))
      g_post_file = MAIN_POST;
    else
      g_post_file = ZOOM_POST;
    finish_making_post(g_name_of_file_ps, g_title, g_make_grid_flag,
		       g_make_label_flag, g_label_i, g_label_j, g_post_adjust,
		       g_mi_flag, g_opt_prob_flag);
  }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  if (g_type_of_file != 2) {
    strcpy(g_name_of_file_img, g_output_name);
    if (g_png_mode)
      strcat(g_name_of_file_img, ".png");
    else if (g_jpg_mode)
      strcat(g_name_of_file_img, ".jpg");
    else
      strcat(g_name_of_file_img, ".gif");
    finish_making_img(g_clear_flag, g_resolution, g_name_of_file_img, g_title,
		      g_make_imgdat_flag, g_make_label_flag, g_label_i, 
		      g_label_j, g_make_grid_flag, g_png_mode, g_jpg_mode, 
		      g_mi_flag, g_opt_prob_flag, 0, NULL, NULL);
  }
#endif
}

int main(int argc, char **argv) {
  printf("%s, %s\n", PROGRAM_NAME, PACKAGE_STRING);
  if ((argc < 2) || (argv[argc - 1][0] == '-')) {
    printf("Usage: %s [options] file\n", PROGRAM_NAME);
    printf("\twhere file.plot is a plot file\n\n");
    printf("Options:\n");
    printf("            Default is to create a PostScript file\n");
    printf("-b          Enable clear background for jpg,png file output.\n");
    printf("-c colors   Specify 2 to 8 colors for dots in plot.\n");
    printf("-d          Create .gifdat file for www zooming\n");
    printf("-f filter   Display only helices of length >= filter.\n");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
    printf("-g           Create only a gif file for output.\n");
#endif
    printf("-go          Grid lines off, default is on\n");
    printf("-i incr      Specify energy increment or probability cutoff.\n");
    printf("-l i j       label i,j and create file.gifeng containing the energy/probability at i,j\n");
#if HAVE_LIBJPEG
    printf("-jpg         Create jpg output\n");
#endif
    printf("-m           Magnify dots in postscript output\n");
    printf("-mi          Treat input data as Mutual Information\n");
    printf("-o name_out  Specify name of output file.\n");
    printf("-p           Treat input as probability (Default is energy)\n");
    printf("-pg          Create postscript and a png/jpg file\n");
    printf("             Use with -png or -jpg\n");
#if HAVE_LIBPNG		
    printf("-png         Create png output\n");
#endif
    printf("-r res       Specify resolution of png/jpg file (50 to 300)\n");
    printf("-t \"TITLE\" Specify TITLE (protect with quotes\n");
    printf("-z l r t b   Specify a zoom region: left, right, top, bottom\n");
    printf("This program reads %s.col to determine colors.\n", PROGRAM_NAME);
    printf("See %s.doc for more information.\n",   PROGRAM_NAME);
    try_exit(4);
  }
  initialize_data1(argc, argv);
  make_output();
  free(g_diag);
  try_exit(0);
  return 0; /* Keep the stupid compiler happy. */
}
