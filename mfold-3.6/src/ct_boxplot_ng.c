/* This program can creates PostScript  files and image files of type
 * gif if libgd is present. jpg images are created if libjpeg is found
 * and png images are created if libpng (or similar) is found
 *
 * Arguments
 *
 *  -a Structure,   Display this structure in the lower triangle
    -b              Enable clear background for Gif output
    -d              Look for ct files in this directory, ends with /
    -g              Create gif output, default is postscript
    -l i j          Label the point i,j in the GIF or zoomed postscript image.
    -m  value       Magnify dot size by this floating point value.
    -o name_out 
                    Use name_out.gif or name_out.ps for output

   -pg              Create both gif and postscript output
   -r resolution    specify gif resolution from 20 to 300
    -t "TITLE"      Specify a title for output. Otherwise, the first
                    line of the first ct file is used.  
   -y s i j         Zoom with scale s about row i and column j
   -w               create .gifdat for web zooming 
   -z left right top bottom        Specifies a zoom region
*/

#define PROGRAM_NAME "ct_boxplot_ng"
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
#include "ct_boxplot_general.h"
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
# include "ct_boxplot_img.h"
#endif
#include "ct_boxplot_ps.h"
#include "ct_boxplot_read_ct.h"
#include "ct_boxplot_setcolor.h"

int g_png_mode;
int g_jpg_mode;

struct helix g_diag[MAXIMUM_HELICES]; 
char g_sequence_name[160];


int g_sequence_name_set;

float g_color[COLOR_TABLE_SIZE][3];

int g_len;  /* rightmost base that is part of a helix */

int *g_gray3_content;  
int g_lower_structure=-2; /* default to show partial and full overlap
			     when right mouse button is clicked,
			     transform coordinates into i,j */ 

/* Variables that could easily be changed for personal preference are
   usually accompanied by $$$.  Colors in particular                
*/ 

int g_last_structure; /* Number of last structure. It is the number of
			 ct files specified + 1.  Colors 2 to g_last
			 structure are used for drawing dots */

/* Size of dots may be multiplied by a factor for the zoom plot
 */

/* Graphic  variables  -------------------------------------------*/

int g_zoom_labels_count ; /* the number of labels used */
int g_zoom_labels_row[50];
int g_zoom_labels_col[50];

char *g_zoom_labels_string[50]; /* value to place at x,y */

/* File Variables _______________________________________________*/
char rec[90];
FILE *fp;   /* for .ct file */
char g_post_filename[120];
int g_ct_files; /* number of ct files used for input */
char filename_ct[MAXIMUM_CT_FILES][255]; /* up to 255 filenames for input */
char file_data_ct[MAXIMUM_CT_FILES][255]; /* first line of ct file */

char filename[120];
char g_img_filename[120]; /* this depends on -o or input file */

int g_diag_count[2*MAXIMUM_SIZE];
int g_diag_start[2*MAXIMUM_SIZE];
int g_helices_in_plot_file; /* number of entries in file */

/* given a row and column, return the color of that position */

void display_all_helices(void) {
  int i;
 for(i=0;i<g_helices_in_plot_file;i++) {
   printf("helix %d row =%d column=%d color=%d length=%d diag is %d or %d\n",
	  i, g_diag[i].row, g_diag[i].column, g_diag[i].color, 
	  g_diag[i].length, g_diag[i].row+g_diag[i].column-1,
	  g_diag[i].diagonal);
 }
 for(i=1;i<=2*g_len;i++) {
   printf("diag %d starts at %d\n",i,g_diag_start[i]);
 }
}

/* Initialization Routines for nongraphics */
void initialize_labels(void) {
  int i;
  g_zoom_labels_count=0; /* initialize to no labels */
  for(i=0;i<50;i++) { /* allocate memory for string */
    g_zoom_labels_string[i]=(char *)malloc(15*sizeof(char));
    if(g_zoom_labels_string[i]==NULL) {
      printf("Error!\tInsufficient memory for zoom labels.\n");
      try_exit(6);
    }
  }
}

void process_input_files(void) {
  g_len=0;
  g_helices_in_plot_file=0;
  read_of_all_ct_files(fp, g_diag_count, g_diag_start, g_ct_files,
		       &g_helices_in_plot_file, &g_len, filename_ct,
		       g_diag, file_data_ct, g_sequence_name, 
		       g_sequence_name_set);
  initialize_colors(g_color, g_ct_files);
}

void display_bad_arguments(int end_flag) {
  printf("\n Bad  arguments:\n");
  printf("\n[ Usage:          ct_boxplot name.ct                       ]");
  printf("\n[             or  ct_boxplot -o name_out name.ct           ]");
  printf("\n[              (.ct suffix is optional)                      ]");
  printf("\n[                                                            ]");
  printf("\n[                Valid arguments                             ]");
  printf("\n[ -a Struc     Structure to display in lower triangle        ]");
  printf("\n[                (-3<=Struc<=number of ct files)             ]");
  printf("\n[                Default is Full and Partial overlap         ]");
  printf("\n[                -3 is Full and partial, -1 is Full overlap, ]");
  printf("\n[                 0 is Partial overlap, -2 is none           ]");
  printf("\n[                For 3 sequences, -11 produces Partial       ]");
  printf("\n[                    Overlap interpreted for content.        ]");
  printf("\n[                    -13 for Partial Overlap interpreted and ]");
  printf("\n[                     Full overlap.                          ]");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("\n[ -b           Enable clear background for png/jpg  output   ]");
#endif
  printf("\n[ -d direct/   Look for ct files in directory direct         ]");
  printf("\n[ -f           Force gray as multicolor                      ]");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("\n[ -g           Create a Gif file, default is ps              ]");
#endif
  printf("\n[ -i           For 3 sequences, gray is shown as 3 colors    ]");
  printf("\n[                to indicate content                         ]");
  printf("\n[ -j           Draw postscript dots at least 1/72 inch       ]");
#if HAVE_LIBJPEG
  printf("\n[ -jpg         Create jpg output                             ]");
#endif
  printf("\n[ -l i j       Label the point i,j in the output             ]");
  printf("\n[ -m value     Magnify dot size by value. (.2 to 100.)       ]");
  printf("\n[ -o name_out  Specify name of output file as                ]");
  printf("\n[              name_out.ps ,name_out.png or name_out.jpg     ]");
  printf("\n[ -pg          Create both postscript and jpg or png  output ]");
  printf("\n                   use -png or -jpg to specify               ]");
#if HAVE_LIBPNG
  printf("\n[ -png         Create png output                             ]");
#endif
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("\n[ -r  RES      Specify  resolution from 20 to 300 for png/jpg]");
#endif
  printf("\n[ -t \"TITLE\"   Specify TITLE as title, (Quotes Required)     ]");
  printf("\n[ -w           Create .gifdat for webzooming                ]");
  printf("\n[ -y s i j     Zoom with scale s about row i, column j       ]");
  printf("\n[ -z l r t b   Specify a zoom region: left,right,top,bottom  ]");
  printf("\n");
  if(end_flag==TRUE) {
    printf("good-bye\n");
    try_exit(7);
    }
}

void initialize_data(int argc, char **argv) {
  int i;
  int img_resolution;
  int zoom_flag;
  int grid_flag;
  int type_of_output_file;
  int lower_structure;
  int post_file;
  int clear_img_flag;
  int dot_zoom_flag, dot_zoom_row = 396, dot_zoom_column = 306;
  int dot_zoom_show;
  int gifdat_flag;
  int big_dot_flag;
  int nongray_for_3_flag;
  int forced_stripe_flag;
  char ct_directory[255];
  char temp_file_name[130];
  float dot_zoom_scale;
  int temp;
  int ct_directory_flag;
  char number_rec[25];
  float dot_magnifier;
  int display_l = 1, display_r = 700, display_t = 1, display_b = 700;
   /* For output,2 is postscript, 1 is gif, 3 is both */
  char argument[255];
  big_dot_flag = FALSE;
  nongray_for_3_flag = FALSE;
  forced_stripe_flag = FALSE;  
  gifdat_flag = FALSE; /* do not create gifdat flag for web zooming */
  dot_zoom_flag = FALSE; /* zoom about a dot, default is FALSE */
  grid_flag = TRUE; /* Default is to have a grid */
  clear_img_flag = FALSE;/* Default is not to use a clear gif */
  g_sequence_name_set = FALSE;
  ct_directory_flag = FALSE;
  zoom_flag = FALSE; /* Default is to display the full image. This zooms
			based on left, right, top,bottom */ 
  img_resolution = 72; /* Default resolution for image */ 
  i = 1;
  lower_structure = -12; /* default is to show both full and partial
			  * overlap in lower triangle. Show interpreted
			  * if 3 sequences */ 
  dot_magnifier = 1.0; /* default is no magnification */
  type_of_output_file = 2; /* default is postscript */
  g_ct_files = 0;
  strcpy(g_sequence_name,"Name not in first ct file?");
  strcpy(g_post_filename,"ct_boxplot_out");  /* set default filenames */
  strcpy(g_img_filename,"ct_boxplot_out");
  initialize_labels();
  g_png_mode=FALSE;
  g_jpg_mode=FALSE;
  if(argc<2) {
    display_bad_arguments(TRUE);
  }
  while(i <= argc-1) {
    if (strcmp(argv[i],"-d")==0) {
      if(i<(argc-1)) { /* make sure there is an argument to use */
	strcpy(ct_directory,argv[i+1]);
	ct_directory_flag = TRUE;
	i++;
      }
    } else if (strcmp(argv[i],"-o")==0) {
      if (i<(argc-1)) { /* make sure there is an argument to use */
	strcpy(g_post_filename,argv[i+1]);
	strcpy(g_img_filename,argv[i+1]);
	i++;
      }
    } else if (strcmp(argv[i],"-b")==0) {
      clear_img_flag=TRUE;
    } else if (strcmp(argv[i],"-png")==0) {
      g_png_mode=TRUE;
      if(type_of_output_file!=3)
	type_of_output_file=1;
    } else if (strcmp(argv[i],"-jpg")==0) {
      g_jpg_mode=TRUE;
      if(type_of_output_file!=3)
	type_of_output_file=1;
    } else if (strcmp(argv[i],"-i")==0) {
      nongray_for_3_flag=TRUE;
    } else if (strcmp(argv[i],"-j")==0) {
      big_dot_flag=TRUE;
    } else if (strcmp(argv[i],"-f")==0) {
      forced_stripe_flag=TRUE;
    } else if (strcmp(argv[i],"-w")==0) {
      gifdat_flag=TRUE;
    } else if ((strcmp(argv[i],"-a")==0) && (i<(argc-1))) {
      lower_structure = atoi(argv[i+1]);
      if(lower_structure>MAXIMUM_CT_FILES) {
	lower_structure=MAXIMUM_CT_FILES;
	lower_structure += 1;
      }
    } else if ((strcmp(argv[i],"-r")==0) && (i<(argc-1))) {
      img_resolution=atoi(argv[i+1]);
      if(img_resolution<20)
	img_resolution=20;
      if(img_resolution>300)
	img_resolution=300;
      i++;
    } else if ((strcmp(argv[i],"-m")==0)&&(i<(argc-1))) { /*set magnification*/
      sscanf(argv[i+1],"%f",&dot_magnifier);
      if(dot_magnifier>100.)
	dot_magnifier=100.0; 
      else if(dot_magnifier<0.2)
	dot_magnifier=0.2;
      i++;  
    } else if (strcmp(argv[i],"-g")==0) {
      if(type_of_output_file!=3)
	type_of_output_file=1;
    } else if (strcmp(argv[i],"-pg")==0) {
      type_of_output_file=3;
    } else if ( (strcmp(argv[i],"-y")==0) && (i<argc-3) ) {
      dot_zoom_flag = TRUE;
      sscanf(argv[i+1],"%f",&dot_zoom_scale);
      if(dot_zoom_scale<.01)
	dot_zoom_scale =.1;
      if(dot_zoom_scale>200.)
	dot_zoom_scale = 200.;
      dot_zoom_row = atoi(argv[i+2]);
      dot_zoom_column = atoi(argv[i+3]);
      if(dot_zoom_row<1)
	dot_zoom_row=1;
      if(dot_zoom_column<1)
	dot_zoom_column=1;
      i+=3;
    } else if ( (strcmp(argv[i],"-z")==0) && (i<argc-4) ) {
      zoom_flag = TRUE;
      display_l = atoi(argv[i+1]);
      display_r = atoi(argv[i+2]);  
      display_t = atoi(argv[i+3]);
      display_b = atoi(argv[i+4]);
      i+=4;
    } else if ( (strcmp(argv[i],"-l")==0) && (i<argc-2) ) {
      g_zoom_labels_count=1;
      g_zoom_labels_row[0]=atoi(argv[i+1]);
      g_zoom_labels_col[0]=atoi(argv[i+2]);
      i+=2;
    } else if ( (strcmp(argv[i],"-t")==0) && (i<argc-1) ) {
      g_sequence_name_set=TRUE;
      strcpy(g_sequence_name,argv[i+1]);
      i++;
    } else {
      strcpy(argument,argv[i]);
      if (argument[0]=='-') {
	printf("Error!\t%s was not recognized.\n", argument);
	printf("Filenames cannot begin with a hyphen.");
	printf("Supply each flag with sufficient parameters\n.");
	display_bad_arguments(FALSE);
      } else {
	strcpy(filename_ct[g_ct_files],argv[i]);
	if(g_ct_files<MAXIMUM_CT_FILES) {
	  g_ct_files++;
	} else {
	  printf("Error!\tToo many ct files.\n");
	  printf("Only %d will be used.\n",MAXIMUM_CT_FILES);
	}
      }
    }
    i++;
  }
  for (i=0; i < g_ct_files; i++) {
    chop_suffix(filename_ct[i], ".ct");
    strcat(filename_ct[i],".ct");
    if(ct_directory_flag) {
      strcpy(temp_file_name, ct_directory);
      strcat(temp_file_name, filename_ct[i]);
      strcpy(filename_ct[i], temp_file_name);
    }
  }
  if(lower_structure<-2) {
    if(g_ct_files!=3)
      lower_structure=-2;
    else {
      if(lower_structure<-10)
        lower_structure=-12;
      else
        lower_structure=-10;
    }
  }
  if(nongray_for_3_flag) {
    if(g_ct_files!=3) {
      printf("Warning!\t-i option ignored.\n");
      printf("This option requires exactly three ct files.\n");
      nongray_for_3_flag = FALSE;
    }
  }
  if(ct_directory_flag) {
    strcpy(temp_file_name,ct_directory);
    strcat(temp_file_name,g_post_filename);
    strcpy(g_post_filename,temp_file_name);
    strcpy(temp_file_name,ct_directory);
    strcat(temp_file_name,g_img_filename);
    strcpy(g_img_filename,temp_file_name);
    printf("Image file is %s and Postscript file is %s.\n",
	   g_img_filename, g_post_filename);
  }
  g_last_structure=g_ct_files+1;
  strcat(g_post_filename,".ps");
  if(g_png_mode)
    strcat(g_img_filename,".png");
  else if (g_jpg_mode)
    strcat(g_img_filename,".jpg");
  else
    strcat(g_img_filename,".gif");
  process_input_files();
  if((nongray_for_3_flag)||lower_structure<=-10) {
    g_gray3_content = 
      store_gray_content(g_gray3_content, g_helices_in_plot_file, g_len,
			 g_diag,g_diag_start);
  }
  if(dot_zoom_flag) {
    dot_zoom_show = (int)(g_len/dot_zoom_scale);
    if(dot_zoom_show<3) 
      dot_zoom_show=3;
    zoom_flag=TRUE;
    display_l=dot_zoom_column-dot_zoom_show/2;
    display_r=dot_zoom_column+dot_zoom_show/2;
    display_t=dot_zoom_row-dot_zoom_show/2;
    display_b=dot_zoom_row+dot_zoom_show/2;
    /* keep zoom region in bounds without altering zoom scale */
    if(display_l<1) {
      display_r = display_r-display_l;
	display_l = 1;
    } else {
      if(display_r>g_len) {
	display_l = display_l - (display_r-g_len);
	display_r = g_len;
      }
    }
    if(display_t<1) {
      display_b = display_b-display_t;
        display_t = 1;
    } else {
      if(display_t>g_len) {
	display_b = display_b - (display_t-g_len);
          display_t = g_len;
      }
    }
  }
  if(zoom_flag) {
    if((display_l>g_len) || (display_l<1))
      display_l=1;
    if((display_r>g_len) || (display_r<1))
      display_r=g_len;
    if((display_t>g_len) || (display_t<1))
      display_t=1;
    if((display_b>g_len) || (display_b<1))
      display_b=g_len;
    /* keep left smaller than right */
    if(display_r<display_l) {
      temp = display_r;
      display_r = display_l;
      display_l = temp;
    }
     /* keep top smaller than bottom */
     if(display_b<display_t) {
       temp = display_b;
       display_b = display_t;
       display_t = temp;
     }
  } else {
    display_l = 1;
     display_r = g_len;
     display_t = 1;
     display_b = g_len;
  }
  if(g_zoom_labels_count==1) {
    /* Force label within bounds and convert to a string */
    if(g_zoom_labels_row[0]<display_t) {
      g_zoom_labels_row[0]=display_t;
    } else {
      if(g_zoom_labels_row[0]>display_b) {
	g_zoom_labels_row[0]=display_b;
      }
    }
    if(g_zoom_labels_col[0]<display_l) {
      g_zoom_labels_col[0]=display_l;
    } else {
      if(g_zoom_labels_col[0]>display_r) {
	g_zoom_labels_col[0]=display_r;
      }
    }
    strcpy(g_zoom_labels_string[0],"\\(");
    sprintf(number_rec,"%d",g_zoom_labels_row[0]);
    strcat(g_zoom_labels_string[0],number_rec);
    strcat(g_zoom_labels_string[0],",");
    sprintf(number_rec,"%d",g_zoom_labels_col[0]);
    strcat(g_zoom_labels_string[0],number_rec);
    strcat(g_zoom_labels_string[0],"\\)");
  }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  if(type_of_output_file!=2)
    finish_making_img(TRUE, g_img_filename, img_resolution,
		      g_zoom_labels_count, g_zoom_labels_row,
		      g_zoom_labels_col, grid_flag, g_last_structure,
		      g_diag, g_color, display_l, display_r, display_t, 
		      display_b, lower_structure, g_sequence_name, 
		      dot_magnifier, g_diag_start, g_diag_count, file_data_ct,
		      clear_img_flag, gifdat_flag, g_len, forced_stripe_flag,
		      nongray_for_3_flag, g_gray3_content, g_png_mode,
		      g_jpg_mode);
#endif
  if(type_of_output_file!=1) {
    if(!zoom_flag)
      post_file=MAIN_POST; /* create post for main */
    else
      post_file=ZOOM_POST; /* create post for zoom */
    finish_making_post(big_dot_flag, post_file, g_post_filename,
		       g_zoom_labels_count, display_l, display_r, display_t,
		       display_b, g_color, g_zoom_labels_row, 
		       g_zoom_labels_col, g_zoom_labels_string, grid_flag,
		       g_last_structure, g_diag_start, g_diag_count, g_diag,
		       g_len, lower_structure, g_sequence_name,
		       file_data_ct, dot_magnifier, forced_stripe_flag,
		       nongray_for_3_flag, g_gray3_content);
  }
}
 
int main(int argc, char **argv)  {
  printf("%s: %s_ng\n", PACKAGE_STRING, PROGRAM_NAME);
  initialize_data(argc,argv);
  try_exit(0);
  return 0; /* Keep the stupid compiler happy. */
}
