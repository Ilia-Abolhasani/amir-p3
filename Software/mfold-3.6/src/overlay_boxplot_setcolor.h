/*
 * Read overlay_boxplot.col to set colors and probabilities for each color 
 */

#define COLOR_FILE "overlay_boxplot.col"
#define TOTAL_COLORS 12
float g_color[TOTAL_COLORS][3];
/* 0 for text and outlines
 * colors 1 for optimal energy
 * color 2-4 for next energy levers
 * color 5 for ct overlap with optimal
 * color 6 for ct overlap with near optimal
 * color 7 for ct 0verlap with non optimal
 * color 8 for background of gif and screen
 * color 9 for label dots
 * color 10 for grid lines 
 */

#define COLOR_TEXT 0
#define COLOR_PLOT_OPTIMAL 1
#define COLOR_OPTIMAL 5
#define COLOR_NEAR_OPTIMAL 6
#define COLOR_NOT_OPTIMAL 7
#define COLOR_BACK 8
#define COLOR_LABEL 9
#define COLOR_GRID 10
#define COLOR_CT_MISSED_OPTIMAL 11

float g_dot_size[4];
/* defines size of dots
 * 0 for optimal gray dot, 1-3 for next gray dots
 * 4 for green optimal dot, 5 for yellow, 6 for red 
 */

FILE *cfp;
char color_filename[80];

int g_number_of_colors;
/* This is always 4 for plot file dots */

/* define color for each range of probabilities
 * 4 color scheme for 0
 * 5 color scheme for 1
 * 6 color scheme for 2
 * 7 color scheme for 3
 * values 0 through 5 are the range of prob. for each color 1 to 7
 * 5 is used only with 7 colors 
 */
float g__prob[4][6];

/* Define and declare
 * #define COLOR_FILE "xxx.col"
 * FILE *cfp;
 * char color_filename[80];
 */

void read_line(char *rec) { 
  if (fgets(rec,90,cfp)==NULL) {
    printf("The %s file is too short.\n",COLOR_FILE);
    try_exit(23);
  }
}

void open_error_color(void) {
  printf("Error!\tThe file %s could not be found.",COLOR_FILE);
  printf("The installation of %s was defective, since the color file, ",
	 PACKAGE_STRING);
  printf("%s, was not placed in %s.\n",COLOR_FILE,PKGDATADIR);
  try_exit(24);
}

void open_color(void) {
  strcpy(color_filename,COLOR_FILE);
  if ( (cfp=fopen(color_filename,"r"))!=NULL ) {
    printf("The color file, %s, was found in the current working directory.\n",
	   COLOR_FILE);
    return; /* COLOR_FILE was found in pwd */
  }
  strcpy(color_filename, PKGDATADIR);
  strcat(color_filename, COLOR_FILE);
  if ((cfp = fopen(color_filename, "rt")) != NULL)
    return; /* COLOR_FILE was found in PKGDATADIR */
  open_error_color(); /* was not found in either directory */
}

void initialize_colors(void) {
  char rec[90];
  int i;
  float prob;
  open_color();
  read_line(rec); /* read first line , ignore it*/
  for(i=0;i<TOTAL_COLORS;i++) { /* read and set colors */
    read_line(rec);
    sscanf(rec,"%f %f %f",&g_color[i][0],&g_color[i][1],&g_color[i][2]);
  }
  read_line(rec); /* read and ignore line */
  read_line(rec); /* read and ignore line */
  /* set dot sizes */
  for(i=0;i<4;i++) {
    read_line(rec);
    sscanf(rec,"%f",&prob);
    g_dot_size[i]=prob;
  }
  read_line(rec); /* read and ignore line */
  read_line(rec);
  sscanf(rec,"%f",&g_close_optimal);
  fclose(cfp);
}

/* set color for probability
 * probabilities have been altered by:
 * energy(int)(-1*(1000.*rec_energy+.5));
 */
int find_color_prob(int energy) {
  int scheme; /* 0 for 4 colors, 3 for 7 colors */
  float prob; /* convert back from energy */
  int current_range;
  scheme = g_number_of_colors - 4;
  prob=((float)energy)/-1000.; /* convert back to float 
				* The following test is for speed
				* It eliminates the need for the bottom line 
				*/
  if(prob<=g__prob[scheme][g_number_of_colors-2])
    return g_number_of_colors;
  for(current_range=0;current_range<=(g_number_of_colors-2);current_range++) {
    if(prob>g__prob[scheme][current_range]) return current_range+1;
  }
  /* the bottom line should never be reached */
  return g_number_of_colors; /* within range of lowest probability */
}
