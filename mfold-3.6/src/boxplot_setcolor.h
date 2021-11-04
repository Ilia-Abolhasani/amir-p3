/* boxplot_setcolor.h
 * Read boxplot.col to set colors and probabilities for each color
 */

#define TOTAL_COLORS 13
#define LAST_PROB_COLOR 8
#define PWD "PWD"
#define COLOR_FILE "boxplot.col"
#define COLOR_BACKGROUND 10
#define COLOR_LABEL 9
#define COLOR_TEXT  0
#define COLOR_GRID  11
#define COLOR_OPTIMAL 1
#define COLOR_COPY  12

float g_color[TOTAL_COLORS][3]; /* 0 for red, 1 for green, 2 for blue */
/* color 0 is for text and outlines
 * colors 1 through 7 are the 7 available colors for dots
 * 8 is for dot labels
 * 9 is background
 */
float g_prob_color[9][3]; /* for 8 color probability */

/* define color for each range of probabilities
 * 4 color scheme for 0
 * 5 color scheme for 1
 * 6 color scheme for 2
 * 7 color scheme for 3
 * 8 color scheme for 4
 * values 0 through 5 are the range of prob. for each color 1 to 7
 * 5 is used only with 7 colors 
 */

float g__prob[5][8];

FILE *cfp;

char color_filename[80];

int color_gl, color_bak, color_b0, color_bk, color_OP, color_re, color_gr,
  color_ye, color_pu, color_br, color_bl, color_fl, color_la, color_co;
const int* const img_color[] = {&color_b0, &color_bk, &color_re, &color_gr,
				&color_ye, &color_pu, &color_br, &color_bl,
				&color_fl, &color_la, &color_bak, &color_gl,
				&color_co};

/* 1-7 are colors,
 *for 2nd  0 is dot magnifier, 1 is miniminum size in points 
 */
float g_ps_dot_min[8][2]; 

int g_number_of_colors; /* assumes values from 4 to 8; default is 4 */
int g_num_color_set = FALSE;

/* The next three functions are from dot_col.h
 * Define and declare
 * #define COLOR_FILE "xxx.col"
 * FILE *cfp;
 * char color_filename[80];
*/

void read_line(char *rec) { 
  if (fgets(rec,90,cfp)==NULL) {
    printf("The %s file is too short.\n",COLOR_FILE);
    try_exit(13);
  }
}

void open_error_color(void) {
  printf("Error!\tThe file %s could not be found.",COLOR_FILE);
  printf("The installation of %s was defective, since the color file, ",
	 PACKAGE_STRING);
  printf("%s, was not placed in %s.\n",COLOR_FILE,PKGDATADIR);
  try_exit(14);
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
  int i,k;
  float prob;
  float size_mag,size_min;
  open_color();
  read_line(rec); /* read first line , ignore it*/
  for(i=0;i<TOTAL_COLORS;i++) { /* read and set colors */
    read_line(rec);
    sscanf(rec,"%f %f %f", &g_color[i][0],&g_color[i][1],&g_color[i][2]);
  }
  /* Set k color probability scheme, k = 4, ... , 8. In each case,
   * skip the first line.
   */
  for (k=4;k<=8;k++) {
    read_line(rec);
    for(i=0;i<=k-2;i++) {
      read_line(rec);
      sscanf(rec,"%f",&prob);
      g__prob[k-4][i]=prob;
    }
  }
  /* set dot magfier and size for postscript */
  read_line(rec);
  read_line(rec);
  for(i=1;i<=7;i++) {
    read_line(rec);
    sscanf(rec,"%f%f",&size_mag,&size_min);
    g_ps_dot_min[i][0]=size_mag;
    g_ps_dot_min[i][1]=size_min;
  }
  /* read 8 probability colors */
  read_line(rec); /* read first line , ignore it*/
  for(i=1;i<=LAST_PROB_COLOR;i++) { /* read and set colors */
    read_line(rec);
    sscanf(rec,"%f %f %f", &g_prob_color[i][0],&g_prob_color[i][1],
	   &g_prob_color[i][2]);
  }
  fclose(cfp);
}

/* set color for probability
 * probabilities have been altered by:
 * energy(int)(-1*(1000000.*rec_energy+.5));
 */
int find_color_prob(int energy) {
  int scheme; /* 0 for 4 colors, 3 for 7 colors */
  float prob; /* convert back from energy */
  int current_range;
  scheme=g_number_of_colors-4;
  /* convert back to float
   * The following test is for speed
   * It eliminates the need for the bottom line 
   */
  prob = ((float)energy)/-1000000.; 
  if(prob<=g__prob[scheme][g_number_of_colors-2])
    return g_number_of_colors;
  for(current_range=0;current_range<=(g_number_of_colors-2);current_range++) {
    if(prob>g__prob[scheme][current_range]) {
      return current_range+1;
    }
  }
  /* the bottom line should never be reached within the 
   * range of lowest probability */
  return g_number_of_colors; 
}

int find_color(int energy, float color_increment, int optimal_energy) {
  int result;
  if(g_prob==TRUE)
    return (find_color_prob(energy));
  if(energy==optimal_energy)
    return (1);
  result=(int) (((float)(energy-optimal_energy))/color_increment + 1.9999);
  return (result);
}
