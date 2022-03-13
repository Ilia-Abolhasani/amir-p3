/* ct_boxplot_setcolor.h */

#define COLOR_FILE "ct_boxplot.col"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 0 for red, 1 for green, 2 for blue 
 * color 0 is overlap
 * color 1 is partial overlap
 * colors 2 through 16 are the 15 available colors for dots
 * 202 is for dot labels
 * 203 is background 
 */

FILE *cfp;
char color_filename[255];

void read_line(char *rec) { 
  if (fgets(rec,90,cfp)==NULL) {
    printf("The %s file is too short.\n",COLOR_FILE);
    try_exit(19);
  }
}

void open_error_color(void) {
  printf("Error!\tThe file %s could not be found.", COLOR_FILE);
  printf("The installation of %s was defective, since the color file, ",
	 PACKAGE_STRING);
  printf("%s, was not placed in %s.\n", COLOR_FILE, PKGDATADIR);
  try_exit(20);
}

void open_color(void) {
  strcpy(color_filename, COLOR_FILE);
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

void initialize_colors(float color_table[][3], int g_ct_files) {
  char rec[90];
  int i, j;
  open_color();
  read_line(rec); /* read first line , ignore it*/
  for(i=0; i<=20; i++) { /* read and set colors */
    read_line(rec);
    if (i <= 16) 
      j = i;
    else
      j = i + 185;
    sscanf(rec, "%f %f %f", &color_table[j][0], &color_table[j][1],
	   &color_table[j][2]);
  }
  for(i=17; i<=201; i++) {
    color_table[i][0] = color_table[16][0];
    color_table[i][1] = color_table[16][1];
    color_table[i][2] = color_table[16][2];
  }
  /* if there are only 2 ct files, copy 5 to 2 */
  if(g_ct_files < 3) {
    color_table[2][0] = color_table[5][0];
    color_table[2][1] = color_table[5][1]; 
    color_table[2][2] = color_table[5][2];
  }
  fclose(cfp);
}
