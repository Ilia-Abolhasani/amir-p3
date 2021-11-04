/* edit colors in main_color_table as you like, rgb format 
 * 00 to FF for each 
 */

/* #define COLOR_BACKGROUND 0 */
#define COLOR_BACKGROUND 1 
#define COLOR_BASES 2
#define COLOR_TEXT 3
#define COLOR_MAIN_CIRCLE 4
#define COLOR_CONNECTING_GC 5
#define COLOR_CONNECTING_AU 6
#define COLOR_CONNECTING_GU 7
#define COLOR_CONNECTING_OTHERS 8
#define COLOR_COPYRIGHT 9
#define COLOR_LABEL_NUMBER 10
#define COLOR_LABEL_LINE 11
#define COLOR_CENTERS    12
#define COLOR_ANN_BLACK_LETTER 13
#define COLOR_ANN_WHITE_LETTER 14

#define MAIN_COLORS 15
char            main_color_table_let[MAIN_COLORS] = " kbtmgauocnlsxz";

char            color_filename[120];
FILE           *COLOR_FP;
int             main_color_table[MAIN_COLORS] = {
  0xffdead,	        /* 0 MZ adds this as alternate background color */
  0xffffff,		/* 1 background */
  0x000000,		/* 2 BASES */
  0x000020,		/* 3 text */
  0x303030,		/* 4 Main Circle */
  0xFF0000,		/* 5 connecting lines GC */
  0x4030A6,		/* 6 connecting lines AU */
  0x107F10,		/* 7 connecting lines GU */
  0xd0ce14,		/* 8 connecting lines other */
  0x000050,		/* 9 for trademark */
  0x000020,		/* 10 for number of label */
  0x242424,		/* 11 for line of label  */
  0x000000,
  0x000000,
  0xFFFFFF};

int img_color[15];
/* const char  *img_color[] = {"b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7",
   "b8", "b9", "b10", "b11", "b12", "b13", "b14"}; */

struct color_table_struct_ps {
  float           red;
  float           green;
  float           blue;
};

struct color_table_struct_ps main_color_table_float[MAIN_COLORS];

int try_open_color(void) {
  if ((COLOR_FP = fopen(color_filename, "r")) == NULL) {
    printf("Could not open color file %s\n", color_filename);
    return FALSE;
  } else
    return TRUE;
}

void set_main_colors(int color_file_exists) {
  int             i, red, green, blue, remainder;
  char            record[120];
  if (color_file_exists)
    color_file_exists = try_open_color();
  if (color_file_exists) {
    for (i = 0; i < MAIN_COLORS; i++) {
      if (fgets(record, 120, COLOR_FP) == NULL) {
	printf("Line %d of color file %s does not exist\n", i, color_filename);
	try_exit(31);
      }
      sscanf(record, "%f %f %f", &main_color_table_float[i].red,
	     &main_color_table_float[i].green, 
	     &main_color_table_float[i].blue);
    }
    fclose(COLOR_FP);
  } else {
    for (i = 0; i < MAIN_COLORS; i++) {
      red = main_color_table[i] / 65536;
      remainder = main_color_table[i] % 65536;
      green = remainder / 256;
      blue = remainder % 256;
      main_color_table_float[i].red = (float) red / 255.;
      main_color_table_float[i].green = (float) green / 255.;
      main_color_table_float[i].blue = (float) blue / 255.;
    }
  }
}
