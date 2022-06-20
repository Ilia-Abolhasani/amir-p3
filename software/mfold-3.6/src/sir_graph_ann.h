#include <stdio.h>
#include <math.h>

#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
#include <gd.h>
#include <gdfontg.h>
#include <stdarg.h>

extern gdFontPtr gdFontGiant;

void gdImageStringCenter(gdImagePtr, gdFontPtr, int, int, char*, int);
void gdImagePrintfCenter(gdImagePtr, gdFontPtr, int, int, int, char*, ...);
#endif

FILE *ctab; /* for color table */
char  ctab_filename[90];
int  *g_ann_to_color; /* ann_to_color[i] holds the color for for base i */ 
int   g_bases_in_file;

/* The 3 values are computed and stored, they are drawn last in order to
 * to be placed on top of all other drawn characters */

int color_table_img[101];

struct color_table_struct_ps color_table_ps[NUM_COLORS + 1];

struct color_table_struct_img {
	int             red;
	int             green;
	int             blue;
};

FILE           *ann_file;	/* used for ann or ss-count file */

/* the extra colors required for the colored bases are created here */

#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
void set_extra_img_colors(gdImagePtr image, int log_flag) {
  int             red, green, blue, i;
  int             remainder;
  if (log_flag) {
    for (i = 0; i <= LOG_NUM_COLORS; i++) {
      red = log_color_table[i] / 65536;
      remainder = log_color_table[i] % 65536;
      green = remainder / 256;
      blue = remainder % 256;
      color_table_img[i] = gdImageColorAllocate(image, red, green, blue);
    }
  } else {
    for (i = 0; i <= NUM_COLORS; i++) {
      red = color_table[i] / 65536;
      remainder = color_table[i] % 65536;
      green = remainder / 256;
      blue = remainder % 256;
      color_table_img[i] = gdImageColorAllocate(image, red, green, blue);
    }
  }
}
#endif

int get_color(double, int, int);

void set_extra_ps_colors(int prob_flag) {
  int             red, green, blue, i, remainder;
  if (prob_flag) {
    for (i = 0; i <= LOG_NUM_COLORS; i++) {
      red = log_color_table[i] / 65536;
      remainder = log_color_table[i] % 65536;
      green = remainder / 256;
      blue = remainder % 256;
      color_table_ps[i].red = (float) red / 256.;
      color_table_ps[i].green = (float) green / 256.;
      color_table_ps[i].blue = (float) blue / 256.;
    }
  } else {
    for (i = 0; i <= NUM_COLORS; i++) {
      red = color_table[i] / 65536;
      remainder = color_table[i] % 65536;
      green = remainder / 256;
      blue = remainder % 256;
      color_table_ps[i].red = (float) red / 255.;
      color_table_ps[i].green = (float) green / 255.;
      color_table_ps[i].blue = (float) blue / 255.;
    }
  }
}

void ann_fix_name(char *filename) {
  int             i;	/* this function strips off the _structure
			 * from the name */
  int             len;
  len = strlen(filename);
  i = len - 1;
  while ((i > 1) && isdigit(filename[i])) i--;
  if ((i > 0) && (filename[i] == '_')) {
    filename[i] = '\0';
    return;
  }
}

int open_ann(char *filename, int g_ann, int g_prob_ann,
	     int no_name_change, char *specific_ann_file) {
  /* if g_ann is true, strip off _structure if it is there.
   * Add  .ann when g_ann is true, else .ss-count
   * if g_prob_ann is true, leave on Underline structure,
   * Then try with it off, add .ss-count
   * if no_name_change is FALSE, add .ann or .ss-count only
   * Take care of name change option first */
  if (no_name_change != TRUE) {
    strcpy(filename, specific_ann_file);
    if (g_ann == TRUE)
      strcat(filename, ".ann");
    else
      strcat(filename, ".ss-count");
    printf("Trying to open  %s\n", filename);
    if ((ann_file = fopen(filename, "r")) == NULL) {
      printf("Could not open file: %s\n", filename);
      return -1;
    }
  } else {
    /* No name change. Use default names */
    if (g_prob_ann != TRUE) {	
      /* This is not a probability annotation. Strip off structure */
      ann_fix_name(filename);
      if (g_ann == TRUE)
	strcat(filename, ".ann");
      else
	strcat(filename, ".ss-count");
      printf("Trying to open %s\n", filename);
      if ((ann_file = fopen(filename, "r")) == NULL) {
	printf("Could not open file: %s\n", filename);
	return -1;
      }
    } else {	/* g_prob_ann is true , try it with
		 * structure on first */
      strcpy(specific_ann_file, filename); /* make a copy */
      strcat(filename, ".ann");
      printf("Trying to open %s\n", filename);
      if ((ann_file = fopen(filename, "r")) != NULL) return 0; /*Open OK*/
      printf("Could not open file: %s\n", filename);
      printf("I will try something else.\n");
      /* strip off structure and try it */
      ann_fix_name(specific_ann_file);
      strcat(specific_ann_file, ".ann");
      printf("Trying to open %s\n", specific_ann_file);
      if ((ann_file = fopen(specific_ann_file, "r")) == NULL) {
	printf("Could not open file: %s\nI give up!\n", specific_ann_file);
	return -1;
      }
    }
  }
  return 0;
}

#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
void finish_col_table_ss_int_img(int maximum_pairings, FILE* ctab) {
  int             width, height, color;
  int             i;
  int             i2;
  int             center1, center2;
  int             stopping_point;
  float           percent_single;
  gdImagePtr      image;
  if (maximum_pairings < 20) {
    stopping_point = maximum_pairings;
    width = 350;
    center1 = 350 / 2;
    center2 = 0 ; /* Added by M. Zuker - was not set before */
    height = 120 + 30 * (maximum_pairings + 1);
  } else {
    stopping_point = maximum_pairings / 2;
    width = 700;
    center1 = 350 / 2;
    center2 = 350 + 350 / 2;
    height = 120 + (30 * stopping_point + 1);
  }
  image = gdImageCreate(width, height); /* declare colors */
  set_extra_img_colors(image, FALSE);
  gdImageStringCenter(image, gdFontGiant, width/2, 30, 
		      "Color Table for SS-Count", 0);
  gdImageStringCenter(image, gdFontGiant, center1 + 12, 65, 
		      "SS-Count   %%   Hex Color", 0);
  for (i = 0; i <= stopping_point; i++) {
    if (maximum_pairings == 0)
      percent_single = 500.0;
    else
      percent_single = 100.0 * (float) (i) / (float) (maximum_pairings);
    gdImagePrintfCenter(image, gdFontGiant, center1-60, 95+30*i, 0, " %d ", i);
    gdImagePrintfCenter(image, gdFontGiant, center1, 95+30*i, 0, " %.1f", 
			percent_single);
    color = (NUM_COLORS *
	     (float) (maximum_pairings - i) / (float) maximum_pairings + .5);
    gdImageFilledRectangle(image,50,83+30*i,95,103+30*i,
			   color_table_img[color]);
    gdImageRectangle(image,95,83+30*i,300,103+30*i, color_table_img[color]);
    gdImagePrintfCenter(image, gdFontGiant, center1 + 80, 95 + 30 * i, 
			color_table_img[color], " %.6x ", color_table[color]);
  }
  if (maximum_pairings != stopping_point) {
    gdImageStringCenter(image, gdFontGiant, center2, 65, 
			" SS-Count  %%  Hex Color ", 0);
    for (i = (stopping_point + 1); i <= maximum_pairings; i++) {
      i2 = i - stopping_point - 1;
      if (maximum_pairings == 0)
	percent_single = 500.0;
      else
	percent_single = 100.0 * (float) (i) / (float) (maximum_pairings);
      gdImagePrintfCenter(image,gdFontGiant,center2-60,95+30*i, 0, " %d ", i);
      gdImagePrintfCenter(image, gdFontGiant,center2-10,95+30*i2,0, 
			  " %.1f ", percent_single);
      color = (NUM_COLORS *
	       (float) (maximum_pairings - i) / (float) maximum_pairings + .5);
      gdImageFilledRectangle(image, 400, 83 + 30*i2, 445, 103 + 30*i2, 
			     color_table_img[color]);
      gdImageRectangle(image, 445, 83 + 30*i2, 650, 103 + 30*i2, 
		       color_table_img[color]);
      gdImagePrintfCenter(image, gdFontGiant, center2 + 80, 95 + 30*i2, 
			  color_table_img[color], " %.6x ",color_table[color]);
    }
  }
  gdImageInterlace(image, 1); /* turn on interlace */
#if defined(HAVE_LIBGD) && defined(HAVE_GDIMAGEGIF)
  gdImageGif(image, ctab);
#endif
}
#endif

void finish_col_table_ann_html(int maximum_pairings, FILE * ctab) {
  int             i, color;
  int             start_color[NUM_COLORS + 1];	/* start_color[i] is the
						 * p-num at which color
						 * i starts */
  int             end_color[NUM_COLORS + 1];	/* end_color[i] is the
						 * p-num at which color
						 * i ends */
  int             i2;
  int             last_color_used;
  int             stopping_point;
  float           percent_start, percent_end;
  /* set up table of where each color starts and ends */
  last_color_used = 0;
  for (i = 0; i <= NUM_COLORS; i++) {
    start_color[i] = 0;
    end_color[i] = 0;
  }
  last_color_used = 0;
  color = 0;
  start_color[0] = 0;
  if (maximum_pairings > NUM_COLORS) {
    for (i = 2; i <= maximum_pairings; i++) {
      color = (int) (NUM_COLORS * ((float) i) / (float) (maximum_pairings));
      if (color > last_color_used) {
	start_color[color] = i;
	end_color[last_color_used] = i - 1;
	last_color_used = color;
      }
    }
    end_color[color] = maximum_pairings;
    end_color[NUM_COLORS - 1] = maximum_pairings;
  }
  fprintf(ctab, "<html><head><title>Color Table: p-num</title></head>\n");
  fprintf(ctab, "<body bgcolor=\"#dfc890\" alink=\"#00ffff\" text=\"#000000\" link=\"#905000\" vlink=\"#905000\">\n");
  fprintf(ctab, "<h1 align=\"center\">Color Table: p-num</h1>\n<table align=\"center\"border=\"8\" width=\"100%%\" ");
  fprintf(ctab, "cellspacing=\"2\" bgcolor=\"#ffffff\">\n");
  if (maximum_pairings > NUM_COLORS) {
    stopping_point = NUM_COLORS / 2 - 1;
    fprintf(ctab, "<tr><th width=\"10%%\">Color</th><th width=\"10%%\">p-num");
    fprintf(ctab, "</th><th width=\"3%%\">%%</th>");
    fprintf(ctab, "<th width=\"12%%\">Hex Color</th>\n");
    fprintf(ctab, "<th width=\"12%%\"><font color=\"#ffffff\">.</th>\n");
    fprintf(ctab, "<th width=\"10%%\">Color</th>\n");
    fprintf(ctab, "<th width=\"10%%\">p-num</th>\n");
    fprintf(ctab, "<th width=\"3%%\">%%</th>\n");
    fprintf(ctab, "<th width=\"12%%\">Hex Color</th></tr>");
    for (i = 0, i2 = stopping_point + 1; i <= stopping_point; i++, i2++) {
      fprintf(ctab, "<tr><td align=\"center\" width=\"10%%\" ");
      fprintf(ctab, "bgcolor=\"#%.6x\"><font color=\"#%.6x\">^</font></td>\n",
	      color_table[i], color_table[i]);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
      if (start_color[i] == end_color[i])
	fprintf(ctab, "%d</td>\n", start_color[i]);
      else
	fprintf(ctab, "%d-%d</td>\n", start_color[i], end_color[i]);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
      percent_start = (100.*(float) start_color[i])/((float) maximum_pairings);
      if (start_color[i] == end_color[i]) {
	fprintf(ctab, "%.1f </td>\n", percent_start);
      } else {
	percent_end = (100.*(float) end_color[i])/((float) maximum_pairings);
	fprintf(ctab, "%.1f-%.1f </td>\n", percent_start, percent_end);
      }
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
      fprintf(ctab, "<font color=\"#%.6x\">%.6x </font></td>", color_table[i], color_table[i]);
      fprintf(ctab, "<td align=\"center\" width=\"15%%\"><font color=\"#ffffff\"> . </td>");
      if (i2 <= last_color_used) {
	fprintf(ctab, "<td align=\"center\" width=\"10%%\" bgcolor=\"#%.6x\"><font color=\"#%.6x\">^ </font> </td>\n", color_table[i2], color_table[i2]);
	fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
	if (start_color[i2] == end_color[i2])
	  fprintf(ctab, "%d</td>\n", start_color[i2]);
	else
	  fprintf(ctab, "%d-%d</td>\n", start_color[i2], end_color[i2]);
	fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
	percent_start = (100. * (float) start_color[i2]) / ((float) maximum_pairings);
	if (start_color[i2] == end_color[i2])
	  fprintf(ctab, "%.1f </td>\n", percent_start);
	else {
	  percent_end = (100. * (float) end_color[i2]) / ((float) maximum_pairings);
	  fprintf(ctab, "%.1f-%.1f </td>\n", percent_start, percent_end);
	}
	fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
	fprintf(ctab, "<font color=\"#%.6x\">%.6x </font></td>", color_table[i2], color_table[i2]);
      } else {
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\"> . </td>");
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\"> . </td>");
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\"> . </td>");
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\"> . </td>");
      }
      fprintf(ctab, " </tr>\n");
    }
  } else {
    fprintf(ctab, "<tr><th width=\"10%%\">Color</th><th width=\"10%%\">p-num</th>\n");
    fprintf(ctab, "<th width=\"3%%\">%%</th><th width=\"12%%\">Hex Color</th></tr>\n");
    for (i = 0; i <= maximum_pairings; i++) {
      if (i > 1) {
	color = (int) (NUM_COLORS * ((float) (i) / (float) (maximum_pairings)));
      } else
	color = 1;
      percent_start = (100. * (float) i) / ((float) maximum_pairings);
      fprintf(ctab, "<tr>");
      fprintf(ctab, "<td align=\"center\" width=\"10%%\" bgcolor=\"#%.6x\"><font color=\"#%.6x\">.</font></td>\n", color_table[color], color_table[color]);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">%d</td>", i);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">%.1f</td>", percent_start);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
      fprintf(ctab, "<font color=\"#%.6x\">%.6x </font></td></tr>\n", color_table[color], color_table[color]);
    }
  }
  fprintf(ctab, "</table><h3>");
  fprintf(ctab, "The color of the <i>i</i><sup>th</sup> base depends on the total number of dots ");
  fprintf(ctab, "in the <i>i</i><sup>th</sup> row and <i>i</i><sup>th</sup> column of the unfiltered energy dot plot. (This is the p-num).</h3>\n");
  fprintf(ctab, "<h3>The pairing at the <i>i</i><sup>th</sup> base is more poorly determined if ");
  fprintf(ctab, "this value is large.<p>The Maximum <i>p-num</i> is %d<p>", maximum_pairings);
  fprintf(ctab, "See <a href=\"http://%s/doc/form1-doc.php#ANN\"> structure annotation:</a></h3><hr></body></html>\n", PACKAGE_URL);
  fclose(ctab);
  return;
}

void finish_col_table_ss_int_html(int maximum_pairings, FILE * ctab) {
  int             i, color = 0, color2 = 0;
  int             i2;
  int             stopping_point;
  float           percent_single, percent_single2;
  if (maximum_pairings < 20) {
    stopping_point = maximum_pairings;
  } else {
    stopping_point = maximum_pairings / 2;
  }
  fprintf(ctab, "<html><head><title>SS-Count Color Table</title></head>\n");
  fprintf(ctab, "<body bgcolor=\"#dfc890\" alink=\"#00FFFF\" text=\"#000000\" link=\"#905000\" vlink=\"#905000\">\n");
  fprintf(ctab, "<h1 align=\"center\">Color Table: ss-count</h1>\n");
  fprintf(ctab, "<table align=\"center\" BORDER=\"8\" width=\"100%%\" CELLSPACING=\"2\" bgcolor=\"#ffffff\"><tr>\n");
  if (maximum_pairings < 20) {
    fprintf(ctab, "<th width=\"10%%\">Color</th>\n");
    fprintf(ctab, "<th width=\"10%%\">ss-count</th>\n");
    fprintf(ctab, "<th width=\"5%%\" >%%</th>\n");
    fprintf(ctab, "<th width=\"10%%\">Hex Color</th>\n");
  } else {
    fprintf(ctab, "<th width=\"10%%\">Color</th>\n");
    fprintf(ctab, "<th width=\"10%%\">ss-count</th>\n");
    fprintf(ctab, "<th width=\"3%%\" >%%</th>\n");
    fprintf(ctab, "<th width=\"12%%\">Hex Color</th>\n");
    fprintf(ctab, "<th width=\"12%%\"><font color=\"#ffffff\">.</th>\n");
    fprintf(ctab, "<th width=\"10%%\">Color</th>\n");
    fprintf(ctab, "<th width=\"10%%\">ss-count</th>\n");
    fprintf(ctab, "<th width=\"3%%\">%%</th>\n");
    fprintf(ctab, "<th width=\"12%%\">Hex Color</th>\n");
  }
  fprintf(ctab, "</tr>");
  if (maximum_pairings < 20) {
    for (i = 0; i <= stopping_point; i++) {
      color = (NUM_COLORS *
	       (float) (maximum_pairings - i) / (float) maximum_pairings + .5);
      if (maximum_pairings == 0)
	percent_single = 500.0;
      else
	percent_single = 100.0 * (float) (i) / (float) (maximum_pairings);
      fprintf(ctab, "<tr>\n");
      fprintf(ctab, "<td align=\"center\" width=\"10%%\" bgcolor=\"#%.6x\"><font color=\"#%.6x\">^</font></td>\n", 
	      color_table[color], color_table[color]);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">%d</td>\n", i);
      fprintf(ctab, "<td align=\"center\" width=\"5%%\">%.1f </td>\n", percent_single);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
      fprintf(ctab, "<font color=\"#%.6x\">%.6x </font></td></tr>\n", color_table[color], color_table[color]);
    }
  } else {
    for (i = 0, i2 = stopping_point + 1; i <= stopping_point; i++, i2++) {
      color = (NUM_COLORS *
	       (float) (maximum_pairings - i) / (float) maximum_pairings + .5);
      color2 = (NUM_COLORS *
		(float) (maximum_pairings - i2) / (float) maximum_pairings + .5);
      percent_single = 100.0 * (float) (i) / (float) (maximum_pairings);
      percent_single2 = 100.0 * (float) (i2) / (float) (maximum_pairings);
      fprintf(ctab, "<tr>\n");
      fprintf(ctab, "<td align=\"center\" width=\"10%%\" bgcolor=\"#%.6x\"><font color=\"#%.6x\">^ </font> </td>\n", color_table[color], color_table[color]);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">%d</td>\n", i);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">%.1f </td>\n", percent_single);
      fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
      fprintf(ctab, "<font color=\"#%.6x\">%.6x </font></td>", color_table[color], color_table[color]);
      fprintf(ctab, "<td align=\"center\" width=\"15%%\"><font color=\"#ffffff\"> . </td>");
      if (i2 <= maximum_pairings) {
	fprintf(ctab, "<td align=\"center\" width=\"10%%\" bgcolor=\"#%.6x\"><font color=\"#%.6x\">.</font></td>\n", color_table[color2], color_table[color2]);
	fprintf(ctab, "<td align=\"center\" width=\"10%%\">%d</td>", i2);
	fprintf(ctab, "<td align=\"center\" width=\"10%%\">%.1f</td>", percent_single2);
	fprintf(ctab, "<td align=\"center\" width=\"10%%\">");
	fprintf(ctab, "<font color=\"#%.6x\">%.6x</font></td>", color_table[color2], color_table[color2]);
      } else {
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\">.</td>");
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\">.</td>");
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\">.</td>");
	fprintf(ctab, "<td align=\"center\" width=\"10%%\"><font color=\"#ffffff\">.</td>");
      }
      fprintf(ctab, " </tr>\n");
    }
  }
  fprintf(ctab, "</table><h3>");
  fprintf(ctab, "The color of the <i>i</i><sup>th</sup> base depends on the number of structures in which it is single stranded.");
  fprintf(ctab, "<p>Total Structures: %d", maximum_pairings);
  fprintf(ctab, "<p>See <a href=\"http://%s/doc/form1-doc.php#ANN\"> structure annotation:</a></h3></body></html>\n", PACKAGE_URL);
  fclose(ctab);
  return;
}

/*
 * The above link is used on html color tables generated during annotation
 */
void make_col_table_ss_int(int maximum_pairings, char *filename, 
			   char color_table_type) {
  strcpy(ctab_filename, filename);
  strcat(ctab_filename, ".col");
  if (color_table_type == 'g') {
    strcat(ctab_filename, ".gif");
  } else
    strcat(ctab_filename, ".php");
  printf("Creating color table %s\n", ctab_filename);
  if (color_table_type == 'g') {
    if ((ctab = fopen(ctab_filename, "wb")) == NULL) {
      printf("Could not open file: %s\n", ctab_filename);
      return;
    }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
    finish_col_table_ss_int_img(maximum_pairings, ctab);
#endif
  } else {
    if ((ctab = fopen(ctab_filename, "r")) != NULL) {
      printf("%s already exists, not created\n", ctab_filename);
      fclose(ctab);
      return;
    }
    if ((ctab = fopen(ctab_filename, "wt")) == NULL) {
      printf("Could not open file: %s\n", ctab_filename);
      return;
    }
    finish_col_table_ss_int_html(maximum_pairings, ctab);
  }
  return;
}

void make_col_table_ann(int maximum_pairings, char *filename, 
			char color_table_type) {
  strcpy(ctab_filename, filename);
  strcat(ctab_filename, ".col");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  if (color_table_type == 'g') {
    strcat(ctab_filename, ".gif");
  }
#endif
  if (color_table_type != 'g') {
    strcat(ctab_filename, ".php");
    if ((ctab = fopen(ctab_filename, "r")) != NULL) {
      printf("%s Already exists, not created\n", ctab_filename);
      fclose(ctab);
      return;
    }
  }
  if ((ctab = fopen(ctab_filename, "w")) == NULL) {
    printf("Could not open file: %s\n", ctab_filename);
    return;
  }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  if (color_table_type == 'g') {
    finish_col_table_ss_int_img(maximum_pairings, ctab);
    printf("color annotation file is %s\n", ctab_filename);
  } 
#endif
  if (color_table_type != 'g') {
    printf("color annotation file is %s\n", ctab_filename);
    finish_col_table_ann_html(maximum_pairings, ctab);
  }
  return;
}

int define_map_ann(char *filename, int g_ann, int g_prob_ann, int
		   no_name_change, char *specific_ann_file, int table_flag,
		   char color_table_type, int sequence_length) {
  /* This function is called for ann or ss-count */
  int        maximum_pairings, base, pairs, i, upper_limit, temp_base;
  float      color;
  char      filename2[120], rec[120];
  strcpy(filename2, filename);
  upper_limit = sequence_length;
  if (open_ann(filename2, g_ann, g_prob_ann, no_name_change, 
	       specific_ann_file) == -1)
    return -1;
  /* Read probability Annotation First */
  if (g_prob_ann == TRUE) {
    maximum_pairings = 1;
    g_bases_in_file = 0;
    while (fgets(rec, 120, ann_file) != NULL) {
      sscanf(rec, "%d%f", &base, &color);
      if ((color > 1.) || (color < 0.)) {
	printf("Error! %f is not a probability.\n", color);
	try_exit(23);
      }
      if (base > upper_limit) {
	printf("Error: Too Many bases\nThere are more bases in the annotation file than in the sequence.\n");
	printf("Sequence length: %d, annotation base was %d\n", 
	       upper_limit, base);
	try_exit(24);
      }
      g_ann_to_color[base] = get_color((double) (color), TRUE, TRUE);
      if (base > g_bases_in_file)
	g_bases_in_file = base;
    }
    if (g_bases_in_file == 0) {
      printf("Error!\tNo data in ann file\n");
      try_exit(25);
    }
  } else {
    if (g_ann == TRUE) { /* ann option: If all 0, avoid 0 division below */
      maximum_pairings = 1;	
      g_bases_in_file = 0;
      temp_base = 0;
      while (fgets(rec, 120, ann_file) != NULL) {
	temp_base++;
	sscanf(rec, "%d%f", &base, &color);
	pairs = (int) (color);
	if ((color < 1.0) && (color > 0.0)) {
	  printf("Error!\tAttempted to treat prob file as a p-num file.\n");
	  try_exit(26);
	}
	if (base > upper_limit) {
	  printf("Error!\tToo many bases.\nSequence length is %d, ",
		 upper_limit);
	  printf("but annotation base was %d.\n", base);
	  try_exit(27);
	}
	if (pairs == 1)
	  pairs = 0;	/* p-num = 0 or 1 are both equally well defined */
	g_ann_to_color[temp_base] = pairs;
	if (pairs > maximum_pairings)
	  maximum_pairings = pairs;
	g_bases_in_file++;
      }
      if (g_bases_in_file == 0) {
	printf("Error!\tNo data in ann file.\n");
	try_exit(28);
      }
      for (i = 1; i <= g_bases_in_file; i++) {
	color = NUM_COLORS*((float)g_ann_to_color[i]/(float)maximum_pairings);
	/* avoid division by zero scaled color based on the possible
	 * number of times a base can pair
	 */
	g_ann_to_color[i] = (int) color;
      }
      if (table_flag) {
	if (color_table_type == 'g') {
	  printf("The -t g option with a p-num file does work.\n");
	  printf("The code to create the .gif format for the ss_count ");
	  printf("color table will be created if necessary.\n");
	} else {
	  if (color_table_type == 'h')
	    make_col_table_ann(maximum_pairings, filename2, 'h');
	}
      }
    } else {
      if (fgets(rec, 120, ann_file) == NULL) {
	printf("Error!\tNo data in ss-count file.\n");
	try_exit(29);
      }
      sscanf(rec, "%d", &maximum_pairings);
      g_bases_in_file = 0;
      while (fgets(rec, 120, ann_file) != NULL) {
	sscanf(rec, "%d%d", &base, &pairs);
	if (base > upper_limit) {
	  printf("Error: Too Many bases\n");
	  printf("Sequence length: %d, annotation base is %d\n", 
		 upper_limit, base);
	  try_exit(30);
	}
	color = NUM_COLORS * (float) (maximum_pairings - pairs) / 
	  (float) maximum_pairings;
	g_ann_to_color[base] = (int) (color + .5);
	/* scaled color based on the possible number
	 * of way a base can be in a base pair 
	 */
	if (base > g_bases_in_file)
	  g_bases_in_file = base;
      }
      if (table_flag) {
	if (color_table_type == 'g')
	  make_col_table_ss_int(maximum_pairings, filename2, 'g');
	else {
	  if (color_table_type == 'h')
	    /*	    make_col_table_ss_int(maximum_pairings, filename2, 'h');*/
	    make_col_table_ss_int(maximum_pairings, filename2, 'g');
	}
      }
    }
  }
  fclose(ann_file);
  return g_bases_in_file;
}
