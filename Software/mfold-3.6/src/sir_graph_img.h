#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
#include "strings_img.h"
#endif

FILE           *bp_fp;

int set_color_base_img(char base1, char base2) {
  int             color;
  /* Zuker adds 2 lines below. May 23, 2003. */
  base1 = toupper(base1);
  base2 = toupper(base2);
  if (base1 == 'T')
    base1 = 'U';
  if (base2 == 'T')
    base2 = 'U';
  if (((base1 == 'G') && (base2 == 'C')) || ((base1 == 'C') && (base2 == 'G')))
    color = COLOR_CONNECTING_GC;
  else if (((base1=='A') && (base2=='U')) || ((base1=='U') && (base2=='A')))
    color = COLOR_CONNECTING_AU;
  else if (((base1=='G') && (base2=='U')) || ((base1=='U') && (base2=='G')))
    color = COLOR_CONNECTING_GU;
  else
    color = COLOR_CONNECTING_OTHERS;
  return color;
}

void finish_img_file(char *img_filename, int png_mode, int jpg_mode, 
		     char *file_type) {
  if (png_mode)
#if HAVE_LIBPNG
    gdImagePng(image, outfp)
#endif
      ;
  else if (jpg_mode)
#if HAVE_LIBJPEG
    gdImageJpeg(image, outfp, -1)
#endif
      ;
  else
#if defined(HAVE_LIBGD) && defined(HAVE_GDIMAGEGIF)
    gdImageGif(image, outfp)
#endif
      ;
  fclose(outfp);
  gdImageDestroy(image);
}

gdImagePtr create_brush(int color) {
  int             red_bk, green_bk, blue_bk;
  int             red, green, blue;
  gdImagePtr brush;
  int color_index;
  red_bk = (int) (255. * main_color_table_float[1].red + .5);
  green_bk = (int) (255. * main_color_table_float[1].green + .5);
  blue_bk = (int) (255. * main_color_table_float[1].blue + .5);
  red = (int) (255. * main_color_table_float[color].red + .5);
  green = (int) (255. * main_color_table_float[color].green + .5);
  blue = (int) (255. * main_color_table_float[color].blue + .5);
  /* create a brush for black lines */
  brush = gdImageCreate(2, 2);
  gdImageColorAllocate(brush, red_bk, green_bk, blue_bk);
  color_index = gdImageColorAllocate(brush, red, green, blue);
  gdImageLine(brush, 0, 0, 1, 0, color_index);
  gdImageLine(brush, 0, 1, 1, 1, color_index);
  return brush;
}

void img_create_colors(int annotation) {
  int             red, green, blue;
  int             j;
  for (j = 0; j < MAIN_COLORS; j++) {
    img_brush[j] = create_brush(j);
  }
  for (j = 1; j < MAIN_COLORS; j++) {
    red = (int) (255. * main_color_table_float[j].red + .5);
    green = (int) (255. * main_color_table_float[j].green + .5);
    blue = (int) (255. * main_color_table_float[j].blue + .5);
    img_color[j] = gdImageColorAllocate(image, red, green, blue);
  }
  if (annotation) {
    if (annotation == PROB)
      set_extra_img_colors(image, TRUE);
    else
      set_extra_img_colors(image, FALSE);
  }
}

void initialize_img(int img_width, int img_height, int img_interlaced, 
		    char *first_line, int annotation) {
  int             center;
  int             vert_location;
  image = gdImageCreate(img_width, img_height);
  img_create_colors(annotation);
  /* interleaved jpg is possible, but not standard */
  if (img_interlaced)	/* perhaps interlaced jpg should be off here */
    gdImageInterlace(image, 1); /* turn on interlace */
  center = img_width / 2 ;
  if (img_width > 100) {
    if (img_height < 500) {
      vert_location = img_height - 10;
      gdImagePrintfCenter(image, gdFontSmall, center, vert_location, 
			  img_color[3], " %s", first_line);
    } else {
      vert_location = img_height - 12;
      gdImagePrintfCenter(image, gdFontGiant, center, vert_location, 
			  img_color[3], " %s", first_line);
    }
  }
}

void adjust_structure_coordinates_img(int img_width,int img_height,
				      int start_base,int end_base,
				      float *loop_center_x,
				      float *loop_center_y,
				      float *loop_center_x_img, 
				      float *loop_center_y_img, 
				      int total_labels, 
				      float *structure_label_x_img,
				      float *structure_label_y_img, 
				      float *structure_label_x,
				      float *structure_label_y,
				      int total_loops, float *base_x_img, 
				      float *base_y_img, float *base_x, 
				      float *base_y, int external_domain_flag,
				      int ex_start_base, int ex_end_base, 
				      float x_shift, float y_shift, 
				      float zoom_factor) {
  
  int             i;
  float           scale;
  float           extra_x_shift;
  float           extra_y_shift;
  scale = img_height / (HT + 72.0);
  x_shift = - x_shift;
  y_shift = - y_shift;
  /* zoom factor <1 makes the structure bigger
   * zoom factor >1 makes the structure smaller */
  zoom_factor = 1. / zoom_factor;
  scale = scale * zoom_factor;
  if (zoom_factor != 1.0) {
    extra_x_shift = img_width / 2. * (zoom_factor - 1.0);
    extra_y_shift = -1.0 * img_height / 2. * (zoom_factor - 1.0);
  } else {
    extra_x_shift = 0.0;
    extra_y_shift = 0.0;
  }
  for (i = start_base; i <= end_base; i++) {
    if (external_domain_flag) {
      if ((i > ex_start_base) && (i < ex_end_base))
	i = ex_end_base;
    }
    base_y_img[i] = img_height - (base_y[i]-y_shift)*scale - extra_y_shift;
    base_x_img[i] = (base_x[i] - x_shift) * scale - extra_x_shift;
  }
  for (i = 1; i <= total_loops; i++) {
    loop_center_y_img[i] =
      img_height - (loop_center_y[i] - y_shift) * scale - extra_y_shift;
    loop_center_x_img[i] = (loop_center_x[i]-x_shift)*scale - extra_x_shift;
  }
  for (i = 1; i <= total_labels; i++) {
    structure_label_y_img[i] = img_height -
      (structure_label_y[i] - y_shift) * scale - extra_y_shift;
    structure_label_x_img[i] = (structure_label_x[i] - x_shift)*scale 
      - extra_x_shift;
  }
}

void display_bases_img(float *structure_label_x_img,
		       float *structure_label_y_img, float scale,
		       int img_width, int img_height, int outline_mode,
		       int total_labels, int *connected, float *base_x_img,
		       float *base_y_img, float font_size, int start_base,
		       int end_base, int *structure_label_value, int lines,
		       char *bases, int loop_labels, int total_loops,
		       float *loop_center_x_img, float *loop_center_y_img,
		       int external_domain_flag, int ex_start_base,
		       int ex_end_base, int history_offset, int annotation,
		       int dots_flag, int bases_flag, int *ann_to_color,
		       float zoom_factor, int is_circular, int *forced_bases,
		       int forced_bases_count, int make_bp_file) {
  int             imod; /* M. Zuker June 7, 2007. imod = i + 1 mod
			* sequence length (0 not used). Introduced for
			* correct plotting of circular molecules
			*/
  int             last_base; /* end_base normally and end_base + 1 when
			      * the molecule is circular */
  int             i, base_1; 
  float           x_shift;
  float           y_shift;
  float           dist;
  float           img_font_size;
  float           angle;
  float           x1, x2, y1, y2, x3, y3, x4, y4;
  float           alpha;
  float           dif_x;
  float           dif_y;
  float           img_scale;
  float           cir_top_edge, cir_bot_edge, cir_left_edge, cir_right_edge;
  gdFontPtr       font;
  int             line_color;
  int             last_line_color = 0;
  int             color, radius;
  int             connected_to;
  int             first_loop;
  int             img_outline_mode;
  int             border = 15;
  int             bp_closeness;
  char            label_content[10];
  float           label_content_offset;
  float           distance_bases;
  
  /* M. Zuker, June 7, 2007. Set last_base */
  if ( is_circular ) {
    last_base = end_base + 1;
  } else {
    last_base = end_base;
  }
  img_scale = scale * img_height / (HT + 72.0);
  img_font_size = font_size;
  img_font_size = font_size * img_height / (HT + 0.0) * .75 * zoom_factor;
  img_outline_mode = FALSE;
  scale = scale * zoom_factor;
  distance_bases = DISTANCE_BASES;
  if (outline_mode)
    img_outline_mode = TRUE;
  else {
    if ((img_font_size < 1.85) && (!annotation))
      img_outline_mode = TRUE;
  }
  /* for image map user must click within 4 character widths + 5 pixels 
   * of the center of the base pair dot */
  if (make_bp_file) {	
    bp_closeness = img_font_size * 4 + 5;
    if (bp_closeness > 72)
      bp_closeness = 72;
    if (bp_closeness < 5)
      bp_closeness = 5;
    bp_closeness = bp_closeness * -1;
    fprintf(bp_fp, "%d \n", bp_closeness);
  }
  if (img_font_size < 1.85) {
    x_shift = 2.0;
    y_shift = 4.0;
    font = gdFontTiny;
  } else if (img_font_size <= 6) {
    font = gdFontTiny;
    x_shift = 2.5;
    y_shift = 4;
  } else if (img_font_size < 7) {
    font = gdFontSmall;
    x_shift = 3;
    y_shift = 6;
  } else if (img_font_size < 8) {
    font = gdFontSmall;
    x_shift = 3;
    y_shift = 6;
  } else if (img_font_size < 9) {
    font = gdFontMediumBold;
    x_shift = 4;
    y_shift = 6;
  } else if (img_font_size < 10) {
    font = gdFontMediumBold;
    x_shift = 4;
    y_shift = 7;
  } else if (img_font_size < 12) {
    font = gdFontLarge;
    x_shift = 5;
    y_shift = 10;
  } else {
    font = gdFontGiant;
    x_shift = 5;
    y_shift = 10;
  }
  if ( (connected[start_base] >= 0) && 
       (connected[start_base] == (connected[start_base - 1] - 1) ) &&
       (connected[start_base + 1] == connected[start_base - 1] - 1) ) {
    first_loop = 2;
  } else
    first_loop = 1;
  /*  draw labels */
  for (i = 1; i <= total_labels; i++) {
    if (!is_circular) {
      if (g_next[g_structure_label_value[i]]==0) { /* 2003-06-21 N Markham */
	strcpy(label_content, "3'");
	label_content_offset = 1.5;
      } else if (g_prev[g_structure_label_value[i]] == 0) { /* 03-06-21 NM */
	strcpy(label_content, "5'");
	label_content_offset = 1.5;
      } else {
	sprintf(label_content, "%d",
		g_history[structure_label_value[i]]); /* 03-06-20 NM */
	label_content_offset =
	  log(structure_label_value[i] + history_offset) / log(10);
      }
    } else {
      sprintf(label_content, "%d",
	      g_history[structure_label_value[i]]); /* 03-06-20 Nick Markham */
      label_content_offset =
	log(structure_label_value[i] + history_offset) / log(10);
    }
    gdImageString(image, font, (int) (structure_label_x_img[i] -
				       x_shift * label_content_offset),
		   (int) (structure_label_y_img[i] - y_shift),
		  (unsigned char*) label_content, img_color[COLOR_LABEL_LINE]);
    x1 = base_x_img[structure_label_value[i]] + 
      (structure_label_x_img[i] - base_x_img[structure_label_value[i]]) / 6;
    y1 = base_y_img[structure_label_value[i]] + 
      (structure_label_y_img[i] - base_y_img[structure_label_value[i]]) / 6;
    x2 = base_x_img[structure_label_value[i]] + 5.*
      (structure_label_x_img[i] - base_x_img[structure_label_value[i]]) / 8;
    y2 = base_y_img[structure_label_value[i]] + 5.*
      (structure_label_y_img[i] - base_y_img[structure_label_value[i]]) / 8;
    gdImageLine(image, x1, y1, x2, y2, img_color[COLOR_LABEL_LINE]);
  }
  /* make connections between bases */
  if (img_font_size > 7.5) {	/* double width line */
    gdImageSetBrush(image, img_brush[COLOR_BASES]);
    line_color = gdBrushed;
  } else
    line_color = img_color[COLOR_BASES];
  for (i = (start_base + 1); i <= last_base; i++) { 
    if (external_domain_flag) {
      if ((i > ex_start_base) && (i <= ex_end_base))
	i = ex_end_base + 1;
    }
    if (i==end_base+1) {
      imod = start_base;
    } else {
      imod = i;
    }
    dif_x = base_x_img[imod] - base_x_img[i - 1];
    dif_y = base_y_img[imod] - base_y_img[i - 1];
    dist = sqrt(dif_x * dif_x + dif_y * dif_y);
    if ( (dist > 5.0*distance_bases*img_scale/8.0) || img_outline_mode ) {
      alpha = atan2(dif_y, dif_x);
      if (img_outline_mode) {
	x1 = base_x_img[i - 1];
	x2 = base_x_img[imod];
	y1 = base_y_img[i - 1];
	y2 = base_y_img[imod];
      } else { 
	x1 = base_x_img[i - 1] + 0.375*distance_bases*img_scale*
	  zoom_factor*cos(alpha);
	x2 = base_x_img[i - 1] + (dist - 0.375*distance_bases * img_scale 
				  * zoom_factor) * cos(alpha);
	y1 = base_y_img[i - 1] + 0.375*distance_bases*
	  img_scale*zoom_factor*sin(alpha);
	y2 = base_y_img[i - 1] + (dist - 0.375*distance_bases * img_scale *
				  zoom_factor) * sin(alpha);
      }
      if (g_next[i - 1] || g_prev[imod]) /* 2003-06-23 Nick Markham */
	gdImageLine(image, x1, y1, x2, y2, line_color);
    }
  } 
  /* draw base_pair lines or dots */
  radius = (int) (img_font_size + .5);
  if (radius < 1)
    radius = 1;
  /* do not fill circles outside or on edge of image
   * radius should have been named diameter */
  cir_top_edge = radius / 2 + 2;
  cir_bot_edge = img_height - (2 + radius / 2);
  cir_left_edge = radius / 2 + 2;
  cir_right_edge = img_width - (2 + radius / 2);
  if (lines) {
    if (img_font_size > 6.5)	/* switch to double width line */
      last_line_color = 0;
    line_color = gdBrushed;
  }
  for (i = start_base; i <= end_base; i++) {
    if (external_domain_flag) {
      if ((i > ex_start_base) && (i < ex_end_base))
	i = ex_end_base;
    }
    connected_to = connected[i];
    if (connected_to > i) { /* connections between base pairs */
      color = set_color_base_img(bases[i], bases[connected_to]);
      if (lines) {
	if (img_font_size > 6.5) {
	  if (color != last_line_color) {
	    gdImageSetBrush(image, img_brush[color]);
	    last_line_color = color;
	  }
	} else
	  line_color = img_color[color];
	if (img_outline_mode) {
	  x1 = base_x_img[i];
	  x2 = base_x_img[connected_to];
	  y1 = base_y_img[i];
	  y2 = base_y_img[connected_to];
	} else {
	  x1 = base_x_img[i] + (base_x_img[connected_to] - base_x_img[i])*.25;
	  x2 = base_x_img[i] + (base_x_img[connected_to] - base_x_img[i])*.75;
	  y1 = base_y_img[i] + (base_y_img[connected_to] - base_y_img[i])*.25;
	  y2 = base_y_img[i] + (base_y_img[connected_to] - base_y_img[i])*.75;
	}
	gdImageLine(image, x1, y1, x2, y2, line_color);
	if (make_bp_file) {
	  if (((y1 > cir_top_edge) && (y1 < cir_bot_edge)) &&
	      ((x1 > cir_left_edge) && (x1 < cir_right_edge)))
	    fprintf(bp_fp, "%d %d %d %d\n", (int) ((x1 + x2) / 2.), 
		    (int) ((y1 + y2) / 2.), i, connected_to);
	}
      } else { /* draw dot for base pair */
	x1 = (base_x_img[i] + base_x_img[connected_to]) / 2.0;
	y1 = (base_y_img[i] + base_y_img[connected_to]) / 2.0;
	gdImageArc(image, x1, y1, radius, radius, 0, 360, img_color[color]);
	/* keep in bounds */
	if (((y1 > cir_top_edge) && (y1 < cir_bot_edge)) &&
	    ((x1 > cir_left_edge) && (x1 < cir_right_edge))) {
	  if (radius > 1)
	    gdImageFill(image, x1, y1, img_color[color]);
	  if (make_bp_file)
	    fprintf(bp_fp, "%d %d %d %d\n", (int) (x1), (int) (y1),
		    i, connected_to);
	}
      }
    }
  }
  if (forced_bases_count > 0) {	/* make diamonds */
    for (i = 0; i < forced_bases_count; i += 2) {
      gdPoint points[4];
      base_1 = forced_bases[i];
      connected_to = connected[base_1];
      x1 = base_x_img[base_1] + 
	(base_x_img[connected_to] - base_x_img[base_1])/4.;
      x2 = base_x_img[base_1] + 
	3.*(base_x_img[connected_to] - base_x_img[base_1])/4.;
      y1 = base_y_img[base_1] + 
	3.*(base_y_img[connected_to] - base_y_img[base_1])/4.;
      y2 = base_y_img[base_1] + 
	3.*(base_y_img[connected_to] - base_y_img[base_1])/4.;
      angle = atan2(y2 - y1, x2 - x1);
      dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2.;
      dist = sqrt(dist * dist + dist * dist);
      x3 = x2 - dist * sin(angle + PI / 4);
      y3 = y2 + dist * cos(angle + PI / 4);
      x4 = x2 - dist * sin(PI / 4 - angle);
      y4 = y2 - dist * cos(PI / 4 - angle);
      color = set_color_base_img(bases[base_1], bases[connected_to]);
      points[0].x = x1; points[0].y = y1;
      points[1].x = x2; points[1].y = y2;
      points[2].x = x3; points[2].y = y3;
      points[3].x = x4; points[3].y = y4;
      gdImageFilledPolygon(image, points, 4, img_color[color]);
    }
  }
  /* draw bases */
  if ((annotation) && (dots_flag)) {
    cir_top_edge = radius + 2;
    cir_bot_edge = img_height - (2 + radius);
    cir_left_edge = radius + 2;
    cir_right_edge = img_width - (2 + radius);
  }
  if (!img_outline_mode) {
    for (i = start_base; i <= end_base; i++) {
      if (external_domain_flag) {
	if ((i > ex_start_base) && (i < ex_end_base))
	  i = ex_end_base;
      }
      if (!annotation)
	gdImageChar(image, font, (int) (base_x_img[i] - x_shift), 
		    (int) (base_y_img[i] - y_shift),  bases[i],
		    img_color[COLOR_BASES]);
      if ((annotation) && (dots_flag)) {
	gdImageArc(image, (int) (base_x_img[i] - .5), (int) base_y_img[i], 
		   radius * 2, radius * 2, 0, 360, 
		   color_table_img[ann_to_color[i]]);
	if (((base_y_img[i] > cir_top_edge) && (base_y_img[i]<cir_bot_edge)) 
	    && ((base_x_img[i] > cir_left_edge) 
		&& (base_x_img[i] < cir_right_edge)))
	  gdImageFillToBorder(image, (int) (base_x_img[i] - .5), 
			      (int) base_y_img[i], 
			      color_table_img[ann_to_color[i]], 
			      color_table_img[ann_to_color[i]]);
      }
      if ((annotation) && (bases_flag)) {
	if (!dots_flag) {
	  gdImageChar(image, font, (int) (base_x_img[i] - x_shift), 
		      (int) (base_y_img[i] - y_shift), bases[i],
		      color_table_img[ann_to_color[i]]); 
	  /* M. Zuker fixes bug. ann_to_color[i] changed to 
	     color_table_img[ann_to_color[i]] */
	} else if (annotation == PROB)
	  border = LOG_WHITE_BLACK_SWITCH;
	else
	  border = WHITE_BLACK_SWITCH;
	if (g_ann_to_color[i] > border) {
	  gdImageChar(image, font, (int) (base_x_img[i] - x_shift), 
		      (int) (base_y_img[i] - y_shift), bases[i],
		      img_color[COLOR_ANN_WHITE_LETTER]);
	} else {
	  gdImageChar(image, font, (int) (base_x_img[i] - x_shift), 
		      (int) (base_y_img[i] - y_shift), bases[i],
		      img_color[COLOR_ANN_BLACK_LETTER]);
	}
      }
    }
  }
  /* loop labels */
  if (loop_labels) {
    for (i = first_loop; i <= total_loops; i++) {
      gdImagePrintf(image, font, 
		    (int) (loop_center_x_img[i] - x_shift * log(i) / log(10)),
		    (int) (loop_center_y_img[i] - y_shift),
		    img_color[COLOR_LABEL_LINE], "%d", i);
    }
  }
}

int             open_output(char *out_filename);
int open_bp_file(char *out_filename) {	/* open specified output file */
  printf("Output: %s\n", out_filename);
  if ((bp_fp = fopen(out_filename, "w")) == NULL) {
    printf("Error! Could not open file base pair file: %s\n", out_filename);
    try_exit(46);
  }
  return (0);
}

void make_structure_img(char *structure_filename,int img_width, int img_height,
			int img_interlaced, char *first_line, int start_base,
			int end_base, float *loop_center_x, 
			float *loop_center_y, int total_labels,
			float *structure_label_x, float *structure_label_y, 
			float *base_x, float *base_y, float scale,
			int outline_mode, int *connected, float font_size, 
			int *structure_label_value, int lines, char *bases, 
			int loop_labels, int total_loops, 
			int external_domain_flag, int ex_start_base, 
			int ex_end_base, int history_offset, int annotation, 
			int *ann_to_color, int dots_flag, int bases_flag, 
			int length, float x_shift, float y_shift, 
			float zoom_factor, int is_circular, int *forced_bases,
			int forced_bases_count, int png_mode,
			int jpg_mode, int make_bp_file) {
  float          *structure_label_x_img;
  float          *structure_label_y_img;
  float          *base_y_img;
  float          *base_x_img;
  float          *loop_center_x_img;
  float          *loop_center_y_img;
  int             size;
  char            file_type[6];
  char            structure_img_filename[120];
  char            bp_filename[120];
  char            action[250];
  int             error;
  size = length + 1;
  if (jpg_mode)		/* interlaced messed up viewing for jpg */
    img_interlaced = FALSE;
  /* allocate memory */
  structure_label_x_img = (float *) malloc((size) * sizeof(float));
  if (structure_label_x_img == NULL)
    report_mem_error("structure_label_x_img");
  structure_label_y_img = (float *) malloc((size) * sizeof(float));
  if (structure_label_y_img == NULL)
    report_mem_error("structure_label_y_img");
  base_y_img = (float *) malloc((size) * sizeof(float));
  if (base_y_img == NULL)
    report_mem_error("base_y_img");
  base_x_img = (float *) malloc((size) * sizeof(float));
  if (base_x_img == NULL)
    report_mem_error("base_x_img");
  loop_center_x_img = (float *) malloc((size) * sizeof(float));
  if (loop_center_x_img == NULL)
    report_mem_error("loop_center_x_img");
  loop_center_y_img = (float *) malloc((size) * sizeof(float));
  if (loop_center_y_img == NULL)
    report_mem_error("loop_center_y_img");
  strcpy(structure_img_filename, structure_filename);
  /* chop .gif, .jpg or .png off the end */
  chop_suffix(structure_img_filename, ".img");
  if (png_mode) {
    strcpy(file_type, "png");
  } else if (jpg_mode) {
    strcpy(file_type, "jpg");
  } else {
    strcpy(file_type, "gif");
  }
  if (make_bp_file) {
    strcpy(bp_filename, structure_img_filename);
    if ((png_mode) || (jpg_mode))
      strcat(bp_filename, ".png2bp");
    else
      strcat(bp_filename, ".gif2bp");
    open_bp_file(bp_filename);
  }
  strcat(structure_img_filename, ".");
  strcat(structure_img_filename, file_type);
  open_output(structure_img_filename);
  initialize_img(img_width, img_height, img_interlaced, first_line, 
		 annotation);
  display_header_message_img(image, img_width, annotation);
  adjust_structure_coordinates_img(img_width, img_height, start_base, end_base,
				   loop_center_x, loop_center_y,
				   loop_center_x_img, loop_center_y_img,
				   total_labels, structure_label_x_img, 
				   structure_label_y_img, structure_label_x,
				   structure_label_y, total_loops, base_x_img,
				   base_y_img, base_x, base_y,
				   external_domain_flag, ex_start_base, 
				   ex_end_base, x_shift, y_shift, zoom_factor);
  display_bases_img(structure_label_x_img, structure_label_y_img, scale,
		    img_width, img_height, outline_mode, total_labels,
		    connected, base_x_img, base_y_img, font_size, start_base,
		    end_base, structure_label_value, lines, bases, loop_labels,
		    total_loops, loop_center_x_img, loop_center_y_img,
		    external_domain_flag, ex_start_base, ex_end_base, 
		    history_offset, annotation, dots_flag, bases_flag,
		    ann_to_color, 1./zoom_factor, is_circular, forced_bases,
		    forced_bases_count, make_bp_file);
  if (make_bp_file) {
    fclose(bp_fp);	/* close file */
    strcpy(action, "sort -n -o ");
    strcat(action, bp_filename);	/* sort the file */
    strcat(action, " ");
    strcat(action, bp_filename);
    strcat(action, "\n");
    error = system(action);
    if (error != 0)
      printf("Error sorting %s from sir_graph action was %s\n",
	     bp_filename, action);
  }
  printf("%s File: %s\n", file_type, structure_img_filename);
  printf("%s dimensions %d x %d\n", file_type, img_width, img_height);
  finish_img_file(structure_img_filename, png_mode, jpg_mode, file_type);
  free(structure_label_x_img);
  free(structure_label_y_img);
  free(base_y_img);
  free(base_x_img);
  free(loop_center_x_img);
  free(loop_center_y_img);
}
