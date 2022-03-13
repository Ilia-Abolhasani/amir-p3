void set_color_ps(int col) {
  fprintf(outfp, " c_%c\n", main_color_table_let[col]);
}

void set_color_ann_ps(int col) {
  fprintf(outfp, " c_%d\n", col);
}

void extra_ps_char(void) {
  /* For ISO Latin encoding beyond \177 Octal */
  fprintf(outfp, "\n/ISOLatin1Encoding where {pop save true}{false} ifelse");
  fprintf(outfp, "\n/reencodeISO {");
  fprintf(outfp, "\n   dup length dict begin");
  fprintf(outfp, "\n        {1 index /FID ne {def}{pop pop} ifelse} forall");
  fprintf(outfp, "\n         /Encoding ISOLatin1Encoding def");
  fprintf(outfp, "\n        currentdict");
  fprintf(outfp, "\n    end");
  fprintf(outfp, "\n} def");
  fprintf(outfp, "\n/findISO {");
  fprintf(outfp, "\n    dup /FontType known {");
  fprintf(outfp, "\n        dup /FontType get 3 ne {");
  fprintf(outfp, "\n              dup /CharStrings known {");
  fprintf(outfp, "\n                 dup /CharStrings get /Thorn known {");
  fprintf(outfp, "\n                     true");
  fprintf(outfp, "\n                }{ false } ifelse");
  fprintf(outfp, "\n            }{ false } ifelse");
  fprintf(outfp, "\n       }{ false } ifelse");
  fprintf(outfp, "\n     }{ false } ifelse");
  fprintf(outfp, "\n} def");
  fprintf(outfp, "\n");
}


void define_colors_ps(int smart_colors, int midpoints, int annotation) {
  int             i;
  fprintf(outfp, "%s c_rad far below defines radius of base pair dots\n", "%%");
  if (annotation)
    fprintf(outfp, 
	    "%s c_drad far below defines radius of annotation dots\n", "%%");
  fprintf(outfp, "%s Define colors for postscript \n", "%%");
  fprintf(outfp, "%s c_k is background, c_t is for general text\n", "%%");
  fprintf(outfp, "%s c_b is bases, c_m is for main circle\n", "%%");
  if (smart_colors) {
    fprintf(outfp, "%s c_g is for connecting arcs for GC, CG \n", "%%");
    fprintf(outfp, "%s c_a is for connecting arcs for  AU,UA \n", "%%");
    fprintf(outfp, "%s c_u is for connecting arcs for  GU,UG \n", "%%");
    fprintf(outfp, "%s c_o is for connecting arcs for other \n", "%%");
  } else {
    fprintf(outfp, 
	    "%s c_g is for connecting lines or dots for GC, CG \n", "%%");
    fprintf(outfp, 
	    "%s c_a is for connecting lines or dots for  AU,UA \n", "%%");
    fprintf(outfp, "%s c_u is for connecting   GU,UG \n", "%%");
    fprintf(outfp, "%s c_o is for connecting arcs for other \n", "%%");
    fprintf(outfp, 
	    "%s in circle graph, c_g may be used for all types \n", "%%");
    fprintf(outfp, "%s of connecting arcs\n", "%%");
  }
  if (midpoints) {
    fprintf(outfp, "%s c_s for centers of arcs and connections \n", "%%");
  }
  fprintf(outfp, "%s c_c for header \n", "%%");
  fprintf(outfp, "%s c_n for number of label \n", "%%");
  fprintf(outfp, "%s c_l for line of label \n", "%%");
  fprintf(outfp, "%s Colors originally set in sir_graph_color.h\n", "%%");
  for (i = 0; i < MAIN_COLORS; i++) {
    fprintf(outfp, "/c_%c { %.3f %.3f %.3f setrgbcolor} def\n", 
	    main_color_table_let[i], main_color_table_float[i].red,
	    main_color_table_float[i].green, main_color_table_float[i].blue);
  }
  extra_ps_char();
  /* define extra colors for ps */
  if (annotation == PROB) {
    set_extra_ps_colors(TRUE);
    for (i = 0; i <= LOG_NUM_COLORS; i++)
      fprintf(outfp, "/c_%d { %.3f %.3f %.3f setrgbcolor} def\n",
	      i, color_table_ps[i].red, color_table_ps[i].green,
	      color_table_ps[i].blue);
  } else {
    if (annotation != NONE) {
      set_extra_ps_colors(FALSE);
      for (i = 0; i <= NUM_COLORS; i++) {
	fprintf(outfp, "/c_%d { %.3f %.3f %.3f setrgbcolor} def\n",
		i, color_table_ps[i].red, color_table_ps[i].green,
		color_table_ps[i].blue);
      }
    }
  }
}


void display_header_message(void){ 
  char time_data[50];
  int cprt = 169 ;
  int length;
  float hor_pos, ver_offset;
  time_t now;
  now = time(NULL);
  strcpy(time_data,"Created ") ;
  strcat(time_data,ctime(&now));
  length = strlen(time_data);
  time_data[length-1] = '\0' ;
  ver_offset = HT - 6.0;  
  hor_pos = WD + 72.0 - 9.0*(float)length/2.0 - 15 ; /* 9.0 is font size */
  /* Try not to make text overwrite image below */
  fprintf(outfp,"/sf 9 def\n");
  fprintf(outfp,"/Helvetica findfont findISO { reencodeISO /Symbol-ISO exch definefont \n");
  fprintf(outfp,"               sf scalefont setfont }\n");
  fprintf(outfp,"               { sf scalefont setfont } ifelse\n");
  fprintf(outfp,"15 %f moveto\n", ver_offset + 66.);
  fprintf(outfp," (Output of %s %c) show\n", PROGRAM_NAME, cprt);
  fprintf(outfp,"15 %f moveto\n", ver_offset + 55);
  fprintf(outfp," (%s) show\n", PACKAGE_STRING);
  fprintf(outfp,"%f %f moveto", hor_pos, ver_offset + 66);
  fprintf(outfp," (%s) show\n", time_data);
}

void finish_ps_file(void) {
  fprintf(outfp, "showpage\n");
  fprintf(outfp, "%sEOF\n", "%%");
  fclose(outfp);
}

void initialize_ps(int smart_colors, int annotation,
		   float x_shift, float y_shift, float zoom_factor) {
  time_t now;
  now = time(NULL);
  fprintf(outfp, "%c!PS-Adobe-3.0 EPSF-2.0 \n", '%');
  fprintf(outfp, "%sBoundingBox: 0 0 %.0f %.0f\n", "%%", WD + 72.0, HT + 72.0);
  fprintf(outfp, " 1.0 1.0 scale\n");
  fprintf(outfp, "%s Created from ct or ss file by sir_graph.c, %s\n", "%%",PACKAGE_STRING);
  fprintf(outfp,"\n%s %s","%%",ctime(&now));
  fprintf(outfp, "%s\n%s\n", "%%", "%%");
  fprintf(outfp, " %.4f %.4f translate\n", - x_shift/zoom_factor + 
	  (WD + 72.0)/2.0, - y_shift/zoom_factor + (HT + 72.0)/2.0);
  fprintf(outfp, "%s     This translate,scale controls zooming \n", "%%");
  fprintf(outfp, " %.5f %.5f scale\n", 1 / zoom_factor, 1 / zoom_factor);
  fprintf(outfp, "%s\n", "%%");
  fprintf(outfp, "%s\n", "%%");
  fprintf(outfp, "%s defaults: 0 0 translate, 1 1 scale for full image\n","%%");
  fprintf(outfp, "%s m n translate, sc sc scale \n", "%%");
  fprintf(outfp, "%s is used to zoom about point (x,y) with scale sc \n", "%%");
  fprintf(outfp, "%s m= -1*(x*scale-%.0f) \n", "%%", (WD + 72.0)/2.0);
  fprintf(outfp, "%s n= -1*(y*scale-%.0f) \n", "%%",  (HT + 72.0)/2.0);
  fprintf(outfp, "%s for (x,y) 0,0 is bottom left, 0<x<%.0f, 0<y<%.0f\n",
	  "%%", WD + 72.0, HT + 72.0);
  fprintf(outfp, "%s\n", "%%");
  fprintf(outfp, "%s\n", "%%");
  fprintf(outfp, "%s To change image size for some viewers ....\n", "%%");
  fprintf(outfp, "%s Either alter Bound ngBox\n", "%%");
  fprintf(outfp, "%s or 1.0 1.0 scale\n", "%%");
  fprintf(outfp, "%s Valid scales are .05 to 1.00\n", "%%");
  fprintf(outfp, "%s ONLY scale and Bound ngBox lines above may be changed.\n", "%%");
  fprintf(outfp, "%s The scale for change of size is\n", "%%");
  fprintf(outfp, "%s before the translate,scale for zooming", "%%");
  fprintf(outfp, "%s\n", "%%");
  fprintf(outfp, "%s\n", "%%");
  fprintf(outfp, "/m {moveto} def\n");
  fprintf(outfp, "/l {lineto} def\n");
  fprintf(outfp, "/s {stroke} def\n");
  fprintf(outfp, "/a {arc} def\n");
  define_colors_ps(smart_colors, FALSE, annotation);
  set_color_ps(COLOR_BACKGROUND);
  fprintf(outfp, "%s Set background color\n", "%%");
  fprintf(outfp, "0 0 m\n");
  fprintf(outfp, "%.0f 0 rlineto\n", WD + 72.0);
  fprintf(outfp, "0 %.0f rlineto\n", HT + 72.0);
  fprintf(outfp, "%.0f 0 rlineto\n", -WD - 72.0);
  fprintf(outfp, "0 %.0f rlineto\n", -HT - 72.0);
  fprintf(outfp, "closepath\n");
  fprintf(outfp, "fill\n");
  set_color_ps(COLOR_COPYRIGHT);
  display_header_message();
  set_color_ps(COLOR_TEXT);
}

void set_color_base_ps(char base1, char base2) {
  /* Zuker adds 2 lines below. May 23, 2003. */
  base1 = toupper(base1);
  base2 = toupper(base2);
  if (base1 == 'T')
    base1 = 'U';
  if (base2 == 'T')
    base2 = 'U';
  if (((base1=='G') && (base2=='C')) || ((base1 == 'C') && (base2 == 'G')))
    set_color_ps(COLOR_CONNECTING_GC);
  else if (((base1=='A') && (base2=='U')) || ((base1 == 'U') && (base2=='A')))
    set_color_ps(COLOR_CONNECTING_AU);
  else if (((base1=='G') && (base2=='U')) || ((base1=='U') && (base2 == 'G')))
    set_color_ps(COLOR_CONNECTING_GU);
  else
    set_color_ps(COLOR_CONNECTING_OTHERS);
}

void display_bases_ps(char *first_line,float line_width, int *connected, 
		      int start_base, int end_base, int total_loops,
		      float *loop_center_x, float *loop_center_y,
		      float font_size, int outline_mode, float scale, 
		      char *bases, float *base_x, float *base_y, int lines, 
		      int total_labels, float *structure_label_x, 
		      float *structure_label_y, int *structure_label_value,
		      int loop_labels, int external_domain_flag,
		      int ex_start_base, int ex_end_base, int history_offset,
		      int annotation, int *ann_to_color, int dots_flag, 
		      int bases_flag, int is_circular, int *forced_pairs,
		      int number_of_forced_pairs) {
  int             imod; /* M. Zuker June 7, 2007. imod = i + 1 mod
			 * sequence length (0 not used). Introduced for
			 * correct plotting of circular molecules
			 */
  int             last_base; /* end_base normally and end_base + 1 when
			      * the molecule is circular */
  int             i;
  int             len;
  float           x_shift;
  float           y_shift;
  float           start;
  float           dist;
  float           x1, x2, y1, y2, x3, y3, x4, y4;
  int             base_1;
  float           angle;
  float           alpha;
  float           dif_x;
  float           dif_y;
  float           ps_font_size;
  char            label_content[10];
  float           label_content_offset;
  int             connected_to;
  int             first_loop = 1;
  float           line_length_adjust;
  float           distance_bases;

  /* M. Zuker, June 7, 2007. Set last_base */
  if ( is_circular ) {
    last_base = end_base + 1;
  } else {
    last_base = end_base;
  }
  distance_bases = DISTANCE_BASES;
  fprintf(outfp, "/Helvetica-Oblique findfont\n");
  fprintf(outfp, "18 scalefont\n");
  fprintf(outfp, "setfont\n");
  len = strlen(first_line);
  start = (WD + 72.0)/2.0 - len * 9.0 / 2.0;
  if (start < 5.)
    start = 5.;
  fprintf(outfp, "%f 15 m\n", start);
  fprintf(outfp, " (%s) show\n", first_line);
  /* display bases */
  fprintf(outfp, " %f setlinewidth\n", line_width);
  /* display loops */
  set_color_ps(COLOR_MAIN_CIRCLE);
  fprintf(outfp, "s\n");
  connected[0] = 0;
  if (loop_labels) {
    if ( (connected[start_base] >= 0) && 
	 (connected[start_base] == (connected[start_base - 1] - 1)) &&
	 (connected[start_base + 1] == connected[start_base - 1] - 1)) {
      first_loop = 2;
    } else
      first_loop = 1;
  }
  ps_font_size = font_size * 1.2;
  if (dots_flag) {
    fprintf(outfp, "%s Adjust the radius annotated base here\n", "%%");
    fprintf(outfp, "/c_drad { %.4f } def \n", font_size * .72);
  }
  if (outline_mode) {
    line_length_adjust = 0;
  } else {
    line_length_adjust = .30;
  }
  x_shift = ps_font_size / 3;
  y_shift = ps_font_size / 3;
  set_color_ps(COLOR_BASES);
  /* make connections between bases */
  set_color_ps(COLOR_BASES);
  fprintf(outfp, " %f setlinewidth\n", line_width * 1.0);
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
    dif_x = base_x[imod] - base_x[i - 1];
    dif_y = base_y[imod] - base_y[i - 1];
    dist = sqrt(dif_x * dif_x + dif_y * dif_y);
    if (dist > ((distance_bases * scale) * (2 * line_length_adjust))) {
      alpha = atan2(dif_y, dif_x);
      x1 = base_x[i-1] + distance_bases*scale*line_length_adjust*cos(alpha);
      x2 = base_x[i-1] + (dist-distance_bases*scale*line_length_adjust)*
	cos(alpha);
      y1 = base_y[i - 1] + distance_bases*scale*line_length_adjust*sin(alpha);
      y2 = base_y[i - 1] + (dist - distance_bases*scale*line_length_adjust)*
	sin(alpha);
      if (g_next[i - 1] || g_prev[imod]) { /* 2003-06-23 Nick Markham */
	fprintf(outfp, "%.2f %.2f m\n", x1, y1);
	fprintf(outfp, "%.2f %.2f l\n", x2, y2);
	fprintf(outfp, "s\n");
      }
    }
  }
  if (!lines) {
    fprintf(outfp, "%s Adjust the radius of base pair dot here\n", "%%");
    fprintf(outfp, "/c_rad { %.4f } def \n", font_size / 3.0);
  }
  if (lines)
    fprintf(outfp, " %f setlinewidth\n", 2 * line_width);
  for (i = start_base; i <= end_base; i++) {
    if (external_domain_flag) {
      if ((i > ex_start_base) && (i <= ex_end_base))
	i = ex_end_base + 1;
    }
    connected_to = connected[i];
    if (connected_to > i) {
      /* base pair connections */
      set_color_base_ps(bases[i], bases[connected_to]);
      if (lines) {
	x1 = base_x[i] + (base_x[connected_to] - base_x[i])*line_length_adjust;
	x2 = base_x[i]+(base_x[connected_to]-base_x[i])*(1-line_length_adjust);
	y1 = base_y[i] + (base_y[connected_to]-base_y[i])*line_length_adjust;
	y2 = base_y[i]+(base_y[connected_to]-base_y[i])*(1-line_length_adjust);
	fprintf(outfp, "%.2f %.2f m\n", x1, y1);
	fprintf(outfp, "%.2f %.2f l\n", x2, y2);
	fprintf(outfp, "s\n");
      } else { /* draw dot for base pair */
	x1 = (base_x[i] + base_x[connected_to]) / 2.0;
	y1 = (base_y[i] + base_y[connected_to]) / 2.0;
	fprintf(outfp, "%.2f %.2f c_rad 0 360 arc closepath fill\n", x1, y1);
      }
    }
  }
  /* draw diamonds for special bases */
  if (number_of_forced_pairs > 0) {	/* make diamonds */
    fprintf(outfp, " %f setlinewidth\n", line_width / 1.5);
    for (i = 0; i < number_of_forced_pairs; i += 2) {
      base_1 = forced_pairs[i];
      connected_to = connected[base_1];
      set_color_ps(COLOR_BASES);
      x1 = base_x[base_1] +
	(base_x[connected_to] - base_x[base_1]) * line_length_adjust / 1.2;
      x2 = base_x[base_1] + (base_x[connected_to] - base_x[base_1]) *
	(1 - line_length_adjust / 1.2);
      y1 = base_y[base_1] + (base_y[connected_to] - base_y[base_1]) *
	line_length_adjust / 1.2;
      y2 = base_y[base_1] + (base_y[connected_to] - base_y[base_1]) *
	(1 - line_length_adjust / 1.2);
      angle = atan2(y2 - y1, x2 - x1);
      dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2.;
      dist = sqrt(dist * dist + dist * dist);
      x3 = x2 - dist * sin(angle + PI / 4);
      y3 = y2 + dist * cos(angle + PI / 4);
      x4 = x2 - dist * sin(PI / 4 - angle);
      y4 = y2 - dist * cos(PI / 4 - angle);
      fprintf(outfp, "%.2f %.2f moveto\n", x1, y1);
      fprintf(outfp, "%.2f %.2f lineto\n", x3, y3);
      fprintf(outfp, "%.2f %.2f lineto\n", x2, y2);
      fprintf(outfp, "%.2f %.2f lineto\n", x4, y4);
      fprintf(outfp, "closepath\n");
      fprintf(outfp, "gsave\n");
      fprintf(outfp, "stroke \n");
      fprintf(outfp, "grestore\n");
      set_color_base_ps(bases[base_1], bases[connected_to]);
      fprintf(outfp, "fill\n");
    }
    fprintf(outfp, " %f setlinewidth\n", 1 * line_width);
  }
  /* display bases */
  if (!outline_mode) {
    fprintf(outfp, "/Helvetica findfont\n");
    fprintf(outfp, "%f scalefont\n", ps_font_size);
    fprintf(outfp, "setfont\n");
    set_color_ps(COLOR_BASES);
    for (i = start_base; i <= end_base; i++) {
      if (external_domain_flag) {
	if ((i > ex_start_base) && (i <= ex_end_base))
	  i = ex_end_base;
      }
      if (annotation) {
	set_color_ann_ps(ann_to_color[i]);
	if (dots_flag) {
	  fprintf(outfp, "%.2f %.2f c_drad 0 360 arc closepath fill \n",
		  base_x[i], base_y[i]);
	  if (bases_flag) {
	    if (annotation == PROB) {
	      if (ann_to_color[i] > LOG_WHITE_BLACK_SWITCH)
		set_color_ps(COLOR_ANN_WHITE_LETTER);
	      else
		set_color_ps(COLOR_ANN_BLACK_LETTER);
	    } else if (ann_to_color[i] > WHITE_BLACK_SWITCH)
	      set_color_ps(COLOR_ANN_WHITE_LETTER);
	    else
	      set_color_ps(COLOR_ANN_BLACK_LETTER);
	  }
	}
      }
      if ((!annotation) || (bases_flag)) {
	fprintf(outfp, "%.2f %.2f m\n", base_x[i]-x_shift, base_y[i]-y_shift);
	fprintf(outfp, "(%c) show\n", bases[i]);
	fprintf(outfp, "stroke\n");
      }
    }
  }
  fprintf(outfp, "/Helvetica findfont\n");
  fprintf(outfp, "%f scalefont\n", font_size);
  fprintf(outfp, "setfont\n");
  /* draw labels */
  if (total_labels > 0) {
    set_color_ps(COLOR_LABEL_LINE);
    fprintf(outfp, " %f setlinewidth\n", line_width / 2.0);
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
	    log(structure_label_value[i] + history_offset ) / log(10); 
	}
      } else { 	/* 2003-06-20 Nick Markham */
	sprintf(label_content, "%d", g_history[structure_label_value[i]]);
	label_content_offset =
	  log(structure_label_value[i] + history_offset) / log(10);
      }
      fprintf(outfp, "%.2f %.2f m\n", structure_label_x[i] -
	      x_shift * label_content_offset, structure_label_y[i] - y_shift);
      fprintf(outfp, "(%s) show\n", label_content);
      x1 = base_x[structure_label_value[i]] +
	(structure_label_x[i] - base_x[structure_label_value[i]]) / 6;
      y1 = base_y[structure_label_value[i]] +
	(structure_label_y[i] - base_y[structure_label_value[i]]) / 6;
      x2 = base_x[structure_label_value[i]] +
	5. * (structure_label_x[i] - base_x[structure_label_value[i]]) / 8;
      y2 = base_y[structure_label_value[i]] +
	5. * (structure_label_y[i] - base_y[structure_label_value[i]]) / 8;
      fprintf(outfp, "%.2f %.2f m\n", x1, y1);
      fprintf(outfp, "%.2f %.2f l\n", x2, y2);
      fprintf(outfp, "s\n");
    }
  }
  if (loop_labels) {
    for (i = first_loop; i <= total_loops; i++) {
      fprintf(outfp, " %.2f %.2f m", loop_center_x[i] - x_shift*log(i)/log(10),
	      loop_center_y[i] - y_shift);
      fprintf(outfp, "(%d) show\n", i);
      fprintf(outfp, " s\n");
    }
  }
}

int             open_output(char *);

void make_structure_ps_real(char *filename, int smart_colors, char *first_line,
			    float line_width, int *connected, int start_base,
			    int end_base, int total_loops,
			    float *loop_center_x, float *loop_center_y,
			    float font_size, int outline_mode, float scale,
			    char *bases, float *base_x, float *base_y,
			    int lines, float *structure_label_x, 
			    float *structure_label_y, 
			    int *structure_label_value, int loop_labels,
			    int total_labels, int external_domain_flag,
			    int ex_start_base, int ex_end_base,
			    int history_offset, int annotation,
			    int *ann_to_color, int dots_flag, int bases_flag,
			    float x_shift, float y_shift, float zoom_factor,
			    int is_circular, int *forced_pairs,
			    int number_of_forced_pairs) {
  char            structure_filename_ps[120];
  strcpy(structure_filename_ps, filename);
  strcat(structure_filename_ps, ".ps");
  open_output(structure_filename_ps);
  initialize_ps(smart_colors, annotation, x_shift, y_shift, zoom_factor);
  display_bases_ps(first_line, line_width, connected, start_base, end_base,
		   total_loops, loop_center_x, loop_center_y, font_size,
		   outline_mode, scale, bases, base_x, base_y, lines,
		   total_labels, structure_label_x, structure_label_y,
		   structure_label_value, loop_labels, external_domain_flag,
		   ex_start_base, ex_end_base, history_offset, annotation,
		   ann_to_color, dots_flag, bases_flag, is_circular,
		   forced_pairs, number_of_forced_pairs);
  finish_ps_file();
}
