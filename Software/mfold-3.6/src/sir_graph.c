/* Michael Zuker, Nick Markham 
   November 25, 2006: -bp option for sir_graph creates 1 or 2 new
   base pairs when 1x1 or 2x2 i-loops are detected, respectively.
   Nick Markham made use of libgd (if available) to
   create images in png, jpg and gif formats. Otherwise,
   these formats cannot be made. He also recreated the
   use of "configure" for the package; first done by
   Alex Yu.
   M. Zuker made many cosmetic changes to make the code
   readable, especially 'call_file_makers' in sir_graph.c
   Nick enabled sir_graph to use historical numbering
   properly.
   Bug created by Nick in converting to libgd corrected
   by M. Zuker on 9/29/06. See Zuker note in
   sir_graph_img.h

   Original sir_graph by Darrin Stewart & Michael Zuker 
   Copyright 1998-2000 
   Washington University 
   Update Nov 13, 2000 sir_graph.c 
   Execute using:
   sir_graph name.ct 
   .ct is optional for the circle graph, bases are
   labeled clockwise with 1 at the bottom 

   * uses external domain from command line 
   * Uses historical numbering 
   * uses counter clockwise drawing 
   * Uses color annotation 
   * non au,gu,gc base pairs made yellow on November 10, 1999 
   * fixed probability annotation to use white letters on black dots 
   * Nov 10, 1999 
   * Nov 18  Loop Stretcher 
   * Nov 19  Natural angles change with no bend on bulge of 1 or 2 
   * Nov 22  Added window  to accept name for ps, gif, or ss file   
   * Nov 23  start on loop naturalizer 
   *
   * Jan 19 2000 Made array size variable to match sequence length with
   * sir_graph_alloc.c
   * Jan 20 2000 Fixed to read GCG format 
   * Jan 20 2000 Made gif output match gif zoom window 
   * Jan 21 2000 Fixed undo to take into account rescale of window 
   * Feb 17 2000 Fixed -ad flag 
   * March 10  flat loop for -f or -fa 
   * It should take a loop and display it as a line with helices coming off it 
   *
   * It works with natural angles, and circle graph angles for clockwise or
   * counter clockwise
   * It uses the exterior loop unless -i specifies bases from some other loop 
   * -fa uses alternating clockwise, counter clockwise to make the helices 
   * alternate up/down from the single strand of the loop 
   * March 15,2000 Added -col file_name option to allow user defined 
   * colors.  sir_graph.col is default colors, sir_graph.1.col 
   * uses a black background.  Annotation colors are still not adjustable 
   * when not used the colors are defined from sir_graph_color.h 
   * March 27, 2000 added -force option 
   * March 31, 2000 added -aj switch 
   * April 10, 2000, fixed natural angle, pairing of first base with   
   * last base problem 
   * June 16, 2000, Added Mouse: single strand linear option 
   * July 20, 2000, Reorganized loop drawing for circle graph angles 
   * Separated circle graph and structure portions of code 
   * Aug 1, 2000 Regularize angles for circle graph angles 
   * Click on loop and get angles regularized          
   * Aug 7, 2000 Fix loop by spreading helices and adjusting single strand
   * bases
   * Sep 14, add -x, -zoom_ps, -zoom options 
   * Last version before creation of sir_graph_core.c 
   */

#define PROGRAM_NAME "sir_graph"
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <time.h>   /* for time within postscript file */ 
#include "util.h"
#define TRUE 1
#define FALSE 0
#include "sir_graph_color.h"

#define FLOAT_MAX 5E32
#define FLOAT_MIN -5E32
#define HT 1000.0
#define WD 720.0
#define DISTANCE_BASES 1
#define DISTANCE_BASE_PAIRS 1
#define PI 3.141592653589793
#define POSTSCRIPT 1
#define GIF        2
#define SS         3
#define JPG        4
#define PNG        5

FILE           *outfp;
FILE           *ctfp;
FILE           *ssfp;		/* for input */
FILE           *outssfp;

void            memory_allocator(int);

/* annotation variables */

#define NONE 0
#define P_NUM 1
#define PROB  2
#define SS_COUNT 3
#include "sir_graph_color_ann.h"
#include "sir_graph_ann.h"

void            report_mem_error(char *);
int* g_history;	//2003-06-20 Nick Markham
int* g_next;	//2003-06-20 Nick Markham
int* g_prev;	//2003-06-20 Nick Markham
int* g_structure_label_value;	//2003-06-21 Nick Markham
#include "sir_graph_ps.h"	/* creates structure output for ps */
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
# include "sir_graph_img.h"	/* creates jpg, png and gif */
#endif

int             g_aj;
int             g_annotation;
int             g_create_ann_table;
char            g_unused[120] = ""; 
char            g_annotation_table_type = 'h';
int             g_counter_clockwise;
int             g_fix_loop_flag; /* spreads out helices around a loop */
int             g_domain_defined;
int             g_domain_row;
int             g_domain_column;
int             g_reg_angle_flag;
int             g_reg_angle;
char            g_annotation_filename[120] = "";
int             g_annotation_bases;	/* true or false */
int             g_annotation_dots;	/* true or false */
int             g_png_mode = FALSE;  /* 2007-01-03 Zuker initializes */
int             g_jpg_mode = FALSE;  /* 2007-01-03 Zuker initializes */
int             g_command_line_zoom_ps_flag; /* for -zoom_ps switch */
int             g_command_line_zoom_ps_x;
int             g_command_line_zoom_ps_y;
float           g_command_line_zoom_ps_s;
int             g_command_line_zoom_flag;
int             g_command_line_zoom_x;
int             g_command_line_zoom_y;
float           g_command_line_zoom_s;
int             g_command_line_x;	/* when true, create list of x,y
					 * locations of base pairs for 
					 * png,jpg images */
float          *g_midpoint_x; /* These 2 used only in traverse_loop */
float          *g_midpoint_y;
/* Below use data  to highlight forced base pairs */
int             g_complot_forced_bases[500];	/* allow at most 500 */
/* Occur as base pairs, 0-1, 2-3, 4-5, not ordered otherwise */
int             g_complot_forced_bases_count;	/* 0 indicates none */
/* k indicates using positions 0 to k-1 of g_complot_forced_bases */
char            g_forced_file[100];
FILE           *forcefp;
int             g_label_forced_row;
int             g_label_forced_column;
float           g_degrees_used;
float           g_rotation_angle;
float           g_close_enough;
float          *g_base_x;
float          *g_base_y;
float          *g_structure_angle;
float          *g_structure_label_x;
float          *g_structure_label_y;
int             g_history_offset;
float          *g_undo1_base_x;
float          *g_undo1_base_y;
float          *g_undo1_structure_angle;
float           g_undo1_scale;
float           g_undo1_ss_distance_basepair;
float           g_undo1_ss_distance_base_adj;
float           g_undo1_font_size;
float           g_undo1_line_width;
float           g_undo2_ss_distance_basepair;
float           g_undo2_ss_distance_base_adj;
float           g_undo2_font_size;
float           g_undo2_line_width;
float          *g_undo2_base_x;
float          *g_undo2_base_y;
float          *g_undo2_structure_angle;
float           g_undo2_scale;
int             g_need_to_store_for_undo = FALSE;
int            *g_base_list;	/* used by traverse_loop and check
				 * multiloop */
/* g_history[i]=j indicates that base i has historical position j */
int            *g_previous_loop;
/* for each base, this stores the previous loop */
int             g_mouse_rotate_helix_set;
int             g_auto_rotate;
int            *g_ss_code;
int             g_img_mode;
int             g_outline_mode;
int             g_interactive_mode;
int             g_img_interlaced;
int             g_img_width;
int             g_img_height;
int             g_ss_mode;
int             g_first_ss_run;
int             g_ss_output_mode;
int             g_ss_backwards;
int             g_external_domain_flag;	/* set by -e or -ew on command line */
int             g_label_frequency;
float           g_helix_base_advance;
float           g_ss_distance_basepair;
float           g_ss_distance_base_adj;
float          *g_base_advance;

/* TRUE for normal, False for small */
float           g_natural;
float           g_small_base_advance;
float          *g_loop_center_x;
float          *g_loop_center_y;
float          *g_undo1_loop_center_x;
float          *g_undo1_loop_center_y;
float          *g_undo2_loop_center_x;
float          *g_undo2_loop_center_y;
float          *g_loop_radius;
int             g_ex_start_base;
int             g_ex_end_base;
char            g_first_line[150];
char           *g_bases;	/* stores base (letter at position i) */
char            g_fig_title[150]; /* M. Zuker, June 26, 2006. 
				   * char  *g_fig_title;      */
float           g_font_size;
float           g_line_width;
int             g_total_loops;
float           g_scale;
int             g_flat;	/* when true, draw exterior loop flat */
int             g_flat_alternate; /* helices alternate between
				   * clockwise and counter clockwise */ 
int             g_flat_edit;	/* true for mouse dragging helix of
				 * external loop */
int             g_flat_fixed_base; /* for above, the j base of the
				    * previous base pair on the loop */  
int             g_lines;
int             g_small_angle;
int             g_total_labels;
int             g_whole_circle;
char            g_ct_filename[120];
char            g_ss_filename[120];
char            g_out_ss_filename[120];
char            g_out_filename[120];
char            g_structure_filename[120];
char            g_structure_img_filename[120];
int             g_bp;  /* M. Zuker, 11/25/07 If TRUE, convert 1x1 and 2x2 
			* interior loops to 1 or 2 base pairs */
int             g_length = 0;
int             g_smart_colors;
int             g_midpoints;
int             g_is_circular;
int             g_angle_list;
int             g_arcs;
int             g_structure;
int             g_start_base = 0; /* M. Zuker, 5/19/07 */
int             g_end_base;
int             g_loop_labels;
int             g_stack_items = 0;	/* points to next available item */
int            *g_connected; /* connected[i]=j indicates base i
			      * is connected to base j */
int            *g_oldconnect; /* As above, but without additional bps
			       * that might be inserted when -bp is used  
			       */
int            *g_stack;
#include "sir_graph_ss.h"	/* read or create ss file */
#include "sir_graph_ct.h"	/* reads ct file */
#include "sir_graph_alloc.h"	/* allocates memory for arrays,pointers */
float           g_diam;

/* each base */

float           g_circle_radius_point;	/* radius of big circle in
					 * points */
void num_string_float(char *text, float value) {
  sprintf(text, "%.4f", value);
}

void num_string_int(char *text, int value) {
  sprintf(text, "%d", value);
}

void    store_for_undo(void);
void    show_both(void) ;

/* M. Zuker, Dec 19, 2006. 
 * Creates error function to replace many similar parts of code where
 * pseudoknots are detected
 */
int error(char *type, int base1, int base2 ) {
  int i,j,count=0;
  if (!strcmp(type, "pk")) {
    for (i=base1; i<=base2; i++) {
      j = g_oldconnect[i];
      if ( (j>0) && ((j < base1) || (j > base2)) ) {
	count++;
	if (count==1) 
	  printf("Error!\tPseudoknot detected.\n");
	printf("The base pair, %d.%d, conflicts with the base pair %d.%d.\n",
	       i,j,base1,base2);
      }
    }
    return (1);
  }
  return (1); /* for safety */
}

void push_stack(int x) {
  if (g_stack_items == g_length) {
    printf("Error!\tstack overflow, enlarge stack size with g_stack from ");
    printf("memory_allocator from memory_allocator()");
    try_exit(1);
  } else {
    g_stack[g_stack_items] = x;
    g_stack_items++;
  }
}

int pop_stack(void) {
  if (g_stack_items < 1) {
    printf("Error!\tCannot pop item off empty stack.\n");
    try_exit(2);
  } else {
    g_stack_items--;
  }
  return (g_stack[g_stack_items]);
}

int open_output(char *out_filename) {
  /* open specified output file */
  printf("Output:  %s\n", out_filename);
  if ((outfp = fopen(out_filename, "wb")) == NULL) {
    printf("Error!\tCould not open file: %s\n", out_filename);
    try_exit(3);
  }
  return (0);
}

void adjust_start_and_end(int *start, int *end, int length, int *connected) {
  int end_base, start_base;
  end_base = *end;
  start_base = *start;
  if (end_base > length)
    end_base = length;
  if (start_base >= 0 && connected[start_base] > end_base) //20/6/03 N Markham
    end_base = connected[start_base];
  else {
    if (start_base >= 0 && (connected[start_base] > 0) &&  //20/6/03 N Markham
	(connected[start_base] < start_base))
      start_base = connected[start_base];
  }
  if ((end_base != length) && (connected[end_base] > 0)
      && (connected[end_base] < start_base))
    start_base = connected[end_base];
  if (start_base < 1) {
    start_base = 1;
    end_base = length;
  }
  if (end_base > length)
    end_base = length;
  *start = start_base;
  *end = end_base;
}

int step_fun(int dis_l, int dis_r) {
  int             step;
  step = (dis_r - dis_l + 1) / 12; /* make 8 steps  (9 tick marks total) */
  if (step < 40) {
    if (step < 3)
      step = 2;
    else if (step > 5) {
      step = (step + 5) / 10;	/* make a large step a
				 * multiple of ten */
      step = step * 10;
    } else if (step == 4)	/* make steps of 4 look like 5 */
      step = 5;
  } else if (step > 90) {
    step = (step + 50) / 100;	/* round to nearest 100 */
    step = step * 100;
  } else {
    step = (step + 25) / 50;	/* round to nearest 50 */
    step = step * 50;
  }
  return step;
}

void draw_label(int label, int circle_radius_point) {
  /* M. Zuker, June 17, 2007. 
   * Add adjust
   */
  float           adjust = 0.0, angle, xpos, ypos;
  angle = g_base_advance[label];
  /* M. Zuker, June 17, 2007. Adjust to 13x17 output 
     Remove following two lines. Don't understand why they were ever there.
     if ((label == g_end_base) || (angle < (PI / -2.0)))
     circle_radius_point = circle_radius_point + 13;
  */
  if ( label!=g_start_base && angle > PI/2.0  ) {
    adjust = 15.;
  }
  set_color_ps(COLOR_LABEL_LINE);
  /* M. Zuker, June 17, 2007. Adjust to 13x17 output 
     xpos = (circle_radius_point * 1.15) * cos(angle);
     ypos = (circle_radius_point * 1.15) * sin(angle);
  */
  xpos = (circle_radius_point * 1.06) * cos(angle);
  ypos = (circle_radius_point * 1.06) * sin(angle);
  fprintf(outfp, "%.2f %.2f m\n", xpos, ypos);
  xpos = (circle_radius_point * 1.01) * cos(angle);
  ypos = (circle_radius_point * 1.01) * sin(angle);
  fprintf(outfp, "%.2f %.2f l\n", xpos, ypos);
  fprintf(outfp, "s\n");
  set_color_ps(COLOR_LABEL_NUMBER);
  /* M. Zuker, June 17, 2007. Adjust to 13x17 output 
     xpos = (circle_radius_point * 1.18) * cos(angle) - 12;
     ypos = (circle_radius_point * 1.18) * sin(angle) - 4;
  */
  xpos = (circle_radius_point * 1.08) * cos(angle) - adjust  ;
  ypos = (circle_radius_point * 1.08) * sin(angle) - 4 ;
  fprintf(outfp, "%.2f %.2f m\n", xpos, ypos);
  fprintf(outfp, "(%d) show\n", g_history[label]); /* 12/12/06 M. Zuker */
}

void swap_angles(float *angle1, float *angle2) {
  float           a1;
  float           a2;
  float           temp;
  a1 = *angle1;
  a2 = *angle2;
  /* keep a1 the biggest */
  if (((a1 * a2) > 0)) {	/* both positive or both negative */
    if (a1 < a2) {
      temp = a1;
      a1 = a2;
      a2 = temp;
      *angle1 = a1;
      *angle2 = a2;
    }
  }
  if ((a1 <= 0) && (a2 > 0)) {	/* draw counter-clockwise from bigger
				 * a1 to smaller a2 */
    temp = a1;	/* -90<=angle<=270 */
    a1 = a2;
    a2 = temp;
    *angle1 = a1;
    *angle2 = a2;
  }
}

void fix_angles(float *angle1, float *angle2) {
  float           a1, a2;
  float           temp;
  a1 = *angle1;
  a2 = *angle2;
  if (g_structure) {	/* the line below may need work */
    /* when drawing of angles cross, add 2 pi and switch them */
    if ((a2 + 2.0 * PI) < a1) {
      a2 = a2 + 2. * PI;
      temp = a2;
      a2 = a1;
      a1 = temp;
    }
    if ((a2 < (3.0 * PI / -2.0)) && (a1 < (3.0 * PI / -2.0))) {
      /* above line is not understood 
       * when both less than -270, add 2 pi and switch them
       * really both on second time around circle */
      a2 = a2 + 2.0 * PI;
      a1 = a1 + 2.0 * PI;
      temp = a2;
      a2 = a1;
      a1 = temp;
    }
  } else {
    while (a1 < (PI / -2.0))
      a1 = a1 + 2.0 * PI;
    while (a2 < (PI / -2.0))
      a2 = a2 + 2.0 * PI;
    swap_angles(&a1, &a2);
  }
  *angle1 = a1;
  *angle2 = a2;
}

void set_midpoints_no_structure(int base_1, int connected_to, float *x,
				float *y) {
  float           angle1, angle2, phi;
  float           midpoint_x, midpoint_y, midpoint_angle;
  float           x1, x2, y1, y2, start_angle, end_angle, dist, X0, Y0;
  float           print_angle=0.0;
  float           gamma, R;
  int             phi_greater_than_pi;
  float           print_angle_perp;
  angle1 = g_base_advance[base_1];
  angle2 = g_base_advance[connected_to];
  fix_angles(&angle1, &angle2);
  phi = (angle1 - angle2);
  while (phi >= (2.0 * PI))
    phi = phi - 2.0 * PI;
  if (((phi - PI) * (phi - PI)) < .00005) {
    /* phi is close to PI; draw line */
    x1 = g_circle_radius_point * cos(angle1);
    y1 = g_circle_radius_point * sin(angle1);
    x2 = g_circle_radius_point * cos(angle2);
    y2 = g_circle_radius_point * sin(angle2);
    midpoint_x = (x1 + x2) / 2.0;
    midpoint_y = (y1 + y2) / 2.0;
    phi_greater_than_pi = FALSE;
    print_angle = angle1 + PI;
  } else {
    if (phi > PI) {	/* Phi is > PI */
      phi = (2 * PI - phi);
      phi = phi / 2.0;
      gamma = angle2 - phi;
      end_angle = angle1 * 180. / PI - 90;
      start_angle = angle2 * 180. / PI + 90;
      phi_greater_than_pi = TRUE;
      print_angle = (gamma + PI / 2.);
    } else {	/* phi is < PI */
      phi = phi / 2.0;
      gamma = angle2 + phi;
      start_angle = angle1 * 180. / PI + 90.;
      end_angle = angle2 * 180. / PI - 90.;
      phi_greater_than_pi = FALSE;
      if (g_angle_list) {
	print_angle = (gamma + PI / 2.0);
	print_angle = print_angle + PI;
	if (print_angle > (2.0 * PI))
	  print_angle = print_angle - 2.0 * PI;
      }
    }
    if (print_angle < 0)
      print_angle = print_angle + 2.0 * PI;
    print_angle_perp = print_angle + PI / 2.0;
    if (print_angle_perp > (2.0 * PI))
      print_angle_perp = print_angle_perp - 2.0 * PI;
    R = tan(phi) * g_circle_radius_point;
    dist = sqrt(R*R + g_circle_radius_point*g_circle_radius_point);
    X0 = dist * cos(gamma);
    Y0 = dist * sin(gamma);
    midpoint_angle = (start_angle + end_angle) / 2.0;
    if (!phi_greater_than_pi)
      midpoint_angle = midpoint_angle + 180;
    midpoint_y = Y0 + sin(midpoint_angle * PI / 180.) * R;
    midpoint_x = X0 + cos(midpoint_angle * PI / 180.) * R;
  }
  *x = midpoint_x;
  *y = midpoint_y;
}

float regularize_angle(float input_angle) {
  /* input in radians, output in radians 
   * For the first base of helix, If this is the same loop as the last
   * helix, if the angles are the same, add 1 degree.
   */
  int             is_neg;	/* is negative flag */
  int             deg_value;
  int             pre_round_deg_value;
  int             remainder;
  /* convert to degrees */
  input_angle = input_angle * 180. / PI;
  /* round to degree based on g_reg_angle */
  is_neg = FALSE;
  /* make positive */
  if (input_angle < 0) {
    is_neg = TRUE;
    input_angle = input_angle * -1.0;
  }
  /* convert to integer */
  pre_round_deg_value = (int) (input_angle + .5);
  /* find remainder when divided by reg_angle */
  remainder = pre_round_deg_value % g_reg_angle;
  deg_value = pre_round_deg_value - remainder;
  if (remainder >= (g_reg_angle / 2.0)) {
    deg_value = deg_value + g_reg_angle;
    /* round up */
  }
  /* otherwise round down */
  if (is_neg)
    deg_value = deg_value * -1;
  return (deg_value * PI / 180.);
}

void get_perpendicular_angle(int base_1, int connected_to,float *helix_angle) {
  /* looks a lot like set_midpoints
   * given bases base_1 and connected_to, return angle perpendicular
   * to circular arc between the two bases
   * when going around a loop, the order must be preserved */
  float           angle1, angle2, phi;
  float           print_angle;
  float           gamma;
  float           print_angle_perp;
  angle1 = g_base_advance[base_1];
  angle2 = g_base_advance[connected_to];
  fix_angles(&angle1, &angle2);
  phi = (angle1 - angle2);
  while (phi >= (2.0 * PI))
    phi = phi - 2.0 * PI;
  if (((phi - PI) * (phi - PI)) < .00005) {
    print_angle = angle1 + PI;
    print_angle_perp = print_angle + PI / 2.0;
  } else {
    if (phi > PI) {
      phi = (2 * PI - phi);
      phi = phi / 2.0;
      gamma = angle2 - phi;
      print_angle = (gamma + PI / 2.);
    } else {
      phi = phi / 2.0;
      gamma = angle2 + phi;
      print_angle = (gamma + PI / 2.0);
      print_angle = print_angle + PI;
    }
  }
  print_angle_perp = print_angle + PI / 2.0;
  while (print_angle_perp < 0)
    print_angle_perp = print_angle_perp + 2.0 * PI;
  while (print_angle_perp > (2.0 * PI))
    print_angle_perp = print_angle_perp - (2.0 * PI);
  *helix_angle = print_angle_perp;
}

void store_base_angles(float angle1, float angle2, int bases_on_list,
		       int *base_list) {
  float           angle_increment;
  int             i, base;
  if (!g_counter_clockwise) {
    if (angle1 < angle2)	/* Aug 9, 1999 */
      angle2 = angle2 - 2.0 * PI;
    while ((angle1 - angle2) > (2.0 * PI + .0001)) {
      angle1 = angle1 - 2.0 * PI;
    }
    if (((angle1 - angle2) * (angle1 - angle2)) < .00001) {
      angle1 = angle2 + 2.0 * PI;
    }
    if ((angle1 < 0) && (angle2 > 0))
      angle1 = angle1 + 2.0 * PI;
  } else {
    while (angle1 > angle2) {
      angle1 = angle1 - 2.0 * PI;
    }
    while ((angle2 - angle1) > (2.0 * PI))
      angle2 = angle2 - 2.0 * PI;
    if (((angle1 - angle2) * (angle1 - angle2)) < .00001) {
      angle2 = angle1 + 2.0 * PI;
    }
  }
  angle_increment = (angle1 - angle2) / (bases_on_list + 3.);
  angle1 = angle1 - angle_increment;
  for (i = 1; i <= bases_on_list; i++) {
    base = base_list[i];
    g_base_x[base] = angle1 - i * angle_increment;
    if (g_label_frequency > 0)
      g_structure_angle[base] = angle1 - i * angle_increment;
  }
}

void draw_basepair(int i, int connected_to, float last_angle_perp, 
		   float current_position_x, float current_position_y, 
		   float distance_base_pairs) {
  int             b1, b2;
  b1 = i;
  b2 = connected_to;
  if (g_counter_clockwise) {
    g_base_x[b2] = current_position_x + distance_base_pairs / 2.0 *
      cos(PI / 2.0 + last_angle_perp);
    g_structure_angle[b2] = PI / 2.0 + last_angle_perp;
    g_structure_angle[b1] = PI / -2.0 + last_angle_perp;
    g_base_y[b2] = current_position_y + distance_base_pairs / 2.0 *
      sin(PI / 2.0 + last_angle_perp);
    g_base_x[b1] = current_position_x + distance_base_pairs / 2.0 *
      cos(last_angle_perp - PI / 2.0);
    g_base_y[b1] = current_position_y + distance_base_pairs / 2.0 *
      sin(last_angle_perp - PI / 2.0);
  } else {
    g_base_x[b1] = current_position_x + distance_base_pairs / 2.0 *
      cos(PI / 2.0 + last_angle_perp);
    g_structure_angle[b1] = PI / 2.0 + last_angle_perp;
    g_structure_angle[b2] = PI / -2.0 + last_angle_perp;
    g_base_y[b1] = current_position_y + distance_base_pairs / 2.0 *
      sin(PI / 2.0 + last_angle_perp);
    g_base_x[b2] = current_position_x + distance_base_pairs / 2.0 *
      cos(last_angle_perp - PI / 2.0);
    g_base_y[b2] = current_position_y + distance_base_pairs / 2.0 *
      sin(last_angle_perp - PI / 2.0);
  }
}

void loop_fix_first_pass(int loop_start,int loop_end,float distance_base_pairs,
			 float distance_bases,int loop_num);

void traverse_loop_structure(int loop_num, int loop_start, int loop_end,
			     int *connected, float angle,
			     float current_position_x,float current_position_y,
			     int extra_base_exists, int largest_base,
			     float distance_base_pairs, float distance_bases,
			     int draw_first_base_pair) {
  /* loop_num is the loop being drawn 
   * loop runs on the base pair loop_start to loop_end 
   * angle is angle of helix this is connected to 
   * from perpendicular 
   * current_position_x,current_position_y are last drawn base pair 
   * largest base is g_length, or the largest base with -iw or -i 
   *
   * extra_base_exists draws trailing single stranded base on first
   * loop of structure
   * This function is for circle graph angles 
   *
   * For a loop beginning with a base pair, draw it if
   * draw_first_base_pair is true.  
   * Normally a loop starts with this base pair already drawn
   * The first loop of a structure is the exception */
  int             bases_on_list;
  int             size;
  int             current_place;
  int             connected_to;
  float           last_helix_angle, helix_angle = 0.0;
  float           base_angle;
  float           placement_angle;
  int             last_base_was_connected;
  int             helices_on_loop;
  int             times_around;
  float           last_regularized_helix_angle=0.0;
  float           first_angle;
  float           base_pair_angle;	/* perpendicular angle of
					 * base_pair on loop */
  times_around = 0;
  current_place = loop_start;
  size = 1;
  if (g_reg_angle_flag) {
    /*    last_regularized_helix_angle = -5000.; May 19,  2007. M. Zuker */
    last_regularized_helix_angle = 0.0;
  }
  last_helix_angle = angle + PI;
  helices_on_loop = 0;
  bases_on_list = 0;
  if (g_angle_list) {
    printf("Traversing loop %d from base %d to base %d\n", loop_num,
	   loop_start, loop_end);
  }
  if (current_place == 0) {
    if ((g_angle_list) && (g_start_base == 1))
      printf("5 prime_end\n");
    current_place = 1;
    if (!g_is_circular)
      size = 2;
  }
  if (connected[current_place] > current_place)
    last_base_was_connected = TRUE;
  else
    last_base_was_connected = FALSE;
  /* last_base_was_connected keeps track if base prior to
   * current_place was paired; Go around loop setting angles 
  */ 
  while (current_place != loop_end) {
    if (current_place > largest_base) {
      if (times_around == 0) {
	times_around++;
	if (g_is_circular)
	  current_place = g_start_base;
	else
	  current_place = g_start_base - 1;
      } else {
	if (error("pk", loop_start, loop_end)==1)
	  try_exit(4);
      }
    }
    connected_to = connected[current_place];
    if ((connected_to > current_place) && (connected_to <= largest_base)) {
      if (g_angle_list)
	printf("%d * %d\n", current_place, connected_to);
      if (current_place > loop_start) {
	push_stack(loop_num);
      }
      get_perpendicular_angle(current_place, connected_to, &helix_angle);
      if (g_reg_angle_flag) {
	helix_angle = regularize_angle(helix_angle);
	if ((helix_angle==last_regularized_helix_angle)&&(helices_on_loop>1)) {
	  /* fix duplicate angles for helices */
	  if (g_counter_clockwise)
	    helix_angle = helix_angle + g_reg_angle * PI / 180.;
	  else
	    helix_angle = helix_angle - g_reg_angle * PI / 180.;
	}
	last_regularized_helix_angle = helix_angle;
      }
      if ((current_place > loop_start) || (draw_first_base_pair))
	g_base_x[current_place] = helix_angle;
      helices_on_loop++;
      if (extra_base_exists) {
	current_place = connected_to;
	push_stack(loop_num);
      }
      if (size > 1) {	/* set angles of bases on list */
	current_place = connected_to;
	last_base_was_connected = TRUE;
	store_base_angles(last_helix_angle, helix_angle, bases_on_list,
			  g_base_list);
	last_helix_angle = helix_angle;
	bases_on_list = 0;
      } else {
	current_place++;
	if (current_place > largest_base) {
	  if (times_around == 0) {
	    times_around++;
	    if (g_is_circular)
	      current_place = g_start_base;
	    else
	      current_place = g_start_base - 1;
	  } else {
	    if (error("pk", loop_start, loop_end)==1)
	      try_exit(5);
	  }
	}
	last_base_was_connected = FALSE;
	if (((loop_num == 1) && (draw_first_base_pair)) &&
	    g_reg_angle_flag) {
	  last_helix_angle = helix_angle; 
	  /* ??????? Added Sep 11, 2000 */
	}
      }
    } else {
      if (g_angle_list)
	printf("%d\n", current_place);
      if (!last_base_was_connected) {
	bases_on_list++;
	g_base_list[bases_on_list] = current_place;
      }
      current_place++;
      if (current_place > largest_base) {
	if (times_around == 0) {
	  times_around++;
	  if (g_is_circular)
	    current_place = g_start_base;
	  else
	    current_place = g_start_base - 1;
	} else {
	  if (error("pk", loop_start, loop_end)==1)
	    try_exit(6);
	}
      }
      last_base_was_connected = FALSE;
    }
    size++;
  }
  if (extra_base_exists) {
    if (g_angle_list)
      printf("%d\n", current_place);
    if (!last_base_was_connected) {
      bases_on_list++;
      g_base_list[bases_on_list] = current_place;
    }
    size++;
  }
  if (bases_on_list > 0) {	/* set radius of bases on list */
    if (((loop_num == 1) && (draw_first_base_pair)) && g_reg_angle_flag) {
      get_perpendicular_angle(loop_start,g_connected[loop_start],&first_angle);
      angle = regularize_angle(first_angle);
      store_base_angles(helix_angle, angle, bases_on_list, g_base_list);
    } else {
      if (g_counter_clockwise) {
	if ((helices_on_loop > 1) || (draw_first_base_pair))
	  store_base_angles(helix_angle, angle + PI, bases_on_list,
			    g_base_list);
	else
	  store_base_angles(angle - PI, angle + PI,
			    bases_on_list, g_base_list);
      } else {
	if ((helices_on_loop > 1) || (draw_first_base_pair))
	  store_base_angles(helix_angle, angle - PI, bases_on_list,
			    g_base_list);
	else
	  store_base_angles(angle + PI, angle - PI,
			    bases_on_list, g_base_list);
      }
    }
  }
  if (g_angle_list) {
    printf("%d\n", loop_end);
    printf("Loop Size: %d\n", size);
  }
  /* set center and radius of the loop */
  g_loop_radius[loop_num] = (float) size / (2.0 * PI) * g_scale;
  if (loop_num > 1) {
    g_loop_center_x[loop_num] = current_position_x +
      g_loop_radius[loop_num] * cos(angle);
    g_loop_center_y[loop_num] = current_position_y +
      g_loop_radius[loop_num] * sin(angle);
  } else {
    g_loop_center_x[1] = 0.0;
    g_loop_center_y[1] = 0.0;
  }
  /* set positions of  bases around loop based on angles and */
  /* number of bases around loop */
  current_place = loop_start;
  size = 1;
  times_around = 0;
  while (current_place != loop_end) {
    connected_to = g_connected[current_place];
    if ((connected_to > current_place) && (connected_to <= g_end_base)) {
      /* draw the base pair here */
      if ((current_place > loop_start) || (draw_first_base_pair)) {
	base_pair_angle = g_base_x[current_place];
	placement_angle = base_pair_angle;
	draw_basepair(current_place, connected_to,base_pair_angle,
		      g_loop_center_x[loop_num] + 
		      (g_loop_radius[loop_num] + .1) * cos(placement_angle),
		      g_loop_center_y[loop_num] +
		      (g_loop_radius[loop_num] + .1) * sin(placement_angle),
		      distance_base_pairs);
      }
      if ((size > 1) || (loop_num == 1))
	current_place = connected_to;
      else {
	current_place++;
	if (current_place > largest_base) {
	  if (times_around == 0) {
	    times_around++;
	    if (g_is_circular)
	      current_place = g_start_base;
	    else
	      current_place = g_start_base - 1;
	  } else {
	    if(error("pk", loop_start, loop_end)==1)
	      try_exit(7);
	  }
	}
      }
    } else {
      if (connected_to < 1) {
	base_angle = g_base_x[current_place];
	g_base_x[current_place] = g_loop_center_x[loop_num] +
	  (g_loop_radius[loop_num] + .1) * cos(base_angle);
	g_base_y[current_place] = g_loop_center_y[loop_num] +
	  (g_loop_radius[loop_num] + .1) * sin(base_angle);
      }
      current_place++;
      if (current_place > largest_base) {
	if (times_around == 0) {
	  times_around++;
	  if (g_is_circular)
	    current_place = g_start_base;
	  else
	    current_place = g_start_base - 1;
	} else {
	  if (error("pk", loop_start, loop_end)==1)
	    try_exit(8);
	}
      }
    }
    size++;
  }
  if (extra_base_exists) {
    connected_to = connected[current_place];
    if (connected_to < 1) {
      base_angle = g_base_x[current_place];
      g_base_x[current_place] = g_loop_center_x[loop_num] +
	(g_loop_radius[loop_num] + .1) * cos(base_angle);
      g_base_y[current_place] = g_loop_center_y[loop_num] +
	(g_loop_radius[loop_num] + .1) * sin(base_angle);
    }
  }
  if (g_fix_loop_flag)
    loop_fix_first_pass(loop_start, loop_end, distance_base_pairs,
			distance_bases, loop_num);
}

void traverse_loop_no_structure(int loop_num, int loop_start, int loop_end,
				int *connected, float x1, float y1,
				int extra_base_exists, int largest_base) {
  /* x1 and y1 are last midpoint calculated
   * angle is angle of helix this is connected to from perpendicular
   * largest base is g_length, or the largest base with -iw or -i */
  int             size;
  float           x, y;
  int             current_place;
  int             connected_to;
  int             i;
  float           sum_x, sum_y;
  int             helices_on_loop;
  int             times_around;
  times_around = 0;
  current_place = loop_start;
  size = 1;
  x = x1;
  y = y1;
  helices_on_loop = 0;
  if (g_angle_list) {
    printf("Traversing loop %d from base %d to base %d\n", loop_num,
	   loop_start, loop_end);
  }
  if (current_place == 0) {
    if ((g_angle_list) && (g_start_base == 1)) {
      printf("5 prime_end\n");
    }
    current_place = 1;
    if (!g_is_circular)
      size = 2;
  }
  while (current_place != loop_end) {
    if (current_place > largest_base) {
      if (times_around == 0) {
	times_around++;
	if (g_is_circular)
	  current_place = g_start_base;
	else
	  current_place = g_start_base - 1;
      } else {
	if (error("pk", loop_start, loop_end)==1)
	  try_exit(9);
      }
    }
    connected_to = connected[current_place];
    if ((connected_to > current_place) && (connected_to <= largest_base)) {
      if (g_angle_list) {
	printf("%d * %d\n", current_place, connected_to);
	if (current_place > loop_start) {
	  push_stack(loop_num);
	}
      }
      if (g_midpoints) {
	helices_on_loop++;
	if (current_place > loop_start) {
	  set_midpoints_no_structure(current_place, connected_to, &x, &y);
	}
	g_midpoint_x[helices_on_loop] = x;
	g_midpoint_y[helices_on_loop] = y;
      }
      if (extra_base_exists) {
	current_place = connected_to;
	push_stack(loop_num);
      }
      if (size > 1) {	/* set radius of bases on list */
	current_place = connected_to;
      } else {
	current_place++;
	if (current_place > largest_base) {
	  if (times_around == 0) {
	    times_around++;
	    if (g_is_circular)
	      current_place = g_start_base;
	    else
	      current_place = g_start_base - 1;
	  } else {
	    if (error("pk", loop_start, loop_end)==1)
	      try_exit(10);
	  }
	}
      }
    } else {
      if (g_angle_list)
	printf("%d\n", current_place);
      current_place++;
      if (current_place > largest_base) {
	if (times_around == 0) {
	  times_around++;
	  if (g_is_circular)
	    current_place = g_start_base;
	  else
	    current_place = g_start_base - 1;
	} else {
	  if (error("pk", loop_start, loop_end)==1)
	    try_exit(11);
	}
      }
    }
    size++;
  }
  if (extra_base_exists) {
    if (g_angle_list)
      printf("%d\n", current_place);
    size++;
  }
  if (g_angle_list) {
    printf("%d\n", loop_end);
    printf("Loop Size: %d\n", size);
  }
  sum_x = 0;
  sum_y = 0;
  for (i = 1; i <= helices_on_loop; i++) {
    sum_x = sum_x + g_midpoint_x[i];
    sum_y = sum_y + g_midpoint_y[i];
  }
  sum_x = sum_x / (float) helices_on_loop;
  sum_y = sum_y / (float) helices_on_loop;
  /* geometric mean of all base pairs around loop */
  set_color_ps(COLOR_CENTERS);
  if ((loop_num > 1) || (helices_on_loop > 1)) {
    for (i = 1; i <= helices_on_loop; i++) {
      fprintf(outfp, "%.2f %.2f m\n", sum_x, sum_y);
      fprintf(outfp, "%.2f %.2f l\n",
	      g_midpoint_x[i], g_midpoint_y[i]);
      fprintf(outfp, " s\n");
    }
  }
  if (!g_smart_colors)
    set_color_ps(COLOR_CONNECTING_GC);
}

int check_multiloop(int start_base, int end_base) {
  int             bases_on_list;
  int             i;
  int             helices_on_loop;
  int             current_base;
  current_base = start_base;
  if (g_connected[current_base] == 0) {
    bases_on_list = 1;
    g_base_list[1] = current_base;
    helices_on_loop = 0;
  } else {
    bases_on_list = 0;
    helices_on_loop = 1;
  }
  current_base++;
  while (current_base != end_base) {
    if (g_connected[current_base] > current_base) {	
      /* start of new helix */
      bases_on_list++;
      g_base_list[bases_on_list] = current_base;
      current_base = g_connected[current_base];
      helices_on_loop++;
    } else {	/* unpaired base */
      if (g_connected[current_base] < 0) {
	bases_on_list++;
	g_base_list[bases_on_list] = current_base;
      }
      current_base++;
      if (current_base > g_length)
	current_base = g_start_base;
    }
  }
  if (helices_on_loop <= 2)
    return 0;
  bases_on_list++;
  g_base_list[bases_on_list] = end_base;
  for (i = 1; i <= bases_on_list; i++) {
    g_base_advance[g_base_list[i]] = 1.0;
  }
  return bases_on_list;
}

float set_circle_angles(float degrees_to_use, int is_circular,
			int whole_circle_flag, int start_base, int end_base,
			int ex_domain_flag, int domain_row, int domain_col,
			float ex_domain_whole_circle_perp_angle) {
  /* start_base, end_base refer to the start and
   * end of the sequence from command line 
   * domain_row, domain_col refer to bounds of interior, exterior domain
   * ex_domain_whole_circle_perp_angle is the angle of the exterior domain
   * The bases must be drawn on each side of it to preserve the angle (MZ ??)
   */
  int             i;
  int             normal_bases;
  int             abnormal_bases;
  int             connected_to;
  float           start_rad;
  float           current_total_angle;
  float           radians_to_use;
  float           base_advance_angle;
  start_rad = 3 * PI / 2.;
  if (is_circular)
    radians_to_use = 2.0*PI - (2.0*PI/(end_base-start_base+1));
  else
    radians_to_use = 2.0 * PI - 4 * PI / 180.;
  if (degrees_to_use > 0)
    radians_to_use = degrees_to_use * PI / 180.;
  printf("Using %.1f degrees of the circle \n",
	 (radians_to_use) * 180. / PI);
  if ((g_external_domain_flag) && (whole_circle_flag)) {
    /*    May18, 2007. M. Zuker. This must be an  error!
	  if (!ex_domain_flag);
    */
    if (!ex_domain_flag) {
      domain_row = g_ex_start_base;
      domain_col = g_ex_end_base;
      ex_domain_flag = TRUE;
      ex_domain_whole_circle_perp_angle = 3. * PI / 2. -
	((g_ex_start_base + g_ex_end_base) / 2.) *
	radians_to_use / (g_end_base - g_start_base);
    }
  }
  if ((ex_domain_flag) && (whole_circle_flag)) {
    base_advance_angle = radians_to_use /
      (float) (g_length - domain_col + domain_row );
    if (g_counter_clockwise) {
      base_advance_angle *= -1;
    }
    while (ex_domain_whole_circle_perp_angle < PI / -2.0)
      ex_domain_whole_circle_perp_angle += 2.0 * PI;
    g_base_advance[domain_row] = ex_domain_whole_circle_perp_angle +
      base_advance_angle / 2;
    for (i = (domain_row - 1); i >= 1; i--) {
      g_base_advance[i] = g_base_advance[i + 1] + base_advance_angle;
    }
    g_base_advance[domain_col] = g_base_advance[domain_row] -
      base_advance_angle;
    for (i = (domain_col + 1); i <= g_length; i++) {
      g_base_advance[i] = g_base_advance[i - 1] - base_advance_angle;
    }
    return (start_rad - radians_to_use);	/* return value is not
						 * used */
  }
  if (!g_small_angle) {
    if (!whole_circle_flag) {
      base_advance_angle = radians_to_use / (float) (g_length - 1);
      if (g_counter_clockwise) {
	base_advance_angle *= -1;
      }
      g_base_advance[1] = start_rad;
      for (i = 2; i <= g_length; i++)
	g_base_advance[i] = start_rad - base_advance_angle * (i - 1);
    } else {
      base_advance_angle = radians_to_use /
	(float) (end_base - start_base);
      if (g_counter_clockwise) {
	base_advance_angle *= -1;
      }
      g_base_advance[start_base] = start_rad;
      for (i = (start_base + 1); i <= end_base; i++)
	g_base_advance[i] = start_rad - base_advance_angle *
	  (i - start_base);
    }
    if (g_counter_clockwise)
      return (start_rad + radians_to_use);
    else
      return (start_rad - radians_to_use);
  } else {
    for (i = 1; i <= g_length; i++)
      g_base_advance[i] = 0.0;
  }
  normal_bases = 0;
  /* place normal advance before and after any helix
   * 1 is big normal advance, 0 is small advance
   * initialize all to zero
   * only bases around a multiloop are assigned 1 
   */
  i = start_base;
  if (g_connected[i] > i) {	/* treat special case when start with
				 * a base pair in a multiloop */ 
    normal_bases = normal_bases + check_multiloop(g_connected[i], end_base);
    if (normal_bases > 0)
      normal_bases += 2;	/* add 2 for i,g_end_base */
  }
  while (i <= end_base) {
    /* find ends of helices */
    if (g_connected[i] > i) {	/* start of helix */
      /* advance to end of helix */
      while (g_connected[i + 1] == (g_connected[i] - 1)) {
	i++;
      }
      /* i is now the last base pair of a helix */
      connected_to = g_connected[i];
      /* connected_to is closing base pair of helix */
      normal_bases = normal_bases + check_multiloop(i, connected_to);
    }
    i++;
  }
  abnormal_bases = end_base - start_base + 1 - normal_bases;
  if (g_base_advance[1] == 0)
    base_advance_angle = radians_to_use /
      (normal_bases + g_helix_base_advance * (abnormal_bases - 1));
  else
    base_advance_angle = radians_to_use /
      (normal_bases - 1 + g_helix_base_advance * abnormal_bases);
  /* normal_base has normal advance, for multiloop
   * abnormal base has advance specified with -z */
  current_total_angle = start_rad;
  if (whole_circle_flag)
    g_base_advance[start_base] = start_rad;
  else
    g_base_advance[1] = start_rad;
  if (whole_circle_flag) {
    for (i = (start_base + 1); i <= end_base; i++) {
      if (g_base_advance[i] == 0.0) {
	current_total_angle = current_total_angle -
	  base_advance_angle * g_helix_base_advance;
      } else {
	current_total_angle = current_total_angle - base_advance_angle;
      }
      g_base_advance[i] = current_total_angle;
    }
  } else {
    for (i = 2; i <= g_length; i++) {
      if (g_base_advance[i] == 0.0) {
	current_total_angle = current_total_angle -
	  base_advance_angle * g_helix_base_advance;
      } else {
	current_total_angle = current_total_angle - base_advance_angle;
      }
      g_base_advance[i] = current_total_angle;
    }
  }
  return (start_rad - radians_to_use);
}

void circle_graph_start(char *first_line) {
  int             len;
  float           start;
  float           small_radius;
  g_circle_radius_point = g_diam / 2 * 72;
  /* display first line from ct file */
  len = strlen(first_line);
  start = (WD + 72.0)/2.0 - 0.5 * len * 13.;
  if (start < 5.)
    start = 5.;
  if (g_midpoints) {
    small_radius = 12.0 / sqrt(g_end_base - g_start_base);
    fprintf(outfp, "%s\n", "%%");
    fprintf(outfp, "%s for -m option, r is radius of dot for base pair on connecting arc\n", "%%");
    fprintf(outfp, "%s try adjusting from %f\n", "%%", small_radius);
    fprintf(outfp, "/r {%f} def\n", small_radius);
    fprintf(outfp, "%s\n", "%%");
  }
  set_color_ps(COLOR_TEXT);
  fprintf(outfp, "/Helvetica findfont\n");
  fprintf(outfp, "18 scalefont\n");
  fprintf(outfp, "setfont\n");
  fprintf(outfp, "%f 50 m\n", start);
  fprintf(outfp, " (%s) show\n", first_line);
}

void draw_bases_around_circle(float last_circle_angle, int is_circular, 
			      char *bases) {
  float           font_size;
  float           font_from_circle;
  int             step, label, i;
  float           angle, extra;
  float           xpos, ypos;
  last_circle_angle = last_circle_angle * 180. / PI;
  set_color_ps(COLOR_MAIN_CIRCLE);
  fprintf(outfp, " 1 setlinewidth\n");
  fprintf(outfp, " %f %f translate\n", (WD + 72.0)/2.0, (HT + 72.0)/2.0);
  fprintf(outfp, " s\n");
  if (is_circular)
    fprintf(outfp, " 0 0 %f -90. 270. a\n", g_circle_radius_point);
  else if (!g_counter_clockwise)
    fprintf(outfp, " 0 0 %f %.2f %.2f a\n", g_circle_radius_point,
	    last_circle_angle, 270.);
  else
    fprintf(outfp, " 0 0 %f %.2f %.2f a\n", g_circle_radius_point,
	    270., last_circle_angle);
  fprintf(outfp, " s\n");
  /* draw bases around circle */
  fprintf(outfp, "newpath\n");
  fprintf(outfp, "%s Draw bases around circle\n", "%%");
  font_size = 12;
  fprintf(outfp, "/Helvetica findfont\n");
  fprintf(outfp, "%f scalefont\n", font_size);
  fprintf(outfp, "setfont\n");
  /* create labels */
  fprintf(outfp, " 1 setlinewidth\n");
  if (g_whole_circle) {
    draw_label(g_start_base, g_circle_radius_point);
    draw_label(g_end_base, g_circle_radius_point);
    step = step_fun(g_start_base, g_end_base);
    
  } else {
    draw_label(1, g_circle_radius_point);
    draw_label(g_length, g_circle_radius_point);
    step = step_fun(1, g_length);
  }
  if (g_whole_circle) {
    for (label = g_start_base - 1 + step; label < g_end_base - 
	   2*step/3;label+=step)
      draw_label(label, g_circle_radius_point);
  } else {
    for (label = step; label < (g_length - 2 * step / 3); label += step)
      draw_label(label, g_circle_radius_point);
  }
  fprintf(outfp, "%s Draw bases around circle\n", "%%");
  if (g_whole_circle) {
    font_size = 15*72./(g_end_base - g_start_base)*g_circle_radius_point/220.;
  } else
    font_size = 15 * 72. / (g_length - step) * g_circle_radius_point / 220.;
  if (font_size > 11)
    font_size = 11;
  fprintf(outfp, "/Helvetica findfont\n");
  fprintf(outfp, "%f scalefont\n", font_size);
  fprintf(outfp, "setfont\n");
     font_from_circle = 4 * font_size;
     fprintf(outfp, " %f setlinewidth\n", font_from_circle / 20);
  if (font_from_circle < 5)
    font_from_circle = 5;
  if (font_from_circle > 20)
    font_from_circle = 20;
  set_color_ps(COLOR_BASES);
  for (i = g_start_base; i <= g_end_base; i++) {
    angle = g_base_advance[i];
    if (angle <= (-PI / 2.0) + .00001)
      extra = 1.;
    else
      extra = 0.;
    xpos = (g_circle_radius_point + font_from_circle * (1 + .6 * extra))
      * cos(angle) - font_size / 3;
    ypos = (g_circle_radius_point + font_from_circle * (1 + .6 * extra))
      * sin(angle) - font_size / 3;
    fprintf(outfp, "%.2f %.2f m\n", xpos, ypos);
    fprintf(outfp, "(%c) show\n", bases[i]);
    xpos = (g_circle_radius_point + .6 * font_from_circle) * cos(angle);
    ypos = (g_circle_radius_point + .6 * font_from_circle) * sin(angle);
    fprintf(outfp, "%.2f %.2f m\n", xpos, ypos);
    xpos = (g_circle_radius_point + .2 * font_from_circle) * cos(angle);
    ypos = (g_circle_radius_point + .2 * font_from_circle) * sin(angle);
    fprintf(outfp, "%.2f %.2f l\n", xpos, ypos);
    fprintf(outfp, "s\n");
  }
}

void draw_circle_graph_arc_no_pi(float start_angle, float end_angle, 
				 int phi_greater_than_pi, float *last_point_x,
				 float *last_point_y, int *point_exists,
				 float phi, float gamma) {
  float           R, dist, X0, Y0, midpoint_angle, midpoint_x, midpoint_y;
  float           advancing_angle, x, y;
  int             j;
  R = tan(phi) * g_circle_radius_point;
  dist = sqrt(R * R + g_circle_radius_point * g_circle_radius_point);
  X0 = dist * cos(gamma);
  Y0 = dist * sin(gamma);
  /* June 17, 2007. M. Zuker adapts to 13x17 
     if (((X0 * X0) < (305. * 305.)) && ((Y0 * Y0) < (400. * 400.))) {
  */
  if (((X0 * X0) < ((WD + 72.0)/2.0 * (WD + 72.0)/2.0)) && 
      ((Y0 * Y0) < ((HT + 72.0)/2.0 * (HT + 72.0)/2.0))) {
    /* center of circle is on page */
    if (g_arcs) {	/* draw the arc */
      fprintf(outfp, "%.2f %.2f %.2f %.2f %.2f a\n",
	      X0, Y0, R, start_angle, end_angle);
      fprintf(outfp, " s\n");
    }
    if (g_midpoints) {	/* draw the midpoint */
      if (g_arcs)
	set_color_ps(COLOR_CENTERS);
      midpoint_angle = (start_angle + end_angle) / 2.0;
      if (!phi_greater_than_pi)
	midpoint_angle = midpoint_angle + 180;
      midpoint_x = X0 + cos(midpoint_angle * PI / 180) * R;
      midpoint_y = Y0 + sin(midpoint_angle * PI / 180) * R;
      fprintf(outfp, "%.2f %.2f r %.2f %.2f a\n",
	      midpoint_x, midpoint_y, 0.0, 360.0);
      fprintf(outfp, " s\n");
      if (*point_exists) {	/* connect to previous
				 * midpoint of helix */
	if (g_smart_colors)
	  set_color_ps(COLOR_CENTERS);
	fprintf(outfp, "%.2f %2.f m\n", *last_point_x, *last_point_y);
	fprintf(outfp, "%.2f %2.f l\n", midpoint_x, midpoint_y);
	fprintf(outfp, " s\n");
	
      } else
	*point_exists = TRUE;
      *last_point_x = midpoint_x;
      *last_point_y = midpoint_y;
      if ((g_arcs) && (!g_smart_colors))
	set_color_ps(COLOR_CONNECTING_GC);
    }
  } else {		/* Draw the circle the hard way. draw 20 lines
			   & convert angle back to radians */
    if (end_angle < 0)
      end_angle = end_angle + 360;
    if (end_angle < start_angle) {
      end_angle = end_angle + 360;
    }
    start_angle = start_angle * PI / 180.;
    end_angle = end_angle * PI / 180.;
    if (g_arcs) {
      advancing_angle = (end_angle - start_angle) / 20.0;
      x = X0 + R * cos(start_angle);
      y = Y0 + R * sin(start_angle);
      fprintf(outfp, "%.2f %.2f m\n", x, y);
      for (j = 1; j <= 20; j++) {
	x = X0 + R * (cos(start_angle + j * advancing_angle));
	y = Y0 + R * (sin(start_angle + j * advancing_angle));
	fprintf(outfp, "%.2f %.2f l\n", x, y);
      }
      fprintf(outfp, "s\n");
    }
    if (g_midpoints) {
      if (g_arcs)
	set_color_ps(COLOR_CENTERS);
      midpoint_angle = (start_angle + end_angle) / 2.0;
      midpoint_x = X0 + cos(midpoint_angle) * R;
      midpoint_y = Y0 + sin(midpoint_angle) * R;
      fprintf(outfp, "%.2f %.2f r %.2f %.2f a\n",
	      midpoint_x, midpoint_y, 0.0, 360.0);
      fprintf(outfp, " s\n");
      if (*point_exists) {
	if (g_smart_colors)
	  set_color_ps(COLOR_CENTERS);
	fprintf(outfp, "%.2f %2.f m\n",
		*last_point_x, *last_point_y);
	fprintf(outfp, "%.2f %2.f l\n", midpoint_x, midpoint_y);
	fprintf(outfp, " s\n");
      } else
	*point_exists = TRUE;
      *last_point_x = midpoint_x;
      *last_point_y = midpoint_y;
      if ((g_arcs) && (!g_smart_colors))
	set_color_ps(COLOR_CONNECTING_GC);
    }
  }
}
 
void draw_circle_graph_arc_near_pi(float angle1, float angle2, 
				   int *point_exists, float *last_point_x,
				   float *last_point_y) {
  float           x, y, x2, y2;
  float           midpoint_x, midpoint_y;
  x = g_circle_radius_point * cos(angle1);
  y = g_circle_radius_point * sin(angle1);
  x2 = g_circle_radius_point * cos(angle2);
  y2 = g_circle_radius_point * sin(angle2);
  if (g_arcs) {
    fprintf(outfp, "%.2f %.2f m\n", x, y);
    fprintf(outfp, "%.2f %.2f l\n", x2, y2);
    fprintf(outfp, "s\n");
  }
  if (g_midpoints) {
    if (g_arcs)
      set_color_ps(COLOR_CENTERS);
    midpoint_x = (x + x2) / 2.0;
    midpoint_y = (y + y2) / 2.0;
    fprintf(outfp, "%.2f %.2f r %.2f %.2f a\n",
	    midpoint_x, midpoint_y, 0.0, 360.0);
    fprintf(outfp, " s\n");
    if (*point_exists) {
      fprintf(outfp, "%.2f %2.f m\n", *last_point_x, *last_point_y);
      fprintf(outfp, "%.2f %2.f l\n", midpoint_x, midpoint_y);
      fprintf(outfp, " s\n");
    } else
      *point_exists = TRUE;
    *last_point_x = midpoint_x;
    *last_point_y = midpoint_y;
    if ((g_arcs) && (!g_smart_colors))
      set_color_ps(COLOR_CONNECTING_GC);
  }
}

void draw_arc(int i, int connected_to, int *midpoint_exists, 
	      float *last_midpoint_x, float *last_midpoint_y, 
	      int helix_length) {
  float    phi, gamma, angle1, angle2, last_angle_perp, print_angle_perp;
  float    print_angle, end_angle, start_angle;
  int      phi_greater_than_pi;
  angle1 = g_base_advance[i];
  angle2 = g_base_advance[connected_to];
  fix_angles(&angle1, &angle2);
  phi = (angle1 - angle2);
  while (phi >= (2.0 * PI))
    phi = phi - 2.0 * PI;
  if (((phi - PI) * (phi - PI)) < .00005) {	/* draw line */
    if (g_angle_list) {
      print_angle = angle1 * 180. / PI;
      print_angle = print_angle + 180.;
      if (print_angle < 0)
	print_angle = print_angle + 360.;
      if (print_angle > 360)
	print_angle = print_angle - 360;
      print_angle_perp = print_angle + 90.;
      if (print_angle_perp > 360)
	print_angle_perp = print_angle_perp - 360;
      if (helix_length == 1)
	last_angle_perp = print_angle_perp * PI / 180.;
      if (g_angle_list) {
	if (helix_length == 1)
	  printf(" Tangent %.2f  Perpendicular %.2f", print_angle,
		 print_angle_perp);
	printf("%d %d \n", i, connected_to);
      }
    }
    draw_circle_graph_arc_near_pi(angle1, angle2, midpoint_exists,
				  last_midpoint_x, last_midpoint_y);
    
  } else {		/* phi is not close to pi */
    if (phi > PI) {	/* phi > pi */
      phi_greater_than_pi = TRUE;
      phi = (2 * PI - phi);
      phi = phi / 2.0;
      gamma = angle2 - phi;
      if (g_angle_list) {
	print_angle = (gamma + PI / 2.0) * 180. / PI;
	if (print_angle < 0)
	  print_angle = print_angle + 360.;
	print_angle_perp = print_angle + 90.;
	if (print_angle_perp > 360)
	  print_angle_perp = print_angle_perp - 360;
	if (helix_length == 1) {
	  last_angle_perp = print_angle_perp;
	  last_angle_perp = last_angle_perp * PI / 180;
	}
	if (g_angle_list) {
	  if (helix_length == 1)
	    printf(" Tangent %.2f  Perpendicular %.2f", print_angle,
		   print_angle_perp);
	  printf("%d %d\n",
		 i, connected_to);
	}
      }
      end_angle = angle1 * 180. / PI - 90;
      start_angle = angle2 * 180. / PI + 90;
    } else {	/* phi < pi */
      phi = phi / 2.0;
      phi_greater_than_pi = FALSE;
      gamma = angle2 + phi;
      if (g_angle_list) {
	print_angle = (gamma + PI / 2.0) * 180. / PI;
	print_angle = print_angle + 180.;
	if (print_angle > 360.)
	  print_angle = print_angle - 360.;
	if (print_angle < 0)
	  print_angle = print_angle + 360.;	/* 360 was 180 */
	print_angle_perp = print_angle + 90;
	if (phi < 0)
	  print_angle_perp = print_angle_perp + 180;
	if (print_angle_perp > 360)
	  print_angle_perp = print_angle_perp - 360;
	if (helix_length == 1) {
	  last_angle_perp = print_angle_perp;
	  last_angle_perp = last_angle_perp * PI / 180.;
	}
	if (g_angle_list) {
	  if (helix_length == 1)
	    printf(" Tangent %.2f  Perpendicular %.2f", print_angle,
		   print_angle_perp);
	  printf("%d %d\n",
		 i, connected_to);
	}
      }
      start_angle = angle1 * 180. / PI + 90.;
      end_angle = angle2 * 180. / PI - 90.;
    }
    draw_circle_graph_arc_no_pi(start_angle, end_angle, phi_greater_than_pi,
				last_midpoint_x, last_midpoint_y, 
				midpoint_exists, phi, gamma);
  }
}

void main_data_ps_no_structure(float degrees_used, int is_circular,
			       char *first_line, char *bases, 
			       int *connected, int start_base_in, 
			       int end_base_in, int whole_circle_flag,
			       float ex_domain_whole_perp_angle) {
  /* This function draws the circle graph in
   * postscript format.
   * - degrees_used is the total number of degrees that are used in the
   *   circle graph
   * - is_circular should be TRUE for circular rna
   * - first line is the information on first line of a ct file
   *   (without the sequence length)
   * - bases is the array from 1 to length of sequence, contains A,C,G,U
   *   or other character for each base.
   * - connected indictes what each base is connected to. 0 means unpaired
   * - start_base_in and end_base_in indicate what included, excluded
   *   domain to  redraw from the screen 
   * - new_perp_angle, returns the angle of the last helix draw    
   * - whole circle flag: true indictes use of the whole circle for the
   *   domain. false indicates use of the whole circle for the entire sequence
   * - ex_domain_whole_perp_angle for excluded domain from screen draw bases
   *   on each side of this angle in order to preserve original angle.
   *   applies only to whole_circle_flag option 
   * - domain degrees is really degrees to use for the defined domain
   */
  int             current_helix;
  int             helix_started;
  int             helix_length;
  int             last_connected_to;
  int             total_loops;
  int             start_base, end_base;
  float           last_circle_angle;
  int             connected_to;
  float           junk_x, junk_y;
  int             i;
  int             midpoint_exists;
  float           last_midpoint_x; /* x coordinate of last midpoint */
  float           last_midpoint_y; /* y coordinate of last midpoint */
  int             helix_loop;
  int             extra_base_flag;
  int             domain_row = -9999, domain_col = -9999;
  start_base = start_base_in;
  end_base = end_base_in;
  midpoint_exists = FALSE;
  g_circle_radius_point = g_diam / 2 * 72;
  /* show first line from ct file */
  circle_graph_start(first_line);
  if (g_external_domain_flag) {
    total_loops = 1;
    domain_row = start_base;
    domain_col = end_base;
    start_base = g_start_base;
    end_base = g_end_base;
    g_small_angle = FALSE;	/* do not use z here */
  } else {		/* for normal execution without domains */
    total_loops = 1;
  }
  /* display first line from ct file
   * draw big circle
   * draw bases around circle */
  last_circle_angle =
    set_circle_angles(degrees_used, is_circular, whole_circle_flag, start_base,
		      end_base, FALSE, start_base_in, end_base_in,
		      ex_domain_whole_perp_angle);
  draw_bases_around_circle(last_circle_angle, is_circular, bases);
  /* create connecting lines */
  fprintf(outfp, "%s Draw connecting lines\n", "%%");
  if (g_arcs)
    set_color_ps(COLOR_CONNECTING_GC);
  else
    set_color_ps(COLOR_CENTERS);
  current_helix = 0;
  helix_started = FALSE;
  helix_length = 0;
  last_connected_to = -10;
  /* find first base_pair */
  i = start_base;
  connected[0] = 0;
  if (connected[i] > i) {
    if ((g_midpoints) || (g_angle_list)) {
      set_midpoints_no_structure(i, connected[i], &junk_x, &junk_y);
      if (g_end_base == connected[i])
	extra_base_flag = FALSE;
      else
	extra_base_flag = TRUE;
      if ((g_angle_list) || (g_midpoints)) {
	traverse_loop_no_structure(total_loops, i, g_end_base, connected,
				   junk_x, junk_y, extra_base_flag, g_length);
      }
      if (connected[i + 1] < 1) { /* when structure begins with a
				   * single base pair, put some extra
				   * on the stack */ 
	push_stack(1);
      }
    }
  } else {
    while (connected[i] < i)
      i++;
    if ((g_angle_list) || (g_midpoints))
      set_midpoints_no_structure(i, connected[i], &junk_x, &junk_y);
    if ((g_angle_list) || (g_midpoints))
      traverse_loop_no_structure(total_loops, i, i-1, connected, junk_x,
				 junk_y, TRUE, end_base);
  }
  total_loops++;
  if (g_external_domain_flag) {
    domain_row = g_ex_start_base;
    domain_col = g_ex_end_base;
  }
  /* Start of main loop to go through bases 
   * though dots of a helix 
   * g_previous_loop is the previous loop for a base pair
   * It is the current loop for a single stranded base */
  while (i <= end_base) {
    if (g_external_domain_flag) {
      if ((i > domain_row) && (i < domain_col)) {
	i = domain_col;
	helix_started = FALSE;
	helix_length = 0;
      }
    }
    connected_to = connected[i];
    if ((connected_to > i) && (connected_to <= g_end_base)) {  
      g_previous_loop[i] = total_loops - 1;
      if (helix_started) {
	if (connected_to == (last_connected_to - 1)) {
	  helix_length++;
	  /* continuing a helix (a new helix with a bulge loop in between) */
	} else {	
	  if (g_angle_list)
	    printf("Helix Length: %d\n", helix_length);
	  if ((g_angle_list) || (g_midpoints)) {
	    traverse_loop_no_structure(total_loops, i - 1, connected[i - 1],
				       connected, last_midpoint_x, 
				       last_midpoint_y, FALSE, g_length);
	  }
	  total_loops++;
	  helix_length = 1;
	  helix_started = TRUE;
	  current_helix++;
	  if (g_angle_list) {
	    helix_loop = pop_stack();
	    if (g_angle_list)
	      printf("Helix: %d of Loop: %d\n", current_helix, helix_loop);
	  }
	  midpoint_exists = FALSE;
	}
      } else {	/* start new helix */
	helix_length = 1;
	current_helix++;
	helix_started = TRUE;
	if (g_angle_list) {
	  helix_loop = pop_stack();
	  if (g_angle_list)
	    printf("Helix: %d of loop: %d\n", current_helix, helix_loop);
	}
	midpoint_exists = FALSE;
      }
      /* draw arc of circle graph */
      if (g_smart_colors)
	set_color_base_ps(bases[i], bases[connected_to]);
      else {
	if (g_annotation)
	  set_color_ann_ps((int) ((g_ann_to_color[i] + 
				   g_ann_to_color[connected_to]-1)/2.));
	/* The -1 above added April 17, 2000 */
      }
      draw_arc(i, connected_to, &midpoint_exists, &last_midpoint_x,
	       &last_midpoint_y, helix_length);
    }
    /* end of if(connected_to > i) single stranded base or column 
     * side of base pair */
    else {		
      if (helix_started) {
	helix_started = FALSE;
	if (g_angle_list)
	  printf("Helix Length: %d\n", helix_length);
	if ((g_angle_list) || (g_midpoints)) {
	  traverse_loop_no_structure(total_loops, i - 1, connected[i - 1],
				     connected, last_midpoint_x, 
				     last_midpoint_y, FALSE, g_length);
	  total_loops++;
	}
	helix_length = 0;
      }
    }
    g_previous_loop[i] = total_loops - 1;
    last_connected_to = connected_to;
    i++;		/* end of while loop */
  }
  if (helix_started) {	/* probably never executed somehow end with helix */
    if (g_angle_list)
      printf("Helix length = %d\n", helix_length);
    if ((g_angle_list) || (g_midpoints)) {
      traverse_loop_no_structure(total_loops, i-1, connected[i-1], connected,
				 last_midpoint_x, last_midpoint_y, FALSE, 
				 g_length);
      total_loops++;
    }
  }
  g_total_loops = total_loops;	/* loop adjust */
}

int main_data_ps_structure(float degrees_used, int is_circular, int *connected,
			   int in_domain_flag, int start_base_in, 
			   int end_base_in, float *new_perp_angle,
			   int whole_circle_flag, int ex_domain_flag, 
			   float ex_domain_whole_perp_angle, 
			   int domain_degrees) {

  /* This function computes the structure for all output formats It
   * does not produce the postscript format circle graph It is similar
   * to main_data_ps_no_structure which does though degrees_used is
   * total degrees of circle within circle graph to use is_circular
   * should be TRUE for circular rna first line is the information on
   * first line of ct file, length has been removed.  bases is the
   * array from 1 to length of sequence, contains A,C,G,U for each
   * base.  connected indictes what each base is connected to
   * in_domain_flag , used to indicate redraw included domain from
   * screen does not apply to included domain from command line
   * g_start_base, g_end_base are used for command line domain
   * start_base_in, end_base_in, indicate what included,excluded
   * domain to redraw from the screen new_perp_angle, returns the
   * angle of the last helix draw whole circle flag , TRUE indicates
   * the use of the whole circle for the domain FALSE indicates "Use
   * the whole circle for the entire sequence ex_domain_flag, used to
   * indicate to redraw excluded domain from screen does not to
   * external domain specified from command line
   * ex_domain_whole_perp_angle for excluded domain from screen draw
   * bases on each side of this angle in order to preserve original
   * angle applies only to whole_circle_flag option domain degrees is
   * really degrees to use for the defined domain
   */
  int   current_helix, helix_started, helix_length, last_connected_to;
  int   total_loops, start_base, end_base, connected_to, i, helix_loop;
  int   extra_loop_flag, extra_base_flag, domain_row = -9999;
  int   domain_col = -9999, i_loop_start, closing_pair;
  float print_angle, print_angle_perp, current_position_x = -9999.;
  float current_position_y = -9999., last_angle_perp = 0.0, distance_bases;
  float distance_base_pairs;

  start_base = start_base_in;
  end_base = end_base_in;
  if (ex_domain_flag)
    degrees_used = (float) domain_degrees;
  if (in_domain_flag) {
    degrees_used = (float) domain_degrees;	/* use full circle */
    g_small_angle = FALSE;	/* do not use z here */
    total_loops = g_previous_loop[start_base];
    distance_bases = DISTANCE_BASES * g_scale;
    distance_base_pairs = DISTANCE_BASE_PAIRS * g_scale;
  } else {
    if ((ex_domain_flag) || (g_external_domain_flag)) {
      total_loops = 1;
      domain_row = start_base;
      domain_col = end_base;
      start_base = g_start_base;
      end_base = g_end_base;
      g_small_angle = FALSE;	/* do not use z here */
      distance_bases = DISTANCE_BASES * g_scale;
      distance_base_pairs = DISTANCE_BASE_PAIRS * g_scale;
    } else {	/* for normal execution without domains */
      total_loops = 1;
      distance_bases = DISTANCE_BASES;
      distance_base_pairs = DISTANCE_BASE_PAIRS;
    }
  }
  set_circle_angles(degrees_used, is_circular, whole_circle_flag, start_base,
		    end_base, ex_domain_flag, start_base_in, end_base_in,
		    ex_domain_whole_perp_angle);
  current_helix = 0;
  helix_started = FALSE;
  helix_length = 0;
  last_connected_to = -10;
  /* find first base_pair */
  i = start_base;
  connected[0] = 0;
  if (connected[i] > i) {
    if (in_domain_flag) {
      extra_loop_flag = FALSE;
      helix_length = 1;
      helix_started = TRUE;
      total_loops = g_previous_loop[i];
      /* set helix_loop to the loop this helix is on */
      /* find start of loop by backing around it */
      last_connected_to = connected[i];
      i_loop_start = i;
      if (i_loop_start > g_start_base)
	i_loop_start--;
      current_helix = g_previous_loop[i];
      while ((i_loop_start >= g_start_base)
	     && (g_connected[i_loop_start] < last_connected_to)) {
	if ((g_connected[i_loop_start] > 0) &&
	    (g_connected[i_loop_start] < i_loop_start))
	  i_loop_start = g_connected[i_loop_start];
	else
	  i_loop_start--;
      }
      helix_loop = g_previous_loop[i_loop_start + 1];
      /* draw base pair for i */
      last_connected_to = connected[i];
      get_perpendicular_angle(i, connected[i],
			      &last_angle_perp);
      *new_perp_angle = last_angle_perp;
      current_position_x = (g_base_x[i] +
			    g_base_x[connected[i]]) / 2.;
      current_position_y = (g_base_y[i] +
			    g_base_y[connected[i]]) / 2.;
      /* draw the base pair */
      draw_basepair(i, last_connected_to, last_angle_perp, current_position_x,
		    current_position_y, distance_base_pairs);
      /* at this point the angle, and the position is known */
      i++;
    } else {
      get_perpendicular_angle(i, connected[i], &last_angle_perp);
      if (g_angle_list)
	printf("Helix: 1 of loop: 1 Tangent = %.2f Perpendicular = %.2f\n",
	       (last_angle_perp - PI / 2) * 180. / PI,
	       last_angle_perp * 180. / PI);
      if (g_end_base == connected[i])
	extra_base_flag = FALSE;
      else
	extra_base_flag = TRUE;
      traverse_loop_structure(total_loops, i, g_end_base, connected,
			      last_angle_perp - PI, 0.0, 0.0, 
			      extra_base_flag, g_length, 
			      distance_base_pairs, distance_bases, TRUE);
      if (g_counter_clockwise)
	last_angle_perp = g_structure_angle[i] + PI / 2;
      else
	last_angle_perp = g_structure_angle[i] - PI / 2;
      current_position_x = (g_base_x[i] + g_base_x[connected[i]])/2.;
      current_position_y = (g_base_y[i] + g_base_y[connected[i]])/2.;
      helix_started = TRUE;
      helix_length = 1;
      current_helix = 1;
      last_connected_to = connected[i];
      if (connected[i + 1] < 1) {	
	/* When a structure begins with a single base
	 * pair, put some extra on the stack */
	push_stack(1);
      }
      g_previous_loop[i] = total_loops;
      i++;
      extra_loop_flag = TRUE;
    }
  } else {
    extra_loop_flag = FALSE;
    i = g_start_base;
    /* find the angle of the last helix */
    closing_pair = g_end_base;
    while (connected[closing_pair] < 1)
      closing_pair--;
    get_perpendicular_angle(connected[closing_pair],
			    closing_pair, &last_angle_perp);
    if (g_reg_angle_flag) {
      last_angle_perp = regularize_angle(last_angle_perp);
    }
    traverse_loop_structure(total_loops, connected[closing_pair],
			    connected[closing_pair] - 1, connected,
			    last_angle_perp - PI, 0.0, 0.0, TRUE, 
			    end_base,distance_base_pairs,distance_bases, TRUE);
    if (!((i == g_start_base) && (g_end_base == connected[i])))
      extra_loop_flag = TRUE;
  }
  total_loops++;
  if ((g_external_domain_flag) && (!ex_domain_flag)) {
    domain_row = g_ex_start_base;
    domain_col = g_ex_end_base;
  }
  /* Start of main loop to go through bases
   * current_position_x, y are positions of dot for base pair
   * last_angle_perpendicular is the angle of the line passing
   * though dots of a helix
   * g_previous_loop is the previous loop for a base pair
   * It is the current loop for a single stranded base
   */
  while (i <= end_base) {
    if ((ex_domain_flag) || (g_external_domain_flag)) {
      if ((i > domain_row) && (i < domain_col)) {
	i = domain_col;
	if (ex_domain_flag)
	  total_loops = g_previous_loop[domain_col] + 1;
	helix_started = FALSE;
	helix_length = 0;
      }
    }
    connected_to = connected[i];
    if ((connected_to > i) && (connected_to <= g_end_base)) {
      /* this base is paired with a larger base */
      if (!ex_domain_flag)
	g_previous_loop[i] = total_loops - 1;
      if (helix_started) {
	if (connected_to == (last_connected_to - 1)) {
	  helix_length++;
	  /* continuing a helix */
	  if (g_counter_clockwise) {
	    current_position_x = current_position_x + distance_bases *
	      cos((double) (last_angle_perp));
	    current_position_y = current_position_y + distance_bases *
	      sin((double) (last_angle_perp));
	  } else {
	    current_position_x = current_position_x + distance_bases *
	      cos((double) last_angle_perp);
	    current_position_y = current_position_y + distance_bases *
	      sin((double) last_angle_perp);
	  }
	  draw_basepair(i, connected_to, last_angle_perp, current_position_x,
			current_position_y, distance_base_pairs);
	} else {	/* this is a new helix with
			 * no single stranded bases
			 * in between it and the previous helix 
			 * and the previous on the 5 prime side
			 */
	  if (g_angle_list)
	    printf("Helix Length: %d\n", helix_length);
	  traverse_loop_structure(total_loops, i-1, connected[i-1], connected,
				  last_angle_perp, current_position_x, 
				  current_position_y, FALSE, g_length, 
				  distance_base_pairs, distance_bases, FALSE);
	  total_loops++;
	  helix_length = 1;
	  helix_started = TRUE;
	  current_helix++;
	  helix_loop = pop_stack();
	  if (g_counter_clockwise)
	    last_angle_perp = g_structure_angle[i] + PI / 2.;
	  else
	    last_angle_perp = g_structure_angle[i] - PI / 2.;
	  current_position_x = (g_base_x[i] + g_base_x[connected[i]]) / 2.;
	  current_position_y = (g_base_y[i] + g_base_y[connected[i]]) / 2.;
	  if (last_angle_perp < PI / -2.0)
	    last_angle_perp += 2.0 * PI;
	  if (g_angle_list)
	    printf("Helix: %d of Loop: %d\n",
		   current_helix, helix_loop);
	}
      } else { /* start new helix */
	helix_length = 1;
	current_helix++;
	helix_started = TRUE;
	helix_loop = pop_stack();
	if (g_angle_list)
	  printf("Helix: %d of loop: %d\n", current_helix, helix_loop);
	if (g_counter_clockwise)
	  last_angle_perp = g_structure_angle[i] + PI / 2.;
	else
	  last_angle_perp = g_structure_angle[i] - PI / 2.;
	current_position_x = (g_base_x[i] + g_base_x[connected[i]]) / 2.;
	current_position_y = (g_base_y[i] + g_base_y[connected[i]]) / 2.;
	if (last_angle_perp < PI / -2.0)
	  last_angle_perp += 2.0 * PI;
      }
      if ((g_angle_list) && (helix_length == 1)) {
	print_angle_perp = last_angle_perp * 180. / PI;
	while (print_angle_perp > 270.)
	  print_angle_perp = print_angle_perp - 360.;
	while (print_angle_perp < -90)
	  print_angle_perp = print_angle_perp + 360;
	print_angle = print_angle_perp - 90.;
	printf("Tangent %.2f  Perpendicular %.2f", 
	       print_angle, print_angle_perp);
      }	/* still in if(connected_to>i) */
    }	/* end of if(connected_to > i) */
    else {	/* single stranded base or column side of base pair */
      if (helix_started) {
	helix_started = FALSE;
	if (g_angle_list)
	  printf("Helix Length: %d\n", helix_length);
	traverse_loop_structure(total_loops, i-1, connected[i-1], connected,
				last_angle_perp, current_position_x,
				current_position_y, FALSE, g_length,
				distance_base_pairs, distance_bases, FALSE);
	total_loops++;
	helix_length = 0;
      }
    }
    if (!ex_domain_flag)
      g_previous_loop[i] = total_loops - 1;
    last_connected_to = connected_to;
    i++;		/* end of while loop */
  }
  if (helix_started) {	/* probably never executed somehow end with helix */
    if (g_angle_list) 
      printf("Weird Helix length = %d\n", helix_length);
    printf("5 traverse loop with total_loops = %d\n", total_loops);
    traverse_loop_structure(total_loops, i - 1, connected[i - 1], connected,
			    last_angle_perp, current_position_x,
			    current_position_y, FALSE, g_length,
			    distance_base_pairs, distance_bases, FALSE);
    total_loops++;
  }
  if ((!in_domain_flag) && (!ex_domain_flag)) {
    g_total_loops = total_loops;	/* loop adjust */
    if (extra_loop_flag)
      g_total_loops--;
  }
  return total_loops;
}

int traverse_loop_natural_count_bases(int loop_num, int i, int connected_to,
				      int first_loop, int *bulge_flag, 
				      int *interior_flag) {
  /* (i,connected_to) is the first base pair of the loop */
  int             current_base;
  int             number_of_bases;
  int             times_around;
  int             bulge;
  int             left_strand, right_strand;
  int             extra_helices;	/* should be 1 for a bulge or
					 * interior loop */
  /* bulge flag set to -2 for left bulge of size 2
   * -1 for left bulge of size 1
   * 0 otherwise
   * 1 for right bulge of size 1
   * 2 for right bulge of size 2
   * hairpin       bulge on left side
   * left   A-C  15
   * G          right
   * 5   A-C  16     
   * interior_flag for interios loops with 1 base more on one side
   * than the other
   * interior_flag set to -1 for extra on left, +1 for right, 0
   * otherwise
   */
  number_of_bases = 2;
  extra_helices = 0;
  current_base = i;
  current_base++;
  bulge = 0;
  left_strand = 0;
  right_strand = 0;
  if ((g_connected[current_base] > 0) && (g_connected[connected_to - 1] < 1)) {
    bulge = 1;	/* bulge possible on right side */
  }
  if ((g_connected[current_base] < 1) && (g_connected[connected_to - 1] > 0)) {
    bulge = -1;	/* bulge possible on left side */
  }
  if (first_loop)
    times_around = 0;
  else
    times_around = 1;
  if (g_angle_list) {
    printf("Traversing loop %d from base %d to %d\n", loop_num,i,connected_to);
    printf("%d\n", i);
  }
  while (current_base != connected_to) {
    if (g_connected[current_base] > current_base) {
      if (g_angle_list)
	printf("%d * %d\n", current_base, g_connected[current_base]);
      current_base = g_connected[current_base];
      number_of_bases++;
      extra_helices++;
    } else {
      if (g_angle_list)
	printf("%d\n", current_base);
      current_base++;
      number_of_bases++;
      if (extra_helices == 0)
	left_strand++;
      else
	right_strand++;
      if (current_base > g_end_base) {
	if (times_around == 0) {
	  current_base = g_start_base;
	  times_around = 1;
	} else {
	  if (error("pk", i, connected_to)==1)
	    try_exit(12);
	}
      }
    }
  }
  if (g_angle_list) {
    printf("%d\n", current_base);
    printf("Loop Size: %d \n", number_of_bases);
  }
  *bulge_flag = 0;
  if ((bulge != 0) && ((!first_loop) && (extra_helices == 1))) {
    if (number_of_bases == 5)
      *bulge_flag = bulge;
    if (number_of_bases == 6)
      *bulge_flag = bulge * 2;
  }
  *interior_flag = 0;
  if ((bulge == 0) && ((!first_loop) && (extra_helices == 1))) {
    right_strand--;
    if ((left_strand - right_strand) == 1) {
      *interior_flag = -1;
    }
    if ((right_strand - left_strand) == 1) {
      *interior_flag = 1;
    }
  }
  if ((first_loop) && (!g_is_circular))
    number_of_bases++;
  return number_of_bases;
}

void traverse_loop_natural_assign_angles(int size, int i, int last_base, 
					 float current_x, float current_y, 
					 int first_loop, 
					 float angle_of_last_basepair, 
					 int loop_num,
					 float distance_base_pairs, 
					 int bulge_flag, int interior_flag) {
  /* assign positions and angles to all bases around the loop
   * size is the number of bases around the loop
   * the loop starts with i and ends with connected_to
   * the connected_to is not always connected to for the first loop
   */
  int             current_base;
  float           x, y, base_angle_increment, current_angle;
  float           base_pair_x, base_pair_y, original_angle;
  int             times_around, connected_to;
  if (g_counter_clockwise)
    angle_of_last_basepair += PI;
  current_base = i;
  x = current_x;
  y = current_y;
  if (angle_of_last_basepair < 0)
    angle_of_last_basepair = angle_of_last_basepair + 2.0 * PI;
  original_angle = angle_of_last_basepair;
  angle_of_last_basepair = angle_of_last_basepair + PI;
  if (angle_of_last_basepair > (3.0 * PI / 2.0))
    angle_of_last_basepair = angle_of_last_basepair - 2.0 * PI;
  if (bulge_flag != 0)
    base_angle_increment = 1.0 * PI / (float) (size - 2);
  else
    base_angle_increment = 2.0 * PI / (float) size;
  if (interior_flag < 0)
    base_angle_increment = 2.0 * PI / (float) (size + 1.);
  if (interior_flag > 0)
    base_angle_increment = 2.0 * PI / (float) (size - 1);
  if (g_counter_clockwise)
    base_angle_increment = base_angle_increment * -1.0;
  g_loop_radius[loop_num] = (float) size / (2.0 * PI)
    * distance_base_pairs;
  g_loop_center_x[loop_num] = x + cos(original_angle)*g_loop_radius[loop_num];
  g_loop_center_y[loop_num] = y + sin(original_angle)*g_loop_radius[loop_num];
  if (g_connected[i] == last_base) { /* the starting base pair of a
				      * helix is assigned before the
				      * call here */ 
    current_angle = angle_of_last_basepair - 1.5 * base_angle_increment;
    current_base++;
  } else {
    current_angle = angle_of_last_basepair;
  }
  if (bulge_flag > 0) {
    if (base_angle_increment > 0)
      current_angle = angle_of_last_basepair-PI+.5*base_angle_increment;
    else
      current_angle = angle_of_last_basepair+PI+.5*base_angle_increment;
  }
  if (first_loop) {
    times_around = 0;
  } else {
    times_around = 1;
  }
  while (current_base != last_base) {	/* not connected base */
    if (g_connected[current_base] <= current_base) {
      g_base_x[current_base] = g_loop_center_x[loop_num] +
	g_loop_radius[loop_num] * cos(current_angle);
      g_base_y[current_base] = g_loop_center_y[loop_num] +
	g_loop_radius[loop_num] * sin(current_angle);
      g_structure_angle[current_base] = current_angle;
      current_base++;
      current_angle = current_angle - base_angle_increment;
      if (current_base > g_end_base) {
	if (times_around == 0) {
	  times_around = 1;
	  current_base = g_start_base;
	} else {
	  current_base = last_base;
	}
      }
    } else {	/* draw base pair */
      connected_to = g_connected[current_base];
      current_angle = current_angle - .5 * base_angle_increment;
      if (g_counter_clockwise) {
	g_structure_angle[current_base] = current_angle - PI / 2.;
	g_structure_angle[connected_to] = current_angle + PI / 2.;
      } else {
	g_structure_angle[current_base] = current_angle + PI / 2.;
	g_structure_angle[connected_to] = current_angle - PI / 2.;
      }
      base_pair_x = g_loop_center_x[loop_num] +
	g_loop_radius[loop_num] * cos(current_angle);
      base_pair_y = g_loop_center_y[loop_num] +
	g_loop_radius[loop_num] * sin(current_angle);
      if (g_counter_clockwise) {
	g_base_x[current_base] = base_pair_x +
	  distance_base_pairs / 2.0 * (cos(PI / -2.0 + current_angle));
	g_base_y[current_base] = base_pair_y +
	  distance_base_pairs / 2.0 * (sin(PI / -2.0 + current_angle));
	g_base_x[connected_to] = base_pair_x + distance_base_pairs / 2.0 *
	  (cos(PI / 2.0 + current_angle));
	g_base_y[connected_to] = base_pair_y + distance_base_pairs / 2.0 *
	  (sin(PI / 2.0 + current_angle));
      } else {
	g_base_x[current_base] = base_pair_x +
	  distance_base_pairs / 2.0 * (cos(PI / 2.0 + current_angle));
	g_base_y[current_base] = base_pair_y +
	  distance_base_pairs / 2.0 * (sin(PI / 2.0 + current_angle));
	g_base_x[connected_to] = base_pair_x + distance_base_pairs / 2.0 *
	  (cos(PI / -2.0 + current_angle));
	g_base_y[connected_to] = base_pair_y + distance_base_pairs / 2.0 *
	  (sin(PI / -2.0 + current_angle));
      }
      if (interior_flag < 0)
	base_angle_increment = 2.0 * PI / (float) (size - 1.);
      if (interior_flag > 0)
	base_angle_increment = 2.0 * PI / (float) (size + 1);
      if ((g_counter_clockwise) && (interior_flag != 0))
	base_angle_increment = base_angle_increment * -1;
      current_angle = current_angle - 1.5 * base_angle_increment;
      current_base = connected_to + 1;
      if (current_base > g_end_base) {
	if (times_around == 0) {
	  times_around = 1;
	  current_base = g_start_base;
	} else {
	  current_base = last_base;
	}
      }
    }
  }
  if (((first_loop) && (g_connected[i] != last_base))
      && (g_connected[last_base] < 1)) { /* draw one last unconnected base */
    g_base_x[last_base] = g_loop_center_x[loop_num] +
      g_loop_radius[loop_num] * cos(current_angle);
    g_base_y[last_base] = g_loop_center_y[loop_num] +
      g_loop_radius[loop_num] * sin(current_angle);
    g_structure_angle[last_base] = current_angle;
  }
}

void regularize_all_angles(void);

void make_natural_structure(int start_base, int end_base, int domain_flag, 
			    int ex_domain_flag) {
  int             total_loops;
  int             i;	/* current base */
  int             connected_to;
  int             helix_length;
  int             loop_bases;
  int             current_helix;
  float           base_average;
  int             size;
  int             domain_row = -9999, domain_col = -9999;
  int             bulge_check;
  /* checks for bulge of size 1 or 2, left or right, 0 for not */
  int             interior_check;
  /* checks for interior loops where 1 side is 1 off from other side
   * -1 for extra base on left, +1 for extra base on right side, 0 for not
   */
  float           base_pair_x;
  float           base_pair_y;
  float           angle;
  float           current_x, current_y;	/* position of last drawn base pair */
  float           distance_base_pairs;
  float           distance_bases;
  if ((!ex_domain_flag) && (!domain_flag))
    distance_base_pairs = DISTANCE_BASE_PAIRS;
  else {
    distance_base_pairs =
      sqrt((g_base_x[start_base] - g_base_x[end_base]) *
	   (g_base_x[start_base] - g_base_x[end_base]) +
	   (g_base_y[start_base] - g_base_y[end_base]) *
	   (g_base_y[start_base] - g_base_y[end_base]));
  }
  if (g_ss_mode)
    distance_bases = distance_base_pairs * g_ss_distance_base_adj /
      g_ss_distance_basepair;
  else
    distance_bases = distance_base_pairs;
  base_average = (distance_bases + distance_base_pairs) / 2.0;
  if (ex_domain_flag) {
    domain_row = start_base;
    domain_col = end_base;
    start_base = g_start_base;
    end_base = g_end_base;
  }
  i = start_base;
  if (domain_flag) {
    total_loops = g_previous_loop[start_base];
  } else {
    total_loops = 1;
    g_loop_center_x[total_loops] = 0.0;
    g_loop_center_y[total_loops] = 0.0;
    g_connected[0] = 0;
  }
  connected_to = g_connected[i];
  if (connected_to > i) {	/* not really any first loop */
    current_helix = total_loops;
    if (!domain_flag) {
      if (!ex_domain_flag)
	g_previous_loop[i] = 1;
      size = 2 + g_end_base - connected_to;
      g_loop_radius[1] = (float) size / (2.0 * PI);
      g_structure_angle[i] = PI;
      g_structure_angle[connected_to] = 0.;
      base_pair_x = g_loop_center_x[total_loops] + 
	g_loop_radius[1] * cos(PI / 2.0);
      base_pair_y = g_loop_center_y[1] +
	g_loop_radius[1] * sin(PI / 2.0);
      g_base_x[i] = base_pair_x +
	distance_base_pairs / 2.0 * (cos(g_structure_angle[i]));
      g_base_y[i] = base_pair_y +
	distance_base_pairs / 2.0 * (sin(g_structure_angle[i]));
      g_base_x[connected_to] = base_pair_x +
	distance_base_pairs / 2.0 *
	(cos(g_structure_angle[connected_to]));
      g_base_y[connected_to] = base_pair_y +
	distance_base_pairs / 2.0 *
	(sin(g_structure_angle[connected_to]));
      current_x = base_pair_x;
      current_y = base_pair_y;
      helix_length = 1;
      /* draw the extra unpaired bases at the end */
      loop_bases = 
	traverse_loop_natural_count_bases(total_loops, g_connected[i], i, 
					  TRUE, &bulge_check, &interior_check);
      /* if statement below added April 10, 2000 */
      /* case with no loop at start did not work */
      if (connected_to < end_base)
	traverse_loop_natural_assign_angles(loop_bases, g_connected[i],
					    i, base_pair_x, base_pair_y, TRUE,
					    g_structure_angle[i]-PI/2.+
					    PI,total_loops,base_average,
					    bulge_check, interior_check);
      if (g_angle_list) {
	printf("Helix: 1  Perpendicular 90 ");
	printf("%d %d\n", start_base, end_base);
      }
    } else {
      current_x = (g_base_x[i] + g_base_x[g_connected[i]]) / 2.;
      current_y = (g_base_y[i] + g_base_y[g_connected[i]]) / 2.;
      helix_length = 1;
    }
  } else {		/* go around the loop counting elements */
    
    if (!ex_domain_flag)
      g_previous_loop[i] = 1;
    loop_bases = 
      traverse_loop_natural_count_bases(1, i, g_end_base, TRUE,
					&bulge_check, &interior_check);
    current_x = 0;
    current_y = 0;
    current_helix = 0;	/* angle argument below is wrong ? */
    traverse_loop_natural_assign_angles(loop_bases, i, g_end_base, current_x,
					current_y, TRUE, PI/2.0 , 1,
					base_average, bulge_check, 
					interior_check);
    /* go around the loop assigning angles */
    helix_length = 0;
  }
  i++;
  while (i <= end_base) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base)) {
	i = g_ex_end_base;
	helix_length = 0;
      }
    }
    if (ex_domain_flag) {
      if ((i > domain_row) && (i < domain_col)) {
	i = domain_col;
	total_loops = g_previous_loop[domain_col];
      }
    }
    connected_to = g_connected[i];
    if (connected_to > i) {
      if (!ex_domain_flag)
	g_previous_loop[i] = total_loops;
      if (helix_length > 0) {
	if (g_connected[i - 1] == (connected_to + 1)) {	/* this is a continued
							 * helix */
	  /* advance and set angle */
	  if (g_angle_list)
	    printf("%d %d\n", i, connected_to);
	  if (g_counter_clockwise)
	    angle = g_structure_angle[i - 1] + PI / 2.0;
	  else
	    angle = g_structure_angle[i - 1] - PI / 2.0;
	  g_structure_angle[i] = g_structure_angle[i - 1];
	  g_structure_angle[connected_to] =
	    g_structure_angle[connected_to + 1];
	  current_x = current_x + distance_bases * cos(angle);
	  current_y = current_y + distance_bases * sin(angle);
	  if (g_counter_clockwise) {
	    g_base_x[i] = current_x+distance_bases/2.*cos(angle-PI/2.0);
	    g_base_y[i] = current_y + distance_base_pairs/
	      2.*sin(angle-PI/ 2.0);
	    g_base_x[connected_to] = current_x + distance_base_pairs / 2.0 *
	      cos(angle + PI / 2.0);
	    g_base_y[connected_to] = current_y + distance_base_pairs / 2.0 *
	      sin(angle + PI / 2.0);
	  } else {
	    g_base_x[i] = current_x + distance_bases / 2.0 *
	      cos(angle + PI / 2.0);
	    g_base_y[i] = current_y + distance_base_pairs / 2.0 *
	      sin(angle + PI / 2.0);
	    g_base_x[connected_to] = current_x + distance_base_pairs / 2.0 *
	      cos(angle - PI / 2.0);
	    g_base_y[connected_to] = current_y + distance_base_pairs / 2.0 *
	      sin(angle - PI / 2.0);
	  }
	  helix_length++;
	} else {	/* start of new helix */
	  if (!ex_domain_flag)
	    g_previous_loop[i] = g_previous_loop[i] + 1;
	  if (g_angle_list)
	    printf("Helix Lngth: %d\n", helix_length);
	  total_loops++;
	  loop_bases =
	    traverse_loop_natural_count_bases(total_loops, i-1, 
					      g_connected[i-1], FALSE, 
					      &bulge_check, &interior_check);
	  traverse_loop_natural_assign_angles(loop_bases, i-1,
					      g_connected[i-1], current_x, 
					      current_y, FALSE,
					      g_structure_angle[i - 1] - 
					      PI / 2., total_loops,
					      base_average, bulge_check, 
					      interior_check);
	  helix_length = 1;
	  current_x = (g_base_x[i] + g_base_x[connected_to]) / 2.0;
	  current_y = (g_base_y[i] + g_base_y[connected_to]) / 2.0;
	  if (g_angle_list) {
	    current_helix++;
	    printf("Helix: %d Perpendicular = %.2f\n", current_helix,
		   (g_structure_angle[i] - PI / 2.0)*180./PI);
	    printf("%d %d\n", i, connected_to);
	  }
	}
      } else {	/* beginning of new helix, it has already been drawn */
	if (g_angle_list) {
	  current_helix++;
	  printf("Helix %d Perpendicular = %.2f\n", current_helix,
		 (g_structure_angle[i] - PI / 2.0) * 180. / PI);
	  printf("%d %d\n", i, connected_to);
	}
	helix_length = 1;
	current_x = (g_base_x[i] + g_base_x[connected_to]) / 2.0;
	current_y = (g_base_y[i] + g_base_y[connected_to]) / 2.0;
      }
    } else {
      if (helix_length > 0) {	/* finish the last helix, traverse the
				 * loop it closed */
	if (g_angle_list)
	  printf("Helix Length: %d\n", helix_length);
	if (!((ex_domain_flag) && (i == domain_col))) {
	  total_loops++;
	  loop_bases =
	    traverse_loop_natural_count_bases(total_loops, i-1, 
					      g_connected[i-1], FALSE,
					      &bulge_check,&interior_check);
	  traverse_loop_natural_assign_angles(loop_bases, i-1, 
					      g_connected[i-1], current_x,
					      current_y, FALSE,
					      g_structure_angle[i-1]-
					      PI/2., total_loops,
					      base_average, bulge_check,
					      interior_check);
	}
	helix_length = 0;
      }
      if (!ex_domain_flag)
	g_previous_loop[i] = total_loops;
    }
    i++;
  }
  /* after while loop, */
  if (helix_length > 0) { /* finish the last helix, traverse the loop
			   * it closed */ 
    if (g_angle_list)
      printf("Helix  Length: %d\n", helix_length);
    total_loops++;
    loop_bases =
      traverse_loop_natural_count_bases(total_loops, i-1, g_connected[i - 1],
					FALSE, &bulge_check, &interior_check);
    traverse_loop_natural_assign_angles(loop_bases, i-1, g_connected[i-1],
					current_x, current_y, FALSE,
					g_structure_angle[i-1] - PI/2.,
					total_loops, base_average, 
					bulge_check, interior_check);
  }
  if ((!domain_flag) && (!ex_domain_flag))
    g_total_loops = total_loops;
  /* M. Zuker enables all angle regularization in natural angles mode.
   * First done for non-interactive  on June 27, 2006, then for all
   * modes on November 27, 2006. The "tiny font" problem was solved
   * on November 27, 2006 by eliminating the call to
   *  scale_structure_to_fit_window within regularize_all_angles
   */
  if (g_reg_angle_flag && !g_flat_alternate)
    regularize_all_angles();
}

int switch_clock(int clock) {
  if (g_flat_alternate) {
    if (clock)
      return FALSE;
    else
      return TRUE;
  } else
    return clock;
}

void rotate_structure(float, float, float, int, int, int, int, int);

void make_flat_structure(int start_base, int end_base, int domain_flag) {
  int             total_loops;
  int             i;	/* current base */
  int             connected_to, j, original_clock;
  float           angle = 0.0;
  float           junk, rot_angle, new_position_x, new_position_y;
  float           current_x, current_y, distance_base_pairs, distance_bases;
  int             whole_flag;
  original_clock = g_counter_clockwise;
  if (g_whole_circle)
    whole_flag = TRUE;
  else
    whole_flag = FALSE;
  if (!domain_flag) {
    distance_base_pairs = DISTANCE_BASE_PAIRS;
  } else {
    distance_base_pairs =
      sqrt((g_base_x[start_base] - g_base_x[end_base]) *
	   (g_base_x[start_base] - g_base_x[end_base]) +
	   (g_base_y[start_base] - g_base_y[end_base]) *
	   (g_base_y[start_base] - g_base_y[end_base]));
    angle = g_structure_angle[start_base];
  }
  if (g_ss_mode)
    distance_bases = distance_base_pairs * g_ss_distance_base_adj /
      g_ss_distance_basepair;
  else
    distance_bases = distance_base_pairs;
  i = start_base;
  total_loops = 1;
  g_loop_center_x[total_loops] = -2.0*distance_base_pairs*cos(angle);
  g_loop_center_y[total_loops] = -2.0*distance_base_pairs*sin(angle);
  g_connected[0] = 0;
  connected_to = g_connected[i];
  g_base_x[i] = 0.0;
  g_base_y[i] = 0.0;
  current_x = 0.0;
  current_y = 0.0;
  angle = 3 * PI / 2;
  g_previous_loop[i] = total_loops;
  if (connected_to == end_base)
    g_connected[i] = 0;
  g_counter_clockwise = switch_clock(g_counter_clockwise);
  if ((connected_to > i) && (connected_to < end_base)) {
    g_counter_clockwise = switch_clock(g_counter_clockwise);
    g_base_x[connected_to] = current_x +
      cos(angle - 3.0 * PI / 2) * distance_base_pairs;
    g_base_y[connected_to] = current_y +
      sin(angle - 3.0 * PI / 2) * distance_base_pairs;
    g_structure_angle[i] = angle - PI / 2;
    printf("angle is %.2f\n", g_structure_angle[i] * 180. / PI);
    g_structure_angle[connected_to] = angle - 3.0 * PI / 2;
    if (g_natural)
      make_natural_structure(i, connected_to, TRUE, FALSE);
    else {
      main_data_ps_structure(1.95 * PI, FALSE, g_connected, TRUE, i, 
			     connected_to, &junk, whole_flag, FALSE, 0.0, 90);
      new_position_x = (g_base_x[i] + g_base_x[connected_to]) / 2;
      new_position_y = (g_base_y[i] + g_base_y[connected_to]) / 2;
      rot_angle = (angle - PI / 2.) - g_structure_angle[i];
      rotate_structure(rot_angle, new_position_x, new_position_y, TRUE, 
		       g_previous_loop[i] + 1, g_previous_loop[connected_to],
		       i, connected_to);
      for (j = i; j <= connected_to; j++)
	g_structure_angle[j] = g_structure_angle[j] +  rot_angle;
    }
    i = connected_to;
    total_loops = g_previous_loop[i];
    current_x = g_base_x[connected_to];
    current_y = g_base_y[connected_to];
  } else {
    g_structure_angle[i] = angle;
  }
  while (i < end_base) {
    i++;
    current_x = current_x + cos(angle - 3.0 * PI / 2.) * distance_bases;
    current_y = current_y + sin(angle - 3.0 * PI / 2.) * distance_bases;
    connected_to = g_connected[i];
    if ((connected_to > i) && (connected_to <= end_base)) {
      g_counter_clockwise = switch_clock(g_counter_clockwise);
      g_base_x[i] = current_x;
      g_base_y[i] = current_y;
      g_previous_loop[i] = total_loops;
      g_base_x[connected_to] = current_x +
	cos(angle - 3.0 * PI / 2) * distance_base_pairs;
      g_base_y[connected_to] = current_y +
	sin(angle - 3.0 * PI / 2) * distance_base_pairs;
      if (g_counter_clockwise) {
	g_structure_angle[i] = angle - PI / 2;
	g_structure_angle[connected_to] = angle - 3.0 * PI / 2;
      } else {
	g_structure_angle[i] = angle - PI / 2;
	g_structure_angle[connected_to] = angle - 3.0 * PI / 2;
      }
      if (g_natural)
	make_natural_structure(i, connected_to, TRUE, FALSE);
      else {
	main_data_ps_structure(1.95 * PI, FALSE,
			       g_connected, TRUE, i, connected_to,
			       &junk, whole_flag, FALSE, 0.0, 90);
	new_position_x = (g_base_x[i] + g_base_x[connected_to]) / 2;
	new_position_y = (g_base_y[i] + g_base_y[connected_to]) / 2;
	rot_angle = (angle - PI / 2.) - g_structure_angle[i];
	
	rotate_structure(rot_angle, new_position_x, new_position_y,
			 TRUE, g_previous_loop[i] + 1,
			 g_previous_loop[connected_to], i, connected_to);
	for (j = i; (j <= connected_to); j++)
	  g_structure_angle[j] = g_structure_angle[j] + rot_angle;
      }
      i = connected_to;
      total_loops = g_previous_loop[i];
      current_x = g_base_x[connected_to];
      current_y = g_base_y[connected_to];
    } else {
      g_base_x[i] = current_x;
      g_base_y[i] = current_y;
      g_structure_angle[i] = angle;
      if (g_connected[i] > g_end_base)
	g_connected[i] = 0;
    }
  }
  if (!domain_flag)
    g_total_loops = total_loops;
  g_loop_center_x[1] = (g_base_x[start_base] + g_base_x[end_base]) / 2.0;
  g_loop_center_y[1] = (g_base_y[start_base] + g_base_y[end_base]) / 2.0;
  g_loop_center_x[1] = g_loop_center_x[1] + 3.0*distance_base_pairs*cos(angle);
  g_loop_center_y[1] = g_loop_center_y[1] + 3.0*distance_base_pairs*sin(angle);
  g_counter_clockwise = original_clock;
}

void make_ps_data(float degrees_used, int is_circular, char *first_line, 
		  char *bases, int *connected) {
  float           junk;
  if (!g_structure) 
    initialize_ps(g_smart_colors, g_annotation, (WD + 72.0)/2.0, 
		  (HT + 72.0)/2.0, 1.0);
  if (g_flat) {
    make_flat_structure(g_start_base, g_end_base, FALSE);
  } else if (g_natural) {
    make_natural_structure(g_start_base, g_end_base, FALSE, FALSE);
  } else if (g_structure) {
    main_data_ps_structure(degrees_used, is_circular, connected, FALSE, 
			   g_start_base, g_end_base, &junk, g_whole_circle, 
			   FALSE, 0.0, 0);
  } else {
    main_data_ps_no_structure(degrees_used, is_circular, first_line, bases, 
			      connected, g_start_base, g_end_base, 
			      g_whole_circle, 0.0);
  }
  if (!g_structure)
    finish_ps_file();
}

void compute_angles_ps(char *out_filename,float degrees_used,int is_circular) {
  if (g_ss_mode) {
    if (g_structure)
      return;	/* nothing to do for this mode */
  }
  if (!g_structure)
    open_output(out_filename);
  make_ps_data(degrees_used, is_circular, g_first_line, g_bases, g_connected);
}

void set_label(int base, int label) {
  float           extra_angle;
  if (g_ss_backwards)
    extra_angle = PI;
  else
    extra_angle = 0.0;
  g_structure_label_x[label] = g_base_x[base] + 1.5 * DISTANCE_BASE_PAIRS * 
    g_scale * cos(g_structure_angle[base] + extra_angle);
  g_structure_label_y[label] = g_base_y[base] + 1.5 * DISTANCE_BASE_PAIRS * 
    g_scale * sin(g_structure_angle[base] + extra_angle);
  g_structure_label_value[label] = base;
}

float find_angle_ss(int b1, int b2) {
  return (float) atan2((g_base_y[b2] - g_base_y[b1]),
		       (g_base_x[b2] - g_base_x[b1])) + PI / 2;
}

float find_angle_ss_connected(int b1) {
  int             b2;
  b2 = g_connected[b1];
  return (float) atan2((g_base_y[b1] - g_base_y[b2]),
		       (g_base_x[b1] - g_base_x[b2]));
}

void set_structure_angles_ss(int start_base, int end_base) {
  int             base;
  float           angle;
  base = start_base;
  if ((base = g_start_base)) {
    if (g_connected[base] > 0) {
      angle = find_angle_ss_connected(base);
      g_structure_angle[base] = angle;
      g_structure_angle[g_connected[base]] = angle - PI;
    } else
      g_structure_angle[base] = find_angle_ss(base, base + 1);
    base++;
  }
  while (base <= end_base) {
    if (g_external_domain_flag) {
      if ((base > g_ex_start_base) && (base < g_ex_end_base))
	base = g_ex_end_base;
    }
    if (g_connected[base] < 1)
      g_structure_angle[base] = find_angle_ss(base - 1, base);
    else {
      if (g_connected[base] > base) {
	angle = find_angle_ss_connected(base);
	g_structure_angle[base] = angle;
	g_structure_angle[g_connected[base]] = angle + PI;
      }
    }
    base++;
  }
}

int set_label_coordinates(void) {
  int number_of_labels, current_place_for_label, i;
  number_of_labels = 0; //2003-06-23 Nick Markham
  for (i = g_start_base; i <= g_end_base; ++i) //2003-06-23 Nick Markham
    if (g_next[i] == 0 || g_prev[i] == 0) // 2003-06-23  Nick Markham 
      set_label(i, ++number_of_labels); // 2003-06-23 Nick Markham
   if (g_label_forced_column > g_end_base)
    g_label_forced_row = 0;
  if (g_label_forced_row > 0) {
    set_label(g_label_forced_row, ++number_of_labels);
    set_label(g_label_forced_column, ++number_of_labels);
  }
  current_place_for_label = g_start_base -
    g_start_base % g_label_frequency + g_label_frequency;
  if (current_place_for_label < (g_start_base + g_label_frequency / 2))
    current_place_for_label++;
    if (g_external_domain_flag)
      if ((current_place_for_label > g_ex_start_base) &&
	  (current_place_for_label < g_ex_end_base)) {
	current_place_for_label = g_ex_end_base;
	number_of_labels++;
	set_label(current_place_for_label, number_of_labels);
	current_place_for_label = current_place_for_label +
	  (g_label_frequency - current_place_for_label % g_label_frequency);
	i = current_place_for_label;
      }
    for (i = current_place_for_label; i <= g_end_base; i+= g_label_frequency) {
      number_of_labels++;
      if (g_external_domain_flag)
	if ((i > g_ex_start_base) && (i < g_ex_end_base)) {
	  i = g_ex_end_base;
	  set_label(i, number_of_labels);
	  number_of_labels++;
	  i = i + (g_label_frequency - i % g_label_frequency);
	  if (i > g_end_base)
	    i = g_end_base;
	}
      set_label(i, number_of_labels);
    }
    g_total_labels = number_of_labels;
    return number_of_labels;
}

void find_structure_width(float *X_width, float *Y_width, float *X_min, 
			  float *X_max, float *Y_min, float *Y_max) {
  float           x_min, x_max, y_min, y_max, x_width, y_width;
  int             i;
  x_min = FLOAT_MAX;
  y_min = FLOAT_MAX;
  x_max = FLOAT_MIN;
  y_max = FLOAT_MIN;	/* find extremes in vertical and horizontal
			 * coordinates */
  /* check bases first */
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag)
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    if (g_base_x[i] > x_max) {
      x_max = g_base_x[i];
    }
    if (g_base_x[i] < x_min) {
      x_min = g_base_x[i];
    }
    if (g_base_y[i] > y_max)
      y_max = g_base_y[i];
    if (g_base_y[i] < y_min)
      y_min = g_base_y[i];
  }
  /* check labels */
  for (i = 1; i <= g_total_labels; i++) {
    if (g_structure_label_x[i] > x_max)
      x_max = g_structure_label_x[i];
    if (g_structure_label_x[i] < x_min)
      x_min = g_structure_label_x[i];
    if (g_structure_label_y[i] > y_max)
      y_max = g_structure_label_y[i];
    if (g_structure_label_y[i] < y_min)
      y_min = g_structure_label_y[i];
  }
  x_width = x_max - x_min;
  y_width = y_max - y_min;	/* 72 points per inch with 13 x 17 */
  *X_width = x_width;
  *Y_width = y_width;
  *X_min = x_min;
  *Y_min = y_min;
  *X_max = x_max;
  *Y_max = y_max;
}

void shift_structure(float hor_shift, float ver_shift, int first_loop, 
		     int last_loop, int start_base, int end_base) {
  int             i;
  for (i = start_base; i <= end_base; i++) {
    g_base_x[i] += hor_shift;
    g_base_y[i] += ver_shift;
  }
  for (i = first_loop; i <= last_loop; i++) {
    g_loop_center_x[i] += hor_shift;
    g_loop_center_y[i] += ver_shift;
  }
}

void rotate_structure(float rotation_angle, float x_center, float y_center, 
		      int labels_too, int first_loop, int last_loop, 
		      int start_base, int end_base) {
  int             i;
  float           x, y;
  float           sin_a, cos_a;
  sin_a = sin(rotation_angle);
  cos_a = cos(rotation_angle);
  /* rotate bases */
  for (i = start_base; i <= end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    x = g_base_x[i];
    y = g_base_y[i];
    g_base_x[i] = cos_a * (x - x_center) - sin_a * (y - y_center) + x_center;
    g_base_y[i] = sin_a * (x - x_center) + cos_a * (y - y_center) + y_center;
  }
  /* rotate labels */
  if (labels_too) {
    for (i = 1; i <= g_total_labels; i++) {
      x = g_structure_label_x[i];
      y = g_structure_label_y[i];
      g_structure_label_x[i]=cos_a*(x-x_center)-sin_a*(y-y_center)+x_center;
      g_structure_label_y[i]=sin_a*(x-x_center)+cos_a*(y-y_center)+y_center;
    }
  }
  for (i = first_loop; i <= last_loop; i++) {
    x = g_loop_center_x[i];
    y = g_loop_center_y[i];
    g_loop_center_x[i] = cos_a*(x-x_center) - sin_a*(y-y_center) + x_center;
    g_loop_center_y[i] = sin_a*(x-x_center) + cos_a*(y-y_center) + y_center;
  }
}

void adjust_structure_coordinates_ss(void) {
  float           x_min, x_max, y_min, y_max;
  float           x_width = WD, y_width = HT;
  float           scale, scale1, scale2;
  float           average;
  float           center_x, center_y;
  float           rotation_angle;
  int             base;
  int             i;
  float           extra_x_translate;
  float           extra_y_translate;
  x_min = FLOAT_MAX;
  y_min = FLOAT_MAX;
  x_max = FLOAT_MIN;
  y_max = FLOAT_MIN;
  set_structure_angles_ss(g_start_base, g_end_base);
  g_scale = g_ss_distance_basepair;
  if (g_label_frequency > 0)
    g_total_labels = set_label_coordinates();
  else
    g_total_labels = 0;
  /* find extremes in vertical and horizontal coordinates */
  if ((g_auto_rotate) || (g_rotation_angle != 0)) {
    find_structure_width(&x_width, &y_width, &x_min, &x_max, &y_min, &y_max);
  }
  rotation_angle = 0;
  if (g_auto_rotate) {
    printf("g_auto_rotate_flag is on\n");
    if (x_width > y_width) {
      printf("Rotating -90 degrees to avoid landscape mode\n");
      rotation_angle = PI / -2.0;
    }
  }
  if (g_rotation_angle != 0)
    rotation_angle = rotation_angle + g_rotation_angle;
  if (rotation_angle != 0) {
    if ((g_auto_rotate) || (g_rotation_angle != 0))
      printf("Structure is rotated %.2f degrees\n",
	     rotation_angle * 180. / PI);
    center_x = x_min + x_width / 2.0;
    center_y = y_min + y_width / 2.0;
    rotate_structure(rotation_angle, center_x, center_y, TRUE,
		     1, g_total_loops, g_start_base, g_end_base);
    if (rotation_angle != 0) {
      for (base = g_start_base; base <= g_end_base; base++)
	g_structure_angle[base] += rotation_angle;
    }
  }
  find_structure_width(&x_width, &y_width, &x_min, &x_max, &y_min, &y_max);
  /* perform rotation */
  scale1 = WD / x_width; /* WD for horizontal width with 1/2 inch
			  * margin */
  scale2 = HT / y_width; /* HT for vertical width with 1/2 inch margin */
  if (scale1 < scale2) {	/* stretching to fill horizontal */
    scale = scale1;
    extra_x_translate = 0;	/* HT is (height in inches - 1)*72 */
    extra_y_translate = (HT - y_width * scale) / 2.;
  }
  /* y_width*scale is height of image */
  else {			/* stretching to fill vertical */
    scale = scale2;
    extra_y_translate = 0;
    /* WD is (width in inches - 1)*72 */
    extra_x_translate = (WD - x_width * scale) / 2.;
  }			/* x_width*scale is width of image */

  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    g_base_x[i] = 36 + extra_x_translate + (g_base_x[i] - x_min) * scale;
    g_base_y[i] = 36 + extra_y_translate + (g_base_y[i] - y_min) * scale;
  }
  for (i = 1; i <= g_total_loops; i++) {
    g_loop_center_x[i] = 36 + extra_x_translate +
      (g_loop_center_x[i] - x_min) * scale;
    g_loop_center_y[i] = 36 + extra_y_translate +
      (g_loop_center_y[i] - y_min) * scale;
    g_loop_radius[i] = g_loop_radius[i] * scale;
  }
  for (i = 1; i <= g_total_labels; i++) {
    g_structure_label_x[i] = 36 + extra_x_translate +
      (g_structure_label_x[i] - x_min) * scale;
    g_structure_label_y[i] = 36 + extra_y_translate +
      (g_structure_label_y[i] - y_min) * scale;
  }
  if (!g_first_ss_run) {
    g_ss_distance_basepair = g_ss_distance_basepair * scale;
    g_ss_distance_base_adj = g_ss_distance_base_adj * scale;
  }
  average = (g_ss_distance_basepair + g_ss_distance_base_adj) / 2.0;
  if (g_first_ss_run)
    g_font_size = average * scale * .5;
  else
    g_font_size = average * scale * .5 * g_scale;
  if (g_font_size > 30.)
    g_font_size = 30.;
  if (g_first_ss_run)
    g_line_width = .05 * scale * g_ss_distance_base_adj * .7;
  else
    g_line_width = .05 * scale * g_ss_distance_base_adj * .7 * g_scale;
  if (g_line_width > 3.0)
    g_line_width = 0.25;
  g_scale = g_scale * scale;
}

void scale_to_fit_window(void) {
  float x_min, x_max, y_min, y_max;
  float x_width = WD, y_width = HT; /* (width in inches - 1)*72 & 
				       (height in inches - 1)*72  */ 
  float scale, scale1, scale2;
  int   i;
  float extra_x_translate;
  float extra_y_translate;
  x_min = FLOAT_MAX;
  y_min = FLOAT_MAX;
  x_max = FLOAT_MIN;
  y_max = FLOAT_MIN;
  find_structure_width(&x_width, &y_width, &x_min, &x_max, &y_min, &y_max);
  scale1 = WD / x_width; /* WD for horizontal width with 1/2 inch margin */
  scale2 = HT / y_width; /* HT for vertical width with 1/2 inch margin */
  if (scale1 < scale2) {	/* stretch to fill horizontal */
    scale = scale1;
    extra_x_translate = 0.0;	/* HT is (height in inches - 1)*72 */
    /* y_width*scale is height of image */
    extra_y_translate = (HT - y_width * scale) / 2.; 
  } else {	    /* stretch to fill vertical */
    scale = scale2;
    extra_y_translate = 0; 
    extra_x_translate = (WD - x_width * scale) / 2.; /* WD = (width in
							inches - 1)*72 */
  }  /* x_width*scale is width of image */
  for (i = g_start_base; i <= g_end_base; i++) {
    g_base_x[i] = 36. + extra_x_translate + (g_base_x[i] - x_min) * scale;
    g_base_y[i] = 36. + extra_y_translate + (g_base_y[i] - y_min) * scale;
  }
  for (i = 1; i <= g_total_loops; i++) {
    g_loop_center_x[i] = 36. + extra_x_translate +
      (g_loop_center_x[i] - x_min) * scale;
    g_loop_center_y[i] = 36. + extra_y_translate +
      (g_loop_center_y[i] - y_min) * scale;
    g_loop_radius[i] = g_loop_radius[i] * scale;
  }
  for (i = 1; i <= g_total_labels; i++) {
    g_structure_label_x[i] = 36. + extra_x_translate +
      (g_structure_label_x[i] - x_min) * scale;
    g_structure_label_y[i] = 36. + extra_y_translate +
      (g_structure_label_y[i] - y_min) * scale;
  }
  g_ss_distance_basepair = g_ss_distance_basepair * scale;
  g_ss_distance_base_adj = g_ss_distance_base_adj * scale;
  g_font_size = g_font_size * scale;
  g_line_width = g_line_width * scale;
  g_scale = g_scale * scale;
}

void adjust_structure_coordinates(void) {
  /* scale image to fit and center it */
  float           x_min, x_max, y_min, y_max;
  float           x_width = WD, y_width = HT;
  float           scale, scale1, scale2;
  float           center_x, center_y;
  float           rotation_angle;
  int             i;
  int             base;
  float           extra_x_translate;
  float           extra_y_translate;
  x_min = FLOAT_MAX;
  y_min = FLOAT_MAX;
  x_max = FLOAT_MIN;
  y_max = FLOAT_MIN;
  if (g_ss_mode) {
    adjust_structure_coordinates_ss();
    return;
  }
  if (g_label_frequency > 0)
    g_total_labels = set_label_coordinates();
  else
    g_total_labels = 0;
  /* find extremes in vertical and horizontal coordinates */
  if ((g_auto_rotate) || (g_rotation_angle != 0)) {
    find_structure_width(&x_width, &y_width, &x_min, &x_max, &y_min, &y_max);
  }
  rotation_angle = 0;
  if (g_auto_rotate) {
    if (x_width > y_width) {
      printf("Rotating -90 degrees to avoid landscape mode\n");
      rotation_angle = PI / -2.0;
    }
  }
  if (g_rotation_angle != 0)
    rotation_angle = rotation_angle + g_rotation_angle;
  if (rotation_angle != 0) {
    if ((g_auto_rotate) || (g_rotation_angle != 0))
      printf("Structure is rotated %.2f degrees\n",
	     rotation_angle * 180. / PI);
    center_x = x_min + x_width / 2.0;
    center_y = y_min + y_width / 2.0;
    rotate_structure(rotation_angle, center_x, center_y, TRUE,
		     1, g_total_loops, g_start_base, g_end_base);
    if (rotation_angle != 0) {
      for (base = g_start_base; base <= g_end_base; base++)
	g_structure_angle[base] += rotation_angle;
    }
  }
  find_structure_width(&x_width, &y_width, &x_min, &x_max, &y_min, &y_max);
    /* perform rotation */
  scale1 = WD / x_width; /* WD for horizontal width with 1/2 inch margin */
  scale2 = HT / y_width; /* HT for vertical width with 1/2 inch margin */ 
  if (scale1 < scale2) { /* stretching to fill horizontal */
    scale = scale1;
    extra_x_translate = 0; 
    extra_y_translate = (HT - y_width * scale) / 2.; /* y_width*scale
						      * is image height */   
  } else { /* stretching to fill vertical */
    scale = scale2;
    extra_y_translate = 0;
    /* WD is (image width in inches - 1)*72 */
    extra_x_translate = (WD - x_width * scale) / 2.;
  } /* x_width*scale is width of image */
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    g_base_x[i] = 36 + extra_x_translate + (g_base_x[i] - x_min) * scale;
    g_base_y[i] = 36 + extra_y_translate + (g_base_y[i] - y_min) * scale;
  }
  for (i = 1; i <= g_total_loops; i++) {
    g_loop_center_x[i] = 36 + extra_x_translate + (g_loop_center_x[i] - 
						   x_min) * scale;
    g_loop_center_y[i] = 36 + extra_y_translate + (g_loop_center_y[i] - 
						   y_min) * scale;
    g_loop_radius[i] = g_loop_radius[i] * scale;
  }
  for (i = 1; i <= g_total_labels; i++) {
    g_structure_label_x[i] = 36 + extra_x_translate +
      (g_structure_label_x[i] - x_min) * scale;
    g_structure_label_y[i] = 36 + extra_y_translate +
      (g_structure_label_y[i] - y_min) * scale;
  }
  g_font_size = DISTANCE_BASE_PAIRS * scale * .5;
  if (g_font_size > 30.)
    g_font_size = 30.;
  g_line_width = .07 * scale;
  if (g_line_width > 1.0)
    g_line_width = 1.0;
  g_scale = g_scale * scale;
}

void make_structure_ps(void) {
  float           scale;
  float           x_center, y_center;
  if (g_command_line_zoom_ps_flag) {
    scale = 1.0 / g_command_line_zoom_ps_s;
    x_center = (float) g_command_line_zoom_ps_x;
    y_center = (float) g_command_line_zoom_ps_y;
  } else {
    scale = 1.0;
    x_center = (WD + 72.0)/2.0;
    y_center = (HT + 72.0)/2.0;
  }
  if (g_command_line_zoom_flag) {
    scale = 1.0 / g_command_line_zoom_s;
    x_center = (float) g_command_line_zoom_x;
    y_center = HT + 72.0 - (float) g_command_line_zoom_y;
    printf("Scale by command_line\n");
  }
  make_structure_ps_real(g_structure_filename, g_smart_colors,
			 g_first_line, g_line_width, g_oldconnect, 
			 g_start_base, g_end_base, g_total_loops,
			 g_loop_center_x, g_loop_center_y, g_font_size,
			 g_outline_mode, g_scale, g_bases, g_base_x,
			 g_base_y, g_lines, g_structure_label_x, 
			 g_structure_label_y, g_structure_label_value,
			 g_loop_labels, g_total_labels,
			 g_external_domain_flag, g_ex_start_base, 
			 g_ex_end_base, g_history_offset, g_annotation,
			 g_ann_to_color, g_annotation_dots, 
			 g_annotation_bases, x_center, y_center, scale, 
			 g_is_circular, g_complot_forced_bases,
			 g_complot_forced_bases_count);
}

void make_structure_img_non_interactive(char *filename, int png_mode,
					int jpg_mode) {
  float           x_shift, y_shift, zoom_factor;
  if (!g_command_line_zoom_flag) {
    x_shift = 0;
    y_shift = 0;
    zoom_factor = 1.0;
  } else { /* convert shifts to (WD + 72) x (HT + 72) */
    x_shift = (WD + 72.0)/2.0 - g_command_line_zoom_x * (WD + 72.0) /
      (float) g_img_width;
    y_shift = -(HT + 72.0)/2.0 + g_command_line_zoom_y * (HT + 72.0) / 
      (float) g_img_height;
    zoom_factor = 1.0 / g_command_line_zoom_s;
  }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  make_structure_img(filename, g_img_width, g_img_height, g_img_interlaced, 
		     g_first_line, g_start_base, g_end_base, g_loop_center_x, 
		     g_loop_center_y, g_total_labels, g_structure_label_x, 
		     g_structure_label_y, g_base_x, g_base_y, g_scale, 
		     g_outline_mode, g_oldconnect, g_font_size, 
		     g_structure_label_value, g_lines,g_bases, 
		     g_loop_labels, g_total_loops, g_external_domain_flag, 
		     g_ex_start_base, g_ex_end_base, g_history_offset, 
		     g_annotation, g_ann_to_color, g_annotation_dots, 
		     g_annotation_bases, g_length, x_shift, y_shift, 
		     zoom_factor, g_is_circular, g_complot_forced_bases,
		     g_complot_forced_bases_count, png_mode, jpg_mode, 
		     g_command_line_x);
#endif
}

void create_ss_output(void) {
  /* May 23, 2007. M. Zuker replaces g_connected with g_oldconnect */
  create_ss_output_real(g_out_ss_filename, g_length, g_oldconnect,
			g_bases, g_base_x, g_base_y, g_ss_code);
}

void store_for_undo(void) {
  int             i;
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    /* for base_x */
    g_undo2_base_x[i] = g_undo1_base_x[i];
    g_undo1_base_x[i] = g_base_x[i];
    /* for base_y */
    g_undo2_base_y[i] = g_undo1_base_y[i];
    g_undo1_base_y[i] = g_base_y[i];
    /* for angle */
    g_undo2_structure_angle[i] = g_undo1_structure_angle[i];
    g_undo1_structure_angle[i] = g_structure_angle[i];
  }
  for (i = 1; i <= g_total_loops; i++) {
    g_undo2_loop_center_x[i] = g_undo1_loop_center_x[i];
    g_undo1_loop_center_x[i] = g_loop_center_x[i];
    g_undo2_loop_center_y[i] = g_undo1_loop_center_y[i];
    g_undo1_loop_center_y[i] = g_loop_center_y[i];
  }
  /* store global scale */
  g_undo2_scale = g_undo1_scale;
  g_undo1_scale = g_scale;
  g_undo2_ss_distance_basepair = g_undo1_ss_distance_basepair;
  g_undo1_ss_distance_basepair = g_ss_distance_basepair;
  g_undo2_ss_distance_base_adj = g_undo1_ss_distance_base_adj;
  g_undo1_ss_distance_base_adj = g_ss_distance_base_adj;
  g_undo2_font_size = g_undo1_font_size;
  g_undo1_font_size = g_font_size;
  g_undo2_line_width = g_undo1_line_width;
  g_undo1_line_width = g_line_width;
}

/* M. Zuker Dec 20, 2006
 * Clean up syntax and allow exterior domain partial
 */
void run_sir_graph_on_domain(int ext_flag, int partial_circle_flag) {
  char            command[200];
  char            number[20];
  int             error_occurred;
  if (!g_domain_defined) {
    printf("Error! No domain defined\n");
    return;
  }
  strcpy(command, "sir_graph ");
  if (g_counter_clockwise)
    strcat(command, "-cc ");
  if (!g_lines)
    strcat(command, "-D ");
  if (g_outline_mode)
    strcat(command, "-outline ");
  if (g_natural)
    strcat(command, "-n ");
  if (g_fix_loop_flag)
    strcat(command, "-fix ");
  if (g_annotation) {
    if (g_annotation == P_NUM)
      strcat(command, "-pnum ");
    if (g_annotation == PROB)
      strcat(command, "-prob ");
    if (g_annotation == SS_COUNT)
      strcat(command, "-ss-count ");
    if (g_annotation_bases && (!g_annotation_dots))
      strcat(command, "-ab ");
    if ((!g_annotation_bases) && g_annotation_dots)
      strcat(command, "-ad ");
    strcat(command, "-af ");
    strcat(command, g_annotation_filename);
  }
  if (ext_flag)
    strcat(command, " -e");
  else
    strcat(command, " -i");
  if (!g_natural && !partial_circle_flag)
    strcat(command, "w ");
  else
    strcat(command, " ");
  num_string_int(number, g_domain_row + g_history_offset);
  strcat(command, number);
  strcat(command, " ");
  num_string_int(number, g_domain_column + g_history_offset);
  strcat(command, number);
  strcat(command, " ");
  strcat(command, g_ct_filename);
  strcat(command, " &"); /* Fork sir_graph */
  printf("Command: %s\n", command);
  error_occurred = system(command);
  if (error_occurred != 0)
    printf("Error trying to execute %s\n", PROGRAM_NAME);
  return;
}

int find_flat_external_loop(int domain_base) {
  /* For the flat external loop, determine if domain_base is part of
     the external loop 
     * return TRUE or FALSE
     * find the previous base pair on the loop and set
     * g_flat_fixed_base to the j coordinate 
     */
	int   i;
	int   start=0; /* first base on loop containing domain_base */ 
	int   found_start;
	int   previous_helix_found;
	i = domain_base;
	if (i == g_start_base)
	  return FALSE;
	previous_helix_found = FALSE;
	found_start = FALSE;
	/* find first base on this loop */
	while (!found_start) {
	  i--;
	  if (g_connected[i] > i) {
	    if ((g_connected[i] > domain_base) && (i > g_start_base))
	      return FALSE;
	    found_start = TRUE;
	    start = i;
	  } else {
	    if (g_connected[i] > 0) {
	      if (!previous_helix_found) {
		g_flat_fixed_base = i; /* leave this base and others
					* to its left fixed for drag */
		previous_helix_found = TRUE;
	      }
	      i = g_connected[i];
	    }
	  }
	  if (i <= g_start_base) {
	    start = g_start_base;
	    found_start = TRUE;
	    {
	      if (!previous_helix_found)
		g_flat_fixed_base = g_start_base; /* leave this base
						     of loop fixed */ 
	    }
	  }
	}
	if (g_previous_loop[start] > 1)
	  return FALSE;
	else
	  return TRUE;
}

void traverse_loop_structure_regularize_loop(int loop_num, int loop_start, 
					     int loop_end, int *connected, 
					     float angle, 
					     float current_position_x, 
					     float current_position_y,
					     int extra_base_exists, 
					     int largest_base,
					     float distance_base_pairs, 
					     int draw_first_base_pair,
					     float angle_loop_rotated) {
  /* loop_num is the loop we are drawing
   * loop runs on the base pair loop_start to loop_end
   * angle is angle of helix this is connected to
   * from perpendicular
   * current_position_x,current_position_y are last drawn base pair
   * largest base is g_length, or the largest base with -iw or -i
   * extra_base_exists draws trailing single stranded base on first
   * loop of structure
   * This function is for circle graph angles
   * For a loop beginning with a base pair, draw it if
   * draw_first_base_pair
   * is true.  Normally a loop starts with this base pair already drawn
   * The first loop of a structure is the exception
   * angle_loop_rotated is the amount this loop has been rotated
   * from its position prior to regularize angles on this loop
   * as angles change on a domain, loops are rotated many times
   */
  int             bases_on_list;
  int             size;
  int             current_place;
  int             connected_to;
  float           last_helix_angle = 0.0 , helix_angle = 0.0;
  float           base_angle;
  float           placement_angle;
  int             last_base_was_connected;
  int             helices_on_loop;
  int             times_around;
  float           last_regularized_helix_angle = 0.0;
  float           base_pair_angle;	/* perpendicular angle of
					 * base_pair on loop */
  times_around = 0;
  /*
    if(angle>2.0*PI)
    angle=angle-2.0*PI; */
  
  current_place = loop_start;
  size = 1;
  if (g_reg_angle_flag)
    /*    last_regularized_helix_angle = -5000.; May 19, 2007. M. Zuker */
    last_regularized_helix_angle = 0.0;
  if (!g_counter_clockwise) {
    last_helix_angle = angle + PI;
    helix_angle = angle + PI;
  } else {
    last_helix_angle = angle + 2.0 * PI;
    helix_angle = angle + 2.0 * PI;
  }
  helices_on_loop = 0;
  bases_on_list = 0;
  
  if (g_angle_list) {
    printf("Traversing loop %d from base %d to base %d\n", loop_num,
	   loop_start, loop_end);
  }
  if (current_place == 0) {
    if ((g_angle_list) && (g_start_base == 1))
      printf("5 prime_end\n");
    current_place = 1;
    if (!g_is_circular)
      size = 2;
  }
  if (loop_start > loop_end)
    size = 2;
  if (connected[current_place] > current_place)
    last_base_was_connected = TRUE;
  else
    last_base_was_connected = FALSE;
  /* last_base_was_connected keeps track if base prior to
   * current_place was paired.  
   * Go around loop setting angles 
   */
  while (current_place != loop_end) {
    if (current_place > largest_base) {
      if (times_around == 0) {
	times_around++;
	if (g_is_circular)
	  current_place = g_start_base;
	else
	  current_place = g_start_base - 1;
      } else {
	if (error("pk", loop_start, loop_end)==1)
	  try_exit(13);
      }
    }
    connected_to = connected[current_place];
    if ((connected_to > current_place) && (connected_to <= largest_base)) {
      if (g_angle_list)
	printf("%d * %d\n", current_place, connected_to);
      if (!g_counter_clockwise)
	helix_angle = g_structure_angle[current_place] - PI / 2;
      else
	helix_angle = g_structure_angle[current_place] + PI / 2;
      if (loop_num > 1)
	helix_angle = helix_angle - angle_loop_rotated;
      helix_angle = regularize_angle(helix_angle);
      if ((helix_angle == last_regularized_helix_angle) && (helices_on_loop > 1)) {
	if (g_counter_clockwise)
	  helix_angle = helix_angle + g_reg_angle * PI / 180.;
	else
	  helix_angle = helix_angle - g_reg_angle * PI / 180.;
      }
      last_regularized_helix_angle = helix_angle;
      if ((current_place > loop_start) || (draw_first_base_pair))
	g_base_x[current_place] = helix_angle;
      helices_on_loop++;
      if (extra_base_exists) {
	current_place = connected_to;
      }
      if (size > 1) {	/* set angles of bases on list */
	current_place = connected_to;
	last_base_was_connected = TRUE;
	store_base_angles(last_helix_angle, helix_angle, bases_on_list,
			  g_base_list);
	last_helix_angle = helix_angle;
	bases_on_list = 0;
      } else {
	current_place++;
	if (current_place > largest_base) {
	  if (times_around == 0) {
	    times_around++;
	    if (g_is_circular)
	      current_place = g_start_base;
	    else
	      current_place = g_start_base - 1;
	  } else {
	    if (error("pk", loop_start, loop_end)==1)
	      try_exit(14);
	  }
	}
	last_base_was_connected = FALSE;
      }
    } else {
      if (g_angle_list)
	printf("%d\n", current_place);
      if (!last_base_was_connected) {
	bases_on_list++;
	g_base_list[bases_on_list] = current_place;
      }
      current_place++;
      if (current_place > largest_base) {
	if (times_around == 0) {
	  times_around++;
	  if (g_is_circular)
	    current_place = g_start_base;
	  else
	    current_place = g_start_base - 1;
	} else {
	  if (error("pk", loop_start, loop_end)==1)
	    try_exit(15);
	}
      }
      last_base_was_connected = FALSE;
    }
    size++;
  }
  if (extra_base_exists) {
    if (g_angle_list)
      printf("%d\n", current_place);
    if (!last_base_was_connected) {
      bases_on_list++;
      g_base_list[bases_on_list] = current_place;
    }
    size++;
  }
  if (bases_on_list > 0) {/* set radius of bases on list */
    if (g_counter_clockwise) {
      if ((helices_on_loop > 1) || (draw_first_base_pair)) {
	store_base_angles(helix_angle, angle + 2.0 * PI, bases_on_list,
			  g_base_list);
      } else {
	store_base_angles(angle, angle + 2.0 * PI, bases_on_list, 
			  g_base_list);
      }
    } else {
      if ((helices_on_loop > 1) || (draw_first_base_pair))
	store_base_angles(helix_angle, angle - PI, bases_on_list,
			  g_base_list);
      else
	store_base_angles(angle + PI, angle - PI, bases_on_list,
			  g_base_list);
    }
  }
  if (g_angle_list) {
    printf("%d\n", loop_end);
    printf("Loop Size: %d\n", size);
  }
  /* set center and radius of the loop */
  g_loop_radius[loop_num] = (float) size / (2.0 * PI) * g_scale;
  if (loop_num > 1) {
    if (g_counter_clockwise)
      angle = angle + PI;
    g_loop_center_x[loop_num] = current_position_x +
      g_loop_radius[loop_num] * cos(angle);
    g_loop_center_y[loop_num] = current_position_y +
      g_loop_radius[loop_num] * sin(angle);
  }
  /* else  keep center in same place for first loop
   * {g_loop_center_x[1]=0.0; g_loop_center_y[1]=0.0; }
   * set positions of  bases around loop based on angles and
   * number of bases around loop */
  current_place = loop_start;
  size = 1;
  times_around = 0;
  while (current_place != loop_end) {
    connected_to = g_connected[current_place];
    if ((connected_to > current_place) && (connected_to <= g_end_base)) {
      /* draw the base pair here */
      if ((current_place > loop_start) || (draw_first_base_pair)) {
	base_pair_angle = g_base_x[current_place];
	placement_angle = base_pair_angle;
	draw_basepair(current_place, connected_to, base_pair_angle, 
		      g_loop_center_x[loop_num] +
		      (g_loop_radius[loop_num] + .1) * cos(placement_angle),
		      g_loop_center_y[loop_num] +
		      (g_loop_radius[loop_num] + .1) * sin(placement_angle),
		      distance_base_pairs);
      }
      if ((size > 1) || (loop_num == 1))
	current_place = connected_to;
      else {
	current_place++;
	if (current_place > largest_base) {
	  if (times_around == 0) {
	    times_around++;
	    if (g_is_circular)
	      current_place = g_start_base;
	    else
	      current_place = g_start_base - 1;
	  } else {
	    if (error("pk", loop_start, loop_end)==1)
	      try_exit(16);
	  }
	}
      }
    } else {
      if (connected_to < 1) {
	base_angle = g_base_x[current_place];
	g_base_x[current_place] = g_loop_center_x[loop_num] +
	  (g_loop_radius[loop_num] + .1) * cos(base_angle);
	g_base_y[current_place] = g_loop_center_y[loop_num] +
	  (g_loop_radius[loop_num] + .1) * sin(base_angle);
      }
      current_place++;
      if (current_place > largest_base) {
	if (times_around == 0) {
	  times_around++;
	  if (g_is_circular)
	    current_place = g_start_base;
	  else
	    current_place = g_start_base - 1;
	} else {
	  if(error("pk", loop_start, loop_end)==1)
	    try_exit(17);
	}
      }
    }
    size++;
  }
  if (extra_base_exists) {
    connected_to = connected[current_place];
    if (connected_to < 1) {
      base_angle = g_base_x[current_place];
      g_base_x[current_place] = g_loop_center_x[loop_num] +
	(g_loop_radius[loop_num] + .1) * cos(base_angle);
      g_base_y[current_place] = g_loop_center_y[loop_num] +
	(g_loop_radius[loop_num] + .1) * sin(base_angle);
    }
  }
}

float find_dif_angles(float old_angle, float new_angle) {
  /* Place both angles within 0 to 2.0*PI */
  float           difference;
  /*
    while(old_angle>=2.0*PI)
    old_angle-=PI;
    while(old_angle<0.)
    old_angle+=PI;
    while(new_angle>=2.0*PI)
    new_angle-=PI;
    while(new_angle<0.)
    new_angle+=PI;*/
  difference = new_angle - old_angle;
  return difference;
}

float angle_dif_loop_fix(float old, float current) {
  /* for clockwise, old should be bigger than current
   * return the difference of the two angles */
  float           dif;
  if (g_counter_clockwise) {
    while ((current - old) > 2.0 * PI)
      current = current - 2.0 * PI;
    while ((current - old) < 0)
      current = current + 2.0 * PI;
    dif = current - old;
  } else {
    while ((old - current) > 2.0 * PI)
      old = old - 2.0 * PI;
    while ((old - current) < 0)
      old = old + 2.0 * PI;
    dif = old - current;
  }
  return dif;
}

float distance(float x1, float x2, float y1, float y2) {
  float           val;
  val = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
  val = sqrt(val);
  return val;
}

float get_angle_dif_1(float angle1, float angle2) {
  float           dif;
  if (g_counter_clockwise) {
    while (angle2 < angle1)
      angle2 += 2.0 * PI;
    while ((angle2 - angle1) >= 2.0 * PI)
      angle2 = angle2 - 2.0 * PI;
    dif = angle2 - angle1;
  } else {
    while (angle1 < angle2)
      angle1 += 2.0 * PI;
    while ((angle1 - angle2) >= 2.0 * PI)
      angle1 = angle1 - 2.0 * PI;
    dif = angle1 - angle2;
  }
  return dif;
}

void fix_loop_border_helices(int domain_start, int domain_end, int *lhp1_i,
			     int *fhp2a_i, int *lhp2a_i, int *fhp2b_i, 
			     int *lhp2b_i, int *fhp3_i, int *lhp3_i) {
  /* find last helix in part 1 lhp1_i
   * find first and last helix in part 2a fhp2a_i lhp2a_i
   * find first and last helix in part2b fhp2b_i lhp2b_i
   * find first and last helix in part3  fhp3_i lhp3_i
   * the parts are as follows, from incoming helix for clockwise
   * clockwise 90 degrees is part 1
   * clockwise 90 to 180 degrees is part2a
   * clockwise 180 to 270 degrees is part2b
   * clockwise 270 to 360 degrees is part3
   * Originally part2a, part2b were intended to have special treatment
   * so far, they might as well have been one part 
   */
  int             base;
  int             connected_to;
  int             last_included_helix_base, current_included_helix_base;
  float           angle_1;
  float           angle_2;
  float           angle_dif;
  last_included_helix_base = -1;
  current_included_helix_base = -1;
  /* start of part 1 */
  angle_1 = g_structure_angle[domain_start];
  base = domain_start;
  if (domain_start == g_start_base) {
    if (g_connected[domain_start] < domain_start) {
      if (g_counter_clockwise)
	angle_1 = angle_1 + PI / 2.;
      else
	angle_1 = angle_1 - PI / 2.;
    } else {
      base = g_connected[domain_start];
      angle_1 = g_structure_angle[g_connected[domain_start]];
      base = domain_start - 1;
    }
  }
  last_included_helix_base = -1;
  current_included_helix_base = -1;
  angle_dif = 0.0;
  while ((angle_dif < (PI) / 2) && (base < domain_end)) {
    base++;
    connected_to = g_connected[base];
    while ((connected_to < base) && (base < domain_end)) {
      base++;
      connected_to = g_connected[base];
    }
    if (base < domain_end) {
      angle_2 = g_structure_angle[connected_to];
      angle_dif = get_angle_dif_1(angle_1, angle_2);
      last_included_helix_base = current_included_helix_base;
      current_included_helix_base = base;
      base = connected_to;
    } else {	/* base is domain_end, set last to current */
      /* set current to -1 so part 3 has no helix */
      last_included_helix_base = current_included_helix_base;
      current_included_helix_base = -1;
    }
  }
  *lhp1_i = last_included_helix_base;
  *fhp2a_i = -1;
  *lhp2a_i = -1;
  if (current_included_helix_base > 0) {
  } else {
    *fhp2b_i = -1;
    *lhp2b_i = -1;
    *fhp3_i = -1;
    *lhp3_i = -1;
    return;
  }
  /* start of part 2 a */
  base = current_included_helix_base;
  if ((angle_dif < PI) && (base < domain_end)) {
    *fhp2a_i = current_included_helix_base;
    last_included_helix_base = current_included_helix_base;
    base = g_connected[base];
    current_included_helix_base = -1;
    while ((angle_dif <= PI) && (base < domain_end)) {
      base++;
      connected_to = g_connected[base];
      while ((connected_to < base) && (base < domain_end)) {
	base++;
	connected_to = g_connected[base];
      }
      if (base < domain_end) {
	angle_2 = g_structure_angle[connected_to];
	angle_dif = get_angle_dif_1(angle_1, angle_2);
	last_included_helix_base = current_included_helix_base;
	current_included_helix_base = base;
	base = connected_to;
      } else {	/* base is domain_end, set last to current */
	/* set current to -1 so part 3 has no helix */
	last_included_helix_base = current_included_helix_base;
	current_included_helix_base = -1;
      }
    }
    *lhp2a_i = last_included_helix_base;
  }
  if (*lhp2a_i < 0) {
    *lhp2a_i = *fhp2a_i;	/* none found beyond first helix */
  }
  if (current_included_helix_base > 0) {
    /* printf("*****  first base pair beyond part 2a is (%d,%d)\n",
     * current_included_helix_base, g_connected[current_included_helix_base]);
     */
  } else {
    *fhp2b_i = -1;
    *lhp2b_i = -1;
    *fhp3_i = -1;
    *lhp3_i = -1;
    return;
  }
  *fhp2b_i = -1;
  *lhp2b_i = -1;
  /* start of part 2b */
  base = current_included_helix_base;
  if ((angle_dif <= (3.0 * PI / 2.)) && (base < domain_end)) {
    *fhp2b_i = current_included_helix_base;
    last_included_helix_base = current_included_helix_base;
    base = g_connected[base];
    current_included_helix_base = -1;
    while ((angle_dif <= (3.0 * PI / 2.0)) && (base < domain_end)) {
      base++;
      connected_to = g_connected[base];
      while ((connected_to < base) && (base < domain_end)) {
	base++;
	connected_to = g_connected[base];
      }
      if (base < domain_end) {
	angle_2 = g_structure_angle[connected_to];
	angle_dif = get_angle_dif_1(angle_1, angle_2);
	last_included_helix_base = current_included_helix_base;
	current_included_helix_base = base;
	base = connected_to;
      } else {	/* base is domain_end, set last to current */
	/* set current to -1 so part 3 has no helix */
	last_included_helix_base = current_included_helix_base;
	current_included_helix_base = -1;
      }
    }
    *lhp2b_i = last_included_helix_base;
  }
  if (*lhp2b_i < 0) {
    *lhp2b_i = *fhp2b_i;	/* none found beyond first helix */
  }
  if (current_included_helix_base > 0) {
    /* printf("***** first base pair beyond part 2b is (%d,%d)\n",
     * current_included_helix_base, g_connected[current_included_helix_base]);
     */
  } else {
    *fhp3_i = -1;
    *lhp3_i = -1;
    return;
  }
  /* start of part 3 */
  base = current_included_helix_base;
  if (base < domain_end) {
    *fhp3_i = current_included_helix_base;
    current_included_helix_base = -1;
    base = g_connected[base];
    while (base < domain_end) {
      base++;
      connected_to = g_connected[base];
      while ((connected_to < base) && (base < domain_end)) {
	base++;
	connected_to = g_connected[base];
      }
      if (base < domain_end) {
	angle_2 = g_structure_angle[connected_to];
	angle_dif = get_angle_dif_1(angle_1, angle_2);
	current_included_helix_base = base;
	base = connected_to;
      }
    }
    *lhp3_i = current_included_helix_base;
  }
  if (*lhp3_i < 0) {
    *lhp3_i = *fhp3_i;	/* none found beyond first helix */
  }
}

void loop_fix_spread_singles(int last_paired_base, int next_paired_base,
			     int number_of_single_bases, float distance_bases, 
			     float distance_base_pairs, float angle_next,
			     float angle_last, float radius, 
			     float loop_center_x, float loop_center_y, 
			     float dif_from_previous_helix_x,
			     float dif_from_previous_helix_y) {
  /* spreads single stranded bases when 2 consecutive helices are OK,
   * but there are too many bases between them */
  int             base;
  float           distance_push;
  int             i;
  int             middle_base;
  float           angle_increment;
  float           pos_x, pos_y;
  float           distance_between;
  float           angular_width;
  /* first spread bases equally between paired bases
   * find angular width of a base pair */
  angular_width = distance_base_pairs / radius;
  angle_increment = get_angle_dif_1(angle_last, angle_next);
  angle_increment = (angle_increment - 2.0 * angular_width) /
    (number_of_single_bases + 1.);
  i = 1;
  /* position the bases between the two helices */
  for (base = (last_paired_base + 1); base < next_paired_base; base++) {
    if (g_counter_clockwise) {
      pos_x = loop_center_x + radius * cos(angle_last + angular_width
					   + i * angle_increment);
      pos_y = loop_center_y + radius * sin(angle_last + angular_width +
					   i * angle_increment);
      g_structure_angle[base] = angle_last + angular_width +
	i * angle_increment;
    } else {
      pos_x = loop_center_x + radius * cos(angle_last - angular_width
					   - i * angle_increment);
      pos_y = loop_center_y + radius * sin(angle_last - angular_width -
					   i * angle_increment);
      g_structure_angle[base] = angle_last - angular_width - i*angle_increment;
    }
    g_base_x[base] = pos_x;
    g_base_y[base] = pos_y;
    i++;
  }
  /* next push bases outward */
  distance_between = angle_increment * radius;
  
  distance_push = sqrt(distance_bases * distance_bases * .85 * .85 -
		       distance_between * distance_between);
  i = 1;
  middle_base = last_paired_base + 1 + number_of_single_bases / 2 - 1 +
    number_of_single_bases % 2;
  for (base = (last_paired_base + 1); base <= middle_base; base++) {
    g_base_x[base] = g_base_x[base] + 
      i * distance_push * cos(g_structure_angle[base]); 
    g_base_y[base] = g_base_y[base] +
      i * distance_push * sin(g_structure_angle[base]);
    i++;
  }
  i = 1;
  for (base = (next_paired_base - 1); base > middle_base; base--) {
    g_base_x[base] = g_base_x[base] + 
      i * distance_push * cos(g_structure_angle[base]);
    g_base_y[base] = g_base_y[base] +
      i * distance_push * sin(g_structure_angle[base]);
    i++;
  }
  for (base = (last_paired_base + 1); base < next_paired_base; base++) {
    g_base_x[base] += dif_from_previous_helix_x;
    g_base_y[base] += dif_from_previous_helix_y;
  }
}

void fix_loop_spread_helices(int loop_num, int domain_start, int domain_end,
			     float distance_bases, float distance_base_pairs,
			     float distance_0, float distance_1, 
			     float distance_2, float distance_2p) {
  int             base;
  int             single_bases;
  float           distance_to_center_last, distance_to_center_current;
  float           loop_center_x, loop_center_y;
  float           x_old_center_of_last, y_old_center_of_last;
  float           angle_old, angle_current;
  /* original location of last helix */
  float           x_old_center_of_current, y_old_center_of_current;
  float           minimum_distance;
  float           angle_dif;
  float           distance_between, distance_change, distance_change_angle;
  float           distance_change_x, distance_change_y;
  float           angle_of_singles;
  float           last_helix_dif_x;
  float           last_helix_dif_y;
  float           push_amount_x, push_amount_y;
  /* original location of current helix */
  float           distance_between_unmoved_helices;
  int             number_of_helices;
  int             connected_to_last;
  int             i;
  int             connected_to_current;
  int             single_base;
  int             single_bases_is_odd;
  int             number_of_single_bases;
  int             stop_left, start_right, middle_base=-5000;
  base = domain_start;
  number_of_helices = 1;
  loop_center_x = g_loop_center_x[loop_num];
  loop_center_y = g_loop_center_y[loop_num];
  x_old_center_of_last =
    (g_base_x[domain_start] + g_base_x[domain_end]) / 2.;
  y_old_center_of_last =
    (g_base_y[domain_start] + g_base_y[domain_end]) / 2.;
  angle_old = g_structure_angle[domain_start];
  if (domain_start != g_start_base) {
    if (g_counter_clockwise)
      angle_old = angle_old - PI / 2.;
    else
      angle_old = angle_old + PI / 2;
  } else {
    if (g_connected[domain_start] > 0) {
      if (g_counter_clockwise) {
	angle_old += PI / 2.;
      } else {
	angle_old -= PI / 2.;
      }
    }
  }
  connected_to_last = domain_start;
  distance_to_center_last = distance(x_old_center_of_last, loop_center_x,
				     y_old_center_of_last, loop_center_y);
  if ((base == g_start_base) && (g_connected[base] > base)) {
    base = g_connected[base];
    connected_to_last = g_connected[g_start_base];
  }
  base++;
  last_helix_dif_x = 0.;
  last_helix_dif_y = 0.;
  while (base < domain_end) {
    connected_to_current = g_connected[base];
    if (connected_to_current > base) {	/* this is a new helix */
      number_of_helices++;
      x_old_center_of_current =
	(g_undo1_base_x[base] +
	 g_undo1_base_x[connected_to_current]) / 2.;
      y_old_center_of_current =
	(g_undo1_base_y[base] +
	 g_undo1_base_y[connected_to_current]) / 2.;
      angle_current = g_structure_angle[base];
      if (g_counter_clockwise)
	angle_current = angle_current + PI / 2.;
      else
	angle_current = angle_current - PI / 2;
      angle_dif = angle_dif_loop_fix(angle_old, angle_current);
      distance_to_center_current = 
	distance(x_old_center_of_current, loop_center_x, 
		 y_old_center_of_current, loop_center_y);
      distance_between = angle_dif *
	(distance_to_center_current + distance_to_center_last) / 2.0;
      if (base == (connected_to_last + 1)) {
	single_bases = 0;
	minimum_distance = distance_0;
      } else {
	if (base == (connected_to_last + 2)) {
	  single_bases = 1;
	  minimum_distance = distance_1;
	} else {
	  if (base == (connected_to_last + 3)) {
	    single_bases = 2;
	    minimum_distance = distance_2;
	  } else {
	    single_bases = 3;
	    minimum_distance = distance_2p;
	  }
	}
      }
      if (minimum_distance > distance_between) {	/* fix the helix */
	distance_change = minimum_distance - distance_between;
	distance_change = distance_change * .90;
	if ((domain_start == g_start_base) &&
	    (connected_to_last == domain_start)) {
	  if (g_connected[domain_start] > domain_start)
	    distance_change_angle = g_structure_angle[domain_start];
	  else {
	    if (g_counter_clockwise)
	      distance_change_angle = g_structure_angle[domain_start] +
		PI / 2.;
	    else
	      distance_change_angle = g_structure_angle[domain_start] -
		PI / 2.;
	  }
	} else {
	  distance_change_angle = g_structure_angle[connected_to_last];
	}
	if (g_counter_clockwise)
	  distance_change_angle += angle_dif / 2;
	else
	  distance_change_angle -= angle_dif / 2;
	distance_change_x = cos(distance_change_angle) * distance_change;
	distance_change_y = sin(distance_change_angle) * distance_change;
	g_base_x[base] += distance_change_x;
	g_base_y[base] += distance_change_y;
	g_base_x[connected_to_current] += distance_change_x;
	g_base_y[connected_to_current] += distance_change_y;
	if (single_bases == 1) {
	  g_base_x[base - 1] += distance_change_x / 2.;
	  g_base_y[base - 1] += distance_change_y / 2.;
	}
	if (single_bases == 2) {
	  g_base_x[base - 2] += distance_change_x / 3.;
	  g_base_y[base - 2] += distance_change_y / 3.;
	  g_base_x[base - 1] += 2. * distance_change_x / 3.;
	  g_base_y[base - 1] += 2. * distance_change_y / 3.;
	}
	if (single_bases == 3) {
	  number_of_single_bases = base - 1 - connected_to_last;
	  g_base_x[connected_to_last + 1] += distance_change_x / 3.;
	  g_base_y[connected_to_last + 1] += distance_change_y / 3.;
	  g_base_x[base - 1] += 2. * distance_change_x / 3.;
	  g_base_y[base - 1] += 2. * distance_change_y / 3.;
	  if (number_of_single_bases == 3) {
	    g_base_x[base - 2] += distance_change_x / 2.;
	    g_base_y[base - 2] += distance_change_y / 2.;
	  } else {
	    if ((number_of_single_bases % 2) == 1) {
	      single_bases_is_odd = TRUE;
	      middle_base = (connected_to_last + base) / 2;
	      stop_left = middle_base - 1;
	      start_right = middle_base + 1;
	      g_base_x[middle_base] += distance_change_x / 2.;
	      g_base_y[middle_base] += distance_change_y / 2.;
	    } else {
	      single_bases_is_odd = FALSE;
	      stop_left = (connected_to_last + base) / 2;
	      start_right = stop_left + 1;
	    }
	    for (single_base = (connected_to_last + 2);
		 single_base <= stop_left; single_base++) {
	      g_base_x[single_base] += distance_change_x / 3.;
	      g_base_y[single_base] += distance_change_y / 3.;
	    }
	    for (single_base = start_right; single_base < (base - 1);
		 single_base++) {
	      g_base_x[single_base] += distance_change_x * 2. / 3.;
	      g_base_y[single_base] += distance_change_y * 2. / 3.;
	    }
	    /* push single bases out from loop
	     */
	    if (g_counter_clockwise)
	      angle_of_singles = angle_current - angle_dif / 2.;
	    else
	      angle_of_singles = angle_current + angle_dif / 2.;
	    /*
	     * the .60 below makes the bases stick out less
	     * the value could be adjusted to 1.0
	     */
	    push_amount_x = 0.60 * distance_bases * cos(angle_of_singles);
	    push_amount_y = 0.60 * distance_bases * sin(angle_of_singles);
	    i = 0;
	    for (single_base = (connected_to_last + 2);
		 single_base <= stop_left; single_base++) {
	      i++;
	      g_base_x[single_base] += (i * push_amount_x);
	      g_base_y[single_base] += (i * push_amount_y);
	    }
	    i = 0;
	    for (single_base = (base - 2);
		 single_base >= start_right;
		 single_base--) {
	      i++;
	      g_base_x[single_base] += (i * push_amount_x);
	      g_base_y[single_base] += (i * push_amount_y);
	    }
	    if (single_bases_is_odd) {
	      i++;
	      g_base_x[middle_base] += (i * push_amount_x);
	      g_base_y[middle_base] += (i * push_amount_y);
	    }
	  }
	  
	}
      /* finished moving the helix if necessary */
      } else {	/* check if the helix is ok, but the bases
		 * are too close together on the loop */
	number_of_single_bases = base - 1 - connected_to_last;
	distance_between_unmoved_helices = distance_between -
	  distance_base_pairs - distance_bases;
	if ((distance_between_unmoved_helices - number_of_single_bases * 
	     distance_bases * .5) < 0) {
	  loop_fix_spread_singles(connected_to_last, base, 
				  number_of_single_bases, distance_bases, 
				  distance_base_pairs, angle_current,
				  angle_old, distance_to_center_current, 
				  loop_center_x, loop_center_y,
				  last_helix_dif_x, last_helix_dif_y);
	}
      }
      g_base_x[base] += last_helix_dif_x;
      g_base_y[base] += last_helix_dif_y;
      g_base_x[connected_to_current] += last_helix_dif_x;
      g_base_y[connected_to_current] += last_helix_dif_y;
      last_helix_dif_x = g_base_x[base] - g_undo1_base_x[base];
      last_helix_dif_y = g_base_y[base] - g_undo1_base_y[base];
      connected_to_last = connected_to_current;
      x_old_center_of_last = x_old_center_of_current;
      y_old_center_of_last = y_old_center_of_current;
      angle_old = angle_current;
      base = connected_to_last + 1;
    } else {
      g_base_x[base] = g_base_x[base] + last_helix_dif_x;
      g_base_y[base] = g_base_y[base] + last_helix_dif_y;
      base++;
    }
  }
  if (loop_num == 1) {
    if (g_connected[domain_end] < 1) {
      g_base_x[domain_end] += last_helix_dif_x;
      g_base_y[domain_end] += last_helix_dif_y;
    }
  }
}

void none_in_3_fix_last_2single(int domain_end, int last_paired_base_j, 
				float incoming_angle, float x0, float y0, 
				float distance_base_pairs) {
  /* if last helix in part 2 was moved at least one
   * distance_base_pairs right, space last 2 single stranded bases
   * after that helix and the closing base equally between closing
   * base and the base previous to the two. 
   */
  int             num_of_single;
  float           x1, y1;
  float           x2, y2;
  float           dif;
  float           x_mid, y_mid;
  int             stationary_base;
  float           distance_moved;
  if (domain_end <= (last_paired_base_j + 1))
    return;		/* no single bases to move */
  if (domain_end - last_paired_base_j > 2) {
    num_of_single = 2;
    stationary_base = domain_end - 3;
  } else {
    num_of_single = 1;
    stationary_base = domain_end - 2;
  }
  /* old position of that base */
  x1 = g_undo1_base_x[stationary_base];
  y1 = g_undo1_base_y[stationary_base];
  /* new positon of that base */
  x2 = g_base_x[stationary_base];
  y2 = g_base_y[stationary_base];
  distance_moved = distance(x1, x2, y1, y2);
  /* 1.4 value below is adjustable
   * Should probably be between .8 and 2.0
   * larger values make fewer changes to the structure 
   */
  if (distance_moved < 1.4 * distance_base_pairs)
    return;
  /* Rotate the points so that incoming helix is parallel to x axis.
   * Find their distance from the x-axis. If have changed by more than
   * 1.4*distance_bases away from the axis then fix the single
   * stranded bases. Find the amount of rotation. Measure difference
   from 0 to incoming angle */ 
  if (g_counter_clockwise)
    incoming_angle = PI - incoming_angle;
  else
    incoming_angle = 0.0 - incoming_angle;
  /* then find how far apart the bases are */
  dif = sin(incoming_angle) * (x1 - x0) + cos(incoming_angle) * (y1 - y0)
    - (sin(incoming_angle) * (x2 - x0) + cos(incoming_angle) * (y2 - y0));
  if ((dif <= 1.4 * distance_base_pairs) && (dif > -1.4*distance_base_pairs)) {
    return;
  }
  /* find middle between stationary base and domain end */
  if (num_of_single == 1) {	/* put the single base half way
				 * between domain_end and domain_end-2
				 */ 
    x_mid = (x2 + g_base_x[domain_end]) / 2.;
    y_mid = (y2 + g_base_y[domain_end]) / 2.;
    g_base_x[domain_end - 1] = x_mid;
    g_base_y[domain_end - 1] = y_mid;
  } else {	/* place 2 single stranded bases equidistant
		 * from domain end and previous base, domain_end-3 */
    x_mid = g_base_x[domain_end];
    y_mid = g_base_y[domain_end];
    g_base_x[domain_end - 2] = x_mid/3. + 2.*g_base_x[domain_end - 3]/3.;
    g_base_y[domain_end - 2] = y_mid/3. + 2.*g_base_y[domain_end - 3]/3.;
    g_base_x[domain_end - 1] = 2.*x_mid/3. + g_base_x[domain_end - 3]/3.;
    g_base_y[domain_end - 1] = 2.*y_mid/3. + g_base_y[domain_end - 3]/3.;
  }
}

void push_out_first_loop(int loop_num, int domain_end, int domain_start,
			 float distance_base_pairs) {
  /* for exterior loop, push last quadrant of loop farther out from loop */
  int             base, connected_to;
  float           percentage, dif_x, dif_y;
  float           angle;
  float           incoming_angle;
  float           angle_dif;
  float           distance_ends;
  /* for exterior loop, push last quadrant of loop farther out from loop
   */
  if (loop_num != 1)
    return;
  /* if first base and last are more than 2 distance_base_pairs apart,
     skip this step 
   */
  distance_ends = distance(g_base_x[domain_start], g_base_x[domain_end],
			   g_base_y[domain_start], g_base_y[domain_end]);
  if (distance_ends > 2.0 * distance_base_pairs)
    return;
  /* if sequence starts with a helix, skip this step */
  if (g_connected[domain_end] == domain_start)
    return;
  base = domain_end;
  if (g_connected[base] > domain_start) {
    if (!g_counter_clockwise)
      angle = g_structure_angle[base] + PI / 2.;
    else
      angle = g_structure_angle[base] - PI / 2.;
  } else {
    angle = g_structure_angle[base];
  }
  incoming_angle = angle;
  angle_dif = 0.0;
  while (angle_dif < (PI / 2.)) {
    percentage = 1.0 - angle_dif / (PI / 2.);
    dif_x = 3.0 * distance_base_pairs * cos(angle) * percentage;
    dif_y = 3.0 * distance_base_pairs * sin(angle) * percentage;
    g_base_x[base] += dif_x;
    g_base_y[base] += dif_y;
    connected_to = g_connected[base];
    if (connected_to >= domain_start) {
      g_base_x[connected_to] += dif_x;
      g_base_y[connected_to] += dif_y;
      base = connected_to - 1;
    } else
      base--;
    if (g_connected[base] > domain_start) {
      if (!g_counter_clockwise)
	angle = g_structure_angle[base] + PI / 2.;
      else
	angle = g_structure_angle[base] - PI / 2.;
    } else {
      angle = g_structure_angle[base];
    }
    angle_dif = angle_dif_loop_fix(angle, incoming_angle);
  }
}

void fix_loop_adjust_pieces(int loop_num, int domain_start, int domain_end,
			    float distance_bases, float distance_base_pairs, 
			    float distance_0, float distance_1, 
			    float distance_2, float distance_2p, int lhp1_i, 
			    int fhp2a_i, int lhp2a_i, int fhp2b_i, int lhp2b_i,
			    int fhp3_i, int lhp3_i) {
  /* find last helix previous to incoming base pair look in part 3
     position it correctly */ 
  float  x_old_center_of_current, y_old_center_of_current;
  float  x_old_center_of_last, y_old_center_of_last;
  float  angle_dif, angle_old, angle_current;
  int    base, connected_to_last;
  float  distance_to_center_current, distance_between, distance_to_center_last;
  int    single_bases;
  float  minimum_distance, loop_center_x, loop_center_y;
  float  distance_change_x, distance_change_y;
  float  other_position_of_last_helix_x;
  float  other_position_of_last_helix_y;
  float  new_position_of_last_helix_x, new_position_of_last_helix_y;
  float  extra_distance, extra_distance_x, extra_distance_y;
  int    number_of_single_bases, single_bases_is_odd;
  int    middle_base = -5000, stop_left, start_right, single_base;
  float  angle_of_singles;
  float  push_amount_x, push_amount_y;
  int    i;
  int    last_base_p2, first_base_p3;
  float  distance_change, distance_between_unmoved_helices, incoming_angle;
  float  original_incoming_angle, distance_change_angle, distance_to_center;
  float  distance_to_move_x, distance_to_move_y;
  float  x3, y3, x2, y2, x0, y0, dif_x, dif_y, dif, total_distance_last_moved;
  float  x_mid, y_mid;
  int    helix_in_part3;
  int    lhp3_j = -5000;
  int    base_m_1;
  if (lhp3_i > -1) {
		lhp3_j = g_connected[lhp3_i];
		helix_in_part3 = TRUE;
		angle_old = g_structure_angle[lhp3_j];
  } else {
    helix_in_part3 = FALSE;
    /* no helix in part3
     *
     * check if last helix in part 2 was moved at least one
     * distance_bases right. if it was, space the last 2 single
     * stranded bases before domain_end equally between  base
     * previous to them and domain_end
     */
    if (lhp2b_i == -1)
      return;	/* no helix in part2b */
    else
      angle_old = g_structure_angle[g_connected[lhp2b_i]];
  }
  /*
   * find out where last helix of part3 belongs relative to the
   * incoming helix
   */
  angle_current = g_structure_angle[domain_start];
  if (domain_start == g_start_base) {
    if (g_connected[domain_start] < domain_start) {
    } else {
      if (g_counter_clockwise)
	angle_current += PI / 2.;
      else
	angle_current -= PI / 2.;
    }
    if (g_counter_clockwise) {
      angle_old -= PI / 2.;
    } else {
      angle_old += PI / 2.;
    }
  } else {
    if (g_counter_clockwise) {
      angle_current -= PI / 2.;
      angle_old -= PI / 2.;
    } else {
      angle_current += PI / 2;
      angle_old += PI / 2.;
    }
  }
  incoming_angle = angle_current - PI;
  x0 = (g_undo1_base_x[domain_start] + g_undo1_base_x[domain_end]) / 2.;
  y0 = (g_undo1_base_y[domain_start] + g_undo1_base_y[domain_end]) / 2.;
  original_incoming_angle = incoming_angle;
  while (incoming_angle >= (3.0 * PI / 2.))
    incoming_angle -= 2.0 * PI;
  while (incoming_angle < (-PI / 2.))
    incoming_angle += 2.0 * PI;
  if (!helix_in_part3) {
    none_in_3_fix_last_2single(domain_end, g_connected[lhp2b_i], 
			       incoming_angle, x0, y0, distance_bases);
    return;
  }
  angle_dif = angle_dif_loop_fix(angle_old, angle_current);
  base = domain_end;
  connected_to_last = lhp3_j;
  x_old_center_of_last =
    (g_undo1_base_x[lhp3_i] + g_undo1_base_x[lhp3_j]) / 2.;
  y_old_center_of_last =
    (g_undo1_base_y[lhp3_i] + g_undo1_base_y[lhp3_j]) / 2.;
  loop_center_x = g_loop_center_x[loop_num];
  loop_center_y = g_loop_center_y[loop_num];
  if (g_connected[domain_start] == domain_end) {
    x_old_center_of_current =
      (g_undo1_base_x[domain_start] +
       g_undo1_base_x[domain_end]) / 2.;
    y_old_center_of_current =
      (g_undo1_base_y[domain_start] + g_undo1_base_y[domain_end]) / 2.;
  } else {
    x_old_center_of_current = g_undo1_base_x[domain_start];
    y_old_center_of_current = g_undo1_base_y[domain_start];
  }
  distance_to_center_current = 
    distance(x_old_center_of_current, loop_center_x, y_old_center_of_current,
	     loop_center_y);
  distance_to_center_last = distance(x_old_center_of_last, loop_center_x,
				     y_old_center_of_last, loop_center_y);
  distance_between = angle_dif *
    (distance_to_center_current + distance_to_center_last) / 2.0;
  single_bases = base - (connected_to_last + 1);
  if (single_bases == 0) {
    minimum_distance = distance_0;
  } else {
    if (single_bases == 1) {
      minimum_distance = distance_1;
    } else {
      if (single_bases == 2) {
	minimum_distance = distance_2;
      } else {
	single_bases = 3;
	minimum_distance = distance_2p;
      }
    }
  }
  other_position_of_last_helix_x = (g_base_x[lhp3_i] + g_base_x[lhp3_j]) / 2.;
  other_position_of_last_helix_y = (g_base_y[lhp3_i] + g_base_y[lhp3_j]) / 2.;
  if ((domain_end == g_end_base) && (g_connected[domain_end] < 1))
    extra_distance = distance_base_pairs;
  else
    extra_distance = 0.0;
  /* extra distance moves a single stranded base at the end away from
     the domain_start base */ 
  if (minimum_distance > distance_between) {	/* fix the helix */
    if (loop_num == 1) {
      if (g_connected[domain_end] < 1)
	g_base_x[domain_end] = g_undo1_base_x[domain_end];
      g_base_y[domain_end] = g_undo1_base_y[domain_end];
    }
    distance_change = minimum_distance - distance_between;
    distance_change = distance_change * .90;
    /* angle_dif/2 below is optional?
     * The PI makes change to avoid incoming helix
     * move the previous instead 
     */
    if (g_counter_clockwise)
      distance_change_angle = g_structure_angle[connected_to_last] - PI;
    else
      distance_change_angle = g_structure_angle[connected_to_last] + PI;
    if (g_counter_clockwise)
      distance_change_angle += angle_dif / 2;
    else
      distance_change_angle -= angle_dif / 2;
    distance_change_x = cos(distance_change_angle) * distance_change;
    distance_change_y = sin(distance_change_angle) * distance_change;
    extra_distance_x = extra_distance * cos(distance_change_angle);
    extra_distance_y = extra_distance * sin(distance_change_angle);
    g_base_x[lhp3_i] = g_undo1_base_x[lhp3_i] + distance_change_x +
      extra_distance_x;
    g_base_y[lhp3_i] = g_undo1_base_y[lhp3_i] + distance_change_y +
      extra_distance_y;
    g_base_x[lhp3_j] = g_undo1_base_x[lhp3_j] + distance_change_x +
      extra_distance_x;
    g_base_y[lhp3_j] = g_undo1_base_y[lhp3_j] + distance_change_y +
      extra_distance_y;
    base = domain_end;
    connected_to_last = lhp3_j;
    if (single_bases > 0) {
      if (single_bases == 1) {
	g_base_x[base - 1] = g_undo1_base_x[base - 1] +
	  distance_change_x / 2. + extra_distance_x;
	g_base_y[base - 1] = g_undo1_base_y[base - 1] +
	  distance_change_y / 2. + extra_distance_y;
      }
      if (single_bases == 2) {
	g_base_x[base - 1] = g_undo1_base_x[base - 1] +
	  distance_change_x / 3. + extra_distance_x;
	g_base_y[base - 1] = g_undo1_base_y[base - 1] +
	  distance_change_y / 3. + extra_distance_y;
	g_base_x[base - 2] = g_undo1_base_x[base - 2] +
	  2. * distance_change_x / 3. + extra_distance_x;
	g_base_y[base - 2] = g_undo1_base_y[base - 2] +
	  2. * distance_change_y / 3. + extra_distance_y;
      }
      if (single_bases == 3) {
	number_of_single_bases = base - 1 - connected_to_last;
	g_base_x[connected_to_last + 1] =
	  g_undo1_base_x[connected_to_last + 1] +
	  2. * distance_change_x / 3. + extra_distance_x;
	g_base_y[connected_to_last + 1] =
	  g_undo1_base_y[connected_to_last + 1] +
	  2. * distance_change_y / 3. + extra_distance_y;
	g_base_x[base - 1] =
	  g_undo1_base_x[base - 1] + distance_change_x / 3.
	  + extra_distance_x;
	g_base_y[base - 1] =
	  g_undo1_base_y[base - 1] + distance_change_y / 3. +
	  extra_distance_y;
	if (number_of_single_bases == 3) {
	  g_base_x[base - 2] = g_undo1_base_x[base - 2] +
	    distance_change_x / 2. + extra_distance_x;
	  g_base_y[base - 2] = g_base_y[base - 2] +
	    distance_change_y / 2. + extra_distance_x;
	} else {
	  if ((number_of_single_bases % 2) == 1) {
	    single_bases_is_odd = TRUE;
	    middle_base = (connected_to_last + base) / 2;
	    stop_left = middle_base - 1;
	    start_right = middle_base + 1;
	    g_base_x[middle_base] =
	      g_undo1_base_x[middle_base] + distance_change_x / 2. +
	      extra_distance_x;
	    g_base_y[middle_base] =
	      g_undo1_base_y[middle_base] +
	      distance_change_y / 2. + extra_distance_y;
	  } else {
	    single_bases_is_odd = FALSE;
	    stop_left = (connected_to_last + base) / 2;
	    start_right = stop_left + 1;
	  }
	  for (single_base = (connected_to_last + 2);
	       single_base <= stop_left; single_base++) {
	    g_base_x[single_base] =
	      g_undo1_base_x[single_base] + 2. *
	      distance_change_x / 3. + extra_distance_x;
	    g_base_y[single_base] =
	      g_undo1_base_y[single_base] + 2. *
	      distance_change_y / 3. + extra_distance_y;
	  }
	  for (single_base = start_right; single_base < (base - 1);
	       single_base++) {
	    g_base_x[single_base] =
	      g_undo1_base_x[single_base] + 1.
	      * distance_change_x / 3. + extra_distance_x;
	    g_base_y[single_base] =
	      g_undo1_base_y[single_base] + 1. *
	      distance_change_y / 3. + extra_distance_y;
	  }
	  /* push single bases out from loop */
	  if (g_counter_clockwise)
	    angle_of_singles = angle_current - angle_dif / 2.;
	  else
	    angle_of_singles = angle_current + angle_dif / 2.;
	  push_amount_x = .6 * distance_bases * cos(angle_of_singles);
	  push_amount_y = .6 * distance_bases * sin(angle_of_singles);
	  i = 0;
	  for (single_base = (connected_to_last + 2);
	       single_base <= stop_left; single_base++) {
	    i++;
	    g_base_x[single_base] += (i * push_amount_x);
	    g_base_y[single_base] += (i * push_amount_y);
	  }
	  i = 0;
	  for (single_base = (base - 2);
	       single_base >= start_right;
	       single_base--) {
	    i++;
	    g_base_x[single_base] += (i * push_amount_x);
	    g_base_y[single_base] += (i * push_amount_y);
	  }
	  if (single_bases_is_odd) {
	    i++;
	    g_base_x[middle_base] += (i * push_amount_x);
	    g_base_y[middle_base] += (i * push_amount_y);
	  }
	}
      }
    }
    if (extra_distance > 0) {
      g_base_x[domain_end] += extra_distance_x;
      g_base_y[domain_end] += extra_distance_y;
    }
  } else { /* Last helix did not initially require fixing. Move it
	      back where it came from */ 
    g_base_x[lhp3_i] = g_undo1_base_x[lhp3_i];
    g_base_y[lhp3_i] = g_undo1_base_y[lhp3_i];
    g_base_x[lhp3_j] = g_undo1_base_x[lhp3_j];
    g_base_y[lhp3_j] = g_undo1_base_y[lhp3_j];
    /* move the single stranded bases back */
    for (base = (lhp3_j + 1); base < domain_end; base++) {
      g_base_x[base] = g_undo1_base_x[base];
      g_base_y[base] = g_undo1_base_y[base];
    }
    if (loop_num == 1) {
      if (g_connected[domain_end] < 1) { /* move the end back too, it
					  * was moved out previously */
	g_base_x[domain_end] = g_undo1_base_x[domain_end];
	g_base_y[domain_end] = g_undo1_base_y[domain_end];
      }
    }
    /* check that bases between incoming helix and last helix are not
       too close together */ 
    number_of_single_bases = domain_end - 1 - lhp3_j;
    loop_center_x = g_loop_center_x[loop_num];
    loop_center_y = g_loop_center_y[loop_num];
    distance_to_center = distance(g_undo1_base_x[domain_end],
				  g_loop_center_x[loop_num],
				  g_undo1_base_y[domain_end],
				  g_loop_center_y[loop_num]);
    distance_between = angle_dif * distance_to_center;
    distance_between_unmoved_helices = distance_between -
      distance_base_pairs - distance_bases;
    if ((distance_between_unmoved_helices - number_of_single_bases *
	 distance_bases * .5) < 0) {
      loop_fix_spread_singles(lhp3_j, domain_end, number_of_single_bases,
			      distance_bases, distance_base_pairs, 
			      angle_current, angle_old, distance_to_center,
			      loop_center_x, loop_center_y, 0., 0.);
    }
  }
  /*
   * at this point the last helix should be positioned correctly
   * relative to the incoming helix
   *
   * Single stranded bases in between should be correct too.
   * Move other bases in 3rd quadrant out to fit last helix.
   * These are the bases between that last helix and part 2
   * new positions */
  new_position_of_last_helix_x =
    (g_base_x[lhp3_i] + g_base_x[lhp3_j]) / 2.;
  new_position_of_last_helix_y =
    (g_base_y[lhp3_i] + g_base_y[lhp3_j]) / 2.;
  /* how can they be fixed */
  distance_to_move_x = new_position_of_last_helix_x -
    other_position_of_last_helix_x;
  distance_to_move_y = new_position_of_last_helix_y -
    other_position_of_last_helix_y;
  base = fhp3_i;
  while (base < lhp3_i) {
    g_base_x[base] += distance_to_move_x;
    g_base_y[base] += distance_to_move_y;
    if (g_connected[base] > base)
      base = g_connected[base];
    else
      base++;
  }
  /*
   * make first base in third quadrant fit last base in second quadrant
   * by moving bases in part2 away from incoming helix
   *
   * find position of first base in part 3 */
  first_base_p3 = fhp3_i;
  x3 = g_base_x[first_base_p3];
  y3 = g_base_y[first_base_p3];
  /* find position of last base in part 2 */
  last_base_p2 = fhp3_i - 1;
  x2 = g_base_x[last_base_p2];
  y2 = g_base_y[last_base_p2];
  /* how far is incoming angle from PI/2. */
  incoming_angle = PI / 2. - incoming_angle;
  /* rotate so that incoming helix at at PI/2. */
  /* then find how far apart the bases are */
  dif = (sin(incoming_angle) * (x2 - x0) + cos(incoming_angle) * (y2 - y0))
    - (sin(incoming_angle) * (x3 - x0) + cos(incoming_angle) * (y3 - y0));
  if (dif < distance_bases) {	/* move outer part loop away from
				 * incoming helix */
    dif = distance_bases - dif;
    dif_y = dif * sin(original_incoming_angle);
    dif_x = dif * cos(original_incoming_angle);
    if (g_connected[lhp1_i] > 0)
      base = g_connected[lhp1_i] + 1;
    else {
      if (fhp2a_i > 0)
	base = fhp2a_i;
      else {
	if (fhp2b_i > 0)
	  base = fhp2b_i;
	else
	  return;
      }
    }
    while (base < first_base_p3) {
      g_base_x[base] += dif_x;
      g_base_y[base] += dif_y;
      if (g_connected[base] > base)
	base = g_connected[base];
      else
	base++;
    }
  }
  /* move single stranded bases between part2b and 3 to look better */
  total_distance_last_moved = sqrt(distance_to_move_x * distance_to_move_x +
				   distance_to_move_y * distance_to_move_y);
  if (total_distance_last_moved <= 1.5 * distance_base_pairs)
    return;
  base = first_base_p3 - 1;
  if (g_connected[base] >= 0) {	/* no single bases in between */
    return;
  }
  if ((fhp2b_i == -1) && (fhp2a_i == -1)) {
    return;		/* no helices here */
  }
  base_m_1 = base - 1;
  /* space the bases equally between neighbors */
  if (g_connected[base_m_1] >= 0) {	/* only 1 single between */
    g_base_x[base] = (g_base_x[base - 1] + g_base_x[base + 1]) / 2.;
    g_base_y[base] = (g_base_y[base - 1] + g_base_y[base + 1]) / 2.;
  } else {		/* 2 singles in between */
    x_mid = (g_base_x[first_base_p3] + g_base_x[first_base_p3 - 3]) / 2.;
    y_mid = (g_base_y[first_base_p3] + g_base_y[first_base_p3 - 3]) / 2.;
    g_base_x[base] = (x_mid + g_base_x[first_base_p3]) / 2;
    g_base_y[base] = (y_mid + g_base_y[first_base_p3]) / 2;
    g_base_x[base_m_1] = (x_mid + g_base_x[base - 2]) / 2;
    g_base_y[base_m_1] = (y_mid + g_base_y[base - 2]) / 2;
  }
}

void fix_bulge(int loop_num, int domain_start, int domain_end,
	       float distance_base_pairs, float distance_bases) {
  int             base;
  int             bulge_left;
  int             size;
  int             b;
  float           angle_b;
  float           radius;
  int             first_single, last_single;
  float           angle_current, incoming_angle, angle_dif;
  float           angle;
  float           angle_repos;
  float           new_mid_y, new_mid_x;
  int             connected_to;
  float           mid_x, mid_y;
  float           angle_change;
  float           loop_center_x, loop_center_y;
  bulge_left = FALSE;
  if (g_connected[domain_start + 1] < domain_start)
    bulge_left = TRUE;  /* for bulge_left base is first paired base on
			 * loop , (not domain start) */
  if (bulge_left) {
    size = 1;
    base = domain_start + 2;
    while (g_connected[base] < base) {
      size++;
      base++;
    }
    angle_current = g_structure_angle[base];
    first_single = domain_start + 1;
    last_single = base - 1;
  } else {
    size = 1;
    base = domain_end - 2;
    while (g_connected[base] < domain_start) {
      size++;
      base--;
    }
    first_single = base + 1;
    last_single = domain_end - 1;
  }
  if (size < 2)
    return;
  /* for bulge right base is last paired base on loop, (not domain end) */
  if (bulge_left)
    incoming_angle = g_structure_angle[domain_start];
  else
    incoming_angle = g_structure_angle[domain_end];
  while (incoming_angle >= (3.0 * PI / 2.))
    incoming_angle -= 2.0 * PI;
  while (incoming_angle < (-PI / 2.))
    incoming_angle += 2.0 * PI;
  angle_current = g_structure_angle[base];
  if (bulge_left) {
    angle_dif = angle_dif_loop_fix(incoming_angle, angle_current);
  } else {
    angle_dif = angle_dif_loop_fix(angle_current, incoming_angle);
  }
  if (g_counter_clockwise) {
    if (bulge_left)
      angle = incoming_angle + (angle_dif / 2);
    else
      angle = angle_current + angle_dif / 2;
  } else {
    if (bulge_left)
      angle = incoming_angle - (angle_dif / 2);
    else
      angle = angle_current - angle_dif / 2;
  }
  /* reposition base pair of helix coming off the loop */
  if (g_counter_clockwise)
    angle_repos = g_structure_angle[domain_start] + PI / 2.;
  else
    angle_repos = g_structure_angle[domain_start] - PI / 2.;
  mid_x = (g_base_x[domain_start] + g_base_x[domain_end]) / 2.;
  mid_y = (g_base_y[domain_start] + g_base_y[domain_end]) / 2.;
  new_mid_x = mid_x + distance_bases * 1.4 * cos(angle_repos);
  new_mid_y = mid_y + distance_bases * 1.4 * sin(angle_repos);
  connected_to = g_connected[base];
  g_base_x[base]=new_mid_x+distance_base_pairs/2.*cos(g_structure_angle[base]);
  g_base_y[base]=new_mid_y+distance_base_pairs/2.*sin(g_structure_angle[base]);
  g_base_x[connected_to] = new_mid_x + distance_base_pairs / 2. *
    cos(g_structure_angle[connected_to]);
  g_base_y[connected_to] = new_mid_y + distance_base_pairs / 2. *
    sin(g_structure_angle[connected_to]);
  /* circumference is size + 3 * distance_bases. radius is
     circumference/(2*pi) */ 
  radius = (size + 3) * distance_bases / (2.0 * PI);
  mid_x = (mid_x + new_mid_x) / 2.;
  mid_y = (mid_y + new_mid_y) / 2.;
  loop_center_x = mid_x + (radius + distance_base_pairs / 2.) * cos(angle);
  loop_center_y = mid_y + (radius + distance_base_pairs / 2.) * sin(angle);
  g_loop_center_x[loop_num] = loop_center_x;
  g_loop_center_y[loop_num] = loop_center_y;
  /* reposition center and single stranded bases */
  angle_change = (2.0 * PI) / (size + 3.);
  angle_change = angle_change * -1.;
  if (g_counter_clockwise)
    angle_b = angle - 2.0 * angle_change - PI;
  else
    angle_b = angle + 2.0 * angle_change - PI;
  for (b = first_single; b <= last_single; b++) {
    g_structure_angle[b] = angle_b;
    g_base_x[b] = loop_center_x + radius * cos(angle_b);
    g_base_y[b] = loop_center_y + radius * sin(angle_b);
    if (g_counter_clockwise)
      angle_b = angle_b - angle_change;
    else
      angle_b = angle_b + angle_change;
  }
}

void loop_fix(int domain_start, int domain_end, float distance_base_pairs, 
	      float distance_bases, int need_to_shift, int loop_num) {
  /* should spread helices around loop and rotate and translate
   * domains connected to them
   * need_to_shift: when true, runs on a predrawn loop
   * shifting must be done on all helices coming off from this loop
   * distance_base_pairs, distance_bases are figured,
   * when need_to_shift is false, this loop is being drawn
   * no shifting is needed, distance_base_pairs, distance_bases are
   * passed in
   *
   * Basic steps (Clockwise)
   *
   * 1. fix_loop_border_helices   determines which helices are in what
   * part of the loop. Part 1 is incoming helix to 90 degrees
   * clockwise. Part 2a is 90 degrees clockwise from part 1 part 2b is
   * 90 degrees clockwise from part 2a part 3 is 90 degrees clockwise
   * from 2b or 90 degrees counter-clockwise from incoming helix. This
   * accounts for all of the loop
   *
   * 2.  fix_loop_spread_helices travels around a loop clockwise moving
   * each helix away from the previous helix. Case 0, no bases between
   * helices. Case 1, 1 base between helices. Case 2, 2 bases between
   * helices. Case 3, At least 3 bases between helices. Helices are
   * moved differently for each case above. For case 3, if there is not
   * room for 3 single stranded bases between helices, the helices are
   * moved to make room for exactly 3 bases.  Extra bases are pushed
   * outward to make a sort of rectangle. If there is room for 3 single
   * stranded bases, the helices are not moved, but the single stranded
   * bases may be moved with loop_fix_spread_singles to make a rounded
   * or pointy looking strand. 3.  fix_loop_adjust_pieces fixes
   * problems created by 2. After number 2, the last helix is
   * positioned correctly relative the previous helix, but it may be a
   * long distance from where it belongs relative to the incoming
   * helix.  This function moves it and single stranded bases between
   * it and the incoming helix where belong relative to the incoming
   * helix.  Other helices and single stranded bases in part 3 are
   * moved to be correct relative to the incoming helix.  At this
   * point, part2b and part 3 no longer match. Part 2 is moved outward
   * from the incoming helix so that there is no overlap on the loop.
   * 4.  On a pre-existing loop shifting is done to move bases
   * connected to base pairs on the loop.  Parts 1 to 3 only move bases
   * on the loop.  Nothing moves the base pair on the incoming loop.
   *
   */
  int             lhp1_i, fhp2a_i, lhp2a_i, fhp2b_i, lhp2b_i, fhp3_i, lhp3_i;
  int             loop, base, start_of_fix, end_of_fix, base_fix;
  int             number_of_helices, bulge_loop_flag;
  float           x_translate, y_translate, y_center_of_rotation;
  float           x_center_of_rotation, old_x_translate, old_y_translate;
  float           angle_translate, angle_old, angle_new;
  float           distance_0, distance_1, distance_2, distance_2p;

  /* distance between 2 helices, with 0, 1, 2 or more than 2 single
     stranded  bases in between
     * avoid dumb case with loop 1 containing only 1 helix and no other
     * bases 
     */
  if ((domain_start == g_start_base) &&
      (g_connected[domain_start] == g_end_base))
    return;
  /* count the helices , stop at 2 */
  if (domain_start == g_start_base) {
    if (g_connected[domain_start] > domain_start) {
      number_of_helices = 1;
      base = g_connected[domain_start] + 1;
    } else {
      base = domain_start + 1;
      number_of_helices = 0;
    }
  } else {
    number_of_helices = 0;
    base = domain_start + 1;
  }
  while ((base < domain_end) && (number_of_helices < 2)) {
    if (g_connected[base] > base) {
      base = g_connected[base] + 1;
      number_of_helices++;
    } else {
      base++;
    }
  }
  if ((number_of_helices < 2) && (domain_start == g_start_base))
    return;
  if (need_to_shift) {
    loop_num = g_previous_loop[domain_start] + 1;
    if (domain_start == g_start_base)
      loop_num--;
  }
  /* find distance base pairs */
  if (need_to_shift) {
    if (loop_num != 1) {
      distance_base_pairs = 
	sqrt((g_base_x[domain_start] - g_base_x[domain_end]) *
	     (g_base_x[domain_start] - g_base_x[domain_end]) +
	     (g_base_y[domain_start] - g_base_y[domain_end]) *
	     (g_base_y[domain_start] - g_base_y[domain_end]));
    } else {
      base = g_start_base;
      while ((g_connected[base] < g_start_base) && (base < g_end_base)) {
	base++;
      }
      if (base < g_end_base) {
	distance_base_pairs = 
	  sqrt((g_base_x[base] - g_base_x[g_connected[base]]) *
	       (g_base_x[base] - g_base_x[g_connected[base]]) +
	       (g_base_y[base] - g_base_y[g_connected[base]]) *
	       (g_base_y[base] - g_base_y[g_connected[base]]));
      } else {
	return;
	/* no base pairs on loop, nothing to do */
      }
    }
  }
  /* determine minimum distances for types of helix adjacent to helix
   * with 0, 1, 2, or 2 or more bases between */
  distance_bases = DISTANCE_BASES * g_scale;
  distance_0 = distance_base_pairs + distance_bases;
  distance_1 = distance_0 + distance_bases;
  distance_2 = distance_1 + distance_bases;
  distance_2p = distance_2 + distance_bases;
  /* get the helices that determine the edges of the quadrants */
  bulge_loop_flag = FALSE;
  if ((loop_num > 1) && (number_of_helices == 1)) {
    if ((g_connected[domain_start + 1] > (domain_start + 1)) ||
	g_connected[domain_end - 1] > domain_start) {
      bulge_loop_flag = TRUE;
    }
  }
  if (bulge_loop_flag)
    fix_bulge(loop_num, domain_start, domain_end, distance_base_pairs,
	      distance_bases);
  if (!bulge_loop_flag)
    fix_loop_border_helices(domain_start, domain_end, &lhp1_i, &fhp2a_i, 
			    &lhp2a_i, &fhp2b_i, &lhp2b_i, &fhp3_i, &lhp3_i);
  if (!bulge_loop_flag)
    fix_loop_spread_helices(loop_num, domain_start, domain_end, distance_bases,
			    distance_base_pairs, distance_0, distance_1, 
			    distance_2, distance_2p);
  /* make adjustments to reconcile last helix with incoming helix */
  if (!bulge_loop_flag)
    fix_loop_adjust_pieces(loop_num, domain_start, domain_end, distance_bases,
			   distance_base_pairs, distance_0, distance_1, 
			   distance_2, distance_2p, lhp1_i, fhp2a_i, lhp2a_i, 
			   fhp2b_i, lhp2b_i, fhp3_i, lhp3_i);
  /* push out last 90 degrees of exterior loop */
  if (!bulge_loop_flag)
    push_out_first_loop(loop_num, domain_end, domain_start, 
			distance_base_pairs);
  /* need_to_shift is true when editing an existing structure
   *
   * Otherwise, return from here and extend helices connected to this loop.
   */
  if (!need_to_shift)
    return;
  base = domain_start + 1;
  while (base < domain_end) {
    if (g_connected[base] > base) {	/* translate */
      start_of_fix = base + 1;
      end_of_fix = g_connected[base] - 1;
      x_translate = -1.0 * (g_base_x[base] + g_base_x[g_connected[base]]) / 2
	+ (g_undo1_base_x[base] +
	   g_undo1_base_x[g_connected[base]]) / 2.;
      y_translate = -1.0 * (g_base_y[base] + g_base_y[g_connected[base]]) / 2
	+ (g_undo1_base_y[base] +
	   g_undo1_base_y[g_connected[base]]) / 2.;
      angle_new = g_structure_angle[base];
      angle_old = g_undo1_structure_angle[base];
      while (angle_new < 0)
	angle_new += 2.0 * PI;
      while (angle_old < 0)
	angle_old += 2.0 * PI;
      angle_translate = angle_new - angle_old;
      x_center_of_rotation = (g_undo1_base_x[base] +
			      g_undo1_base_x[g_connected[base]]) / 2.;
      y_center_of_rotation = (g_undo1_base_y[base] +
			      g_undo1_base_y[g_connected[base]]) / 2.;
      
      /* rotate */
      for (base_fix = start_of_fix; base_fix <= end_of_fix; base_fix++) {
	old_x_translate = g_undo1_base_x[base_fix] -
	  x_center_of_rotation;
	old_y_translate = g_undo1_base_y[base_fix] -
	  y_center_of_rotation;
	g_base_x[base_fix] = x_center_of_rotation +
	  cos(angle_translate) * old_x_translate +
	  -1.0 * sin(angle_translate) * old_y_translate;
	g_base_y[base_fix] = y_center_of_rotation +
	  1.0 * sin(angle_translate) * old_x_translate +
	  cos(angle_translate) * old_y_translate;
      }
      for (loop = (g_previous_loop[base] + 1);
	   loop <= g_previous_loop[g_connected[base]]; loop++) {
	old_x_translate = g_loop_center_x[loop] - x_center_of_rotation;
	old_y_translate = g_loop_center_y[loop] - y_center_of_rotation;
	g_loop_center_x[loop] = x_center_of_rotation +
	  cos(angle_translate) * old_x_translate +
	  -1.0 * sin(angle_translate) * old_y_translate;
	g_loop_center_y[loop] = y_center_of_rotation +
	  1.0 * sin(angle_translate) * old_x_translate +
	  cos(angle_translate) * old_y_translate;
      }
      /* translate */
      for (base_fix = start_of_fix; base_fix <= end_of_fix; base_fix++) {
	g_base_x[base_fix] = g_base_x[base_fix] - x_translate;
	g_base_y[base_fix] = g_base_y[base_fix] - y_translate;
	g_structure_angle[base_fix] = g_structure_angle[base_fix] 
	  + angle_translate;
      }
      for (loop = (g_previous_loop[base] + 1);
	   loop <= g_previous_loop[g_connected[base]]; loop++) {
	g_loop_center_x[loop] = g_loop_center_x[loop] - x_translate;
	g_loop_center_y[loop] = g_loop_center_y[loop] - y_translate;
      }
      base = g_connected[base];
    }
    base++;
  }
}

void loop_fix_first_pass(int domain_start, int domain_end, 
			 float distance_base_pairs, float distance_bases, 
			 int loop_num) {
  /* This loop_fix calles regular loop_fix after storing the positions
   * of the bases around the loop. The next loop_fix needs the
   * original positions
   */ 
  int             base;
  float          *temp_x, *temp_y;
  int             bases_to_store;
  if (loop_num == 1) {
    domain_start = g_start_base;
    domain_end = g_end_base;
  }
  if (domain_end < domain_start) { /* cannot handle this case yet */
    printf("Hit invalid loop %d with start = %d, end = %d\n", loop_num,
	   domain_start, domain_end);
    return;
  }
  bases_to_store = domain_end - domain_start + 1;
  temp_x = (float *) malloc((bases_to_store) * sizeof(float));
  temp_y = (float *) malloc((bases_to_store) * sizeof(float));
  /* note that too many bases are stored, but the code is simple */
  if (temp_x == NULL)
    report_mem_error("temp_x from loop_fix_first_pass");
  if (temp_y == NULL)
    report_mem_error("temp_y from loop_fix_first_pass");
  /* store the original undo values for x and y */
  /* temp should not be needed in a non-graphical mode */
  temp_x[domain_start - domain_start] = g_undo1_base_x[domain_start];
  temp_y[domain_start - domain_start] = g_undo1_base_y[domain_start];
  temp_x[domain_end - domain_start] = g_undo1_base_x[domain_end];
  temp_y[domain_end - domain_start] = g_undo1_base_y[domain_end];
  /* loop_fix needs base_x,base_y values in undo1 to keep track of how
   * bases are moved. This is crude here, but they are automatically
   * there when loop_fix is called interactively to fix a single loop
   * with the mouse 
   */ 
  g_undo1_base_x[domain_start] = g_base_x[domain_start];
  g_undo1_base_y[domain_start] = g_base_y[domain_start];
  g_undo1_base_x[domain_end] = g_base_x[domain_end];
  g_undo1_base_y[domain_end] = g_base_y[domain_end];
  base = domain_start + 1;
  while (base < domain_end) {
    temp_x[base - domain_start] = g_undo1_base_x[base];
    temp_y[base - domain_start] = g_undo1_base_y[base];
    g_undo1_base_x[base] = g_base_x[base];
    g_undo1_base_y[base] = g_base_y[base];
    if (g_connected[base] > base) {
      base = g_connected[base];
    } else
      base++;
  }
  loop_fix(domain_start, domain_end, distance_base_pairs, distance_bases, 
	   FALSE, loop_num);
  /* put the original undo values back */
  base = domain_start + 1;
  while (base < domain_end) {
    g_undo1_base_x[base] = temp_x[base - domain_start];
    g_undo1_base_y[base] = temp_y[base - domain_start];
    if (g_connected[base] > base) {
      base = g_connected[base];
    } else
      base++;
  }
  g_undo1_base_x[domain_start] = temp_x[domain_start - domain_start];
  g_undo1_base_y[domain_start] = temp_y[domain_start - domain_start];
  g_undo1_base_x[domain_end] = temp_x[domain_end - domain_start];
  g_undo1_base_y[domain_end] = temp_y[domain_end - domain_start];
  free(temp_x);
  free(temp_y);
}

void set_undo_data(void) {
  int             i;
  /* copy base_x, base_y structure_angle to two previous states */
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    g_undo1_base_x[i] = g_base_x[i];
    g_undo2_base_x[i] = g_base_x[i];
    g_undo1_base_y[i] = g_base_y[i];
    g_undo2_base_y[i] = g_base_y[i];
    g_undo1_structure_angle[i] = g_structure_angle[i];
    g_undo2_structure_angle[i] = g_structure_angle[i];
  }
  for (i = 1; i <= g_total_loops; i++) {
    g_undo1_loop_center_x[i] = g_loop_center_x[i];
    g_undo2_loop_center_x[i] = g_loop_center_x[i];
    g_undo1_loop_center_y[i] = g_loop_center_y[i];
    g_undo2_loop_center_y[i] = g_loop_center_y[i];
  }
  g_undo1_scale = g_scale;
  g_undo2_scale = g_scale;
  g_undo1_ss_distance_basepair = g_ss_distance_basepair;
  g_undo2_ss_distance_basepair = g_ss_distance_basepair;
  g_undo1_ss_distance_base_adj = g_ss_distance_base_adj;
  g_undo2_ss_distance_base_adj = g_ss_distance_base_adj;
  g_undo1_font_size = g_font_size;
  g_undo2_font_size = g_font_size;
  g_undo1_line_width = g_line_width;
  g_undo2_line_width = g_line_width;
}

void define_forced_bases(void) {
  char            rec[120];
  char            junk;
  int             i;
  int             base_1, base_2, helix_length, ignored_bp;
  printf("Using forced base file %s\n", g_forced_file);
  g_complot_forced_bases_count = 0;
  if ((forcefp = fopen(g_forced_file, "r")) == NULL) {
    printf("Error, Could not open %s\n", g_forced_file);
    return;
  }
  ignored_bp = 0;
  while (fgets(rec, 98, forcefp) != NULL) {
    sscanf(rec, "%c%d%d%d", &junk, &base_1, &base_2, &helix_length);
    /* check in bounds */
    base_1 = base_1 - g_history_offset - 1;
    base_2 = base_2 - g_history_offset + 1;
    for (i = 0; i < helix_length; i++) {
      base_1++;
      base_2--;
      if (((base_1 >= g_start_base) && (base_1 <= g_end_base)) &&
	  ((base_2 >= g_start_base) && (base_2 <= g_end_base)) &&
	  (g_connected[base_1] == base_2)) {
	if (base_1 < base_2) {
	  g_complot_forced_bases[g_complot_forced_bases_count] = base_1;
	  g_complot_forced_bases[g_complot_forced_bases_count + 1] = base_2;
	} else { 
	  g_complot_forced_bases[g_complot_forced_bases_count + 1] = base_1;
	  g_complot_forced_bases[g_complot_forced_bases_count + 2] = base_2;
	}
	if (g_complot_forced_bases_count < 498)
	  g_complot_forced_bases_count += 2;
	else {
	  printf("More than allowed 498 bases in forced file\n");
	  printf("Extras ignored\n");
	}
      } else {
	printf("forced file data row %s does not match input data\n", rec);
	ignored_bp++;
      }
    }
  }
  printf("%s annotated %d forced base pairs and ignored %d.\n", PROGRAM_NAME,
	 g_complot_forced_bases_count/2, ignored_bp);
  fclose(forcefp);
}

int check_ss_filetype(char *filename) {
  int             len;
  len = strlen(filename);
  if ((filename[len-1]=='s')&&(filename[len-2]=='s')&&(filename[len-3]=='.')) {
    printf("Input file treated as ss file\n");
    return TRUE;
  }
  return FALSE;
}

void display_flags(void) {
  printf("Usage: %s [options] infile[.ct] | infile.ss\nOptions:\n",
	 PROGRAM_NAME);
  printf("-v, -V\t(display version information)\n");
  printf("-h,\t(display this information)\n");
  printf("-a,\t(Output all tangent & perpendicular angles in degrees)\n");
  printf("-ab\t(Annotate bases for -pnum, -prob or -ss-count)\n"),
  printf("-ad\t(Annotate dots for -pnum, -prob or -ss-count)\n");
  printf("-af <af_file>\t(Specify annotation file: af_file.ann\n");
  printf("\t\t Default is infile.ann for -pnum,\n");
  printf("\t\t infile.ss-count for -ss-count and\n");
  printf("\t\t infile.ann or infile_k.ann for -prob, k = 1, 2, ... )\n");
  printf("-aj\t(Show sequence name and free energy in title bar)\n");
  printf("-ar\t(AutoRotate structure for best fit.)\n");
  printf("-bp\t(Display 1x1 & 2x2 i-loops as if mismatched bases pair)\n");
  printf("-c \t(Circular RNA or DNA. Default is linear)\n");
  printf("-cc\t(Draw counter-clockwise. Default is clockwise)\n");
  printf("-cg\t(Draw a circle graph. Output file is infile_cir.ps)\n");
  printf("-col <color_file>\t(User supplied color file.)\n");
  printf("-col_ann <color_ann_file>\t(As above for annotation colors)\n");
  printf("-d <diam>\t(Diameter of circle when -cg is set, max = 6in)\n");
  printf("-D\t(Use dots, not lines for structure base pairs)\n");
  printf("-e <i, j>\t(Draw circle graph or structure for domain excluded by base pair i,j)\n");
  printf("-ew <i, j>\t(As above, but use entire circle for the domain in untangle mode.)\n");
  printf("-f\t(Draw a flat exterior loop\n");
  printf("\tMay be used with -i or -iw to make any loop flat)\n");
  printf("-fa\t(Same as -f, but helices alternate between\n");
  printf("\tclockwise and counter clockwise (up & down).\n");
  printf("-fix\t(In untangle mode, loops are \"fixed\" by spreading.)\n");
  printf("-force <force_file>\t(Use diamonds to annotate base pairs from force_file\n");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  printf("-g wid\t(Create structure plot in gif format.\n");
  printf("\t wid = width of image in pts. 50 <= wid <= 3200.\n");
  printf("\t height = wid*17./13. Default: %dx%d)\n", 
	 (int) (WD + 72.0), (int) (HT + 72.0));
#endif
  printf("-i <i, j>\t(Draw circle graph or structure for included domain: i...j)\n");
  printf("-iw <i, j>\t(As above, but use entire circle for the domain in untangle mode.)\n");
#if HAVE_LIBJPEG
  printf("-jpg <wid>\t(Create structure plot in jpg, [jpeg] format.\n");
  printf("\t wid = width of image in pts. 50 <= wid <= 3200.\n");
  printf("\t height = wid*17./13. Default: %dx%d)\n",
	 (int) (WD + 72.0), (int) (HT + 72.0));
#endif
  printf("-lab n\t(Label every nth base. 0 for none. Default varies with sequence length)\n");
  printf("-lij <i j>\t(Force labels on bases i & j, paired or not\n");
  printf("-loop\t(Label all structure loops)\n");
  printf("-m\t(Show arc midpoints in circle graph)\n");
  printf("-mo\t(Show arc midpoints only in circle graph)\n");
  printf("-n\t(Generate structure based on natural angles)\n");
  printf("-ni\t(Create a non-interlaced image for -png)\n");
  printf("-o <out_file>\t(Specify output file prefix. Default is infile. Output files are out_file_cir.ps or out_file.suffix, where suffix = ps, gif, jpg or png as appropriate.)");
  printf("-outline\t(Do not show bases in structure)\n");
  printf("-p\t(PostScript output: suffix = ps)\n");
#if HAVE_LIBPNG
  printf("-png [<width>]\t(Create png output of structure\n");
  printf("\t\t height = width*17./13.,\n");
  printf("\t\t Default: %d pts, 50 < width < 3200)\n", (int) (WD + 72.0));
#endif
  printf("-pnum\t(Use p-num annotation  from .ann file)\n");
  printf("-prob\t(Use probability annotation from .ann file)\n");
  printf("-r <deg>\t(Degrees of circle for untangle mode\n");
  printf("\t\t Default: 356)\n");
  printf("-reg <deg>\t(Round all helix angles to integral multples of 'deg' ");
  printf("\t\t around loops. Must be integral. 1 <= deg <= 90\n");
  printf("-rot <deg>\t(Rotate structure deg degrees counterclockwise)\n");
  printf("-s\t(Color base pair connecting lines or dots as follows:\n");
  printf("\t GC / CG: red; AU / UA: blue; GU / UG: green; other: yellow)\n");
  printf("-ss\t(Create infile.new.ss in ss format.\n");
  printf("\t Does not work in circle graph mode.)\n");
  printf("-ss-count\t(Use ss-count annotation from .ss-count file)\n");
  printf("-t <\"TITLE\">\t(Specify a title. Protect with quotes\n");
  printf("\t\t Default title: ct file header or ss file name.)\n");
  printf("-tab_ann_html\t(Create html format table for annotation\n");
  printf("\t\t Works with -pnum, -prob & -ss-count\n");
  printf("\t\t Filename is infile.'annotation-type'.col.html)\n");
  printf("-x\t(Create infile.png2bp of base pair locations ");
  printf("for png, jpg or gif images: (i,j) base pair <--> (x,y) in pts.)\n");
  printf("-z <ang>\t(Circle graph angle for bases not in multiloop\n");
  printf("\t\t 0 <= ang <= 2.0. Default: 1.0\n");
  printf("\t\t Use 0.1 to advance non multiloop bases 10%s\n", "%");
  printf("\t\t of the advance given to multiloop bases.)\n");
  printf("-zoom <s, x, y>\t(Zooms about point x,y at magnification s\n");
  printf("\t\t x,y is a png/jpg/gif coordinate, (0,0) is top left\n");
  printf("\t\t s = 1.0 is normal size, s = 5.0 magnifies 5 times\n");
  printf("\t\t s = 0.5 makes image half as big as normal.\n");
  printf("\t\t This option works in non-interactive mode only.)\n");
  printf("-zoom_ps <s, x, y,>\t(Zooms about point x,y at magnification s.\n");
  printf("\t\t x,y is a PostScript coordinate, (0,0) is bottom left.\n");
  printf("\t\t For PostScript ouput in non-interactive mode only.)\n");
  /*  printf("See sir_graph.doc for more information\n"); */
  try_exit(0);
}

void select_annotation_type_ng(int choice) {
  /* same function as select_annotation_type except this one does not
   * redraw the interactive windows 
   */ 
  int  length = 0;
  if (choice == NONE) {
    g_annotation = NONE;
    return;
  }
  if (choice == g_annotation)
    return;
  else
    g_annotation = choice;
  if (choice == P_NUM) {
    length = define_map_ann(g_annotation_filename, TRUE, FALSE, TRUE, g_unused,
			    g_create_ann_table, g_annotation_table_type,
			    g_length);
    set_extra_ps_colors(FALSE);
  }
  if (choice == SS_COUNT) {
    length = define_map_ann(g_annotation_filename, FALSE, FALSE, TRUE, 
			    g_unused, g_create_ann_table, 
			    g_annotation_table_type, g_length);
    set_extra_ps_colors(FALSE);
  }
  if (choice == PROB) {
    length = define_map_ann(g_annotation_filename, TRUE, TRUE, TRUE, g_unused,
			    g_create_ann_table, g_annotation_table_type,
			    g_length);
    set_extra_ps_colors(TRUE);
  }
  if (length != g_length) {
    g_annotation = NONE;
    if (length < 1)
      printf("Warning!\tNonsense length %d from annotation file.\n", length);
    else
      printf("Warning!\tLength %d from annotation file does not match %d\n",
	     length, g_length);
    printf("%s continues without annotation.\n", PROGRAM_NAME);
  }
}

void regularize_all_angles(void) {
  int       loop_num;
  float     current_x, current_y;
  float     distance_base_pairs;
  float     last_angle_perp = 0.0;
  int       base;
  float     original_angle, angle_change;
  float     angle_adjust = 0.0; /* M. Zuker: angle_adjust was not
				 * initialized in the original!
				 * Caused this feature to fail under Mac OS X 
				 */
  float     counter_clockwise_angle;
  int       extra_base_exists;
  int       i;
  int       domain_start, domain_end;
  int       helix_started;
  float     distance_bases;
  if (g_counter_clockwise)
    counter_clockwise_angle = -PI;
  else
    counter_clockwise_angle = 0.0;
  domain_start = g_start_base;
  domain_end = g_end_base;
  loop_num = 1;
  while (g_connected[domain_start] < domain_start)
    domain_start++;
  domain_end = g_connected[domain_start];
  current_x = (g_base_x[domain_start] + g_base_x[domain_end]) / 2;
  current_y = (g_base_y[domain_start] + g_base_y[domain_end]) / 2;
  /* find distance base pairs */
  last_angle_perp = g_structure_angle[domain_start] + PI / 2;
  /*
   * This is the first loop and the base is paired, rotate to give the
   * first helix the correct angle
   */
  if (!g_counter_clockwise)
    original_angle = g_structure_angle[domain_start] - PI / 2.;
  else
    original_angle = g_structure_angle[domain_start] + PI / 2.;
  while (original_angle > (2.0 * PI))
    original_angle -= (2.0 * PI);
  while (original_angle < 0)
    original_angle += (2.0 * PI);
  angle_change = regularize_angle(original_angle) - original_angle;
  rotate_structure(angle_change, (WD + 72.0)/2.0, (HT + 72.0)/2.0, TRUE, 1, 
		   g_total_loops, g_start_base, g_end_base);
  for (i = g_start_base; i <= g_end_base; i++) {
    g_structure_angle[i] = g_structure_angle[i] + angle_change;
  }
  /*   scale_to_fit_window(); M. Zuker removes this 27/11/06
   *   Caused tiny font problem when -f and -reg were used together
   */
  distance_base_pairs = 
    sqrt((g_base_x[domain_start] - g_base_x[domain_end]) *
	 (g_base_x[domain_start] - g_base_x[domain_end]) +
	 (g_base_y[domain_start] - g_base_y[domain_end]) *
	 (g_base_y[domain_start] - g_base_y[domain_end]));
  distance_bases = DISTANCE_BASES * g_scale;
  store_for_undo();
  if (g_end_base == g_connected[domain_start])
    extra_base_exists = FALSE;
  else
    extra_base_exists = TRUE;
  domain_start = g_start_base ;
  domain_end = g_end_base;
  /* for first helix, let domain start after it, finish before it */
  while (g_connected[domain_start] < domain_start)
    domain_start++;
  domain_end = domain_start - 1;
  current_x = 0.0;	/* not used though */
  current_y = 0.0;
  /* the traverse loop below should fix all helices on the first loop */
  /*  if ( (loop_num>1) || (!g_flat && !g_flat_alternate)) */
  /*  if ( (loop_num>1) || !g_flat ) */
  /*  if ( !g_flat )  */
  /* M. Zuker. There has been a bug since I allowed automatic all
  angle regularize with the natural angle option. The first helix does
  not fit smoothly into the exterior loop.
  */
  if ( !g_flat )
    traverse_loop_structure_regularize_loop(loop_num, domain_start, domain_end,
					    g_connected, last_angle_perp + 
					    angle_change, current_x, current_y,
					    extra_base_exists, g_end_base,
					    distance_base_pairs, FALSE, 0.0);
  base = g_start_base;
  domain_end = g_end_base;
  if (g_connected[g_start_base] > g_start_base) {
    helix_started = TRUE;
    current_x = (g_base_x[base] + g_base_x[g_connected[base]]) / 2.;
    current_y = (g_base_y[base] + g_base_y[g_connected[base]]) / 2.;
    last_angle_perp = g_structure_angle[base] - PI / 2.0
      + counter_clockwise_angle;
  } else
    helix_started = FALSE;
  base++;
  while (base < domain_end) {
    if (g_connected[base] < base) {
      if (helix_started) {
	loop_num++;
	traverse_loop_structure_regularize_loop(loop_num, base - 1, 
						g_connected[base-1], 
						g_connected, last_angle_perp + 
						counter_clockwise_angle, 
						current_x, current_y, FALSE, 
						g_end_base,distance_base_pairs,
						FALSE, angle_adjust);
	helix_started = FALSE;
      }
    } else {
      if (g_connected[base] == (g_connected[base - 1] - 1)) {	
	/* continue helix draw next base pair */
	if (g_counter_clockwise) {
	  current_x = current_x + distance_bases *
	    cos((double) (last_angle_perp));
	  current_y = current_y + distance_bases *
	    sin((double) (last_angle_perp));
	} else {
	  current_x = current_x + distance_bases *
	    cos((double) last_angle_perp);
	  current_y = current_y + distance_bases *
	    sin((double) last_angle_perp);
	}
	draw_basepair(base, g_connected[base], last_angle_perp, current_x,
		      current_y, distance_base_pairs);
      } else {/* first base pair of helix */
	if (helix_started == TRUE) {
	  loop_num++;
	  traverse_loop_structure_regularize_loop(loop_num, base - 1, 
						  g_connected[base - 1], 
						  g_connected, last_angle_perp 
						  + counter_clockwise_angle, 
						  current_x, current_y, FALSE,
						  g_end_base,
						  distance_base_pairs, FALSE, 
						  angle_adjust);
	}
	helix_started = TRUE;
	current_x = (g_base_x[base] + g_base_x[g_connected[base]]) / 2.0;
	current_y = (g_base_y[base] + g_base_y[g_connected[base]]) / 2.0;
	if (g_counter_clockwise)
	  last_angle_perp = g_structure_angle[base] + PI / 2.0;
	else
	  last_angle_perp = g_structure_angle[base] - PI / 2.0;
	if (g_counter_clockwise)
	  angle_adjust = 
	    find_dif_angles(last_angle_perp, g_undo1_structure_angle[base] + 
			    PI / 2.0);
	else
	  angle_adjust = 
	    find_dif_angles(last_angle_perp, g_undo1_structure_angle[base] - 
			    PI / 2.0);
      }
    }
    base++;
  }
  set_label_coordinates();
}

/* On the Mac, libglut does not exist separately. Search instead for
   the function glutInit. Under Linux, glut.h should be in
   /usr/include/GL and /usr/include is automatically searched.
*/

#ifdef HAVE_GLUTINIT
#if HAVE_LIBGLUT
# include <GL/glut.h>
#else
# include <glut.h>
#endif
#endif

char      g_input_main_message[80];
char      g_input_default_message[80];
char      g_input_name[80];
char      g_input_abort_message[80];
char      g_input_true_input[80];
char      g_name_string[140]; /* Jan23-07 MZ Copy of name_string to be
				 used when input window is created */
int       g_window_1; /* for zoom window */
int       g_window_nozoom; /* for unzoomed window */
int       g_current_window; /* M. Zuker April 15, 2009. */
int       g_window_text;
int       g_window_input;
int       g_window_input_exists = FALSE; /* Jan23-07 MZ */
int       g_window_1_width ;
int       g_window_1_height ;
int       g_window_nozoom_width ;
int       g_window_nozoom_height ;
int       g_img_is_main = FALSE ; /*true for image selected from main window*/
int       g_input_type;
int       g_post_is_zoom;
int       g_clicked_row;
int       g_clicked_column;
/* set by mouse menu, single strand stretch a single stranded region */
int       g_mouse_strand_stretch;	
/* set by mouse menu, fix loop by spreading helices and single strand bases */
int       g_mouse_loop_fix; 
/* set by mouse menu, single strand linear. From a base in a single
 * stranded region, space other bases linear from the base to helix on
 * both sides */ 
int       g_mouse_strand_linear; 
int       g_mouse_strand_linear_base; /* the base to be	dragged */
float     g_mouse_strand_linear_x;
float     g_mouse_strand_linear_y;
int       g_mouse_strand_linear_loop;
/* set by g_strand_stretch and g_strand linear
 * when a single stranded region has been selected by left mouse button */
int       g_strand_set;
float     g_mouse_strand_stretch_first_x;
float     g_mouse_strand_stretch_first_y;
float     g_mouse_strand_stretch_center_x;
float     g_mouse_strand_stretch_center_y;
/* set by mouse menu loop regularize.
 * When true a clicked loop is drawn with helix angles regularized */
int       g_mouse_loop_regularize; 
/* set by mouse menu Loop Natural.
 * When true a clicked loop is drawn with natural angles */
int       g_mouse_loop_natural;	
/* set by mouse menu: edit when true a helix can be rotated */
int       g_mouse_rotate_helix;	
float     g_mouse_rotate_helix_center_x;
float     g_mouse_rotate_helix_center_y;
float     g_mouse_rotate_helix_last_x;
float     g_mouse_rotate_helix_last_y;
/* Set by mouse edit. When true, a helix or base can be dragged */
int       g_mouse_translate;	
int       g_mouse_translate_pos_set;
int       g_mouse_translate_base;
int       g_mouse_loop_stretch = TRUE;
int       g_mouse_loop_stretch_pos_set = FALSE;
int       g_mouse_loop_stretch_pos_x;
int       g_mouse_loop_stretch_pos_y;
int       g_mouse_loop_stretch_loop;
int       g_mouse_loop_stretch_domain_row;
int       g_mouse_loop_stretch_domain_column;
int       g_mouse_loop_stretch_domain_row1;
int       g_mouse_loop_stretch_domain_column1;
int       g_mouse_loop_stretch_domain_row2;
int       g_mouse_loop_stretch_domain_column2;
float     g_mouse_loop_stretch_center_x;
float     g_mouse_loop_stretch_center_y;
float     g_mouse_loop_stretch_r;
float     g_mouse_loop_stretch_start_r;
/* set by mouse function zoom.
 * when true, zooms on a clicked point in the zoom window */
int       g_mouse_display_zoom_mode;
float     g_display_window_zoom_factor;
float     g_display_window_x_shift;
float     g_display_window_y_shift;

int call_file_makers(int, char **) ; 

int main(int argc, char **argv) {
  
  if (argc == 1) {
    printf("Error. Type '%s -h' for help.\n", PROGRAM_NAME);
    try_exit(18);
  } else if (argc == 2) {
    if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "-V") == 0) {
      printf("%s:\t%s\n", PROGRAM_NAME, PACKAGE_STRING);
      try_exit(0);
    } else if (strcmp(argv[1], "-h") == 0 ) {
      display_flags();
    }
  }
  if ( argv[argc - 1][0] == '-') {
    printf("Invalid file name. Enter '%s -h' for help.\n", PROGRAM_NAME);
    try_exit(19);
  }
  g_start_base = -1;
  g_bp = FALSE;
  g_is_circular = FALSE;
  g_angle_list = FALSE;
  g_interactive_mode = TRUE;
  g_midpoints = FALSE;
  g_structure = TRUE;
  g_natural = FALSE;
  g_post_is_zoom = FALSE;
  g_mouse_loop_stretch = FALSE;
  g_loop_labels = FALSE;
  g_auto_rotate = FALSE;
  g_rotation_angle = 0;
  g_degrees_used = 356.0;
  g_whole_circle = FALSE;
  g_arcs = TRUE;
  g_lines = TRUE;
  g_outline_mode = FALSE;
  g_label_frequency = -50;
  g_img_interlaced = TRUE;
  g_small_angle = FALSE;
  g_ss_mode = FALSE;
  g_helix_base_advance = 1.0;
  g_img_mode = FALSE;
  g_img_is_main = TRUE;
  g_img_width = (int) (WD + 72.0);
  g_img_height = (int) ( ((float) g_img_width) * 17.0/13.0 );
  g_create_ann_table = FALSE;
  g_ss_output_mode = FALSE;
  g_ss_backwards = FALSE;
  g_scale = 1.0;
  g_label_forced_row = 0;
  g_label_forced_column = 0;
  g_external_domain_flag = FALSE;
  g_annotation = NONE;
  g_annotation_bases = TRUE;
  g_annotation_dots = TRUE;
  g_counter_clockwise = FALSE;
  g_aj = FALSE;
  g_command_line_x = FALSE;
  g_command_line_zoom_flag = FALSE;
  g_command_line_zoom_ps_flag = FALSE;
  call_file_makers(argc, argv);
  try_exit(20);
  return 0;	/* this line is never hit */
}

/* End of main */
#if HAVE_GLUTINIT
void inputSetWindow(void);

void make_structure_ps_interactive_start(void) {
  g_input_type = POSTSCRIPT;
  inputSetWindow();
  strcpy(g_input_name, g_structure_filename);
  strcat(g_input_name, ".ps");
  g_input_true_input[0] = '\0';
  glutShowWindow();
  glutSwapBuffers();
  glFlush();
  glutPostRedisplay();
}

void make_structure_ps_interactive_finish(void) {
  char            filename[100];
  float           x_shift, y_shift, zoom_factor;
  glutDestroyWindow(g_window_input);
  g_window_input_exists = FALSE;
  if (strlen(g_input_true_input)==1 && g_input_true_input[0]=='/')
    return;
  if (strlen(g_input_true_input) == 0)
    strcpy(filename, g_structure_filename);
  else
    strcpy(filename, g_input_true_input);
  /* chop .ps off the end */
  chop_suffix(filename, ".ps");
  if (g_post_is_zoom) {
    x_shift = (WD + 72.0)/2.0 - g_display_window_x_shift;
    y_shift = (HT + 72.0)/2.0 - g_display_window_y_shift;
    zoom_factor = g_display_window_zoom_factor;
  } else {
    x_shift = (WD + 72.0)/2.0;
    y_shift = (HT + 72.0)/2.0;
    zoom_factor = 1.0;
  }
  make_structure_ps_real(filename, g_smart_colors, g_first_line, 
			 g_line_width, g_oldconnect, g_start_base, 
			 g_end_base, g_total_loops, g_loop_center_x, 
			 g_loop_center_y, g_font_size, g_outline_mode, 
			 g_scale, g_bases, g_base_x, g_base_y, 
			 g_lines, g_structure_label_x, 
			 g_structure_label_y, g_structure_label_value,
			 g_loop_labels, g_total_labels,
			 g_external_domain_flag, g_ex_start_base, 
			 g_ex_end_base, g_history_offset, g_annotation,
			 g_ann_to_color, g_annotation_dots, 
			 g_annotation_bases, x_shift, y_shift, zoom_factor, 
			 g_is_circular, g_complot_forced_bases,
			 g_complot_forced_bases_count);
}

void            make_structure_img_interactive_finish(int);
void            create_output_ss_interactive_finish(void);

void keyboard_input(unsigned char key, int x, int y) {
  int             string_len1;	/* accepts keys from text window */
  x = y;
  y = x;		/* get rid of warning message */
  if (key == 13) {	/* on <RETURN> mark end of text */
    if (g_input_type == POSTSCRIPT)
      make_structure_ps_interactive_finish();
    if (((g_input_type == GIF) || (g_input_type == PNG)) ||
	g_input_type == JPG)
      make_structure_img_interactive_finish(g_input_type);
    if (g_input_type == SS)
      create_output_ss_interactive_finish();
    return;
  }
  string_len1 = (int) strlen(g_input_true_input);
  if (string_len1 > 0)
    if ((key == 8)) {	/* process backspace or delete */
      g_input_true_input[string_len1 - 1] = '\0';
      inputSetWindow();
      glutPostRedisplay();
      return;
    }
  if (key != ' ') {
    g_input_true_input[string_len1] = key;	/* add key to end of input */
    g_input_true_input[string_len1 + 1] = '\0';
  }
  inputSetWindow();	/* show the text */
  glutPopWindow();
  glutPostRedisplay();
}

/* ************************************************************ *
 * Start of interactive window functions                        *
 * ************************************************************ *
 */

void sniff_for_opengl_errors(void) {
  GLenum          errCode;
  const GLubyte  *errString;
  const char     *my_string;
  if ((errCode = glGetError()) != GL_NO_ERROR) {
    errString = gluErrorString(errCode);
    my_string = (char *) errString;
    fprintf(stderr, "OpenGL Error: %s\n", my_string);
  }
}

void fixcolor_plot(int color) {
  /* for the screen */
  if ((color >= 0) && (color < MAIN_COLORS))
    glColor3f(main_color_table_float[color].red,
	      main_color_table_float[color].green,
	      main_color_table_float[color].blue);
}

void fixcolor_plot_ann(int color) {
  if ((color >= 0) && (color <= NUM_COLORS)) {
    glColor3f(color_table_ps[color].red,
	      color_table_ps[color].green,
	      color_table_ps[color].blue);
  }
}

void set_color_base_window(char base1, char base2) {
  /* Zuker adds 2 lines below. May 23, 2003. */
  base1 = toupper(base1);
  base2 = toupper(base2);
  if (base1 == 'T') base1 = 'U';
  if (base2 == 'T') base2 = 'U';
  if (((base1=='G') && (base2 == 'C')) || ((base1 == 'C') && (base2 == 'G'))) 
    fixcolor_plot(COLOR_CONNECTING_GC);
  else if (((base1=='A') && (base2=='U')) || ((base1=='U') && (base2 == 'A')))
    fixcolor_plot(COLOR_CONNECTING_AU);
  else if (((base1=='G') && (base2=='U')) || ((base1=='U') && (base2 == 'G')))
    fixcolor_plot(COLOR_CONNECTING_GU);
  else
    fixcolor_plot(COLOR_CONNECTING_OTHERS);
}

void window_draw_int(float x, float y, int number, int test_3) {
  int             i;
  int             len;
  float           window_font_size;
  char            number_string[12];
  sprintf(number_string, "%d", number);
  if (test_3)
    strcat(number_string, " 3'");
  len = strlen(number_string);
  glPushMatrix();
  window_font_size = (WD + 72.0) * g_font_size/112320.0 ;
  glTranslatef(x, y, 0.0);
  glScalef(window_font_size, window_font_size, window_font_size);
  for (i = 0; i < len; i++)
    glutStrokeCharacter(GLUT_STROKE_ROMAN, number_string[i]); glPopMatrix();
}

void window_draw_char(float x, float y, char *number_string) {
  int             i;
  int             len;
  float           window_font_size;
  len = strlen(number_string);
  glPushMatrix();
  window_font_size = (WD + 72.0) * g_font_size/112320.0;
  glTranslatef(x, y, 0.0);
  glScalef(window_font_size, window_font_size, window_font_size);
  for (i = 0; i < len; i++)
    glutStrokeCharacter(GLUT_STROKE_ROMAN, number_string[i]);
  glPopMatrix();
}

void output_text(int x, int y, char *string) {
  int             len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
  }
}

void zoom_entry(int state) {
  if (state == GLUT_ENTERED) {
    g_current_window = g_window_1;
    g_img_is_main = FALSE;
  }
}

void nozoom_entry(int state) {
  if (state == GLUT_ENTERED) {
    g_current_window = g_window_nozoom;
    g_img_is_main = TRUE;
  }
}

void window_output_text(float x, float y, char *text) {
  int             len, i;
  printf("Printing %s at %f,%f\n", text, x, y);
  glRasterPos2f(x + 150, y + 150);
  len = strlen(text);
  for (i = 0; i < len; i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, text[i]);
}

void show_text(void) {
  glutSetWindow(g_window_text);
  glutPostRedisplay();
}

void show_zoom(void) {
  glutSetWindow(g_window_1);
  glutPostRedisplay();
}

void show_nozoom(void) {
  glutSetWindow(g_window_nozoom);
  glutPostRedisplay();
}

void show_both(void) {
  show_zoom();
  show_nozoom();
}

void display_window_input(void) {
  char            text[100];
  int             i;
  glClearColor(main_color_table_float[COLOR_BACKGROUND].red,
	       main_color_table_float[COLOR_BACKGROUND].green,
	       main_color_table_float[COLOR_BACKGROUND].blue, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  fixcolor_plot(COLOR_BASES);
  output_text(30, 80, g_input_main_message);
  strcpy(text, g_input_default_message);
  strcat(text, g_input_name);
  output_text(30, 55, text);
  output_text(15, 10, ">");
  output_text(30, 30, g_input_abort_message);
  for (i = 0; i < strlen(g_input_true_input); i++)
    text[i] = g_input_true_input[i];
  text[i] = '_';
  text[i + 1] = '\0';
  output_text(30, 10, text);
  glutSwapBuffers();
  glFlush();
  glutPostRedisplay();
  /*  sniff_for_opengl_errors(); */
}

void display_window_text(void) {
  char            string[120];
  char            number[10];
  glClearColor(main_color_table_float[COLOR_BACKGROUND].red,
	       main_color_table_float[COLOR_BACKGROUND].green,
	       main_color_table_float[COLOR_BACKGROUND].blue, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  fixcolor_plot(COLOR_BASES);
  if (g_domain_defined) {
    strcpy(string, "Base Pair: (");
    num_string_int(number, g_history[g_clicked_row]); //2003-06-20 N Markham
    strcat(string, number);
    strcat(string, ",");
    num_string_int(number, g_history[g_clicked_column]); //2003-06-20 N Markham
    strcat(string, number);
    strcat(string, ")");
  } else {
    strcpy(string, "Base:    ");
    if (g_clicked_row > 0) {
      num_string_int(number, g_history[g_clicked_row]); //2003-06-20 N Markham
      strcat(string, number);
    }
  }
  output_text(5, 82, string);
  if (!g_domain_defined) {
    strcpy(string, "Domain undefined");
  } else {
    strcpy(string, "Domain:    (");
    num_string_int(number, g_history[g_domain_row]); //2003-06-20 N Markham
    strcat(string, number);
    strcat(string, ",");
    num_string_int(number, g_history[g_domain_column]); //2003-06-20 N Markham
    strcat(string, number);
    strcat(string, ")");
  }
  output_text(5, 60, string);
  if (g_mouse_display_zoom_mode) {
    strcpy(string, "Mouse: Zooming");
  }
  else if (g_mouse_loop_stretch) {
      strcpy(string, "Mouse: Loop Stretch");
  }
  else if (g_mouse_loop_regularize) {
    strcpy(string, "Mouse: Loop Regularize");
  }
  else if (g_mouse_loop_natural) {
    strcpy(string, "Mouse: Loop Natural");
  }
  else if (g_mouse_strand_stretch) {
    strcpy(string, "Mouse: Single Strand Stretch");
  }
  else if (g_mouse_strand_linear) {
    strcpy(string, "Mouse: Single Strand Linear");
  }
  else if (g_mouse_loop_fix) {
    strcpy(string, "Mouse: Loop Fix");
  }
  else {
    strcpy(string, "Mouse: Editing");
  }
  output_text(5, 38, string);
  if (!g_counter_clockwise) {
    strcpy(string, "Draw: clockwise");
  } else { 
    strcpy(string, "Draw: counter-clockwise");
  }
  if (g_degrees_used < 1) {
    strcpy(number, " 356");
  } else {
    num_string_int(number, (int) (g_degrees_used));
    strcat(string, " ");
    strcat(string, number);
    strcat(string, " deg");
  }
  if (g_reg_angle_flag) {
    strcat(string, " Reg= ");
    num_string_int(number, g_reg_angle);
    strcat(string, number);
  }
  if (g_fix_loop_flag) {
    strcat(string, " Fix");
  }
  output_text(5, 16, string);
  glutSwapBuffers();
  sniff_for_opengl_errors();
}

void display_window_1(void) {
  float           dif_x, dif_y, alpha;
  int             imod; /* M. Zuker June 7, 2007. imod = i + 1 mod
			* sequence length (0 not used). Introduced for
			* correct plotting of circular molecules
			*/
  int             last_base; /* end_base normally and end_base + 1 when
			      * the molecule is circular */
  int             i;
  float           distance_bases;
  int             first_loop;
  float           dist, x1, y1, x2, y2, x3, y3, x4, y4;
  float           x_shift, y_shift;
  float           radius;
  float           window_font_size;
  float           up_vector_x;
  float           up_vector_y;
  float           dot_size = -9999.;
  float           oct_size = -9999.;
  int             base_1;
  double          angle;
  int             connected_to;
  int             window_outline_mode;
  int             test_3;
  int             test_5;

  /* M. Zuker, June 7, 2007. Set last_base */
  if ( g_is_circular ) {
    last_base = g_end_base + 1;
  } else {
    last_base = g_end_base;
  }
  distance_bases = DISTANCE_BASES;
  glPushMatrix();
  if ( (g_outline_mode) || ( (g_annotation == NONE) &&
       ( (g_font_size / g_display_window_zoom_factor) < 2.0) ) ) {
    window_outline_mode = TRUE;
  } else
    window_outline_mode = FALSE;
  glClearColor(main_color_table_float[COLOR_BACKGROUND].red,
	       main_color_table_float[COLOR_BACKGROUND].green,
	       main_color_table_float[COLOR_BACKGROUND].blue, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  up_vector_x = 0.;
  up_vector_y = 1.0;
  gluLookAt((WD + 72.0)/2.0 - g_display_window_x_shift,
	    (HT + 72.0)/2.0 - g_display_window_y_shift, 2.0,
	    (WD + 72.0)/2.0 - g_display_window_x_shift,
	    (HT + 72.0)/2.0 - g_display_window_y_shift,
	    0.0, up_vector_x, up_vector_y, 0.0);
  /*
   * camera x,y,z object x,y,z which way is up vector
   */
  fixcolor_plot(COLOR_BASES);
  glLineWidth(g_line_width * 1.3/g_display_window_zoom_factor *
	      g_window_1_width/(WD + 72.0) );
  glDisable(GL_LINE_SMOOTH);
  /* draw labels first */
  x_shift = .4 * g_font_size;
  y_shift = .4 * g_font_size;
  if (g_total_labels > 0) {
    fixcolor_plot(COLOR_LABEL_LINE);
    for (i = 1; i <= g_total_labels; i++) {
      if ((!g_external_domain_flag) || 
	  (((g_structure_label_value[i] >= g_start_base) &&
	    (g_structure_label_value[i] <= g_ex_start_base)) ||
	   ((g_structure_label_value[i] >= g_ex_end_base)))) {
	test_3 = FALSE;
	test_5 = FALSE;
	if (!g_is_circular) {
	  /* 2003-06-21 N. Markham adds next & prev test */
	  if (g_next[g_structure_label_value[i]] == 0) { 
	    test_3 = TRUE;
	  } else if (g_prev[g_structure_label_value[i]] == 0) { 
	    test_5 = TRUE;
	  }
	}
	if (test_5) {
	  window_draw_char(g_structure_label_x[i] - x_shift * 2,
			   g_structure_label_y[i] - 2 * y_shift, "5'");
	} else if (test_3) {
	  window_draw_char(g_structure_label_x[i] - x_shift * 2,
			   g_structure_label_y[i] - y_shift, "3'");
	} else {
	  /* g_history added 2003-06-20 Nick Markham */
	  window_draw_int(g_structure_label_x[i] - x_shift *
			  log(g_structure_label_value[i] + 
			      g_history_offset)/log(10),
			  g_structure_label_y[i] - y_shift,
			  g_history[g_structure_label_value[i]], FALSE); 
	}
	x1 = g_base_x[g_structure_label_value[i]] + (g_structure_label_x[i] -
	   g_base_x[g_structure_label_value[i]]) / 6 ;
	y1 = g_base_y[g_structure_label_value[i]] + (g_structure_label_y[i] -
	   g_base_y[g_structure_label_value[i]]) / 6;
	x2 = g_base_x[g_structure_label_value[i]] + 
	  5.*(g_structure_label_x[i] - g_base_x[g_structure_label_value[i]])/8;
	y2 = g_base_y[g_structure_label_value[i]] + 
	  5.*(g_structure_label_y[i] - g_base_y[g_structure_label_value[i]])/8;
	glBegin(GL_LINES);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glEnd();
      }
    }
  }
  /* draw loop numbers 2nd 
   * loop labels */
  if (g_loop_labels) {
    if ((g_oldconnect[g_start_base] >= 0) &&
	(g_oldconnect[g_start_base] ==
	 (g_oldconnect[g_start_base - 1] - 1)) &&
	(g_oldconnect[g_start_base + 1] ==
	 g_oldconnect[g_start_base - 1] - 1)) {
      first_loop = 2;
    } else
      first_loop = 1;
    for (i = first_loop; i <= g_total_loops; i++) {
      window_draw_int(g_loop_center_x[i] - x_shift * log(i) / log(10),
		      g_loop_center_y[i] - y_shift, i, FALSE);
      
    }
  }
  glDisable(GL_LINE_SMOOTH);
  /* make connections between bases 3rd */
  fixcolor_plot(COLOR_BASES);
  glBegin(GL_LINES);
  for (i = (g_start_base + 1); i <= last_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i <= g_ex_end_base))
	i = g_ex_end_base + 1;
    }
    if (i==g_end_base+1) {
      imod = g_start_base;
    } else {
      imod = i;
    }
    dif_x = g_base_x[imod] - g_base_x[i - 1];
    dif_y = g_base_y[imod] - g_base_y[i - 1];
    dist = sqrt(dif_x * dif_x + dif_y * dif_y);
    if ( dist>5.0*distance_bases*g_scale/8.0 || window_outline_mode) {
      alpha = (float) atan2(dif_y, dif_x);
      if (window_outline_mode) {
	x1 = g_base_x[i - 1];
	x2 = g_base_x[imod];
	y1 = g_base_y[i - 1];
	y2 = g_base_y[imod];
      } else {
	x1 = g_base_x[i - 1] + distance_bases*3.*g_scale/8.*cos(alpha);
	x2 = g_base_x[i - 1] +
	  (dist - distance_bases*3.*g_scale/8.)*cos(alpha);
	y1 = g_base_y[i - 1] + distance_bases*3.*g_scale/8.*sin(alpha);
      y2 = g_base_y[i - 1] +
	(dist - distance_bases*3.*g_scale/8.)*sin(alpha);
      }
      
      /* draw connecting line */
      if (g_next[i - 1] || g_prev[imod]) { //2003-06-23 Nick Markham
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
      }
    }
  }
  glEnd();
  /* draw base_pair lines or dots 4th */
  if (!g_lines) {
    radius = g_window_1_width /(WD + 72.0) *
      g_font_size / 2.0 * 1. / g_display_window_zoom_factor;
    /* diameter in this case */
    if (radius < 1)
      radius = 1;	/* switch back to 1 later */
    glPointSize(radius);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
  } else {
    glLineWidth(g_line_width*2.0/g_display_window_zoom_factor *
		g_window_1_width /(WD + 72.0) );
    glBegin(GL_LINES);
  }
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    connected_to = g_oldconnect[i];
    if (connected_to > i) {	/* connections between base pairs */
      set_color_base_window(g_bases[i], g_bases[connected_to]);
      if (g_lines) {
	if (window_outline_mode) {
	  x1 = g_base_x[i];
	  x2 = g_base_x[connected_to];
	  y1 = g_base_y[i];
	  y2 = g_base_y[connected_to];
	} else {
	  x1 = g_base_x[i] + (g_base_x[connected_to] - g_base_x[i]) * .30;
	  x2 = g_base_x[i] + (g_base_x[connected_to] - g_base_x[i]) * .70;
	  y1 = g_base_y[i] + (g_base_y[connected_to] - g_base_y[i]) * .30;
	  y2 = g_base_y[i] + (g_base_y[connected_to] - g_base_y[i]) * .70;
	}
	glVertex3f(x1, y1, 0.0);
	glVertex3f(x2, y2, 0.0);
      } else {/* draw dot for base pair */
	x1 = (g_base_x[i] + g_base_x[connected_to]) / 2.0;
	y1 = (g_base_y[i] + g_base_y[connected_to]) / 2.0;
	glVertex3f(x1, y1, 0.0);
      }
    }
  }
  glEnd();
  /* draw diamond for forced base pairs */
  if (g_complot_forced_bases_count > 0) {	/* make diamonds */
    glLineWidth(1.0);
    for (i = 0; i < g_complot_forced_bases_count; i += 2) {
      base_1 = g_complot_forced_bases[i];
      connected_to = g_oldconnect[base_1];
      set_color_base_window(g_bases[base_1],
			    g_bases[connected_to]);
      glBegin(GL_POLYGON);
      x1 = g_base_x[base_1] +
	(g_base_x[connected_to] - g_base_x[base_1]) * .25;
      x2 = g_base_x[base_1] +
	(g_base_x[connected_to] - g_base_x[base_1]) * .75;
      y1 = g_base_y[base_1] +
	(g_base_y[connected_to] - g_base_y[base_1]) * .25;
      y2 = g_base_y[base_1] +
	(g_base_y[connected_to] - g_base_y[base_1]) * .75;
      angle = atan2(y2 - y1, x2 - x1);
      dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2.;
      dist = sqrt(dist * dist + dist * dist);
      glVertex3f(x1, y1, 0.0);
      x3 = x2 - dist * sin(angle + PI / 4);
      y3 = y2 + dist * cos(angle + PI / 4);
      glVertex3f(x3, y3, 0.0);
      glVertex3f(x2, y2, 0.0);
      x4 = x2 - dist * sin(PI / 4 - angle);
      y4 = y2 - dist * cos(PI / 4 - angle);
      glVertex3f(x4, y4, 0.0);
      glEnd();
    }
  }
  /* draw bases 5th */
  if (!window_outline_mode) {
    glDisable(GL_LINE_SMOOTH);
    glPushMatrix();
    glLineWidth(g_line_width * 1.0 / g_display_window_zoom_factor *
		g_window_1_width /(WD + 72.0) );
    window_font_size = (WD + 72.0) * g_font_size/112320.0;
    fixcolor_plot(COLOR_BASES);
    x_shift = 0.4 * g_font_size;
    y_shift = 0.4 * g_font_size;
    if (g_annotation_dots) {
      dot_size = 1.8 * g_window_1_width / (WD + 72.0) *
	g_font_size / 2.0 * 1. / g_display_window_zoom_factor;
      oct_size = dot_size * g_display_window_zoom_factor * .95;
      glPointSize(dot_size);
    }
    for (i = g_start_base; i <= g_end_base; i++) {
      if (g_external_domain_flag) {
	if ((i > g_ex_start_base) && (i < g_ex_end_base))
	  i = g_ex_end_base;
      }
      if (g_annotation != NONE) {
	fixcolor_plot_ann(g_ann_to_color[i]);
	if (g_annotation_dots) {
	  if (dot_size > 2.0) {
	    glBegin(GL_POLYGON);
	    glVertex2f(g_base_x[i] + oct_size, g_base_y[i]);
	    glVertex2f(g_base_x[i] + .70 * oct_size, g_base_y[i] +
		       .70 * oct_size);
	    glVertex2f(g_base_x[i], g_base_y[i] + oct_size);
	    glVertex2f(g_base_x[i] - .70 * oct_size, g_base_y[i] +
		       .70 * oct_size);
	    glVertex2f(g_base_x[i] - oct_size, g_base_y[i]);
	    glVertex2f(g_base_x[i] - .70 * oct_size, g_base_y[i] -
		       .70 * oct_size);
	    glVertex2f(g_base_x[i], g_base_y[i] - oct_size);
	    glVertex2f(g_base_x[i] + .70 * oct_size, g_base_y[i] -
		       .70 * oct_size);
	    glEnd();
	  } else {
	    glBegin(GL_POINTS);
	    glVertex3f(g_base_x[i], g_base_y[i], 0.0);
	    glEnd();
	  }
	  if (g_annotation_bases) {
	    if (g_annotation == PROB) {
	      if (g_ann_to_color[i] > LOG_WHITE_BLACK_SWITCH)
		fixcolor_plot(COLOR_ANN_WHITE_LETTER);
	      else
		fixcolor_plot(COLOR_ANN_BLACK_LETTER);
	    } else {
	      if (g_ann_to_color[i] > WHITE_BLACK_SWITCH)
		fixcolor_plot(COLOR_ANN_WHITE_LETTER);
	      else
		fixcolor_plot(COLOR_ANN_BLACK_LETTER);
	    }
	  }
	}
      }
      if ((g_annotation == NONE) || ((g_annotation != NONE)
				     && (g_annotation_bases))) {
	glPushMatrix();
	glTranslatef(g_base_x[i] - x_shift, g_base_y[i] - y_shift, 0.0);
	glScalef(window_font_size, window_font_size, window_font_size);
	glutStrokeCharacter(GLUT_STROKE_ROMAN, g_bases[i]);
	glPopMatrix();
      }
    }
    glPopMatrix();
  }
  glPopMatrix();
  glutSwapBuffers();
  sniff_for_opengl_errors();
}

void display_window_nozoom(void) {
  float           dif_x, dif_y, alpha;
  int             imod; /* M. Zuker June 7, 2007. imod = i + 1 mod
			 * sequence length (0 not used). Introduced for
			 * correct plotting of circular molecules
			 */
  int             last_base; /* end_base normally and end_base + 1 when
			      * the molecule is circular */
  int             i;
  float           distance_bases;
  int             first_loop, base_1;
  float           dist, x1, y1, x2, y2, x3, y3, x4, y4, angle;
  float           x_shift, y_shift;
  float           radius;
  float           window_font_size;
  float           up_vector_x;
  float           up_vector_y;
  float           dot_size=-9999.;
  int             connected_to;
  int             window_outline_mode;
  int             test_3;
  int             test_5;
  float           box_x_center, box_y_center, box_width, box_height;

  /* M. Zuker, June 7, 2007. Set last_base */
  if ( g_is_circular ) {
    last_base = g_end_base + 1;
  } else {
    last_base = g_end_base;
  } 
  distance_bases = DISTANCE_BASES;
  glPushMatrix();
  if ((g_outline_mode) || ((g_annotation == NONE) && 
				((g_font_size) < 2.0))) {
    window_outline_mode = TRUE;
  } else
    window_outline_mode = FALSE;
  glClearColor(main_color_table_float[COLOR_BACKGROUND].red,
	       main_color_table_float[COLOR_BACKGROUND].green,
	       main_color_table_float[COLOR_BACKGROUND].blue, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  up_vector_x = 0.;
  up_vector_y = 1.0;
  gluLookAt((double) (WD + 72.0)/2.0, (double) (HT + 72.0)/2.0, 2.0, 
	    (double) (WD + 72.0)/2.0, (double) (HT + 72.0)/2.0, 0.0, 
	    up_vector_x, up_vector_y, 0.0); 
  /*
   * camera x,y,z object x,y,z which way is up vector
   */
  fixcolor_plot(COLOR_BASES);
  glLineWidth(g_line_width * 1.3 * g_window_nozoom_width / (WD + 72.));
  glDisable(GL_LINE_SMOOTH);
  /* draw box around zoom window region */ 
  box_x_center = (WD + 72.0)/2.0 - g_display_window_x_shift;
  box_y_center = (HT + 72.0)/2.0 - g_display_window_y_shift;
  box_width = (WD + 72.0) * g_display_window_zoom_factor;
  box_height = (HT + 72.0) * g_display_window_zoom_factor;
  glBegin(GL_LINE_LOOP);
  glVertex2f(box_x_center - box_width / 2., box_y_center - box_height / 2.);
  glVertex2f(box_x_center - box_width / 2., box_y_center + box_height / 2.);
  glVertex2f(box_x_center + box_width / 2., box_y_center + box_height / 2.);
  glVertex2f(box_x_center + box_width / 2., box_y_center - box_height / 2.);
  glEnd();
  /* draw labels first */
  x_shift = .4 * g_font_size;
  y_shift = .4 * g_font_size;
  if (g_total_labels > 0) {
    fixcolor_plot(COLOR_LABEL_LINE);
    for (i = 1; i <= g_total_labels; i++) {
      if ((!g_external_domain_flag) ||
	  (((g_structure_label_value[i] >= g_start_base) &&
	    (g_structure_label_value[i] <= g_ex_start_base)) ||
	   ((g_structure_label_value[i] >= g_ex_end_base)))) {
	test_3 = FALSE;
	test_5 = FALSE;
	if (!g_is_circular) {
	  if (i == g_total_labels) {
	    test_3 = TRUE;
	  } else {
	    if (i == 1)
	      test_5 = TRUE;
	  }
	}
	if (test_5) {
	  window_draw_char(g_structure_label_x[i] - x_shift * 2,
			   g_structure_label_y[i] - 2 * y_shift, "5'");
	} else {
	  if (test_3) {
	    window_draw_char(g_structure_label_x[i] - x_shift * 2,
			     g_structure_label_y[i] - y_shift, "3'");
	  } else {
	    window_draw_int(g_structure_label_x[i] - x_shift *
			    log(g_structure_label_value[i] + g_history_offset)
			    /log(10), g_structure_label_y[i] - y_shift,
			    g_history[g_structure_label_value[i]], FALSE); 
	    /* 2003-06-20 Nick Markham */
	  }
	}
	x1 = g_base_x[g_structure_label_value[i]] + 
	  (g_structure_label_x[i] - g_base_x[g_structure_label_value[i]])/6.;
	y1 = g_base_y[g_structure_label_value[i]] +
	  (g_structure_label_y[i] - g_base_y[g_structure_label_value[i]])/6.;
	x2 = g_base_x[g_structure_label_value[i]] +
	  5.*(g_structure_label_x[i]-g_base_x[g_structure_label_value[i]])/8.;
	y2 = g_base_y[g_structure_label_value[i]] + 
	  5.*(g_structure_label_y[i]-g_base_y[g_structure_label_value[i]])/8.;
	glBegin(GL_LINES);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glEnd();
      }
    }
  }
  /* draw loop numbers 2nd */
  /* loop labels */
  if (g_loop_labels) {
    if ((g_oldconnect[g_start_base] >= 0) &&
	(g_oldconnect[g_start_base] ==
	 (g_oldconnect[g_start_base - 1] - 1)) &&
	(g_oldconnect[g_start_base + 1] ==
	 g_oldconnect[g_start_base - 1] - 1)) {
      first_loop = 2;
    } else
      first_loop = 1;
    for (i = first_loop; i <= g_total_loops; i++) {
      window_draw_int(g_loop_center_x[i] - x_shift * log(i) / log(10),
		      g_loop_center_y[i] - y_shift, i, FALSE);
      
    }
  }
  /* glEnable(GL_LINE_SMOOTH); */
  glLineWidth(g_line_width * 0.5 * g_window_nozoom_width/(WD + 72.0));
  /* make connections between bases 3rd */
  fixcolor_plot(COLOR_BASES);
  glBegin(GL_LINES);
  for (i = (g_start_base + 1); i <= last_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i <= g_ex_end_base))
	i = g_ex_end_base + 1;
    }
    if (i==g_end_base+1) {
      imod = g_start_base;
    } else {
      imod = i;
    }
    dif_x = g_base_x[imod] - g_base_x[i - 1];
    dif_y = g_base_y[imod] - g_base_y[i - 1];
    dist = sqrt(dif_x * dif_x + dif_y * dif_y);
    if ((dist>5.0*distance_bases*g_scale/8.) || window_outline_mode) {
	  alpha = (float) atan2(dif_y, dif_x);
      if (window_outline_mode) {
	x1 = g_base_x[i - 1];
	x2 = g_base_x[imod];
	y1 = g_base_y[i - 1];
	y2 = g_base_y[imod]; 
      } else {
	x1 = g_base_x[i - 1] + distance_bases*3.*g_scale/8.*cos(alpha);
	x2 = g_base_x[i - 1] + (dist-distance_bases*3.*g_scale/8.)*cos(alpha);
	y1 = g_base_y[i - 1] + distance_bases*3.*g_scale/8.*sin(alpha);
	y2 = g_base_y[i - 1] + (dist-distance_bases*3.*g_scale/8.)*sin(alpha);
      }
      /* draw connecting line */
      if (g_next[i - 1] || g_prev[imod]) { //2003-06-23 Nick Markham
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
      }
    }
  }
  glEnd();
  /* draw base_pair lines or dots 4th */
  if (!g_lines) {
    radius = g_window_nozoom_width / (WD + 72.0) * g_font_size / 2.0;
    /* diameter in this case */
    if (radius < 1)
      radius = 1;	/* switch back to 1 later */
    glPointSize(radius);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
  } else {
    glLineWidth(g_line_width * 2.0 * g_window_nozoom_width / (WD + 72.0) );
    glBegin(GL_LINES);
  }
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    connected_to = g_oldconnect[i];
    if (connected_to > i) {
      /* connections between base pairs */
      set_color_base_window(g_bases[i], g_bases[connected_to]);
      if (g_lines) {
	if (window_outline_mode) {
	  x1 = g_base_x[i];
	  x2 = g_base_x[connected_to];
	  y1 = g_base_y[i];
	  y2 = g_base_y[connected_to];
	} else {
	  x1 = g_base_x[i] + (g_base_x[connected_to] - g_base_x[i]) * .30;
	  x2 = g_base_x[i] + (g_base_x[connected_to] - g_base_x[i]) * .70;
	  y1 = g_base_y[i] + (g_base_y[connected_to] - g_base_y[i]) * .30;
	  y2 = g_base_y[i] + (g_base_y[connected_to] - g_base_y[i]) * .70;
	}
	glVertex3f(x1, y1, 0.0);
	glVertex3f(x2, y2, 0.0);
      } else {/* draw dot for base pair */
	x1 = (g_base_x[i] + g_base_x[connected_to]) / 2.0;
	y1 = (g_base_y[i] + g_base_y[connected_to]) / 2.0;
	glVertex3f(x1, y1, 0.0);
      }
    }
  }
  glEnd();
  /* draw diamond for forced base pairs */
  if (g_complot_forced_bases_count > 0) {	/* make diamonds */
    glLineWidth(1.0);
    for (i = 0; i < g_complot_forced_bases_count; i += 2) {
      base_1 = g_complot_forced_bases[i];
      connected_to = g_oldconnect[base_1];
      set_color_base_window(g_bases[base_1],
			    g_bases[connected_to]);
      glBegin(GL_POLYGON);
      x1 = g_base_x[base_1] +
	(g_base_x[connected_to] - g_base_x[base_1]) * .25;
      x2 = g_base_x[base_1] +
	(g_base_x[connected_to] - g_base_x[base_1]) * .75;
      y1 = g_base_y[base_1] +
	(g_base_y[connected_to] - g_base_y[base_1]) * .25;
      y2 = g_base_y[base_1] +
	(g_base_y[connected_to] - g_base_y[base_1]) * .75;
      angle = atan2(y2 - y1, x2 - x1);
      dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2.;
      dist = sqrt(dist * dist + dist * dist);
      glVertex3f(x1, y1, 0.0);
      x3 = x2 - dist * sin(angle + PI / 4);
      y3 = y2 + dist * cos(angle + PI / 4);
      glVertex3f(x3, y3, 0.0);
      glVertex3f(x2, y2, 0.0);
      x4 = x2 - dist * sin(PI / 4 - angle);
      y4 = y2 - dist * cos(PI / 4 - angle);
      glVertex3f(x4, y4, 0.0);
      glEnd();
    }
  }
  /* draw bases 5th */
  if (!window_outline_mode) {
    glDisable(GL_LINE_SMOOTH);
    glPushMatrix();
    glLineWidth(g_line_width * g_window_nozoom_width / (WD + 72.0));
    window_font_size = (WD + 72.0) * g_font_size/112320.0;
    fixcolor_plot(COLOR_BASES);
    x_shift = .4 * g_font_size;
    y_shift = .4 * g_font_size;
    if (g_annotation_dots) {
      dot_size = x_shift * 2.0 * g_window_nozoom_width / (WD + 72.0);
      if (dot_size < 1.0)
	dot_size = 1.0;
      glPointSize(dot_size * 1.0);
    }
    for (i = g_start_base; i <= g_end_base; i++) {
      if (g_external_domain_flag) {
	if ((i > g_ex_start_base) && (i < g_ex_end_base))
	  i = g_ex_end_base;
      }
      if (g_annotation != NONE) {
	fixcolor_plot_ann(g_ann_to_color[i]);
	if (g_annotation_dots) {
	  if (dot_size > 2.25) {
	    glBegin(GL_POLYGON);
	    glVertex2f(g_base_x[i] + dot_size, g_base_y[i]);
	    glVertex2f(g_base_x[i] + .70 * dot_size, g_base_y[i] +
		       .70 * dot_size);
	    glVertex2f(g_base_x[i], g_base_y[i] + dot_size);
	    glVertex2f(g_base_x[i] - .70 * dot_size, g_base_y[i] +
		       .70 * dot_size);
	    glVertex2f(g_base_x[i] - dot_size, g_base_y[i]);
	    glVertex2f(g_base_x[i] - .70 * dot_size, g_base_y[i] -
		       .70 * dot_size);
	    glVertex2f(g_base_x[i], g_base_y[i] - dot_size);
	    glVertex2f(g_base_x[i] + .70 * dot_size, g_base_y[i] -
		       .70 * dot_size);
	    glEnd();
	    if (g_annotation_bases) {
	      if (g_annotation == PROB) {
		if (g_ann_to_color[i] > LOG_WHITE_BLACK_SWITCH)
		  fixcolor_plot(COLOR_ANN_WHITE_LETTER);
		else
		  fixcolor_plot(COLOR_ANN_BLACK_LETTER);
	      } else {
		if (g_ann_to_color[i] > WHITE_BLACK_SWITCH)
		  fixcolor_plot(COLOR_ANN_WHITE_LETTER);
		else
		  fixcolor_plot(COLOR_ANN_BLACK_LETTER);
	      }
	    }
	  } else {
	    glBegin(GL_POINTS);
	    glVertex2f(g_base_x[i], g_base_y[i]);
	    glEnd();
	  }
	}
      }
      if ((g_annotation==NONE)||((g_annotation!=NONE)&&(g_annotation_bases))) {
	glPushMatrix();
	glTranslatef(g_base_x[i] - x_shift, g_base_y[i] - y_shift, 0.0);
	glScalef(window_font_size, window_font_size, window_font_size);
	glutStrokeCharacter(GLUT_STROKE_ROMAN, g_bases[i]);
	glPopMatrix();
      }
    }
    glPopMatrix();
  }
  glPopMatrix();
  glutSwapBuffers();
  sniff_for_opengl_errors();
}

void make_structure_img1(char *filename, int png_mode, int jpg_mode) {
  float           x_shift, y_shift, zoom_factor;
  if (g_img_is_main) {
    x_shift = 0;
    y_shift = 0;
    zoom_factor = 1.0;
  } else {
    x_shift = g_display_window_x_shift;
    y_shift = g_display_window_y_shift;
    zoom_factor = g_display_window_zoom_factor;
  }
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  make_structure_img(filename, g_img_width, g_img_height, g_img_interlaced,
		     g_first_line, g_start_base, g_end_base, g_loop_center_x,
		     g_loop_center_y, g_total_labels, g_structure_label_x, 
		     g_structure_label_y, g_base_x, g_base_y, g_scale,
		     g_outline_mode, g_oldconnect, g_font_size, 
		     g_structure_label_value, g_lines, g_bases,
		     g_loop_labels, g_total_loops, g_external_domain_flag, 
		     g_ex_start_base, g_ex_end_base, g_history_offset, 
		     g_annotation, g_ann_to_color, g_annotation_dots, 
		     g_annotation_bases, g_length, x_shift, y_shift, 
		     zoom_factor, g_is_circular, g_complot_forced_bases,
		     g_complot_forced_bases_count, png_mode, jpg_mode, 
		     g_command_line_x);
#endif
}

void reshape_window_nozoom(int w, int h) {
  if (w > (int) WD)
    w = (int) WD;
  if (h == g_window_nozoom_height)
    h = (int) (((float) w) * (HT + 72.0)/(WD + 72.0) + 0.5);
  else
    w = (int) (((float) h) * (WD + 72.0)/(HT + 72.0) + 0.5);
  g_window_nozoom_width = w;
  g_window_nozoom_height = h;
  glutReshapeWindow(w, h);
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	/* adjust last values below for z distance shown */
  /*  glFrustum(-(WD + 72.0)/2.0, (WD + 72.0)/2.0, -(HT + 72.0)/2.0, 
      (HT + 72.0)/2.0, 50, 200); */
  glFrustum(-(WD + 72.0)/2.0, (WD + 72.0)/2.0, -(HT + 72.0)/2.0, 
	    (HT + 72.0)/2.0, 1.9999, 2.0001);
  glMatrixMode(GL_MODELVIEW);
  glutPostRedisplay();
}

void reshape_window(int w, int h) {
  if (h == g_window_1_height) /* for new width, set height */ 
    h = (int) (((float) w) * (HT + 72.0)/(WD + 72.0) + 0.5); 
  else /* for new height, set width */
    w = (int) (((float) h) * (WD + 72.0)/(HT + 72.0) + 0.5); 
  g_window_1_width = w;
  g_window_1_height = h;
  glutReshapeWindow(w, h);
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	/* adjust last values below for z distance
			 * shown */
  glFrustum(-(WD/2.0 + 36.) * g_display_window_zoom_factor,
	    (WD/2.0 + 36.)* g_display_window_zoom_factor,
	    -(HT/2.0 + 36.) * g_display_window_zoom_factor,
	    (HT/2.0 + 36.) * g_display_window_zoom_factor, 1.9999, 2.0001);
  glMatrixMode(GL_MODELVIEW);
  glutPostRedisplay();
}

void window_create_img_res(int p) {
  g_png_mode = FALSE;
  g_jpg_mode = FALSE;
  if (p > 10) {
    p -= 10;
    g_img_is_main = FALSE;
  } else
    g_img_is_main = TRUE;
  if (p == 1) {
    g_img_width = 300;
  } else if (p == 2) {
    g_img_width = 510;
  } else if (p == 3) {
    g_img_width = 612;
  } else if (p == 4) {
    g_img_width = 765;
  } else if (p == 5) {
    g_img_width = 900;
  } else if (p == 6) {
    g_img_width = 1210;
  }
  g_img_height = (int) ( ((float) (g_img_width))*17.0/13.0);
}

void window_create_gif_res(int p) {
  window_create_img_res(p);
  g_input_type = GIF;
  strcpy(g_input_name, g_structure_filename);
  strcat(g_input_name, ".gif");
  g_input_true_input[0] = '\0';
  inputSetWindow();
  glutShowWindow();
}

void window_create_jpg_res(int p) {
  window_create_img_res(p);
  g_input_type = JPG;
  strcpy(g_input_name, g_structure_filename);
  strcat(g_input_name, ".jpg");
  g_input_true_input[0] = '\0';
  inputSetWindow();
  glutShowWindow();
}

void window_create_png_res(int p) {
  window_create_img_res(p);
  g_input_type = PNG;
  strcpy(g_input_name, g_structure_filename);
  strcat(g_input_name, ".png");
  g_input_true_input[0] = '\0';
  inputSetWindow();
  glutShowWindow();
}

void make_structure_img_interactive_finish(int type) {
  char            filename[100];
  if (type == PNG)
    g_png_mode = TRUE;
  if (type == JPG)
    g_jpg_mode = TRUE;
  if (strlen(g_input_true_input) == 1) {
    if (g_input_true_input[0] == '/') {
      glutDestroyWindow(g_window_input);
      g_window_input_exists = FALSE;
      return;
    }
  }
  if (strlen(g_input_true_input) == 0)
    strcpy(filename, g_structure_filename);
  else
    strcpy(filename, g_input_true_input);
  if (g_img_is_main) {
    scale_to_fit_window();
    glutSetWindow(g_window_1);
    g_display_window_zoom_factor = 1.0;
    reshape_window(g_window_1_width, g_window_1_height);
  }
  make_structure_img1(filename, g_png_mode, g_jpg_mode);
  glutDestroyWindow(g_window_input);
  g_window_input_exists=FALSE;
}

void window_set_label(int step) {
  if (step == 0)
    g_label_frequency = g_length + 10;
  else
    g_label_frequency = step;
  set_label_coordinates();
  scale_to_fit_window(); 
  show_both();
}

void create_ss_output_interactive_start(void) {
  g_input_type = SS;
  inputSetWindow();
  strcpy(g_input_name, g_out_ss_filename);
  strcat(g_input_name, ".new.ss");
  g_input_true_input[0] = '\0';
  glutShowWindow();
}

void create_output_ss_interactive_finish(void) {
  char            filename[100];
  if (strlen(g_input_true_input) == 1) {
    if (g_input_true_input[0] == '/') {
      glutDestroyWindow(g_window_input);
      g_window_input_exists = FALSE;
      return;
    }
  }
  if (strlen(g_input_true_input) == 0) {
    strcpy(filename, g_out_ss_filename);
    strcat(filename, ".new");
  } else
    strcpy(filename, g_input_true_input);
  chop_suffix(filename, ".ss");
  /* May 23, 2007. M. Zuker replaces g_connected with g_oldconnect */
  create_ss_output_real(filename, g_length, g_oldconnect, g_bases, g_base_x, 
			g_base_y, g_ss_code);
  glutDestroyWindow(g_window_input);
  g_window_input_exists = FALSE;
}

void mouse_function_submenu_choice(int choice) {
  g_mouse_display_zoom_mode = FALSE;
  g_mouse_rotate_helix = FALSE;
  g_mouse_translate = FALSE;
  g_mouse_loop_stretch_pos_set = FALSE;
  g_mouse_loop_stretch = FALSE;
  g_mouse_loop_natural = FALSE;
  g_mouse_loop_regularize = FALSE;
  g_mouse_strand_stretch = FALSE;
  g_mouse_strand_linear = FALSE;
  g_mouse_loop_fix = FALSE;
  if (choice == 1)
    g_mouse_display_zoom_mode = TRUE;
  else if (choice == 4) {
    g_mouse_translate = TRUE;
    g_mouse_rotate_helix = TRUE;
    g_mouse_rotate_helix_set = FALSE;
    g_mouse_translate_pos_set = FALSE;
  }
  else if (choice == 9) {
    g_mouse_loop_stretch = TRUE;
  }
  else if (choice == 15) {
    g_mouse_loop_natural = TRUE;
  }
  else if (choice == 16) {
    g_mouse_loop_regularize = TRUE;
  }
  else if (choice == 20) {
    g_mouse_strand_stretch = TRUE;
    g_strand_set = FALSE;
  }
  else if (choice == 25) {
    g_mouse_strand_linear = TRUE;
    g_strand_set = FALSE;
  }
  else if (choice == 35) {
    g_mouse_loop_fix = TRUE;
  }
  show_text();
}

void execute_undo(void) {
  int             i;
  float           temp;
  /* copy base_x, base_y structure_angle to two previous states */
  for (i = g_start_base; i <= g_end_base; i++) {
    if (g_external_domain_flag) {
      if ((i > g_ex_start_base) && (i < g_ex_end_base))
	i = g_ex_end_base;
    }
    /* for base_x */
    temp = g_base_x[i];
    g_base_x[i] = g_undo1_base_x[i];
    g_undo1_base_x[i] = g_undo2_base_x[i];
    g_undo2_base_x[i] = temp;
    temp = g_base_y[i];
    g_base_y[i] = g_undo1_base_y[i];
    g_undo1_base_y[i] = g_undo2_base_y[i];
    g_undo2_base_y[i] = temp;
    temp = g_structure_angle[i];
    g_structure_angle[i] = g_undo1_structure_angle[i];
    g_undo1_structure_angle[i] = g_undo2_structure_angle[i];
    g_undo2_structure_angle[i] = temp;
  }
  for (i = 1; i <= g_total_loops; i++) {
    /* for loop at _x */
    temp = g_loop_center_x[i];
    g_loop_center_x[i] = g_undo1_loop_center_x[i];
    g_undo1_loop_center_x[i] = g_undo2_loop_center_x[i];
    g_undo2_loop_center_x[i] = temp;
    /* for loop at y */
    temp = g_loop_center_y[i];
    g_loop_center_y[i] = g_undo1_loop_center_y[i];
    g_undo1_loop_center_y[i] = g_undo2_loop_center_y[i];
    g_undo2_loop_center_y[i] = temp;
  }
  temp = g_scale;
  g_scale = g_undo1_scale;
  g_undo1_scale = g_undo2_scale;
  g_undo2_scale = temp;
  temp = g_ss_distance_basepair;
  g_ss_distance_basepair = g_undo1_ss_distance_basepair;
  g_undo1_ss_distance_basepair = g_undo2_ss_distance_basepair;
  g_undo2_ss_distance_basepair = temp;
  temp = g_ss_distance_base_adj;
  g_ss_distance_base_adj = g_undo1_ss_distance_base_adj;
  g_undo1_ss_distance_base_adj = g_undo2_ss_distance_base_adj;
  g_undo2_ss_distance_base_adj = temp;
  temp = g_font_size;
  g_font_size = g_undo1_font_size;
  g_undo1_font_size = g_undo2_font_size;
  g_undo2_font_size = temp;
  temp = g_line_width;
  g_line_width = g_undo1_line_width;
  g_undo1_line_width = g_undo2_line_width;
  g_undo2_line_width = temp;
  set_label_coordinates();
}

/* M. Zuker Dec 20,2006
 * Add option 96 and improve syntax
 * Note that whole/partial is ignored in natural
 * angle mode (g_natural is TRUE)
 */
void window_run_sir_graph(int option) {
  if (option == 99) {	
    /* run sir_graph on domain, whole */
    run_sir_graph_on_domain(FALSE, FALSE);
  } else if (option == 97) {	
    /* run sir_graph on domain, partial */
    run_sir_graph_on_domain(FALSE, TRUE);
  } else if (option == 98) {
    /* run sir_graph on excluded domain */
    run_sir_graph_on_domain(TRUE, FALSE);
  } else {
    /* run sir_graph on excluded domain */
    run_sir_graph_on_domain(TRUE, TRUE);
  }
}

void            regularize_all_angles(void); 

void No_Action_Menu(int option) {
  if (option < 0)
    return;
  if (option == 1)
    try_exit(0);
  if ((option == 3) || (option == 35)) {
    if (option == 35)
      g_post_is_zoom = TRUE;
    else {
      g_post_is_zoom = FALSE;
      scale_to_fit_window(); 
      adjust_structure_coordinates(); 
    }
    glutSetWindow(g_window_1);
    reshape_window(g_window_1_width, g_window_1_height);
    make_structure_ps_interactive_start();
    return;
  }
  if (option == 68) {
    if (g_fix_loop_flag == TRUE)
      g_fix_loop_flag = FALSE;
    else
      g_fix_loop_flag = TRUE;
    show_text();
  }
  if (option == 42) {
    if (g_counter_clockwise == TRUE) {
      g_counter_clockwise = FALSE;
      printf("Drawing clockwise.\n");
    } else {
      g_counter_clockwise = TRUE;
      printf("Drawing counterclockwise.\n");
    }
    show_text();
    store_for_undo();
    /* The excluded code immediately redraws the entire structure in
     * the opposite mode.
     g_ss_backwards = FALSE;
     g_small_angle = FALSE;
     g_ss_mode = FALSE;
     g_scale = 1.0;
     compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
     adjust_structure_coordinates();
     show_both();
    */
    return;
  }
  if (option == 4) {
    if (g_outline_mode == TRUE)
      g_outline_mode = FALSE;
    else
      g_outline_mode = TRUE;
    show_both();
    return;
  }
  if (option == 5) {
    if (g_lines == TRUE)
      g_lines = FALSE;
    else
      g_lines = TRUE;
    show_both();
    return;
  }
  if (option == 7) {
    if (g_loop_labels == TRUE) {
      g_loop_labels = FALSE;
    } else
      g_loop_labels = TRUE;
    show_both();
    return;
  }
  if (option == 6) { /* redraw entire structure with natural angles */
    store_for_undo();
    g_ss_backwards = FALSE;
    g_natural = TRUE;
    g_small_angle = FALSE;	/* turn off z value stuff */
    g_ss_mode = FALSE;
    g_scale = 1.0;
    compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
    adjust_structure_coordinates();
    show_both();
    return;
  }
  if (option == 72) { /* redraw entire structure with all angles regularized */
    if (!g_reg_angle_flag) {
      printf("The Regularize angle must be set first.\n");
      return;
    }
    store_for_undo();
    regularize_all_angles();
    show_both();
    return;
  }
  if (option == 61) {	/* redraw entire structure with circle graph
			 * angles */
    store_for_undo();
    g_ss_backwards = FALSE;
    g_small_angle = FALSE;
    g_natural = FALSE;
    g_ss_mode = FALSE;
    g_scale = 1.0;
    compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
    adjust_structure_coordinates();
    show_both();
    return;
  }
  if (option == 19) {	/* make the structure fit in the window */
    scale_to_fit_window();
    g_display_window_zoom_factor = 1.0;
    glutSetWindow(g_window_1);
    reshape_window(g_window_1_width, g_window_1_height);
    show_nozoom();
    return;
  }
  if (option == 25) {
    scale_to_fit_window();
    show_both();
    create_ss_output_interactive_start();
    return;
  }
  if (option == 26) {	/* undo option */
    printf("Undo\n");
    execute_undo();
    show_both();
    return;
  }
}

/* Rewritten by M. Zuker, Dec 19, 2006 */
void window_set_z(int z_val) {
  store_for_undo();
  g_ss_mode = FALSE;
  g_ss_backwards = FALSE;
 if (z_val == 2) {
    g_helix_base_advance = 2;
  } else if (z_val == 3) {
    g_helix_base_advance = .5;
  }  else if (z_val == 4) {
    g_helix_base_advance = .1;
  }  else if (z_val == 5) {
    g_helix_base_advance = .05;
  }  else if (z_val == 6) {
    g_helix_base_advance = .01;
  } else if (z_val == 7) {
    g_helix_base_advance = .005;
  }  else if (z_val == 8) {
    g_helix_base_advance = .0001;
  } else if (z_val == 0) {
    g_small_angle = FALSE;
    g_small_angle = TRUE;
    g_natural = FALSE;
    printf("z value turned off\n");
  } else {
    printf("Error! z_val = %d, which is invalid.\n",z_val);
  }
 if (z_val != 0) 
   printf("z value set to %f\n", g_helix_base_advance);
 g_scale = 1.0;
 compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
 adjust_structure_coordinates();
 show_both();
}

void draw_flat(int choice) {
  if (choice==80) {
    if (g_flat==TRUE || g_flat_alternate==TRUE) {
      g_flat = FALSE;
      g_flat_alternate = FALSE;
      printf("Flat mode turned off.\n");
      store_for_undo();
      g_ss_backwards = FALSE;
      g_small_angle = FALSE;
      g_ss_mode = FALSE;
      g_scale = 1.0;
      compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
      adjust_structure_coordinates();
      show_both();
    }
  } else if (choice==81) {
    if (g_flat==FALSE) {
      g_flat = TRUE;
      g_flat_alternate = FALSE;
      printf("Flat mode turned on.\n"); 
      store_for_undo();
      g_ss_backwards = FALSE;
      g_small_angle = FALSE;
      g_ss_mode = FALSE;
      g_scale = 1.0;
      compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
      adjust_structure_coordinates();
      show_both();
    }
    g_flat_alternate = FALSE;
  } else if (choice==82) {
    if (g_flat_alternate==FALSE) {
      printf("Flat alternate mode turned on.\n");
      g_flat = TRUE;
      g_flat_alternate = TRUE;
      store_for_undo();
      g_ss_backwards = FALSE;
      g_small_angle = FALSE;
      g_ss_mode = FALSE;
      g_scale = 1.0;
      compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
      adjust_structure_coordinates();
      show_both();
    }
  }
} 

void window_set_deg(int deg) {
  store_for_undo();
  g_natural = FALSE;
  g_ss_mode = FALSE;
  g_degrees_used = (float) deg;
  g_scale = 1.0;
  compute_angles_ps(g_out_filename,
		    g_degrees_used, g_is_circular);
  adjust_structure_coordinates();
  show_both();
  glutSetWindow(g_window_text);
  glutPostRedisplay();
}
void window_set_reg(int deg) {
  if (deg == 0) {
    g_reg_angle_flag = FALSE;
    g_reg_angle = 0.;
    printf("Regularize_angle_turned_off\n");
  } else {
    g_reg_angle = (float) deg;
    g_reg_angle_flag = TRUE;
    printf("Regularize Angles to %d degrees\n", deg);
  }
  glutSetWindow(g_window_text);
  glutPostRedisplay();
}

void window_set_rot(int deg) {
  int             i;
  store_for_undo();
  rotate_structure(deg * PI / 180., (WD + 72.0)/2.0, (HT + 72.0)/2.0,
		   TRUE, 1, g_total_loops, g_start_base, g_end_base);
  printf("rotating %.2f degrees\n", (float) deg);
  for (i = g_start_base; i <= g_end_base; i++)
    g_structure_angle[i] = g_structure_angle[i] + deg * PI / 180.;
  scale_to_fit_window();
  /* g_display_window_zoom_factor = 1.0; */
  show_both();
}

void c_or_cc(int choice) {
  return;
}

void window_redraw_ex_domain(int choice) {     /* draw the excluded domain */
  float           original_angle, new_angle;
  int             i;
  float           circle_graph_perp;
  int             domain_row, domain_column;
  float           last_position_x, last_position_y;
  float           new_position_x, new_position_y;
  if (!g_domain_defined) {
    printf("No domain defined\n");
    return;
  }
  store_for_undo();
  if (choice == 1) {
    printf("Selected natural_angles for excluded domain (%d,%d)\n",
	   g_domain_row, g_domain_column);
    original_angle = g_structure_angle[g_domain_row];
    last_position_x =
      (g_base_x[g_domain_row] + g_base_x[g_domain_column]) / 2.;
    last_position_y =
      (g_base_y[g_domain_row] + g_base_y[g_domain_column]) / 2.;
    
    make_natural_structure(g_domain_row, g_domain_column, FALSE, TRUE);
    new_position_x =
      (g_base_x[g_domain_row] + g_base_x[g_domain_column]) / 2.;
    new_position_y =
      (g_base_y[g_domain_row] + g_base_y[g_domain_column]) / 2.;
    new_angle = g_structure_angle[g_domain_row];
    shift_structure(-1.0 * (last_position_x - new_position_x),
		    -1.0 * (last_position_y - new_position_y),
		    g_previous_loop[g_domain_row] + 1,
		    g_previous_loop[g_domain_column], g_domain_row + 1,
		    g_domain_column - 1);
    rotate_structure(new_angle-original_angle, new_position_x, new_position_y,
		     FALSE,
		     g_previous_loop[g_domain_row] + 1,
		     g_previous_loop[g_domain_column],
		     g_domain_row + 1, g_domain_column - 1);
    for (i = (g_domain_row + 1); (i < g_domain_column); i++)
      g_structure_angle[i] = g_structure_angle[i] +
	(new_angle - original_angle);
    set_label_coordinates();
    rotate_structure(original_angle - new_angle, (WD + 72.0)/2.0, 
		     (HT + 72.0)/2.0, TRUE, 1, g_total_loops, 
		     g_start_base, g_end_base);
    for (i = g_start_base; i <= g_end_base; i++)
      g_structure_angle[i] = g_structure_angle[i] +
	(original_angle - new_angle);
    scale_to_fit_window();
    show_both();
    return;
  }
  if (choice > 1) {	/* normal circle graph structure on this
			 * domain */
    i = g_domain_row;
    domain_row = i;
    domain_column = g_connected[domain_row];
    last_position_x =
      (g_base_x[g_domain_row] + g_base_x[g_domain_column]) / 2.;
    last_position_y =
      (g_base_y[g_domain_row] + g_base_y[g_domain_column]) / 2.;
    original_angle = g_structure_angle[g_domain_row];
    if (choice == 2) {
      printf("Selected Circle Graph for excluded Domain (%d,%d)\n",
	     domain_row, domain_column);
      main_data_ps_structure(g_degrees_used,g_is_circular,g_oldconnect, FALSE,
			     domain_row, domain_column, &circle_graph_perp, 
			     FALSE, TRUE, 0.0, 0);
    } else {
      printf("Selected Full Circle Graph for excluded Domain (%d,%d)\n",
	     domain_row, domain_column);
      main_data_ps_structure(g_degrees_used,g_is_circular, g_oldconnect, FALSE,
			     domain_row, domain_column, &circle_graph_perp, 
			     TRUE, TRUE, g_structure_angle[domain_row] 
			     - PI / 2., choice);
    }
    new_position_x =
      (g_base_x[g_domain_row] + g_base_x[g_domain_column]) / 2.;
    new_position_y =
      (g_base_y[g_domain_row] + g_base_y[g_domain_column]) / 2.;
    new_angle = g_structure_angle[g_domain_row];
    shift_structure(-1.0 * (last_position_x - new_position_x),
		    -1.0 * (last_position_y - new_position_y),
		    g_previous_loop[g_domain_row] + 1,
		    g_previous_loop[g_domain_column], g_domain_row + 1,
		    g_domain_column - 1);
    rotate_structure(new_angle-original_angle, new_position_x, new_position_y,
		     FALSE,
		     g_previous_loop[g_domain_row] + 1,
		     g_previous_loop[g_domain_column],
		     g_domain_row + 1, g_domain_column - 1);
    for (i = (g_domain_row + 1); (i < g_domain_column); i++)
      g_structure_angle[i] = g_structure_angle[i] +
	(new_angle - original_angle);
    set_label_coordinates();
    rotate_structure(original_angle - new_angle, (WD + 72.0)/2.0, 
		     (HT + 72.0)/2.0, TRUE, 1, g_total_loops, 
		     g_start_base, g_end_base);
    for (i = g_start_base; (i <= g_end_base); i++)
      g_structure_angle[i] = g_structure_angle[i] +
	(original_angle - new_angle);
    scale_to_fit_window();
    show_both();
    return;
  }
}

void window_redraw_domain(int choice) {
  float           original_angle, new_angle;
  int             i;
  float           circle_graph_perp;
  int             last_loop_drawn;
  int             domain_row, domain_column;
  if (!g_domain_defined) {
    printf("No domain defined\n");
    return;
  }
  store_for_undo();
  if (choice == 1) {
    printf("Selected natural_angles for domain (%d,%d)\n",
	   g_domain_row, g_domain_column);
    make_natural_structure(g_domain_row, g_domain_column, TRUE, FALSE);
    set_label_coordinates();
    show_both();
    return;
  }
  if (choice == 2) {	/* normal circle graph structure on this
			 * domain */
    i = g_domain_row;
    while (g_connected[i - 1] == (g_connected[i] + 1))
      i--;
    domain_row = i;
    if (i < g_start_base)
      i = g_start_base;
    domain_column = g_connected[domain_row];
    original_angle = g_structure_angle[domain_row];
    printf("Selected Normal Circle Graph for Domain (%d,%d)\n",
	   domain_row, domain_column);
    last_loop_drawn = 
      main_data_ps_structure(g_degrees_used, g_is_circular,
			     g_oldconnect, TRUE, domain_row, domain_column,
			     &circle_graph_perp, FALSE, FALSE, 0.0, 0);
    if (g_counter_clockwise)
      circle_graph_perp -= PI / 2;
    else
      circle_graph_perp += PI / 2;
    new_angle = original_angle - circle_graph_perp;
    rotate_structure(new_angle, (g_undo1_base_x[domain_row] +
		      g_undo1_base_x[g_connected[domain_row]]) / 2.,
		     (g_undo1_base_y[domain_row] +
		      g_undo1_base_y[g_connected[domain_row]]) / 2.,
		     FALSE, g_previous_loop[domain_row] + 1,
		     last_loop_drawn - 1, domain_row, domain_column);
    g_structure_angle[g_connected[i]] = circle_graph_perp + PI;
    for (i = (domain_row); (i <= domain_column); i++)
      g_structure_angle[i] = g_structure_angle[i] + (new_angle);
    set_label_coordinates();
    show_both();
    return;
  }
  /* normal circle graph except use full circle on this domain */
  i = g_domain_row;
  while (g_connected[i - 1] == (g_connected[i] + 1))
    i--;
  if (i < g_start_base)
    i = g_start_base;
  domain_row = i;
  domain_column = g_connected[domain_row];
  original_angle = g_structure_angle[domain_row];
  printf("Full Circle Graph for Domain (%d,%d)\n",
	 domain_row, domain_column);
  last_loop_drawn = 
    main_data_ps_structure(g_degrees_used, g_is_circular,
			   g_oldconnect, TRUE, domain_row, domain_column, 
			   &circle_graph_perp, TRUE, FALSE, 0.0, choice);
  if (g_counter_clockwise)
    circle_graph_perp -= PI / 2;
  else
    circle_graph_perp += PI / 2;
  new_angle = original_angle - circle_graph_perp;
  rotate_structure(new_angle,
		   (g_base_x[domain_row] +
		    g_base_x[g_connected[domain_row]]) / 2.,
		   (g_base_y[domain_row] +
		    g_base_y[g_connected[domain_row]]) / 2.,
		   FALSE,
		   g_previous_loop[domain_row] + 1,
		   last_loop_drawn - 1,
		   domain_row, domain_column);
  g_structure_angle[g_connected[i]] = circle_graph_perp + PI;
  for (i = domain_row; i <= domain_column; i++)
    g_structure_angle[i] = g_structure_angle[i] + (new_angle);
  set_label_coordinates();
  show_both();
  return;
}

void select_annotation_type(int choice) {
  int  length = -5000;
  if (choice == NONE) {
    g_annotation = NONE;
    show_both();
    return;
  }
  if (choice == g_annotation)
    return;
  else
    g_annotation = choice;
  if (choice == P_NUM) {
    length = 
      define_map_ann(g_annotation_filename, TRUE, FALSE, TRUE, 
		     g_unused, g_create_ann_table, 
		     g_annotation_table_type, g_length);
    set_extra_ps_colors(FALSE);
  }
  if (choice == SS_COUNT) {
    length = 
      define_map_ann(g_annotation_filename, FALSE, FALSE, TRUE,
		     g_unused, g_create_ann_table, 
		     g_annotation_table_type, g_length);
    set_extra_ps_colors(FALSE);
  }
  if (choice == PROB) {
    length = 
      define_map_ann(g_annotation_filename, TRUE, TRUE, TRUE,
		     g_unused, g_create_ann_table, 
		     g_annotation_table_type, g_length);
    set_extra_ps_colors(TRUE);
  }
  if (length != g_length) {
    if (length < 1)
      g_annotation = NONE;
    else
      printf("Warning!\tLength %d from annotation file does not match %d.\n",
	     length, g_length);
  }
  show_both();
}

void select_annotation_option(int choice) {
  g_annotation_bases = TRUE;
  g_annotation_dots = TRUE;
  if (choice == 0) {
    g_annotation_dots = FALSE;
	} else {
	  if (choice == 1) {
	    g_annotation_bases = FALSE;
	  }
	}
  show_both();
}

void make_menu(int show_mouse) {
  int             gif_submenu, jpg_submenu, png_submenu, run_submenu;
  int             label_submenu, menu_shift;
  int             z_submenu;
  int             flat_submenu;
  int             deg_submenu;
  int             reg_submenu;
  int             rot_submenu;
  int             rot_c_submenu; /* M. Zuker Dec 20, 2006 */
  int             rot_cc_submenu; /* M. Zuker Dec 20, 2006 */
  int             mouse_function_submenu = -5000;
  int             redraw_domain_submenu;
  int             redraw_ex_domain_submenu;
  int             annotation_submenu;
  int             annotation_option_submenu;

  if (g_current_window == g_window_nozoom) {
    menu_shift = 0;
  } else {
    menu_shift = 10;
  }
  run_submenu = glutCreateMenu(window_run_sir_graph);
  glutAddMenuEntry("Included Domain, whole", 99);
  glutAddMenuEntry("Included Domain, partial", 97);
  glutAddMenuEntry("Excluded Domain, whole", 98);
  glutAddMenuEntry("Excluded Domain, partial", 96); /* M. Zuker Dec
						       20, 2006. */ 
  flat_submenu = glutCreateMenu(draw_flat);
  glutAddMenuEntry("Default", 80);
  glutAddMenuEntry("Flat", 81);
  glutAddMenuEntry("Flat Alt", 82);
  
  annotation_submenu = glutCreateMenu(select_annotation_type);
  glutAddMenuEntry("Off", 0);
  glutAddMenuEntry("p-num", 1);
  glutAddMenuEntry("probability", 2);
  glutAddMenuEntry("ss-count", 3);
  annotation_option_submenu = glutCreateMenu(select_annotation_option);
  glutAddMenuEntry("Bases", 0);
  glutAddMenuEntry("Dot", 1);
  glutAddMenuEntry("Both", 2);
  redraw_ex_domain_submenu = glutCreateMenu(window_redraw_ex_domain);
  glutAddMenuEntry("Natural Angles", 1);
  glutAddMenuEntry("Circle Graph", 2);
  glutAddMenuEntry("Circle Graph 356", 356);
  glutAddMenuEntry("Circle Graph 270", 270);
  glutAddMenuEntry("Circle Graph 180", 180);
  glutAddMenuEntry("Circle Graph 120", 120);
  glutAddMenuEntry("Circle Graph 90", 90);
  glutAddMenuEntry("Circle Graph 60", 60);
  glutAddMenuEntry("Circle Graph 45", 45);
  redraw_domain_submenu = glutCreateMenu(window_redraw_domain);
  glutAddMenuEntry("Natural Angles", 1);
  glutAddMenuEntry("Circle Graph", 2);
  glutAddMenuEntry("Circle Graph 356", 356);
  glutAddMenuEntry("Circle Graph 270", 270);
  glutAddMenuEntry("Circle Graph 180", 180);
  glutAddMenuEntry("Circle Graph 120", 120);
  glutAddMenuEntry("Circle Graph 90", 90);
  glutAddMenuEntry("Circle Graph 60", 60);
  glutAddMenuEntry("Circle Graph 45", 45);
  /* M. Zuker Dec 20, 2006.  begin */
  rot_submenu = glutCreateMenu(c_or_cc);
  rot_c_submenu = glutCreateMenu(window_set_rot);
  glutAddMenuEntry("None", 0);
  glutAddMenuEntry(" 1", -1);
  glutAddMenuEntry(" 5", -5);
  glutAddMenuEntry(" 10", -10);
  glutAddMenuEntry(" 30", -30);
  glutAddMenuEntry(" 45", -45);
  glutAddMenuEntry(" 90", -90);
  glutAddMenuEntry(" 60", -60);
  glutAddMenuEntry(" 120", -120);
  glutAddMenuEntry(" 135", -135);
  glutAddMenuEntry(" 180", 180);
  glutSetMenu(rot_submenu);
  glutAddSubMenu("Clockwise", rot_c_submenu);
  rot_cc_submenu = glutCreateMenu(window_set_rot);
  glutAddMenuEntry(" 1", 1);
  glutAddMenuEntry(" 5", 5);
  glutAddMenuEntry(" 10", 10);
  glutAddMenuEntry(" 30", 30);
  glutAddMenuEntry(" 45", 45);
  glutAddMenuEntry(" 90", 90);
  glutAddMenuEntry(" 60", 60);
  glutAddMenuEntry(" 120", 120);
  glutAddMenuEntry(" 135", 135);
  glutSetMenu(rot_submenu);
  glutAddSubMenu("Counterclockwise", rot_cc_submenu);
  /* M. Zuker Dec 20, 2006.  end */
  deg_submenu = glutCreateMenu(window_set_deg);
  glutAddMenuEntry("356", 356);
  glutAddMenuEntry("270", 270);
  glutAddMenuEntry("180", 180);
  glutAddMenuEntry("120", 120);
  glutAddMenuEntry("90", 90);
  glutAddMenuEntry("60", 60);
  glutAddMenuEntry("45", 45);
  glutAddMenuEntry("20", 20);
  reg_submenu = glutCreateMenu(window_set_reg);
  glutAddMenuEntry("Off", 0);
  glutAddMenuEntry("90", 90);
  glutAddMenuEntry("60", 60);
  glutAddMenuEntry("45", 45);
  glutAddMenuEntry("40", 40);
  glutAddMenuEntry("36", 36);
  glutAddMenuEntry("30", 30);
  glutAddMenuEntry("24", 24);
  glutAddMenuEntry("20", 20);
  glutAddMenuEntry("18", 18);
  glutAddMenuEntry("10", 10);
  glutAddMenuEntry("5", 5);
  z_submenu = glutCreateMenu(window_set_z);
  glutAddMenuEntry("Off", 0);
  /*  glutAddMenuEntry("1.0", 1); */
  glutAddMenuEntry("2.0", 2);
  glutAddMenuEntry("0.5", 3);
  glutAddMenuEntry("0.1", 4);
  glutAddMenuEntry("0.05", 5);
  glutAddMenuEntry("0.01", 6);
  glutAddMenuEntry("0.005", 7);
  glutAddMenuEntry("0.0001", 8);
  label_submenu = glutCreateMenu(window_set_label);
  glutAddMenuEntry("None", 0);
  glutAddMenuEntry("2", 2);
  glutAddMenuEntry("5", 5);
  glutAddMenuEntry("10", 10);
  glutAddMenuEntry("20", 20);
  glutAddMenuEntry("25", 25);
  glutAddMenuEntry("40", 40);
  glutAddMenuEntry("50", 50);
  glutAddMenuEntry("100", 100);
  gif_submenu = glutCreateMenu(window_create_gif_res);
  glutAddMenuEntry("Dimensions:  Default", menu_shift);
  glutAddMenuEntry("300 x 392", menu_shift + 1);
  glutAddMenuEntry("510 x 666", menu_shift + 2);
  glutAddMenuEntry("612 x 800 ", menu_shift + 3);
  glutAddMenuEntry("765 x 1000 ", menu_shift + 4);
  glutAddMenuEntry("900 x 1176 ", menu_shift + 5);
  glutAddMenuEntry("1210 x 1582", menu_shift + 6);
  jpg_submenu = glutCreateMenu(window_create_jpg_res);
  glutAddMenuEntry("Dimensions:  Default", menu_shift);
  glutAddMenuEntry("300 x 392", menu_shift + 1);
  glutAddMenuEntry("510 x 666", menu_shift + 2);
  glutAddMenuEntry("612 x 800 ", menu_shift + 3);
  glutAddMenuEntry("765 x 1000 ", menu_shift + 4);
  glutAddMenuEntry("900 x 1176 ", menu_shift + 5);
  glutAddMenuEntry("1210 x 1582", menu_shift + 6);
  png_submenu = glutCreateMenu(window_create_png_res);
  glutAddMenuEntry("Dimensions:  Default", menu_shift);
  glutAddMenuEntry("300 x 392", menu_shift + 1);
  glutAddMenuEntry("510 x 666", menu_shift + 2);
  glutAddMenuEntry("612 x 800 ", menu_shift + 3);
  glutAddMenuEntry("765 x 1000 ", menu_shift + 4);
  glutAddMenuEntry("900 x 1176 ", menu_shift + 5);
  glutAddMenuEntry("1210 x 1582", menu_shift + 6);
  mouse_function_submenu = glutCreateMenu(mouse_function_submenu_choice);
  glutAddMenuEntry("Zoom", 1);
  glutAddMenuEntry("Edit", 4);
  glutAddMenuEntry("Loop Stretch", 9);
  glutAddMenuEntry("Loop Natural", 15);
  glutAddMenuEntry("Loop Regularize", 16);
  glutAddMenuEntry("Loop Fix", 35);
  glutAddMenuEntry("Single Strand Stretch", 20);
  glutAddMenuEntry("Single Strand linear", 25);
  /* Begin */
  glutCreateMenu(No_Action_Menu);
  glutAddMenuEntry("Undo", 26);
  glutAddSubMenu("Mouse Function:", mouse_function_submenu);
  glutAddMenuEntry("Scale to fit window", 19);
  glutAddSubMenu("Rotate Structure", rot_submenu);
  glutAddSubMenu("Base Labels", label_submenu);
  glutAddSubMenu("Set Regularize angle", reg_submenu);
  glutAddMenuEntry("Regularize All Angles", 72);
  glutAddMenuEntry("Exit", 1);
  glutAddMenuEntry("________Toggle Options________", -1);
  glutAddMenuEntry("Base Pair: lines/dots", 5);
  glutAddMenuEntry("Outline mode: On/Off", 4);
  glutAddMenuEntry("Loop labels: On/Off", 7);
  glutAddMenuEntry("Draw Clockwise/Counter-Clockwise", 42);
  glutAddMenuEntry("Fix Loops On/Off", 68);
  glutAddMenuEntry("______Annotation Options______", -1);
  glutAddSubMenu("Type", annotation_submenu);
  glutAddSubMenu("Mode", annotation_option_submenu);
  glutAddMenuEntry("_________Draw Options_________", -1);
  glutAddMenuEntry("Draw with Natural Angles", 6);
  glutAddMenuEntry("Draw with Circle Graph Angles", 61);
  glutAddSubMenu("Flat options", flat_submenu);
  glutAddSubMenu("Draw with Z value Circle Graph Angles", z_submenu);
  glutAddSubMenu("Use degrees of Circle for above two", deg_submenu);
  glutAddSubMenu("Draw Included domain with", redraw_domain_submenu);
  glutAddSubMenu("Draw Excluded domain with", redraw_ex_domain_submenu);
  glutAddMenuEntry("________Output________", -4);
  glutAddMenuEntry("Create Postscript", 35);
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  glutAddSubMenu("Create gif", gif_submenu);
#endif
#if HAVE_LIBPNG
  glutAddSubMenu("Create png", png_submenu);
#endif
#if HAVE_LIBJPEG
  glutAddSubMenu("Create jpg", jpg_submenu);
#endif
  glutAddMenuEntry("Create ss output", 25);
  glutAddSubMenu("Run sir_graph", run_submenu);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void gfxinit_nontext(void) {
  /* initalize graphics for each window */
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
  glutReshapeFunc(reshape_window);
  make_menu(TRUE);
}

void gfxinit_nontext_nozoom(void) {	
  /* initalize graphics for each window */
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
  glutReshapeFunc(reshape_window_nozoom);
  make_menu(FALSE);
}

void input_reshape(int w, int h) {
  if (w > 700) {
    w = 700; 
    glutReshapeWindow(w, h);
  }
  if (h > 100) {
    h = 100;
    glutReshapeWindow(w, h);
  }
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 600.0, 0.0, 100.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
}

void text_reshape(int w, int h) {
  if (w > 500) {
    w = 500;
    glutReshapeWindow(w, h);
  }
  if (h > 100) {
    h = 100;
    glutReshapeWindow(w, h);
  }
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 300.0, 0.0, 100.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
}

void gfxinit_text(void) {
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
  glutReshapeFunc(text_reshape);
  make_menu(TRUE);
}

static void Key2(int key, int x, int y) {
  float           x_change;
  float           y_change;
  switch (key) {
  case GLUT_KEY_LEFT:
    x_change = 1.0;
    y_change = 0;
    g_display_window_x_shift -= 10 * g_display_window_zoom_factor * x_change;
    g_display_window_y_shift -= 10 * g_display_window_zoom_factor * y_change;
    break;
  case GLUT_KEY_RIGHT:
    x_change = 1.0;
    y_change = 0.0;
    g_display_window_x_shift += 10 * g_display_window_zoom_factor * x_change;
    g_display_window_y_shift += 10 * g_display_window_zoom_factor * y_change;
    break;
  case GLUT_KEY_UP:
    x_change = 0.0;
    y_change = 1.0;
    g_display_window_y_shift += 10 * g_display_window_zoom_factor * y_change;
    g_display_window_x_shift += 10 * g_display_window_zoom_factor * x_change;
    break;
  case GLUT_KEY_DOWN:
    x_change = 0.0;
    y_change = 1.0;
    g_display_window_x_shift -= 10 * g_display_window_zoom_factor * x_change;
    g_display_window_y_shift -= 10 * g_display_window_zoom_factor * y_change;
    break;
  case GLUT_KEY_PAGE_UP: {
    g_display_window_zoom_factor = g_display_window_zoom_factor * .9;
    if (g_display_window_zoom_factor < .004)
      g_display_window_zoom_factor = .004;
    reshape_window(g_window_1_width, g_window_1_height);
  }
    break;
  case GLUT_KEY_PAGE_DOWN: {
    g_display_window_zoom_factor =
      g_display_window_zoom_factor / .9;
    if (g_display_window_zoom_factor > 10000.)
      g_display_window_zoom_factor = 10000;
    reshape_window(g_window_1_width, g_window_1_height);
  }
    break;
  default:
    return;
  }
  show_both();
}

static void Key2_nozoom(int key, int x, int y) {
  float           x_change;
  float           y_change;
  switch (key) {
  case GLUT_KEY_LEFT:
    x_change = 1.0;
    y_change = 0;
    g_display_window_x_shift += 10 * g_display_window_zoom_factor * x_change;
    g_display_window_y_shift += 10 * g_display_window_zoom_factor * y_change;
    break;
  case GLUT_KEY_RIGHT:
    x_change = 1.0;
    y_change = 0.0;
    g_display_window_x_shift -= 10 * g_display_window_zoom_factor * x_change;
    g_display_window_y_shift -= 10 * g_display_window_zoom_factor * y_change;
    break;
  case GLUT_KEY_UP:
    x_change = 0.0;
    y_change = 1.0;
    g_display_window_y_shift -= 10 * g_display_window_zoom_factor * y_change;
    g_display_window_x_shift -= 10 * g_display_window_zoom_factor * x_change;
    break;
  case GLUT_KEY_DOWN:
    x_change = 0.0;
    y_change = 1.0;
    g_display_window_x_shift += 10 * g_display_window_zoom_factor * x_change;
    g_display_window_y_shift += 10 * g_display_window_zoom_factor * y_change;
    break;
  case GLUT_KEY_PAGE_UP: {
    g_display_window_zoom_factor =
      g_display_window_zoom_factor * .9;
    if (g_display_window_zoom_factor < .004)
      g_display_window_zoom_factor = .004;
    glutSetWindow(g_window_1);
    reshape_window(g_window_1_width, g_window_1_height);
  }
    break;
  case GLUT_KEY_PAGE_DOWN: {
    g_display_window_zoom_factor =
      g_display_window_zoom_factor / .9;
    if (g_display_window_zoom_factor > 10000.)
      g_display_window_zoom_factor = 10000;
    glutSetWindow(g_window_1);
    reshape_window(g_window_1_width, g_window_1_height);
  }
    break;
  default:
    return;
  }
  show_both();
}

void mouse_zoom(int button, int state, int x, int y) {
  float           center_x, center_y;
  float           amount_to_shift;
  float           angle;
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      center_x = g_window_1_width / 2.0;
      center_y = g_window_1_height / 2.0;
      /* end of rotation changes */
      y = g_window_1_height - 1 - y;
      
      amount_to_shift = sqrt((x - center_x) * (x - center_x) +
			     (y - center_y) * (y - center_y));
      amount_to_shift = amount_to_shift * g_display_window_zoom_factor;
      /* Dec 1, 2006. M. Zuker */
      if (x==center_x && y==center_y) { 
	angle = 1.0; /* Has no effect */
      } else {
	angle = (float) atan2(y - center_y, x - center_x);
      }
      g_display_window_x_shift = g_display_window_x_shift -
	amount_to_shift * cos(angle) * (WD + 72.0)/g_window_1_width;
      g_display_window_y_shift = g_display_window_y_shift -
	amount_to_shift*sin(angle)*(HT + 72.0)/g_window_1_height;
      g_display_window_zoom_factor = 0.6 * g_display_window_zoom_factor;
      if (g_display_window_zoom_factor < .004)
	g_display_window_zoom_factor = .004;
      reshape_window(g_window_1_width, g_window_1_height);
      glutSetWindow(g_window_nozoom);
      glutPostRedisplay();
    }
  }
  if (button == GLUT_MIDDLE_BUTTON) {
    if (state == GLUT_DOWN) {
      center_x = g_window_1_width / 2.0;
      center_y = g_window_1_height / 2.0;
      g_display_window_zoom_factor =
	g_display_window_zoom_factor / .6;
      if (g_display_window_zoom_factor > 10000.)
	g_display_window_zoom_factor = 10000;
      reshape_window(g_window_1_width, g_window_1_height);
      glutSetWindow(g_window_nozoom);
      glutPostRedisplay();
    }
  }
}

void mouse_nozoom_zoom(int button, int state, int x, int y) {
  float           center_x, center_y;
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      center_x = g_window_nozoom_width / 2.0;
      center_y = g_window_nozoom_height / 2.0;
      /* end of rotation changes */
      y = g_window_nozoom_height - 1 - y;
      g_display_window_x_shift = -(WD + 72) * (x - center_x) / 
	g_window_nozoom_width;
      g_display_window_y_shift = -(HT + 72.0) * (y - center_y) /
	g_window_nozoom_height;
    }
  }
}

float distance_base(int base, float x, float y) {
  return ((x - g_base_x[base]) * (x - g_base_x[base]) +
	  (y - g_base_y[base]) * (y - g_base_y[base]));
}

float find_base(float x_scaled, float y_scaled) {
  float           closeness, dist;
  int             best_base;
  int             base;
  closeness = distance_base(g_start_base, x_scaled, y_scaled);
  best_base = g_start_base;
  for (base = (g_start_base + 1); base <= g_end_base; base++) {
    dist = distance_base(base, x_scaled, y_scaled);
    if (dist < closeness) {
      closeness = dist;
      best_base = base;
    }
  }
  if (g_oldconnect[best_base] > 0) {
    if (g_oldconnect[best_base] < best_base)
      best_base = g_oldconnect[best_base];
  }
  if (g_oldconnect[best_base] > 0) {
    g_domain_defined = TRUE;
    if (best_base < g_oldconnect[best_base]) {
      g_domain_row = best_base;
      g_domain_column = g_oldconnect[best_base];
    } else {
      g_domain_column = best_base;
      g_domain_row = g_oldconnect[best_base];
    }
  } else {
    g_domain_defined = FALSE;
    show_text();
    g_domain_row = best_base;
  }
  g_clicked_row = g_domain_row;
  g_clicked_column = g_domain_column;
  show_text();
  return sqrt(closeness);
}

void rotate_helix_set_start(float x, float y) {
  g_mouse_rotate_helix_last_x = x;
  g_mouse_rotate_helix_last_y = y;
  g_mouse_rotate_helix_center_x =
    (g_base_x[g_domain_column] + g_base_x[g_domain_row])/2.;
  g_mouse_rotate_helix_center_y =
    (g_base_y[g_domain_column] + g_base_y[g_domain_row])/2.;
	g_mouse_rotate_helix_set = TRUE;
}

void translate_base_or_domain_set_start(float x, float y) {
  int             i;
  g_mouse_rotate_helix_last_x = x;
  g_mouse_rotate_helix_last_y = y;
  if ((g_connected[g_domain_row] > g_domain_row) &&
      (g_connected[g_domain_row] <= g_end_base)) {
    g_mouse_translate_base = FALSE;
  } else {
    g_mouse_translate_base = TRUE;
    g_flat_edit = FALSE;
  }
  if (!g_mouse_translate_base) {
    i = g_domain_row;
    while (g_connected[i - 1] == (g_connected[i] + 1))
      i--;
    g_domain_row = i;
    g_domain_column = g_connected[g_domain_row];
    if (g_flat) {
      g_flat_edit = find_flat_external_loop(g_domain_row);
    } else
      g_flat_edit = FALSE;
  }
  g_mouse_translate_pos_set = TRUE;
}


void adjust_radius_of_loop(float new_radius, float old_radius,
			   float center_x, float center_y, int start_base,
			   int end_base, int start_base2, int end_base2) {
  int             i;
  int             base;
  float           x_change, y_change;
  float           change;
  float           angle;
  float           pos_x, pos_y;
  int             loop;
  int             first_loop, end_loop;
  int             exterior_loop;
  exterior_loop = FALSE;
  if ((start_base == start_base2) && (end_base == end_base2))
    exterior_loop = TRUE;
  change = new_radius - old_radius;
  i = start_base + 1;
  first_loop = g_previous_loop[i] + 1;
  end_loop = g_previous_loop[i];
  if (exterior_loop) {
    i = g_start_base;
    end_base = g_end_base + 1;
  }
  while (i < end_base) {
    if (g_connected[i] < 1) {
      angle = (float) atan2(g_base_y[i] - center_y, g_base_x[i] - center_x);
      x_change = change * cos(angle);
      y_change = change * sin(angle);
      
      g_base_x[i] = g_base_x[i] + x_change;
      g_base_y[i] = g_base_y[i] + y_change;
      i++;
    } else {
      pos_x = (g_base_x[i] + g_base_x[g_connected[i]]) / 2.;
      pos_y = (g_base_y[i] + g_base_y[g_connected[i]]) / 2.;
      angle = (float) atan2(pos_y - center_y, pos_x - center_x);
      x_change = change * cos(angle);
      y_change = change * sin(angle);
      
      for (base = i; base <= g_connected[i]; base++) {
	g_base_x[base] = g_base_x[base] + x_change;
	g_base_y[base] = g_base_y[base] + y_change;
      }
      end_loop = g_previous_loop[g_connected[i]];
      for (loop = (g_previous_loop[i] + 1);
	   loop <= g_previous_loop[g_connected[i]]; loop++) {
	g_loop_center_x[loop] += x_change;
	g_loop_center_y[loop] += y_change;
      }
      if (g_connected[i] > i)
	i = g_connected[i];
      i++;
    }
  }
  /* shift excluded domain */
  if (!exterior_loop) {
    i = start_base;
    pos_x = (g_base_x[i] + g_base_x[g_connected[i]]) / 2.;
    pos_y = (g_base_y[i] + g_base_y[g_connected[i]]) / 2.;
    angle = (float) atan2(pos_y - center_y, pos_x - center_x);
    x_change = change * cos(angle);
    y_change = change * sin(angle);
    for (base = g_start_base; base <= g_end_base; base++) {
      if (base == (start_base + 1))
	base = end_base;
      
      g_base_x[base] = g_base_x[base] + x_change;
      g_base_y[base] = g_base_y[base] + y_change;
    }
    for (loop = 1; loop <= g_total_loops; loop++) {
      g_loop_center_x[loop] = g_loop_center_x[loop] + x_change;
      g_loop_center_y[loop] = g_loop_center_y[loop] + y_change;
      if (loop == (first_loop - 2))
	loop = end_loop;
    }
  }
}

void redo_loop_regularize(int domain_start, int domain_end) {
  float   current_x, current_y, distance_base_pairs, last_angle, angle_new;
  float   angle_old, old_x_translate, old_y_translate, x_translate;
  float   y_translate, angle_translate, x_center_of_rotation;
  float   y_center_of_rotation, original_angle;
  float   angle_change = 0.0; /* M. Zuker found this uninitialized! */  
  int     loop_num, base, base_fix, start_of_fix, end_of_fix, loop, first_loop;
  int     extra_base_exists, i;
  loop_num = g_previous_loop[domain_start] + 1;
  if ((domain_start == g_start_base) &&
      (g_connected[domain_start] < domain_start))
    loop_num--;
  if (domain_start == g_start_base) {
    first_loop = TRUE;
    while (g_connected[domain_start] < domain_start)
      domain_start++;
    domain_end = g_connected[domain_start];
    current_x = (g_base_x[domain_start] + g_base_x[domain_end]) / 2;
    current_y = (g_base_y[domain_start] + g_base_y[domain_end]) / 2;
    loop_num = 1;
  } else {
    first_loop = FALSE;
    current_x = (g_base_x[domain_start] + g_base_x[domain_end]) / 2;
    current_y = (g_base_y[domain_start] + g_base_y[domain_end]) / 2;
  }
  /* find distance base pairs */
  if (!first_loop) {
    distance_base_pairs = 
      sqrt((g_base_x[domain_start] - g_base_x[domain_end]) *
	   (g_base_x[domain_start] - g_base_x[domain_end]) +
	   (g_base_y[domain_start] - g_base_y[domain_end]) *
	   (g_base_y[domain_start] - g_base_y[domain_end]));
  } else {
    base = g_start_base;
    while ((g_connected[base] < g_start_base) && (base < g_end_base)) {
      base++;
    }
    if (base < g_end_base) {
      distance_base_pairs = 
	sqrt((g_base_x[base] - g_base_x[g_connected[base]]) *
	     (g_base_x[base] - g_base_x[g_connected[base]]) +
	     (g_base_y[base] - g_base_y[g_connected[base]]) *
	     (g_base_y[base] - g_base_y[g_connected[base]]));
    } else
      distance_base_pairs = DISTANCE_BASE_PAIRS * g_scale;
  }
  last_angle = g_structure_angle[domain_start];
  if (!first_loop)
    last_angle -= PI / 2.;
  else {
    last_angle = g_structure_angle[domain_start] + PI / 2;
  }
  
  store_for_undo();
  if (first_loop) { /* if this is the first loop and the first base is
		     * paired, rotate to give the first helix the
		     * correct angle */ 
    if (!g_counter_clockwise)
      original_angle = g_structure_angle[domain_start] - PI / 2.;
    else
      original_angle = g_structure_angle[domain_start] + PI / 2.;
    while (original_angle > (2.0 * PI))
      original_angle -= (2.0 * PI);
    while (original_angle < 0)
      original_angle += (2.0 * PI);
    angle_change = regularize_angle(original_angle) - original_angle;
    rotate_structure(angle_change, (WD + 72.0)/2.0, (HT + 72.0)/2.0, 
		     TRUE, 1, g_total_loops, g_start_base, g_end_base);
    for (i = g_start_base; i <= g_end_base; i++)
      g_structure_angle[i] = g_structure_angle[i] + angle_change;
    scale_to_fit_window();
    store_for_undo();
    if (g_end_base == g_connected[domain_start])
      extra_base_exists = FALSE;
    else
      extra_base_exists = TRUE;
    domain_start = g_start_base;
    domain_end = g_end_base;
  } else
    extra_base_exists = FALSE;
  /*
   * if this is the first loop and the first base is paired, rotate to
   * give the first helix the correct angle 
   */
  
  if (first_loop) { /* find first helix, let domain start after it,
		       finish before it */ 
    i = domain_start;
    while (g_connected[i] < i)
      i++;
    domain_start = i;
    domain_end = i - 1;
    traverse_loop_structure_regularize_loop(loop_num, domain_start, domain_end,
					    g_connected, last_angle + 
					    angle_change, current_x, current_y,
					    extra_base_exists, g_end_base, 
					    distance_base_pairs, FALSE, 0.0);
  } else {
    traverse_loop_structure_regularize_loop(loop_num, domain_start, domain_end,
					    g_connected, last_angle, current_x,
					    current_y, extra_base_exists, 
					    domain_end, distance_base_pairs, 
					    FALSE, 0.0);
  }
  /* fix structure beyond loop; rotate and shift helices coming off
     loop to match the loop */ 
  if (first_loop) {
    base = g_start_base;
    domain_end = g_end_base;
  } else
    base = domain_start + 1;
  while (base < domain_end) {
    if (g_connected[base] > base) { /* translate */
      start_of_fix = base + 1;
      end_of_fix = g_connected[base] - 1;
      x_translate = -1.0 * (g_base_x[base] + g_base_x[g_connected[base]]) / 2
	+ (g_undo1_base_x[base] + g_undo1_base_x[g_connected[base]]) / 2.;
      y_translate = -1.0 * (g_base_y[base] + g_base_y[g_connected[base]]) / 2
	+ (g_undo1_base_y[base] + g_undo1_base_y[g_connected[base]]) / 2.;
      angle_new = g_structure_angle[base];
      angle_old = g_undo1_structure_angle[base];
      while (angle_new < 0)
	angle_new += 2.0 * PI;
      while (angle_old < 0)
	angle_old += 2.0 * PI;
      angle_translate = angle_new - angle_old;
      x_center_of_rotation = (g_undo1_base_x[base] +
			      g_undo1_base_x[g_connected[base]]) / 2.;
      y_center_of_rotation = (g_undo1_base_y[base] +
			      g_undo1_base_y[g_connected[base]]) / 2.;
      
      /* rotate */
      for (base_fix = start_of_fix; base_fix <= end_of_fix; base_fix++) {
	old_x_translate = g_undo1_base_x[base_fix] - x_center_of_rotation;
	old_y_translate = g_undo1_base_y[base_fix] - y_center_of_rotation;
	g_base_x[base_fix] = x_center_of_rotation + 
	  cos(angle_translate) * old_x_translate -
	  sin(angle_translate) * old_y_translate;
	g_base_y[base_fix] = y_center_of_rotation +
	  sin(angle_translate) * old_x_translate +
	  cos(angle_translate) * old_y_translate;
      }
      for (loop = (g_previous_loop[base] + 1);
	   loop <= g_previous_loop[g_connected[base]]; loop++) {
	old_x_translate = g_loop_center_x[loop] - x_center_of_rotation;
	old_y_translate = g_loop_center_y[loop] - y_center_of_rotation;
	g_loop_center_x[loop] = x_center_of_rotation +
	  cos(angle_translate) * old_x_translate -
	  sin(angle_translate) * old_y_translate;
	g_loop_center_y[loop] = y_center_of_rotation +
	  sin(angle_translate) * old_x_translate +
	  cos(angle_translate) * old_y_translate;
      }
      /* translate */
      for (base_fix = start_of_fix; base_fix <= end_of_fix; base_fix++) {
	g_base_x[base_fix] = g_base_x[base_fix] - x_translate;
	g_base_y[base_fix] = g_base_y[base_fix] - y_translate;
	g_structure_angle[base_fix] = g_structure_angle[base_fix]
	  + angle_translate;
      }
      for (loop = (g_previous_loop[base] + 1);
	   loop <= g_previous_loop[g_connected[base]]; loop++) {
	g_loop_center_x[loop] = g_loop_center_x[loop] - x_translate;
	g_loop_center_y[loop] = g_loop_center_y[loop] - y_translate;
      }
      base = g_connected[base];
    }
    base++;
  }
  set_label_coordinates();
  show_both();
}

void redo_loop_natural( int domain_start, int domain_end) {
  int             bulge_flag, interior_flag;
  int             loop_size;
  int             loop_num;
  float           current_x, current_y;
  float           distance_base_pairs;
  float           last_angle;
  float           angle_new, angle_old;
  int             base;
  int             base_fix;
  int             start_of_fix, end_of_fix;
  int             loop;
  int             first_loop;
  float           old_x_translate, old_y_translate;
  float           x_translate, y_translate, angle_translate;
  float           x_center_of_rotation;
  float           y_center_of_rotation;
  loop_num = g_previous_loop[domain_start] + 1;
  if (domain_start == g_start_base) {
    first_loop = TRUE;
    while (g_connected[domain_start] < domain_start)
      domain_start++;
    domain_end = g_connected[domain_start];
    current_x = (g_base_x[domain_start] + g_base_x[domain_end]) / 2;
    current_y = (g_base_y[domain_start] + g_base_y[domain_end]) / 2;
    loop_num = 1;
  } else {
    first_loop = FALSE;
    current_x = (g_base_x[domain_start] + g_base_x[domain_end]) / 2;
    current_y = (g_base_y[domain_start] + g_base_y[domain_end]) / 2;
  }
  /* find distance base pairs */
  if (!first_loop) {
    distance_base_pairs = 
      sqrt((g_base_x[domain_start] - g_base_x[domain_end]) *
	   (g_base_x[domain_start] - g_base_x[domain_end]) +
	   (g_base_y[domain_start] - g_base_y[domain_end]) *
	   (g_base_y[domain_start] - g_base_y[domain_end]));
  } else {
    base = g_start_base;
    while ((g_connected[base] < g_start_base) && (base < g_end_base)) {
      base++;
    }
    if (base < g_end_base) {
      distance_base_pairs = 
	sqrt((g_base_x[base] - g_base_x[g_connected[base]]) *
	     (g_base_x[base] - g_base_x[g_connected[base]]) +
	     (g_base_y[base] - g_base_y[g_connected[base]]) *
	     (g_base_y[base] - g_base_y[g_connected[base]]));
    } else
      distance_base_pairs = DISTANCE_BASE_PAIRS * g_scale;
  }
  last_angle = g_structure_angle[domain_start];
  if (!first_loop)
    last_angle -= PI / 2.;
  else {
    last_angle = g_structure_angle[domain_start] + PI / 2;
    bulge_flag = FALSE;
    interior_flag = FALSE;
  }
  
  if (first_loop) {
    loop_size = 
      traverse_loop_natural_count_bases(loop_num, domain_end, domain_start,
					first_loop,&bulge_flag,&interior_flag);
    traverse_loop_natural_assign_angles(loop_size, domain_end, domain_start, 
					current_x, current_y, first_loop, 
					last_angle, loop_num, 
					distance_base_pairs, bulge_flag, 
					interior_flag);
    domain_start = g_start_base;
    domain_end = g_end_base;
  } else {
    loop_size = 
      traverse_loop_natural_count_bases(loop_num, domain_start, domain_end,
					first_loop,&bulge_flag,&interior_flag);
    traverse_loop_natural_assign_angles(loop_size, domain_start, domain_end, 
					current_x, current_y, first_loop,
					last_angle, loop_num, 
					distance_base_pairs, bulge_flag, 
					interior_flag);
  }
  /* fix structure beyond loop */
  if (first_loop) {
    base = g_start_base;
    domain_end = g_end_base;
  } else
    base = domain_start + 1;
  while (base < domain_end) {
    if (g_connected[base] > base) {	/* translate */
      start_of_fix = base + 1;
      end_of_fix = g_connected[base] - 1;
      x_translate = -1.0 * (g_base_x[base] + g_base_x[g_connected[base]]) / 2
	+ (g_undo1_base_x[base] +
	   g_undo1_base_x[g_connected[base]]) / 2.;
      y_translate = -1.0 * (g_base_y[base] + g_base_y[g_connected[base]]) / 2
	+ (g_undo1_base_y[base] +
	   g_undo1_base_y[g_connected[base]]) / 2.;
      angle_new = g_structure_angle[base];
      angle_old = g_undo1_structure_angle[base];
      while (angle_new < 0)
	angle_new += 2.0 * PI;
      while (angle_old < 0)
	angle_old += 2.0 * PI;
      angle_translate = angle_new - angle_old;
      x_center_of_rotation = (g_undo1_base_x[base] +
			      g_undo1_base_x[g_connected[base]]) / 2.;
      y_center_of_rotation = (g_undo1_base_y[base] +
			      g_undo1_base_y[g_connected[base]]) / 2.;
      
      /* rotate */
      for (base_fix = start_of_fix; base_fix <= end_of_fix; base_fix++) {
	old_x_translate = g_undo1_base_x[base_fix] -
	  x_center_of_rotation;
	old_y_translate = g_undo1_base_y[base_fix] -
	  y_center_of_rotation;
	g_base_x[base_fix] = x_center_of_rotation +
	  cos(angle_translate) * old_x_translate -
	  sin(angle_translate) * old_y_translate;
	g_base_y[base_fix] = y_center_of_rotation +
	  sin(angle_translate) * old_x_translate +
	  cos(angle_translate) * old_y_translate;
	
      }
      for (loop = (g_previous_loop[base] + 1);
	   loop <= g_previous_loop[g_connected[base]]; loop++) {
	old_x_translate = g_loop_center_x[loop] - x_center_of_rotation;
	old_y_translate = g_loop_center_y[loop] - y_center_of_rotation;
	g_loop_center_x[loop] = x_center_of_rotation +
	  cos(angle_translate) * old_x_translate - 
	  sin(angle_translate) * old_y_translate;
	g_loop_center_y[loop] = y_center_of_rotation +
	  sin(angle_translate) * old_x_translate +
	  cos(angle_translate) * old_y_translate;
      }
      /* translate */
      for (base_fix = start_of_fix; base_fix <= end_of_fix; base_fix++) {
	g_base_x[base_fix] = g_base_x[base_fix] - x_translate;
	g_base_y[base_fix] = g_base_y[base_fix] - y_translate;
	g_structure_angle[base_fix] = g_structure_angle[base_fix] 
	  + angle_translate;
      }
      for (loop = (g_previous_loop[base] + 1);
	   loop <= g_previous_loop[g_connected[base]]; loop++) {
	g_loop_center_x[loop] = g_loop_center_x[loop] - x_translate;
	g_loop_center_y[loop] = g_loop_center_y[loop] - y_translate;
      }
      base = g_connected[base];
    }
    base++;
  }
  set_label_coordinates();
  show_both();
}

void strand_set_start(float x, float y) {
  /* x and y are real coordinates of clicked point */
  int             base;
  int             current_loop;
  /* find center of loop */
  
  base = g_mouse_loop_stretch_domain_row1;
  while ((g_oldconnect[base] < g_mouse_loop_stretch_domain_column1) &&
	 (base >= g_start_base)) {
    base--;
  }
  if (base < g_start_base)
    current_loop = 1;
  else
    current_loop = g_previous_loop[base] + 1;
  g_mouse_strand_stretch_first_x = x;
  g_mouse_strand_stretch_first_y = y;
  g_mouse_strand_stretch_center_x = g_loop_center_x[current_loop];
  g_mouse_strand_stretch_center_y = g_loop_center_y[current_loop];
}

void mouse_strand_set_start_linear(void) {
  int             base;
  int             current_loop;
  /* find loop */
  base = g_mouse_loop_stretch_domain_row1;
  while ((g_oldconnect[base] < g_mouse_loop_stretch_domain_column1) &&
	 (base >= g_start_base)) {
    base--;
  }
  if (base < g_start_base)
    current_loop = 1;
  else
    current_loop = g_previous_loop[base] + 1;
  g_mouse_strand_linear_loop = current_loop;
}

void mouse_stretch_set_start(float x, float y) {
  /* x and y are real coordinates of clicked point */
  g_mouse_strand_stretch_first_x = x;
  g_mouse_strand_stretch_first_y = y;
}

void mouse_nozoom(int button, int state, int x1, int y1) {
  mouse_nozoom_zoom(button, state, x1, y1);
  show_both();
  return;
}

void mouse(int button, int state, int x1, int y1) {
  float           x, y;
  float           distance;
  float           base_pair_x, base_pair_y;
  int             i;
  int             current_loop;
  int             domain_row, domain_column;
  if (g_mouse_display_zoom_mode) {
    mouse_zoom(button, state, x1, y1);
    return;
  }
  if ((g_mouse_translate) || (g_mouse_rotate_helix)) {
    if (button == GLUT_LEFT_BUTTON) {
      g_need_to_store_for_undo = TRUE;
      if (state == GLUT_DOWN) {
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor +
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor +
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	if (distance < (1.5 * g_scale)) { /* adjust for sensitivity to drag */
	  g_close_enough = TRUE;
	  if (!g_domain_defined) {
	    g_mouse_translate_pos_set = TRUE;
	    translate_base_or_domain_set_start(x, y);
	  } else {
	    i = g_domain_row;
	    while (g_oldconnect[i - 1] == (g_oldconnect[i] + 1))
	      i--;
	    domain_row = i;
	    domain_column = g_oldconnect[domain_row];
	    if (domain_row == g_domain_row) {	
	      /* This is the start of a helix. Translate */
	      translate_base_or_domain_set_start(x, y);
	    } else {	
	      /* This is not the start of a helix. Rotate */
	      g_domain_row = domain_row;
	      g_domain_column = domain_column;
	      rotate_helix_set_start(x, y);
	    }
	  }
	  
	} else {
	  g_close_enough = FALSE;
	}
      } else {
	g_need_to_store_for_undo = FALSE;
	if ((g_mouse_rotate_helix) && (g_mouse_rotate_helix_set == TRUE)) {
	  g_mouse_rotate_helix_set = FALSE;
	  
	}
	if (g_mouse_translate) {
	  g_mouse_translate_pos_set = FALSE;
	}
      }
    }
  }
  /* end of translate or rotate helix */
  
  if (g_mouse_loop_stretch) {
    if (button == GLUT_LEFT_BUTTON) {
      g_need_to_store_for_undo = TRUE;
      if (state == GLUT_DOWN) {
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor +
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor +
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	/* go back to start of domain */
	if (g_domain_defined) {	
	  i = g_domain_row;
	  /* go to the first base pair of the helix if this is not the
	     last base pair of the helix 
	  */
	  if (g_oldconnect[i + 1] == (g_oldconnect[i] - 1)) {
	    while (g_oldconnect[i - 1] == (g_oldconnect[i] + 1)) {
	      i--;
	    }
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row2 column2 to first base pair of a helix on this loop */
	    g_mouse_loop_stretch_domain_column2 = g_domain_column;
	    g_mouse_loop_stretch_domain_row2 = i;
	    /* set row1, column1 to last base pair on first helix on loop */
	    if (i > g_start_base)
	      i--;
	    while ((g_oldconnect[i] < g_domain_column)
		   && (i > g_start_base)) {
	      if (g_oldconnect[i] > 0)
		i = g_oldconnect[i];
	      i--;
	    }
	    if (i < g_start_base)
	      i = g_start_base;
	    g_mouse_loop_stretch_domain_row1 = i;
	    if (i > g_start_base)
	      g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	    else
	      g_mouse_loop_stretch_domain_column1 = g_end_base;
	  } else {
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1,column 1 to last base pair on first helix on loop */
	    g_mouse_loop_stretch_domain_row1 = i;
	    g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	    /* set row2,column2 to some helix on the loop */
	    i++;
	    while ((g_oldconnect[i] < i) && (i < g_domain_column)) {
	      i++;
	    }
	    if (i < g_domain_column) {
	      g_mouse_loop_stretch_domain_row2 = i;
	      g_mouse_loop_stretch_domain_column2 = g_oldconnect[i];
	    } else {
	      g_mouse_loop_stretch_domain_row2 = i;
	      g_mouse_loop_stretch_domain_column2 = g_oldconnect[i];
	    }
	  }
	}
	/* proceed only if this base is not paired */
	if (!g_domain_defined) {
	  g_mouse_loop_stretch_pos_x = x;
	  g_mouse_loop_stretch_pos_y = y;
	  g_mouse_loop_stretch_pos_set = TRUE;
	  /* find start and end of single stranded part
	   *
	   * go back to a paired base */
	  i = g_domain_row;
	  while ((g_oldconnect[i] < 1) && (i >= g_start_base)) {
	    i--;
	  }
	  i++;
	  g_mouse_loop_stretch_domain_row = i;
	  /* g_mouse_loop_stretch_domain_row currently set single stranded
	   *
	   * base in loop just after a helix or g_start_base
	   */
	  i = g_domain_row;
	  while ((g_oldconnect[i] < 1) && (i <= (g_end_base - 1)))
	    i++;
	  i--;
	  /* g_mouse_loop_stretch_domain_column currently is set to
	   * single stranded base in loop just before a helix, or
	   * g_end_base 
	   */
	  g_mouse_loop_stretch_domain_column = i;
	  if ((g_mouse_loop_stretch_domain_row - 1) <
	      g_oldconnect[g_mouse_loop_stretch_domain_row - 1]) {
	    g_mouse_loop_stretch_domain_row1 = 
	      g_mouse_loop_stretch_domain_row - 1;
	    g_mouse_loop_stretch_domain_column1 =
	      g_oldconnect[g_mouse_loop_stretch_domain_row - 1];
	    g_mouse_loop_stretch_domain_row2 =
	      g_mouse_loop_stretch_domain_column + 1;
	    g_mouse_loop_stretch_domain_column2 =
	      g_oldconnect[g_mouse_loop_stretch_domain_column + 1];
	  } else {
	    g_mouse_loop_stretch_domain_row1 =
	      g_oldconnect[g_mouse_loop_stretch_domain_column + 1];
	    g_mouse_loop_stretch_domain_column1 =
	      g_mouse_loop_stretch_domain_column + 1;
	    g_mouse_loop_stretch_domain_row2 =
	      g_oldconnect[g_mouse_loop_stretch_domain_row - 1];
	    g_mouse_loop_stretch_domain_column2 =
	      g_mouse_loop_stretch_domain_row - 1;
	  }
	  i = g_mouse_loop_stretch_domain_row2 - 1;
	  if (i < g_start_base)
	    i = g_start_base;
	  if (g_mouse_loop_stretch_domain_column2 > g_end_base)
	    g_mouse_loop_stretch_domain_column2 = g_end_base;
	  while ((g_oldconnect[i] < g_mouse_loop_stretch_domain_column2) && 
		 (i > g_start_base)) {	/* go to starting bp of loop */
	    if (g_oldconnect[i] > 0)
	      i = g_oldconnect[i];
	    i--;
	  }
	  g_mouse_loop_stretch_domain_row1 = i;
	  g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	  if (g_mouse_loop_stretch_domain_column1 < 1)
	    g_mouse_loop_stretch_domain_column = g_end_base;
	}
	if (g_mouse_loop_stretch_domain_column1 < 1)
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	if (g_mouse_loop_stretch_domain_column2 == 0) {
	  g_mouse_loop_stretch_domain_row2 = g_start_base;
	  g_mouse_loop_stretch_domain_column2 = g_end_base;
	}
	if ((g_mouse_loop_stretch_domain_row1 == g_start_base) &&
	    (g_mouse_loop_stretch_domain_column == g_end_base)) {
	  g_mouse_loop_stretch_domain_row2 = g_start_base;
	  g_mouse_loop_stretch_domain_column2 = g_end_base;
	}
	if (g_mouse_loop_stretch_domain_row1 < g_start_base)
	  g_mouse_loop_stretch_domain_row1 = g_start_base;
	if (g_mouse_loop_stretch_domain_row1 == g_start_base) {
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	  g_mouse_loop_stretch_domain_row2 = g_start_base;
	  g_mouse_loop_stretch_domain_column2 = g_end_base;
	  current_loop = 1;
	} else {
	  current_loop =
	    g_previous_loop[g_mouse_loop_stretch_domain_row1] + 1;
	}
	/* find center of loop based on these coordinates */
	g_mouse_loop_stretch_center_x = g_loop_center_x[current_loop];
	g_mouse_loop_stretch_center_y = g_loop_center_y[current_loop];
	base_pair_x = (g_base_x[g_mouse_loop_stretch_domain_row1] + 
		       g_base_x[g_mouse_loop_stretch_domain_column1]) / 2.;
	base_pair_y = (g_base_y[g_mouse_loop_stretch_domain_row1] + 
		       g_base_y[g_mouse_loop_stretch_domain_column1]) / 2.;
	g_mouse_loop_stretch_r = 
	  sqrt((g_mouse_loop_stretch_center_x - base_pair_x) *
	       (g_mouse_loop_stretch_center_x - base_pair_x) +
	       (g_mouse_loop_stretch_center_y - base_pair_y) *
	       (g_mouse_loop_stretch_center_y - base_pair_y));
	g_mouse_loop_stretch_start_r =
	  sqrt((g_mouse_loop_stretch_center_x - x) *
	       (g_mouse_loop_stretch_center_x - x) +
	       (g_mouse_loop_stretch_center_y - y) *
	       (g_mouse_loop_stretch_center_y - y));
	g_mouse_loop_stretch_pos_set = TRUE;
      } else {
	g_need_to_store_for_undo = FALSE;
	g_mouse_loop_stretch_pos_set = FALSE;
      }
    }
  }
  /* end of loop stretch 
   *
   * start of Loop Natural */
  if (g_mouse_loop_natural) {
    if (button == GLUT_LEFT_BUTTON) {
      if (state == GLUT_DOWN) {
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor +
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor +
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	/* go back to start of domain */
	if (g_domain_defined) {	
	  i = g_domain_row;
	  /* go to first base pair of helix if this is not the last
	   * base pair of a helix 
	   */
	  if (g_oldconnect[i + 1] == (g_oldconnect[i] - 1)) {
	    while (g_oldconnect[i - 1] == (g_oldconnect[i] + 1))
	      i--;
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1, column1 to last base pair on first helix on loop */
	    if (i > g_start_base)
	      i--;
	    while ((g_oldconnect[i] < g_domain_column) && (i > g_start_base)) {
	      if (g_oldconnect[i] > 0)
		i = g_oldconnect[i];
	      i--;
	    }
	    if (i < g_start_base)
	      i = g_start_base;
	    g_mouse_loop_stretch_domain_row1 = i;
	    if (i > g_start_base)
	      g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	    else
	      g_mouse_loop_stretch_domain_column1 = g_end_base;
	  } else {
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1,column 1 to last base pair on first helix on loop */
	    g_mouse_loop_stretch_domain_row1 = i;
	    g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	  }
	}
	/* proceed only if this base is not paired */
	if (!g_domain_defined) {
	  g_mouse_loop_stretch_pos_x = x;
	  g_mouse_loop_stretch_pos_y = y;
	  g_mouse_loop_stretch_pos_set = TRUE;
	  /* find first last base pair of first helix on this loop;
	   * go back to a paired base 
	   */
	  i = g_domain_row;
	  while ((g_oldconnect[i] < g_domain_row) && (i >= g_start_base))
	    i--;
	  g_mouse_loop_stretch_domain_row1 = i;
	  g_mouse_loop_stretch_domain_column1 =
	    g_oldconnect[g_mouse_loop_stretch_domain_row1];
	}
	if (g_mouse_loop_stretch_domain_column1 < 1)
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	if (g_mouse_loop_stretch_domain_row1 < g_start_base)
	  g_mouse_loop_stretch_domain_row1 = g_start_base;
	if (g_mouse_loop_stretch_domain_row1 == g_start_base) {
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	}
	store_for_undo();
	redo_loop_natural(g_mouse_loop_stretch_domain_row1,
			  g_mouse_loop_stretch_domain_column1);
      }
    }
  }
  /* end of mouse loop natural
   * start of Loop Regularize 
   */
  if (g_mouse_loop_regularize) {
    if (button == GLUT_LEFT_BUTTON) {
      if (state == GLUT_DOWN) {
	if (!g_reg_angle_flag) {
	  printf("Ignored!\tYou must first set a regularize angle.\n");
	  return;
	}
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor +
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor +
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	/* go back to start of domain */
	if (g_domain_defined) {	
	  i = g_domain_row;
	  /* Go to the first base pair of the helix if this is not the
	   * last base pair of the helix.
	   */
	  if (g_oldconnect[i + 1] == (g_oldconnect[i] - 1)) {
	    while (g_oldconnect[i - 1] == (g_oldconnect[i] + 1)) {
	      i--;
	    }
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1, column1 to last base pair on first helix on loop */ 
	    if (i > g_start_base)
	      i--;
	    while ((g_oldconnect[i] < g_domain_column)
		   && (i > g_start_base)) {
	      if (g_oldconnect[i] > 0)
		i = g_oldconnect[i];
	      i--;
	    }
	    if (i < g_start_base)
	      i = g_start_base;
	    g_mouse_loop_stretch_domain_row1 = i;
	    if (i > g_start_base)
	      g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	    else
	      g_mouse_loop_stretch_domain_column1 = g_end_base;
	  } else {
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1,column 1 to last base pair on first helix on loop */
	    g_mouse_loop_stretch_domain_row1 = i;
	    g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	  }
	}
	/* proceed only if this base is not paired */
	if (!g_domain_defined) {
	  g_mouse_loop_stretch_pos_x = x;
	  g_mouse_loop_stretch_pos_y = y;
	  g_mouse_loop_stretch_pos_set = TRUE;
	  /* find first last base pair of first helix on this loop
	   * go back to a paired base 
	   */
	  i = g_domain_row;
	  while ((g_oldconnect[i] < g_domain_row) && (i >= g_start_base))
	    i--;
	  g_mouse_loop_stretch_domain_row1 = i;
	  g_mouse_loop_stretch_domain_column1 =
	    g_oldconnect[g_mouse_loop_stretch_domain_row1];
	}
	if (g_mouse_loop_stretch_domain_column1 < 1)
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	if (g_mouse_loop_stretch_domain_row1 < g_start_base)
	  g_mouse_loop_stretch_domain_row1 = g_start_base;
	if (g_mouse_loop_stretch_domain_row1 == g_start_base) {
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	}
	redo_loop_regularize(g_mouse_loop_stretch_domain_row1,
			     g_mouse_loop_stretch_domain_column1);
      }
    }
  }
  /* end of mouse loop Regularize
   * start of Loop Fix 
   */
  if (g_mouse_loop_fix) {
    if (button == GLUT_LEFT_BUTTON) {
      if (state == GLUT_DOWN) {
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor + 
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	/* go back to start of domain */
	if (g_domain_defined) {	
	  i = g_domain_row;
	  /* go to first base pair of helix if this is not the last
	     base pair of a helix 
	  */
	  if (g_oldconnect[i + 1] == (g_oldconnect[i] - 1)) {
	    while (g_oldconnect[i - 1] == (g_oldconnect[i] + 1))
	      i--;
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1, column1 to last base pair on first helix on loop */
	    if (i > g_start_base)
	      i--;
	    while ((g_oldconnect[i] < g_domain_column) && (i > g_start_base)) {
	      if (g_oldconnect[i] > 0)
		i = g_oldconnect[i];
	      i--;
	    }
	    if (i < g_start_base)
	      i = g_start_base;
	    g_mouse_loop_stretch_domain_row1 = i;
	    if (i > g_start_base)
	      g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	    else
	      g_mouse_loop_stretch_domain_column1 = g_end_base;
	  } else {
	    g_domain_row = i;
	    g_domain_column = g_oldconnect[i];
	    /* set row1,column 1 to last base pair on first helix on loop */
	    g_mouse_loop_stretch_domain_row1 = i;
	    g_mouse_loop_stretch_domain_column1 = g_oldconnect[i];
	  }
	}
	/* proceed only if this base is not paired */
	if (!g_domain_defined) {
	  g_mouse_loop_stretch_pos_x = x;
	  g_mouse_loop_stretch_pos_y = y;
	  g_mouse_loop_stretch_pos_set = TRUE;
	  /* find first last base pair of first helix on this loop 
	   * go back to a paired base 
	   */
	  i = g_domain_row;
	  while ((g_oldconnect[i] < g_domain_row) && (i >= g_start_base))
	    i--;
	  g_mouse_loop_stretch_domain_row1 = i;
	  g_mouse_loop_stretch_domain_column1 =
	    g_oldconnect[g_mouse_loop_stretch_domain_row1];
	}
	if (g_mouse_loop_stretch_domain_column1 < 1)
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	if (g_mouse_loop_stretch_domain_row1 < g_start_base)
	  g_mouse_loop_stretch_domain_row1 = g_start_base;
	if (g_mouse_loop_stretch_domain_row1 == g_start_base) {
	  g_mouse_loop_stretch_domain_column1 = g_end_base;
	}
	store_for_undo();
	loop_fix(g_mouse_loop_stretch_domain_row1,
		 g_mouse_loop_stretch_domain_column1, 0.0, 0.0, TRUE, 0);
	set_label_coordinates();
      }
      show_both();
    }
  }
  /* end of mouse loop Fix */
  if (g_mouse_strand_stretch) {
    if (button == GLUT_LEFT_BUTTON) {
      g_need_to_store_for_undo = TRUE;
      if (state == GLUT_DOWN) {
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0)* g_display_window_zoom_factor + 
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0)*g_display_window_zoom_factor + 
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	if (!g_domain_defined) {
	  g_mouse_loop_stretch_pos_x = x;
	  g_mouse_loop_stretch_pos_y = y;
	  g_mouse_loop_stretch_pos_set = TRUE;
	  /* find first last base pair of first helix on this loop 
	   * go back to a paired base 
	   */
	  i = g_domain_row;
	  while ((g_oldconnect[i] < 1) && (i >= g_start_base))
	    i--;
	  i++;
	  g_mouse_loop_stretch_domain_row1 = i;
	  i = g_domain_row;
	  while ((g_oldconnect[i] < 1) && (i <= g_end_base))
	    i++;
	  i--;
	  g_mouse_loop_stretch_domain_column1 = i;
	  strand_set_start(x, y);
	  g_strand_set = TRUE;
	} else {
	  g_strand_set = FALSE;
	}
      } else { /* button is up */
	g_strand_set = FALSE;
      }
    }
  }
  if (g_mouse_strand_linear) {
    if (button == GLUT_LEFT_BUTTON) {
      g_need_to_store_for_undo = TRUE;
      if (state == GLUT_DOWN) {
	y = (float) (g_window_1_height - 1 - y1);
	x = (float) x1;
	y = y * (HT + 72.0) / g_window_1_height;
	x = x * (WD + 72.0) / g_window_1_width;
	y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
	  (HT + 72.0)/2.0 - g_display_window_y_shift;
	x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor + 
	  (WD + 72.0)/2.0 - g_display_window_x_shift;
	distance = find_base(x, y);
	if (!g_domain_defined) {
	  g_mouse_loop_stretch_pos_x = x;
	  g_mouse_loop_stretch_pos_y = y;
	  g_mouse_loop_stretch_pos_set = TRUE;
	  /* find first last base pair of first helix on this loop 
	   * go back to a paired base 
	   */
	  i = g_domain_row;
	  g_mouse_strand_linear_base = i;
	  while ((g_oldconnect[i] < 1) && (i >= g_start_base))
	    i--;
	  i++;
	  g_mouse_loop_stretch_domain_row1 = i;
	  i = g_domain_row;
	  /* go forward to a paired base */
	  while ((g_oldconnect[i] < 1) && (i <= g_end_base))
	    i++;
	  i--;
	  g_mouse_loop_stretch_domain_column1 = i;
	  mouse_strand_set_start_linear();
	  g_strand_set = TRUE;
	} else
	  g_strand_set = FALSE;
      } else { /* button is up */
	g_strand_set = FALSE;
      }
    }
  }
}

void mouse_motion(int x1, int y1) {
  float           angle1, angle2, angle;
  float           x, y;
  float           x_dif, y_dif;
  float           new_radius;
  float           base_angle;
  int             i;
  int             base;
  int             extra_base;
  float           dist_x, dist_y;
  float           base_offset;
  float           linear_x_adjust;
  float           linear_y_adjust;
  int             linear_stationary_base;
  float           new_distance, old_distance, change_distance;
  if ( !g_close_enough && !g_mouse_loop_stretch && 
       !g_mouse_strand_stretch && !g_mouse_strand_linear )
    return;
  if (g_need_to_store_for_undo) {
    g_need_to_store_for_undo = FALSE;
    store_for_undo();
  }
  if ((g_mouse_rotate_helix) && (g_mouse_rotate_helix_set == TRUE)) {
    y = (float) (g_window_1_height - 1 - y1);
    x = (float) x1;
    y = y * (HT + 72.0) / g_window_1_height;
    x = x * (WD + 72.0) / g_window_1_width;
    y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
      (HT + 72.0)/2.0 - g_display_window_y_shift;
    x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor + 
      (WD + 72.0)/2.0 - g_display_window_x_shift;
    if ((x == g_mouse_rotate_helix_last_x) && 
	(y == g_mouse_rotate_helix_last_y)) {
      printf("No rotation performed\n");
      return;
    }
    angle1 = (float) atan2(g_mouse_rotate_helix_last_y -
			   g_mouse_rotate_helix_center_y,
			   g_mouse_rotate_helix_last_x -
			   g_mouse_rotate_helix_center_x);
    angle2 = (float) atan2(y - g_mouse_rotate_helix_center_y,
			   x - g_mouse_rotate_helix_center_x);
    angle = angle2 - angle1;
    rotate_structure(angle, g_mouse_rotate_helix_center_x,
		     g_mouse_rotate_helix_center_y, FALSE, 
		     g_previous_loop[g_domain_row] + 1,
		     g_previous_loop[g_domain_column],
		     g_domain_row, g_domain_column);
    for (i = (g_domain_row); (i <= g_domain_column); i++)
      g_structure_angle[i] = g_structure_angle[i] + angle;
    set_label_coordinates();
    g_mouse_rotate_helix_last_x = x;
    g_mouse_rotate_helix_last_y = y;
    show_both();
  }
  if ((g_mouse_translate == TRUE) && (g_mouse_translate_pos_set == TRUE)) {
    y = (float) (g_window_1_height - 1 - y1);
    x = (float) x1;
    y = y * (HT + 72.0) / g_window_1_height;
    x = x * (WD + 72.0) / g_window_1_width;
    y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
      (HT + 72.0)/2.0 - g_display_window_y_shift;
    x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor + 
      (WD + 72.0)/2.0 - g_display_window_x_shift;
    if ( (x == g_mouse_rotate_helix_last_x) && 
	 (y == g_mouse_rotate_helix_last_y)) {
      printf("No translation performed\n");
      return;
    }
    if (g_mouse_translate_base) { /* translate only a base */
      g_base_x[g_domain_row] = g_base_x[g_domain_row] +
	x - g_mouse_rotate_helix_last_x;
      g_base_y[g_domain_row] = g_base_y[g_domain_row] +
	y - g_mouse_rotate_helix_last_y;
      /* fix structure angle on this base */
      if (g_domain_row != g_end_base)
	g_structure_angle[g_domain_row] =
	  find_angle_ss(g_domain_row, g_domain_row + 1);
      else
	g_structure_angle[g_domain_row] =
	  find_angle_ss(g_domain_row - 1, g_domain_row);
      if ((g_ss_backwards) || (g_counter_clockwise))
	g_structure_angle[g_domain_row] += PI;
      /* fix structure angle on previous base */
      if (g_domain_row < g_end_base) {
	extra_base = g_domain_row - 1;
	if (g_oldconnect[extra_base] < 1) {
	  g_structure_angle[extra_base]=find_angle_ss(extra_base,extra_base+1);
	  if ((g_ss_backwards) || (g_counter_clockwise))
	    g_structure_angle[extra_base] += PI;
	}
      }
    } else {	/* translating domain */
      x_dif = x - g_mouse_rotate_helix_last_x;
      y_dif = y - g_mouse_rotate_helix_last_y;
      for (base = g_domain_row; base <= g_domain_column; base++) {
	g_base_x[base] += x_dif;
	g_base_y[base] += y_dif;
      }
      for (i = (g_previous_loop[g_domain_row] + 1);
	   i <= g_previous_loop[g_domain_column]; i++) {
	g_loop_center_x[i] += x_dif;
	g_loop_center_y[i] += y_dif;
      }
      /* for flat edit, stretch out the single stranded bases */
      if (g_flat_edit) {
	if (g_domain_row - g_flat_fixed_base > 1) {
	  dist_x = (g_base_x[g_domain_row] - g_base_x[g_flat_fixed_base]) /
	    (g_domain_row - g_flat_fixed_base);
	  dist_y = (g_base_y[g_domain_row] - g_base_y[g_flat_fixed_base]) /
	    (g_domain_row - g_flat_fixed_base);
	  for (i = (g_flat_fixed_base + 1); i < g_domain_row; i++) {
	    g_base_x[i] = g_base_x[g_flat_fixed_base] +
	      (float) (i - g_flat_fixed_base) * dist_x;
	    g_base_y[i] = g_base_y[g_flat_fixed_base] +
	      (float) (i - g_flat_fixed_base) * dist_y;
	  }
	}
	for (i = (g_domain_column + 1); i <= g_end_base; i++) {
	  g_base_x[i] += x_dif;
	  g_base_y[i] += y_dif;
	}
	for (i = (g_previous_loop[g_domain_column] + 1); i <= g_total_loops;
	     i++) {
	  g_loop_center_x[i] += x_dif;
	  g_loop_center_y[i] += y_dif;
	}
      }
    }
    set_label_coordinates();
    g_mouse_rotate_helix_last_x = x;
    g_mouse_rotate_helix_last_y = y;
    show_both();
  }
  if ((g_mouse_loop_stretch) && (g_mouse_loop_stretch_pos_set)) {
    y = (float) (g_window_1_height - 1 - y1);
    x = (float) x1;
    y = y * (HT + 72.0) / g_window_1_height;
    x = x * (WD + 72.0) / g_window_1_width;
    y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
      (HT + 72.0)/2.0 - g_display_window_y_shift;
    x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor + 
      (WD + 72.0)/2.0 - g_display_window_x_shift;
    new_radius = sqrt((g_mouse_loop_stretch_center_x - x) *
		      (g_mouse_loop_stretch_center_x - x) +
		      (g_mouse_loop_stretch_center_y - y) *
		      (g_mouse_loop_stretch_center_y - y));
    /* travel around shifting each base to account for new radius */
    adjust_radius_of_loop(new_radius, g_mouse_loop_stretch_start_r,
			  g_mouse_loop_stretch_center_x,
			  g_mouse_loop_stretch_center_y,
			  g_mouse_loop_stretch_domain_row1,
			  g_mouse_loop_stretch_domain_column1,
			  g_mouse_loop_stretch_domain_row2,
			  g_mouse_loop_stretch_domain_column2);
    g_mouse_loop_stretch_start_r = new_radius;
    set_label_coordinates();
    show_both();
  }
  if ((g_mouse_strand_stretch) && (g_strand_set)) {
    y = (float) (g_window_1_height - 1 - y1);
    x = (float) x1;
    y = y * (HT + 72.0) / g_window_1_height;
    x = x * (WD + 72.0) / g_window_1_width;
    y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
      (HT + 72.0)/2.0 - g_display_window_y_shift;
    x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor + 
      (WD + 72.0)/2.0 - g_display_window_x_shift;
    old_distance = 
      sqrt((g_mouse_strand_stretch_first_x - g_mouse_strand_stretch_center_x) *
	   (g_mouse_strand_stretch_first_x - g_mouse_strand_stretch_center_x) +
	   (g_mouse_strand_stretch_first_y - g_mouse_strand_stretch_center_y) *
	   (g_mouse_strand_stretch_first_y - g_mouse_strand_stretch_center_y));
    new_distance = sqrt((x - g_mouse_strand_stretch_center_x) *
			(x - g_mouse_strand_stretch_center_x) +
			(y - g_mouse_strand_stretch_center_y) *
			(y - g_mouse_strand_stretch_center_y));
    change_distance = new_distance - old_distance;
    for (base = g_mouse_loop_stretch_domain_row1;
	 base <= g_mouse_loop_stretch_domain_column1; base++) {
      base_offset = (g_mouse_loop_stretch_domain_column1 +
		     g_mouse_loop_stretch_domain_row1) / 2.;
      if (base >= base_offset)
	base_offset = base - base_offset;
      else
	base_offset = base_offset - base;
      base_angle = (float) 
	atan2((g_undo1_base_y[base] - g_mouse_strand_stretch_center_y),
	      (g_undo1_base_x[base] - g_mouse_strand_stretch_center_x));
      base_offset *= .92;
      g_base_x[base] = g_undo1_base_x[base] + 
	cos(base_angle)*change_distance*
	(1 - base_offset*2./(g_mouse_loop_stretch_domain_column1 -
			     g_mouse_loop_stretch_domain_row1 + 1));
      g_base_y[base] = g_undo1_base_y[base] + sin(base_angle)*change_distance* 
	(1 - base_offset*2./(g_mouse_loop_stretch_domain_column1 -
			     g_mouse_loop_stretch_domain_row1 + 1));;
    }
    set_label_coordinates();
    show_both();
  }
  if ((g_mouse_strand_linear) && (g_strand_set)) {
    y = (float) (g_window_1_height - 1 - y1);
    x = (float) x1;
    y = y * (HT + 72.0) / g_window_1_height;
    x = x * (WD + 72.0) / g_window_1_width;
    y = (y - (HT + 72.0)/2.0) * g_display_window_zoom_factor + 
      (HT + 72.0)/2.0 - g_display_window_y_shift;
    x = (x - (WD + 72.0)/2.0) * g_display_window_zoom_factor +
      (WD + 72.0)/2.0 - g_display_window_x_shift;
    /* move the initial base */
    base = g_mouse_strand_linear_base;
    g_base_x[base] = x;
    g_base_y[base] = y;
    /* adjust bases prior to g_mouse_strand_linear_base */
    linear_stationary_base = g_mouse_loop_stretch_domain_row1 - 1;
    if (g_mouse_loop_stretch_domain_row1 == g_start_base)
      linear_stationary_base++;
    linear_x_adjust = (x - g_base_x[linear_stationary_base]) /
      (base - linear_stationary_base);
    linear_y_adjust = (y - g_base_y[linear_stationary_base]) /
      (base - linear_stationary_base);
    for (base = (linear_stationary_base + 1); base<g_mouse_strand_linear_base; 
	 base++) {
      g_base_x[base] = g_base_x[linear_stationary_base] +
	(base - (linear_stationary_base)) * linear_x_adjust;
      g_base_y[base] = g_base_y[linear_stationary_base] +
	(base - (linear_stationary_base)) * linear_y_adjust;
    }
    /* adjust bases after to g_mouse_strand_linear_base */
    linear_stationary_base = g_mouse_loop_stretch_domain_column1 + 1;
    if (g_mouse_loop_stretch_domain_column1 == g_end_base)
      linear_stationary_base--;
    base = g_mouse_strand_linear_base;
    linear_x_adjust = (g_base_x[linear_stationary_base] - x) /
      (linear_stationary_base - base);
    linear_y_adjust = (g_base_y[linear_stationary_base] - y) /
      (linear_stationary_base - base);
    for (base = (g_mouse_strand_linear_base + 1); base<linear_stationary_base;
	 base++) {
      g_base_x[base] = x + (base - g_mouse_strand_linear_base) *
	linear_x_adjust;
      g_base_y[base] = y + (base - g_mouse_strand_linear_base) *
	linear_y_adjust;
    }
    for (base = g_mouse_loop_stretch_domain_row1;
	 base <= g_mouse_loop_stretch_domain_column1; base++) {
      if (base != g_end_base)
	g_structure_angle[base] = find_angle_ss(base, base + 1);
      else
	g_structure_angle[base] = find_angle_ss(base - 1, base);
      if ((g_ss_backwards) || (g_counter_clockwise))
	g_structure_angle[base] += PI;
    }
    set_label_coordinates();
    show_both();
  }
}

void adjust_name_line(char *old_name) {
  /* take energy value from start of name and place it on the end
     searches for dG = value    */
  char  new_name[140] = "";
  int   len, starting_pos, starting_pos_name, bracket_pair_exists, i;
  /* find dg = */
  len = strlen(old_name);
  i = 0;
  starting_pos = -1;
  while ((i < len - 3) && (starting_pos < 0)) {
    if (((old_name[i] == 'd') && (old_name[i + 1] == 'G')) && 
	((old_name[i + 2] == ' ') && (old_name[i + 3] == '='))) {
      starting_pos = i + 4;
    }
    i++;
  }
  if (starting_pos < 0)
    return;
  /* advance through number */
  i = starting_pos + 2;	/* advance until a number is found */
  starting_pos_name = -1;
  while ((starting_pos_name < 0) && (i < len)) {
    if (isdigit(old_name[i]))
      starting_pos_name = i;
    i++;
  }
  if (starting_pos_name < 0)
    return;
  /* find a space */
  i = starting_pos_name + 1;
  starting_pos_name = -1;
  while ((starting_pos_name < 0) && (i < len)) {
    if (isspace(old_name[i]))
      starting_pos_name = i;
    i++;
  }
  if (starting_pos_name < 0)
    return;
  /* if there is a [ ] pair, put its contents at the end */
  bracket_pair_exists = FALSE;
  i = starting_pos_name + 1;
  while ((bracket_pair_exists == FALSE) && (i < len)) {
    if (old_name[i] == ']') {
      while ((bracket_pair_exists == FALSE) && (i < len)) {
	if (old_name[i] == ']') {
	  bracket_pair_exists = TRUE;
	  starting_pos_name = i + 1;
	}
	i++;
      }
      i = len;
    }
    i++;
  }
  /* starting_pos_name gives the start of the name */
  for (i = starting_pos_name; i < len; i++)
    new_name[i - starting_pos_name] = old_name[i];
  new_name[len-starting_pos_name] = ' ';
  for (i = 0; i < starting_pos_name; i++)
    new_name[len-starting_pos_name + i + 1] = old_name[i];
  if (len>=140) len = 139;
  new_name[len+1] = '\0';
  strcpy(old_name, new_name);
  /* returns a better old_name */
}
/* Don't open input window until necessary */
void inputSetWindow(void) {
  if (g_window_input_exists==FALSE) {
    g_window_input_exists = TRUE;
    glutInitWindowSize(600, 100);
    glutInitWindowPosition(6, 350); 
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    g_window_input = glutCreateWindow(g_name_string);
    glutDisplayFunc(display_window_input);
    glutSetCursor(GLUT_CURSOR_HELP);
    glutReshapeFunc(input_reshape);
    glutKeyboardFunc(keyboard_input);
  }
  glutSetWindow(g_window_input);
}
#endif
/* M. Zuker, June 12, 2007. No need for sir_graph_ng. Use sir_graph
   for all.
   Need function below as a dummy even if no glut 
   If glut is not present, the call will return 0 and exit if
   interactive mode is attempted.
*/
int enable_window(int argc, char **argv) {
#if HAVE_GLUTINIT
  char            name_string[140];
  char            fixed_name[140];
  glutInit(&argc, argv);
  set_undo_data();
  g_need_to_store_for_undo = FALSE;
  g_window_1_width = (WD + 72.0);
  g_window_1_height = (HT + 72.0);
  g_mouse_display_zoom_mode = FALSE;
  g_mouse_translate = TRUE;
  g_mouse_rotate_helix = TRUE;
  g_mouse_loop_natural = FALSE;
  g_mouse_loop_regularize = FALSE;
  g_mouse_strand_stretch = FALSE;
  g_mouse_strand_linear = FALSE;
  g_mouse_loop_fix = FALSE;
  g_strand_set = FALSE;
  g_mouse_rotate_helix_set = FALSE;
  g_mouse_translate_pos_set = FALSE;
  g_display_window_x_shift = 0.0;
  g_display_window_y_shift = 0.0;
  g_display_window_zoom_factor = 1.0;	/* try .8 */
  g_domain_row = -1;
  g_clicked_row = 0;
  g_close_enough = FALSE;
  g_domain_column = -1;
  g_domain_defined = FALSE;
  strcpy(g_input_main_message, "Enter name for output file");
  strcpy(g_input_default_message, "<ENTER> for default: ");
  strcpy(g_input_name, "undefined");
  strcpy(g_input_abort_message, "/ to abort");
    /* start of zoom window */
  glutInitWindowSize((WD + 72.0), (HT + 72.0));
  glutInitWindowPosition(175, 225); 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  strcpy(fixed_name, g_first_line);
  if (g_aj)
    adjust_name_line(fixed_name);
  strcpy(name_string, PROGRAM_NAME": Edit Window ");
  strcat(name_string, fixed_name);
  g_window_1 = glutCreateWindow(name_string);
  glutEntryFunc(zoom_entry);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_motion);
  glutSpecialFunc(Key2);
  gfxinit_nontext();
  glutDisplayFunc(display_window_1);
  
  /* start of main window */
  glutInitWindowSize(153, 198);
  glutInitWindowPosition(14, 225); 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  strcpy(name_string, "Main ");
  strcat(name_string, fixed_name);
  g_window_nozoom = glutCreateWindow(name_string);
  glutEntryFunc(nozoom_entry);
  glutMouseFunc(mouse_nozoom);
  glutSpecialFunc(Key2_nozoom);
  gfxinit_nontext_nozoom();
  glutDisplayFunc(display_window_nozoom);
  /* start of text window */  
  glutInitWindowSize(450, 120);
  glutInitWindowPosition(6, 85);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  strcpy(fixed_name, g_first_line);
  if (g_aj)
    adjust_name_line(fixed_name);
  strcpy(name_string, "Text ");
  strcat(name_string, fixed_name);
  g_window_text = glutCreateWindow(name_string);
  glutDisplayFunc(display_window_text);
  gfxinit_text();
  /* M. Zuker, 2006: start of input window used to be here 
   * It has been moved and the window is created as needed.
   */
  glutMainLoop();
#endif
  printf("%s cannot be used interactively because glut was not found.\n", 
	 PROGRAM_NAME);
  printf("Re-run with a file output flag (-p -png -jpg or -g)\n");
  try_exit(21);
  return 0;  /* Dummy line to keep compiler happy. */
}

/* ******************************************************** 
 * end of interactive window functions 
 ********************************************************* */

int call_file_makers(int argc, char **argv) {
  int i;
  int temp;
  int need_to_adjust_his_domain;
  int color_file_exists;
  int label_range;
  float degrees_used;
  g_smart_colors = FALSE;
  g_interactive_mode = TRUE;
  g_flat = FALSE;
  color_file_exists = FALSE;
  g_color_ann_file_exists = FALSE;
  g_flat_alternate = FALSE;
  g_png_mode = FALSE;
  g_jpg_mode = FALSE;
  g_reg_angle_flag = FALSE;
  g_fix_loop_flag = FALSE;
  g_degrees_used = -1;
  g_diam = 11.0;
  g_complot_forced_bases_count = 0;
  need_to_adjust_his_domain = FALSE;
  strcpy(g_out_filename, argv[argc - 1]);
  strcpy(g_out_ss_filename, g_out_filename);
  chop_suffix(g_out_ss_filename, ".ss");
  g_ss_mode = check_ss_filetype(g_out_filename);
  g_annotation_filename[0] = '\0';
  if (g_ss_mode) {
    strcpy(g_ss_filename, g_out_filename);
    chop_suffix(g_out_filename, ".ss");
    strcpy(g_annotation_filename, g_out_filename);
  } else {
    chop_suffix(g_out_filename, ".ct"); /* removes suffix */
    strcpy(g_annotation_filename, g_out_filename);
    strcpy(g_ct_filename, g_out_filename);
  }
  i = 1;
  argc--;
  while (i <= (argc - 1)) { 
    /* specify name of output file */
    if (strcmp(argv[i], "-o") == 0) {
      strcpy(g_out_filename, argv[i + 1]);
      strcpy(g_out_ss_filename, g_out_filename);
      i += 2;
    } else if (strcmp(argv[i], "-af") == 0) {
      strcpy(g_annotation_filename, argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-t") == 0) {
      strcpy(g_fig_title, argv[i + 1]);
      g_fig_title[strlen(g_fig_title)] = '\0';
      i += 2;
    } else if (strcmp(argv[i], "-force") == 0) {
      strcpy(g_forced_file, argv[i + 1]);
      g_complot_forced_bases_count = 1;
      /* indicates that file exists */
      i += 2;
    } else if (strcmp(argv[i], "-bp") == 0) {
      g_bp = TRUE;
      i++;
    } else if (strcmp(argv[i], "-c") == 0) {
      g_is_circular = TRUE;
      i++;
    } else if (strcmp(argv[i], "-n") == 0) {
      g_natural = TRUE;
      i++;
    } else if (strcmp(argv[i], "-s") == 0) {
      g_smart_colors = TRUE;
      i++;
    } else if (strcmp(argv[i], "-tab_ann_html") == 0) {
      g_create_ann_table = TRUE;
      i++;
    } else if (strcmp(argv[i], "-p") == 0) {
      g_structure = TRUE;
      g_interactive_mode = FALSE;
      i++;
    } else if (strcmp(argv[i], "-a") == 0) {
      g_angle_list = TRUE;
      i++;
    } else if (strcmp(argv[i], "-m") == 0) {
      g_midpoints = TRUE;
      i++;
    } else if (strcmp(argv[i], "-D") == 0) {
      g_lines = FALSE;
      i++;
    } else if (strcmp(argv[i], "-ar") == 0) {
      g_auto_rotate = TRUE;
      i++;
    } else if (strcmp(argv[i], "-ni") == 0) {
      g_img_interlaced = FALSE;
      i++;
    } else if (strcmp(argv[i], "-fix") == 0) {
      g_fix_loop_flag = TRUE;
      i++;
    } else if (strcmp(argv[i], "-aj") == 0 ) {
      g_aj = TRUE;
      i++;
    } else if (strcmp(argv[i], "-f") == 0 ) {
      g_flat = TRUE;
      g_structure = TRUE;
      printf("flat flag detected\n");
      i++;
    } else if (strcmp(argv[i], "-fa") == 0 ) {
      g_flat = TRUE;
      g_structure = TRUE;
      g_flat_alternate = TRUE;
      printf("flat flag detected with alternating\n");
      i++;
    } else if (strcmp(argv[i], "-outline") == 0 ) {
      g_outline_mode = TRUE;
      i++;
    } else if (strcmp(argv[i], "-cg") == 0 ) {
      g_interactive_mode = FALSE;
      g_structure = FALSE;
      i++;
    } else if (strcmp(argv[i], "-mo") == 0 ) {
      g_midpoints = TRUE;
      g_arcs = FALSE;
      i++;
    } else if (strcmp(argv[i], "-loop") == 0 ) {
      g_loop_labels = TRUE;
      i++;
    } else if (strcmp(argv[i], "-ss") == 0 ) {
      g_ss_output_mode = TRUE;
      g_structure = TRUE;
      g_interactive_mode = FALSE;
      i++;
    } else if (strcmp(argv[i], "-z") == 0 ) {
      g_small_angle = TRUE;
      sscanf(argv[i + 1], "%f", &g_helix_base_advance);
      i += 2;
    } else if (strcmp(argv[i], "-col") == 0 ) {
      if ((i + 1) <= (argc - 1)) {
	color_file_exists = TRUE;
	strcpy(color_filename, argv[i + 1]);
	i += 2;
      }
    } else if (strcmp(argv[i], "-col_ann") == 0 ) {
      if ((i + 1) <= (argc - 1)) {
	g_color_ann_file_exists = TRUE;
	strcpy(color_ann_filename, argv[i + 1]);
	printf("Using annotation colors from file: %s\n", color_ann_filename);
	i += 2;
      }
    } else if (strcmp(argv[i], "-g") == 0 ) {
      g_img_mode = TRUE;
      g_structure = TRUE;
      g_interactive_mode = FALSE;
      if ((i < (argc - 1)) && (isdigit(argv[i + 1][0]))) {
	sscanf(argv[i + 1], "%d", &g_img_width);
	i += 2;
	if (g_img_width < 50)
	  g_img_width = 50;
	if (g_img_width > 3200)
	  g_img_width = 3200;
	g_img_height = (int) ( ((float) g_img_width) * 17.0 /13.0);
      } else {
	i++;
	g_img_width = (int) (WD + 72.0);
	g_img_height = (int) (HT + 72.0);
      }
    } else if (strcmp(argv[i], "-png") == 0 ) {
      g_png_mode = TRUE;
      g_img_mode = TRUE;
      g_structure = TRUE;
      g_interactive_mode = FALSE;
      if ((i < (argc - 1)) && (isdigit(argv[i + 1][0]))) {
	sscanf(argv[i + 1], "%d", &g_img_width);
	i += 2;
	if (g_img_width < 50)
	  g_img_width = 50;
	if (g_img_width > 3200)
	  g_img_width = 3200;
	g_img_height = (int) ( ((float) g_img_width) * 17.0/13.0);
      } else {
	i++;
	g_img_width = (int) (WD + 72.0);
	g_img_height = (int) (HT + 72.0);
      }
    } else if (strcmp(argv[i], "-jpg") == 0 ) {
      g_jpg_mode = TRUE;
      g_img_mode = TRUE;
      g_structure = TRUE;
      g_interactive_mode = FALSE;
      if ((i < (argc - 1)) && (isdigit(argv[i + 1][0]))) {
	sscanf(argv[i + 1], "%d", &g_img_width);
	i += 2;
	if (g_img_width < 50)
	  g_img_width = 50;
	if (g_img_width > 3200)
	  g_img_width = 3200;
	g_img_height = (int) ( ((float) g_img_width) * 17.0/13.0);
      } else {
	i++;
	g_img_width = (int) (WD + 72.0);
	g_img_height = (int) (HT + 72.0);
      }
    } else if (strcmp(argv[i], "-lab") == 0 ) {
      sscanf(argv[i + 1], "%d", &g_label_frequency);
      if (g_label_frequency < 0)
	g_label_frequency = 0;
      if (g_label_frequency > 1000)
	g_label_frequency = 1000;
      i += 2;
      /* Check to adjust radians between each base */
    } else if (strcmp(argv[i], "-r") == 0 ) {
      sscanf(argv[i + 1], "%f", &degrees_used);
      if (degrees_used < 0.) {
	degrees_used = -1.0;
	printf("Warning, value for -r was ignored\n");
      } else
	g_degrees_used = degrees_used;
      i += 2;
    } else if (strcmp(argv[i], "-pnum") == 0 ) { /* Check p-num */
      g_annotation = P_NUM;
      i++;
    } else if (strcmp(argv[i], "-prob") == 0 ) {  /* Check prob */
      g_annotation = PROB;
      i++;
    } else if (strcmp(argv[i], "-ss-count") == 0 ) { /* Check ss-count */
      g_annotation = SS_COUNT;
      i++;
    } else if (strcmp(argv[i], "-ab") == 0 ) {   /* Check annotate bases */
      g_annotation_bases = TRUE;
      g_annotation_dots = FALSE;
      i++;
    } else if (strcmp(argv[i], "-ad") == 0 ) {  /* Check annotate dots */
      g_annotation_dots = TRUE;
      g_annotation_bases = FALSE;
      i++;
    } else if (strcmp(argv[i], "-cc") == 0 ) {  /* Check counter-clockwise */
      g_counter_clockwise = TRUE;
      i++;
    } else if (strcmp(argv[i], "-x") == 0 ) {
      g_command_line_x = TRUE;
      i++;
    } else if (strcmp(argv[i], "-rot") == 0 ) { /* Check to rotate structure */
      sscanf(argv[i + 1], "%f", &g_rotation_angle);
      printf("Rotation angle: %.2f degrees\n", g_rotation_angle);
      g_rotation_angle = g_rotation_angle * PI / 180.;
      i += 2;
    } else if (strcmp(argv[i], "-zoom_ps") == 0 ) {  /* Check zoom_ps */
      g_command_line_zoom_ps_flag = TRUE;
      sscanf(argv[i + 1], "%f", &g_command_line_zoom_ps_s);
      sscanf(argv[i + 2], "%d", &g_command_line_zoom_ps_x);
      sscanf(argv[i + 3], "%d", &g_command_line_zoom_ps_y);
      if (g_command_line_zoom_ps_s < .000001) {
	g_command_line_zoom_ps_s = 1.0;
	printf("zoom s of zoom_ps must be greater than 0\n");
      }
      if (g_command_line_zoom_ps_x < 0)
	g_command_line_zoom_ps_x = 0;
      if (g_command_line_zoom_ps_y < 0)
	g_command_line_zoom_ps_y = 0;
      i = i + 4;
    } else if (strcmp(argv[i], "-zoom") == 0 ) {   /* Check zoom */
      g_command_line_zoom_flag = TRUE;
      sscanf(argv[i + 1], "%f", &g_command_line_zoom_s);
      sscanf(argv[i + 2], "%d", &g_command_line_zoom_x);
      sscanf(argv[i + 3], "%d", &g_command_line_zoom_y);
      if (g_command_line_zoom_s < .000001) {
	g_command_line_zoom_s = 1.0;
	printf("zoom s of zoom must be greater than 0\n");
      }
      if (g_command_line_zoom_x < 0)
	g_command_line_zoom_x = 0;
      if (g_command_line_zoom_y < 0)
	g_command_line_zoom_y = 0;
      i = i + 4;
    } else if (strcmp(argv[i], "-e") == 0 ) { /* Check excluded range */
      if (i < argc - 2) {
	sscanf(argv[i + 1], "%d", &g_ex_start_base);
	sscanf(argv[i + 2], "%d", &g_ex_end_base);
	if (g_ex_start_base < 1)
	  g_ex_start_base = 1;
	/* if (g_ex_end_base>MAXLEN) g_ex_end_base=MAXLEN; */
	i = i + 2;
	g_external_domain_flag = TRUE;
      } else
	printf("Insufficient arguments for -e\n");
      i++;
    } else if (strcmp(argv[i], "-ew") == 0 ) { /* Check excluded whole range */
      if (i < argc - 2) {
	sscanf(argv[i + 1], "%d", &g_ex_start_base);
	sscanf(argv[i + 2], "%d", &g_ex_end_base);
	if (g_ex_start_base < 1)
	  g_ex_start_base = 1;
	g_whole_circle = TRUE;
	i += 2;
	g_external_domain_flag = TRUE;
      } else
	printf("Insufficient arguments for -ew\n");
      i++;
    } else if (strcmp(argv[i], "-i") == 0 ) { /* Check included range */
      if (i < argc - 2) {
	sscanf(argv[i + 1], "%d", &g_start_base);
	sscanf(argv[i + 2], "%d", &g_end_base);
	if (g_start_base < 1)
	  g_start_base = 1;
	need_to_adjust_his_domain = TRUE;
	i += 2;
      } else
	printf("Insufficient arguments for -i\n");
      i++;
    } else if (strcmp(argv[i], "-reg") == 0 ) { /* Check to regularize
						   angles of helices
						   about loops */ 
      sscanf(argv[i + 1], "%d", &g_reg_angle);
      g_reg_angle_flag = TRUE;
      printf("Regularize angle is %d degrees\n", g_reg_angle);
      if (g_reg_angle < 1) {
	g_reg_angle = 1;
	printf("Regularize angle raised to 1 degree\n");
      } else if (g_reg_angle > 90) {
	g_reg_angle = 90;
	printf("Regularize angle lowered to 90 degrees\n");
      } else {
	printf("Using regularize angle of %d\n", g_reg_angle);
      }
      i += 2;
    } else if (strcmp(argv[i], "-iw") == 0 ) { /* Check included whole range */
      if (i < argc - 2) {
	sscanf(argv[i + 1], "%d", &g_start_base);
	sscanf(argv[i + 2], "%d", &g_end_base);
	if (g_start_base < 1)
	  g_start_base = 1;
	g_whole_circle = TRUE;
	need_to_adjust_his_domain = TRUE;
	i += 2;
      } else
	printf("Insufficient arguments for -iw\n");
      i++;
    } else if (strcmp(argv[i], "-lij") == 0 ) { /* Check force labels on
						   i and j */ 
      if (i < argc - 2) {
	sscanf(argv[i + 1], "%d", &g_label_forced_row);
	sscanf(argv[i + 2], "%d", &g_label_forced_column);
	if (g_label_forced_row < 1)
	  g_label_forced_row = 0;
	if (g_label_forced_row > g_label_forced_column) {
	  temp = g_label_forced_column;
	  g_label_forced_column = g_label_forced_row;
	  g_label_forced_row = temp;
	}
	i += 2;
      } else
	printf("Insufficient arguments for -lij.\tIgnored.\n");
      i++;
    } else if (strcmp(argv[i], "-d") == 0 ) { /* Check to adjust angles
						 between bases n and 1 */ 
      sscanf(argv[i + 1], "%f", &g_diam);
      if (g_diam > 11.0) {
	g_diam = 11.0;
	printf("Size lowered to 11 inches\n");
      }
      if (g_diam < 2.) {
	g_diam = 2.;
	printf("Size raised to 2 inches\n");
      }
      i += 2;
    } else {
      printf("Warning!\tFlag %s was not recognized.\n", argv[i]);
      i++;
    }
  } /* End of command line parsing */
  if (g_img_mode)
    g_structure = TRUE;
  if (g_natural) {
    g_structure = TRUE;
    g_small_angle = FALSE;
    g_whole_circle = FALSE;
    /*    g_reg_angle_flag = FALSE; M. Zuker removes. 11/27/06 */
  }
  if (g_interactive_mode) {
    g_structure = TRUE;
    g_ss_output_mode = FALSE;
    g_command_line_zoom_flag = FALSE;
    g_command_line_zoom_ps_flag = FALSE;
  }
  if (g_ss_mode)
    g_structure = TRUE;
  if (g_structure) {
    g_midpoints = FALSE;
    g_smart_colors = FALSE;
    g_arcs = FALSE;
  } else {
    g_label_frequency = 0;
    g_reg_angle_flag = FALSE;
  }
  chop_suffix(g_out_filename, ".ps" );
  strcpy(g_structure_filename, g_out_filename);
  strcat(g_out_filename, "_cir.ps");
  set_main_colors(color_file_exists);
  if (g_color_ann_file_exists)
    set_annotation_colors_from_file();
  g_first_ss_run = TRUE;
  if (g_ss_mode) {
    general_read_ss_file(g_ss_filename,g_first_line,&g_length,&g_start_base,
			 &g_end_base,&g_ss_distance_basepair,
			 &g_ss_distance_base_adj,&g_ss_backwards,
			 &g_total_loops);
    g_history_offset = 0;
  } else {
    general_read_ct_file(g_ct_filename, &g_length, g_first_line);
    g_history_offset = g_history[1] - 1;
    if (need_to_adjust_his_domain) {
      g_start_base -= g_history_offset;
      g_end_base -= g_history_offset;
      if (g_label_forced_row > 0) {
	g_label_forced_row -= g_history_offset;
	g_label_forced_column -= g_history_offset;
      }
      if (g_external_domain_flag) {
	g_ex_end_base -= g_history_offset;
	g_ex_start_base -= g_history_offset;
      }
    }
    adjust_start_and_end(&g_start_base, &g_end_base, g_length, g_connected);
  }
  if (g_label_frequency == 0) /* turns off labels */
    g_label_frequency = g_length + 10;
  if (g_external_domain_flag) {
    if (g_connected[g_ex_start_base] != g_ex_end_base) {
      printf("Invalid external domain, Drawing entire structure\n");
      g_external_domain_flag = FALSE;
    }
    if (g_label_frequency < 0)
      g_label_frequency = 50;
  }
  if (g_label_frequency < 0) {
    label_range = g_end_base - g_start_base + 1;
    g_label_frequency = 50;
    if (label_range < 500)
      g_label_frequency = 20;
    if (label_range < 100)
      g_label_frequency = 10;
    if (label_range < 20)
      g_label_frequency = 5;
  }
  if (g_complot_forced_bases_count > 0)
    define_forced_bases();
  temp = g_annotation;
  g_annotation = NONE;
  if (temp != NONE)
    select_annotation_type_ng(temp);
  compute_angles_ps(g_out_filename, g_degrees_used, g_is_circular);
  if (!g_structure)
    try_exit(0);
  adjust_structure_coordinates(); /* set for (WD + 72.0) x (HT + 72.0) */
  g_auto_rotate = FALSE;
  g_rotation_angle = 0;
  if ((!g_interactive_mode)) {
    if (g_img_mode) {
      make_structure_img_non_interactive(g_structure_filename,
					 g_png_mode, g_jpg_mode);
    } else {
      make_structure_ps();
    }
    if (g_ss_output_mode)
      create_ss_output();
  } else {
    return (enable_window(argc, argv));
  }
  /* The program cannot reach this line in glut mode */
  try_exit(0);
  return 0; /* Keep the stupid compiler happy. */
}
