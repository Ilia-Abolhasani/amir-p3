/* Include file for sir_graph
 * Read ss files and create ss output
 */

int open_out_ss(char *filename) {
  if ((outssfp = fopen(filename, "w")) == NULL) {
    printf("Error!\nCould not open output file: %s \n", filename);
    return TRUE;
  }
  return FALSE;
}

void create_ss_output_real(char *filename, int length, int *connected, 
			   char *bases, float *base_x, float *base_y,
			   int *ss_code) {
  int             i;
  int             connected_to;
  float ssx=0.0, ssy=0.0, lastx=0.0, lasty=0.0; /* M. Zuker Dec 18, 2006. */
  int sscode = 0; /* M. Zuker Dec 18, 2006. */
  char base = '\1'; /* M. Zuker Dec 18, 2006. 
		     * ^A is used to denote a base that will appear as a blank.
		     * The null character is already used to denote 5' and 3'
		     */
  chop_suffix(filename, ".ss");
  chop_suffix(filename, ".ps");
#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
  chop_suffix(filename, ".img");
#endif
  strcat(filename, ".ss");
  printf("ss file output is %s\n", filename);
  if (open_out_ss(filename))
    return;
  /* M. Zuker Dec 18, 2006.
   * Find first x,y pair that is not 0,0
   */
  i = 1;
  while (i <= length && base_x[i]==0.0 && base_y[i]==0.0) i++;
  if (i<length) {
    lastx = base_x[i];
    lasty = base_y[i];
  }
  for (i = 1; i <= length; i++) {
    if (connected[i] > 0)
      connected_to = connected[i];
    else 
      connected_to = 0;
    /* M. Zuker Dec 18, 2006.
     * When -iw or -ew are used, x and y are both zero. The resulting .ss file
     * cannot be read by sir_graph(_ng) unless columns 5 and 6 are
     * also set to 0 
     */
    if (base_x[i]==0.0 && base_y[i]==0.0) {
      base = '\1';
      ssx = lastx;
      ssy = lasty;
      sscode = 0;
      connected_to = 0;
    } else {
      base = bases[i];
      ssx = base_x[i];
      ssy = base_y[i];
      lastx = ssx;
      lasty = ssy;
      sscode = ss_code[i];
    }
    fprintf(outssfp, "%d\t%c\t%.2f\t%.2f\t%d\t%d\n", i, base, ssx, ssy, 
	    sscode, connected_to);
  }
  fclose(outssfp);
}

/* input of ss file */
  
void open_ss(char *filename) {
  if ((ssfp = fopen(filename, "r")) == NULL) {
    printf("Error!\tCould not open input file: %s\n", filename);
    try_exit(47);
  }
}

/* returns TRUE when the nth line does not begin with n */
int process_ss_record(char *record, int length, int *connected, int *ss_code,
		      float *base_x,float *base_y, char *bases, int *history) {
  int             row, connected_to, code;
  char            base;
  float           x, y;
  sscanf(record,"%d %c %f %f %d %d", &row,&base,&x,&y,&code,&connected_to);
  history[row] = row;
  if (length != row) {
    printf("Extra information in ss file lines from %d onward ignored\n",
	   length);
    printf("First ignored line is %s\n", record);
    return TRUE;
  }
  if (connected_to == 0)
    connected[row] = -1;
  else if (connected_to > row) {
    connected[row] = connected_to;
    connected[connected_to] = row;
  } 
  else if (connected[connected_to] != row && connected[connected_to] > 0) {
    printf("Error!\tBase %d pairs with %d,", row, connected_to);
    printf(" but base %d pairs with %d\n", connected[connected_to],
	   connected_to);
    try_exit(48);
  }
  ss_code[length] = code;
  base_x[length] = x;
  base_y[length] = y;
  bases[length] = base;
  return FALSE;
}

void adjust_start_and_end(int *, int *, int, int *);

void set_ss_distance_basepair(int length, int *connected,
			      float *ss_distance_basepair, 
			      float *ss_distance_base_adj, float *base_x, 
			      float *base_y, int *ss_backwards, int *ss_code) {
  int             i;
  for (i = 2; i <= length; i++) {
    if (connected[i] > 0) {
      *ss_distance_basepair =
	sqrt((base_x[i] - base_x[connected[i]]) *
	     (base_x[i] - base_x[connected[i]]) +
	     (base_y[i] - base_y[connected[i]]) *
	     (base_y[i] - base_y[connected[i]]));
      if (connected[i + 1] > 0)
	*ss_distance_base_adj = sqrt((base_x[i] - base_x[i + 1]) *
				     (base_x[i] - base_x[i + 1]) +
				     (base_y[i] - base_y[i + 1]) *
				     (base_y[i] - base_y[i + 1]));
      else
	*ss_distance_base_adj = *ss_distance_basepair;
      printf("Base pair distance = %.2f, base adj distance = %.2f\n",
	     *ss_distance_basepair, *ss_distance_base_adj);
      if (*ss_distance_base_adj < (*ss_distance_basepair / 4.0))
	*ss_distance_base_adj = *ss_distance_basepair / 4.0;
      if (ss_code[i] < 0)
	*ss_backwards = TRUE;
      return;
    }
  }
  *ss_distance_basepair = sqrt((base_x[1] - base_x[connected[2]]) *
			       (base_x[1] - base_x[connected[2]]) +
			       (base_y[1] - base_y[connected[2]]) *
			       (base_y[1] - base_y[connected[2]]));
  *ss_distance_base_adj = *ss_distance_basepair;
  return;
}

void set_limits(int current_base, float *left_x, float *right_x,
		float *top_y, float *bottom_y, float *base_x, float *base_y) {
  float           x, y;
  x = base_x[current_base];
  y = base_y[current_base];
  if (x < *left_x)
    *left_x = x;
  else if (x > *right_x)
    *right_x = x;
  if (y < *bottom_y)
    *bottom_y = y;
  else if (y > *top_y)
    *top_y = y;
}

void traverse_ss_loop(int start_base, int end_base, int current_loop,
		      float *base_x, float *base_y, int gl_start_base, 
		      float *loop_center_x, float *loop_center_y,
		      int *connected) {
  float           left_x, right_x, top_y, bottom_y;
  int             current_base;
  /* find extreme right, left, top,bottom values
   * place loop center at middle of these */
  left_x = base_x[start_base];
  right_x = base_x[start_base];
  top_y = base_y[start_base];
  bottom_y = base_y[start_base];
  if (start_base == gl_start_base) {
    if (connected[start_base] > start_base)
      current_base = connected[start_base];
    else
      current_base = start_base + 1;
  } else
    current_base = start_base + 1;
  set_limits(current_base,&left_x,&right_x,&top_y,&bottom_y, base_x, base_y);
  while (current_base < end_base) {
    if (connected[current_base] > current_base)
      current_base = connected[current_base];
    else
      current_base++;
    set_limits(current_base,&left_x,&right_x,&top_y,&bottom_y,base_x,base_y);
  }
  loop_center_x[current_loop] = (left_x + right_x) / 2.;
  loop_center_y[current_loop] = (top_y + bottom_y) / 2.;
}

void set_ss_loops(int start_base, int end_base, int *connected, int length,
		  int *total_loops, float *base_x, float *base_y,
		  float *loop_center_x, float *loop_center_y, 
		  int *previous_loop) {
  int             current_loop;
  int             current_base;
  current_base = start_base;
  current_loop = 1;
  traverse_ss_loop(start_base, end_base, current_loop, base_x, base_y,
		   start_base, loop_center_x, loop_center_y, connected);
  current_base = start_base + 1;
  connected[0] = 0;
  connected[length + 1] = 0;
  while (current_base <= end_base) {	/* If find last helix of base
					 * pair, traverse loop */
    previous_loop[current_base] = current_loop;
    if ((connected[current_base] > current_base) &&
	(connected[current_base + 1] != (connected[current_base] - 1))) {
      current_loop++;
      traverse_ss_loop(current_base, connected[current_base],
		       current_loop, base_x, base_y,
		       start_base, loop_center_x, loop_center_y, connected);
    }
    current_base++;
  }
  *total_loops = current_loop;
}

void read_ss_file(char *first_line, int *length, char *ss_filename, 
		  int *connected, int *ss_code, float *base_x, float *base_y, 
		  char *bases,int *start,int *end,float *ss_distance_basepair,
		  float *ss_distance_base_adj, int *ss_backwards, 
		  float *loop_center_x, float *loop_center_y, 
		  int *total_loops, int *previous_loop, int *history) {
  char            record[120];
  int             sequence_length;
  int             hit_extra_junk;
  int             tit_len; /* M. Zuker, June 27, 2006. */
  
  if (fgets(record, 120, ssfp) == NULL) {
    printf("Error!\tss file is completely empty\n");
    try_exit(49);
  } else {
    strcpy(first_line, "ss file: ");
    strcat(first_line, ss_filename);
  }
  tit_len = strlen(g_fig_title);
  if (tit_len > 0) { strcpy(first_line, g_fig_title); }
  sequence_length = 1;
  hit_extra_junk = process_ss_record(record, sequence_length, connected,
				     ss_code, base_x, base_y, bases, history);
  if (hit_extra_junk) {
    printf("Did not find any valid lines in file\n");
    try_exit(50);
  }
  while ((fgets(record, 120, ssfp) != NULL) && (!hit_extra_junk)) {
    sequence_length++;
    hit_extra_junk = 
      process_ss_record(record, sequence_length, connected,
			ss_code, base_x, base_y, bases, history);
    if (hit_extra_junk)
      sequence_length--;
  }
  *length = sequence_length;
  /* M. Zuker, May 19, 2007. Fix bug. When base 1 pairs with the last
  base, display with either -f or -fa causes these base to be drawn as
  unpaired. The simple fix is to set the global flat variables to
  false. 
  */
  if (connected[1]==sequence_length) {
    g_flat = FALSE;
    g_flat_alternate = FALSE;
  }
  adjust_start_and_end(start, end, *length, connected);
  set_ss_distance_basepair(*length, connected, ss_distance_basepair, 
			   ss_distance_base_adj, base_x, base_y, ss_backwards,
			   ss_code);
  set_ss_loops(*start, *end, connected, *length, total_loops, base_x, base_y, 
	       loop_center_x, loop_center_y, previous_loop);
}

void general_read_ss_file(char *ss_filename, char *first_line, int *length, 
			  int *start, int *end, float *ss_distance_basepair, 
			  float *ss_distance_base_adj,
			  int *ss_backwards, int *total_loops) {
  int             i;
  int             found_length;
  int            *connected;
  int            *ss_code;
  float          *base_x;
  float          *base_y;
  char           *bases;
  float          *loop_center_x;
  float          *loop_center_y;
  int            *previous_loop;
  int            *history;
  char            record[120];
  open_ss(ss_filename);
  found_length = 0;
  while (fgets(record, 120, ssfp) != NULL) {
    found_length++;
  }
  fclose(ssfp);
  printf("Lines in ss file: %d\n", found_length);
  /* memory is allocated assuming no extra lines at end of ss file
   * This could be done better */
  memory_allocator(found_length);
  connected = g_oldconnect;
  connected[0] = 0;                         /* 2003-06-24 Nick Markham */
  ss_code = g_ss_code;
  base_x = g_base_x;
  base_y = g_base_y;
  bases = g_bases;
  loop_center_x = g_loop_center_x;
  loop_center_y = g_loop_center_y;
  previous_loop = g_previous_loop;
  history = g_history;
  open_ss(ss_filename);
  read_ss_file(first_line,length,ss_filename,connected,ss_code,base_x,base_y,
	       bases,start,end,ss_distance_basepair,ss_distance_base_adj,
	       ss_backwards,loop_center_x,loop_center_y,total_loops, 
	       previous_loop, history);
  printf("Sequence length: %d\n", *length);
  for (i = 1; i < *length; ++i) {      /* 2003-06-24 Nick Markham */
    g_prev[i] = i - 1;         /* 2003-06-24 Nick Markham */
    g_next[i] = i + 1;         /* 2003-06-24 Nick Markham */
    g_connected[i] = g_oldconnect[i]; /* 2007-05-23 M. Zuker */
  }                                 /* 2003-06-24 Nick Markham */
  /* M. Zuker, June 7, 2007.
   * if the molecule is circular because the -c flag is used on the
   * command line, make sure previous and next are set correctly.
   * This adds to what Nick Markham did on 2003-06-24
   * In addition, make sure flat is turned off in such cases, although
   * I'm not sure it is really needed for ss file input.
   */
  if (g_is_circular) {
    g_prev[1] = *length;
    g_next[*length] = 1;
    g_flat = FALSE;
    g_flat_alternate = FALSE;
  } else {
    g_prev[*length] = *length - 1; /* 2003-06-24 Nick Markham */
    g_next[*length] = 0;           /* 2003-06-24 Nick Markham */
  }
  for (i = 1; i <= *total_loops; ++i) /* 2003-06-24 Nick Markham */
    g_loop_radius[i] = 1;        /* 2003-06-24 Nick Markham */
  
  if (ctfp) /* 2003-06-20 Nick Markham */
    fclose(ctfp);
}
