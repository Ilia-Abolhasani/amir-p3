
void open_ct(char *ct_filename) {
  char            long_name[120];
  if ((ctfp = fopen(ct_filename, "r")) == NULL) {
    strcpy(long_name, ct_filename);
    strcat(long_name, ".ct");
    if ((ctfp = fopen(long_name, "r")) == NULL) {
      printf("Error!\tCould not open input files %s or %s\n", 
	     ct_filename, long_name);
      try_exit(36);
    } else
      printf("Input File: %s\n", long_name);
  } else
    printf("Input File: %s\n", ct_filename);
}

void fix_first_line(char *first_line) {
  int             len;
  int             i, j;
  int             found;
  /* remove trailing spaces */
  len = strlen(first_line);
  found = FALSE;
  for (i = (len - 2); (i >= 1) && (found == FALSE); i--) {
    if (first_line[i] != ' ') {
      first_line[i + 1] = '\0';
      found = TRUE;
    }
  }

  /* Replace all control characters (tabs) by spaces and compress
     multiple spaces */ 
  len = strlen(first_line);
  if (first_line[0] < 32)
    first_line[0] = 32;
  for (i = 1; i < len ; i++) {
    if (first_line[i] < 32)
      first_line[i] = 32;
    if (first_line[i]==32 && first_line[i-1]==32) {
      for (j=i ; j < len ; j++ ) {
	first_line[j-1] = first_line[j] ;
      }
      len-- ;
      first_line[len] = '\0' ;
    }
  }

  /* find first non-numeric, non space */
  len = strlen(first_line);
  for (i = 0; i < len; i++) {
    if (isalpha(first_line[i])) {  /* i is first character to save  */
      for (j = i; j <= len; j++) {
	first_line[j - i] = first_line[j];
      }
      return;
    }
  }
}

void remove_junk_gcg(char *record, char *first_line) {
  /* find Check: from end of string */
  int             i, j, length, found;
  int             end, start;
  length = strlen(record);
  i = length - 1;
  found = FALSE;
  /* find end */
  
  while ((i > 12) && (found == FALSE)) {
    if ((record[i] == ':') &&
	(record[i - 1] == 'k') &&
	(record[i - 2] == 'c') &&
	(record[i - 3] == 'e') &&
	(record[i - 4] == 'h') &&
	(record[i - 5] == 'C'))
      found = TRUE;
    i--;
  }
  if (found == TRUE)
    end = i - 5;	/* i is space in front of C */
  else {
    printf("Error, gcg format file does not have CHECK in first line.\n");
    printf("The first line is %s\n", first_line);
    try_exit(37);
  }
  /* find start
   * find first space after first alpha, since MFOLD is first parameter */
  i = 4;
  found = FALSE;
  
  while ((i < end) && (found == FALSE)) {
    if ((isalpha(record[i])) && (!isalpha(record[i + 1])))
      found = TRUE;
    i++;
  }
  if (found == TRUE)
    start = i + 4;	/* i is space in front of of: */
  else {
    printf("Error! gcg format first line is weird\n");
    printf("The first line is %s\n", first_line);
    try_exit(38);
  }
  j = 0;
  for (i = start; i <= end; i++) {
    first_line[j] = record[i];
    j++;
  }
  first_line[j] = '\0';
}

int check_gcg(int *length, char *record, char *first_line) {
  /* first line should start with MFOLD
   * second line has length: (the length) and ends with .. */
  char            mfold_string[90];
  char            energy_string[90];
  char            energy_s_value[90];
  int             i;
  i = sscanf(record, "%s", mfold_string); /* assign length */
  if (i != 1)
    return FALSE;
  for (i = 0; i < strlen(mfold_string); i++)
    mfold_string[i] = toupper(mfold_string[i]);
  if (strcmp(mfold_string, "MFOLD"))
    return FALSE;
  /* does contain MFOLD */
  remove_junk_gcg(record, first_line);
  if (fgets(record, 120, ctfp) == NULL) {
    printf("Error! gcg file contains only first line\n");
    try_exit(39);
  }
  if (strlen(record) < 5)
    if (fgets(record, 120, ctfp) == NULL) {
      printf("Error! gcg  file contains only first line?\n");
      try_exit(40);
    }
  sscanf(record, "%s%d%s%s", mfold_string, length, energy_string, 
	 energy_s_value);
  strcpy(energy_string, "dG= ");
  strcat(energy_string, energy_s_value);
  strcat(energy_string, first_line);
  strcpy(first_line, energy_string);
  printf("GCG format assumed\n");
  return TRUE;
}

void read_ct_file(int *len, char *first_line) {
  char  record[120], base, *bases;
  int   file_length, row, previous, next, column, his, *connected;
  int  *oldconnect, *history, gcg_file_flag, sequence_length, data_items_found;
  /* diag is the diagonal that each row,col corresponds to */
  if (fgets(record, 120, ctfp) == NULL) {
    printf("Error! CT file is completely empty\n");
    try_exit(41);
  } else { /* read again for near empty first line */
    if (strlen(record) < 4) {
      if (fgets(record, 120, ctfp) == NULL) {
	printf("Error!  CT file is completely empty\n");
	try_exit(42);
      }
    }
    gcg_file_flag = check_gcg(&file_length, record, first_line);
    if (!gcg_file_flag) {
      sscanf(record, "%d", &file_length); /* assign length */
    }
    if (file_length > 0) {
      memory_allocator(file_length);
      bases = g_bases;
      connected = g_connected;
      oldconnect = g_oldconnect;
      history = g_history;
    } else {
      printf("The first record of ct file does not include length\n");
      printf("The first record is:\t%s\n", record);
      try_exit(43);
    }
    if (!gcg_file_flag) {
      strcpy(first_line, record);
      fix_first_line(first_line);
    }
  }
  /* Read the ct file */
  sequence_length = 0;
  while (fgets(record, 120, ctfp) != NULL) {
    data_items_found = sscanf(record, "%d %c %d %d %d %d", &row, &base, 
			      &previous, &next, &column, &his);
    if (data_items_found == 6) { /* skip blank or odd lines */
      history[row] = his;
      if (row != (sequence_length + 1)) {
	printf("Error in ct file, Format is wrong\n");
	printf("Record %d of the file contains %s, ",
	       sequence_length + 1, record);
	printf("but it should start with %d\n", sequence_length + 1);
	try_exit(44);
      }
      g_next[row] = next;     /* 2003-06-21 Nick Markham */
      g_prev[row] = previous; /* 2003-06-21 Nick Markham */
      sequence_length++;
      /*	Zuker comments out: May 22,2003. */
      /*	bases[sequence_length] = toupper(base);*/
      bases[sequence_length] = base;
      if (column == 0) {
	connected[sequence_length] = -1;
	oldconnect[sequence_length] = -1;
      } else {
	connected[sequence_length] = column;
	oldconnect[sequence_length] = column;
	/* if (column<row) insert block:
	 * November 25, 2006 insert by M. Zuker
	 * 1x1 and 2x2 interior loops are converted to single or tandem
	 * mismatches, respectively. 
	 * Enabled with -bp flag (g_bp is TRUE)
	 */
	if (column < row && g_bp) {
	  if (connected[column] != sequence_length) {
	    printf("Error!\tBbase %d is paired with %d, but", 
		   sequence_length, column);
	    printf(" base %d is paired with %d\n", column, connected[column]);
	  }
	  if (connected[sequence_length-1] == -1 ) {
	    if (connected[sequence_length-2] == column+2) {
	      connected[sequence_length-1] = column + 1;
	      connected[column + 1] = sequence_length-1;
	    } 
	    else if (connected[sequence_length-2] == -1 && 
		     connected[sequence_length-3] == column+3) {
	      connected[sequence_length-1] = column + 1;
	      connected[sequence_length-2] = column + 2;
	      connected[column+1] = sequence_length-1;
	      connected[column+2] = sequence_length-2;
	    }
	  }
	}
      }
    }
  }
  /* end of while loop */
  /* M. Zuker, May 19, 2007.
   * Fix bug. If base 1 pairs with the last base, then they are
   * displayed as unpaired when the -f or -fa flag is used. Check for
   * this and set the global flat variables to FALSE
   */
  if (connected[1]==sequence_length) {
    g_flat = FALSE;
    g_flat_alternate = FALSE;
  }
  /* M. Zuker, June 7, 2007.
   * if prev[1] equals sequences_length and next[sequence_length]
   * equals 1, then the molecule is circular. Set g_is_circular to
   * TRUE and make sure neither flat option is used.
   * On the other hand, if the molecular is circular because the -c
   * flag is used on the command line, make sure previous and next are
   * set correctly.
   */
  if (g_prev[1]==sequence_length && g_next[sequence_length]==1) {
    g_is_circular = TRUE;
    g_flat = FALSE;
    g_flat_alternate = FALSE;
  } else if (g_is_circular) {
    g_prev[1] = sequence_length;
    g_next[sequence_length] = 1;
    g_flat = FALSE;
    g_flat_alternate = FALSE;
  }    
  if (file_length != sequence_length) {
    printf("Error! There are %d proper records, but the first line\n",
	   sequence_length);
    printf("indicates that there should be %d records.\n",
	   file_length);
    try_exit(45);
  }
  *len = sequence_length;
}

void traverse_ss_loop_set_code(int start_base, int end_base, int gl_start_base,
			       int *connected, int *ss_code) {
  int             current_base;
  /* find extreme right, left, top,bottom values 
   * place loop center at middle of these */
  if (start_base == gl_start_base) {
    if (connected[start_base] > start_base) {
      current_base = start_base;
      if (connected[current_base + 1] == (connected[current_base] - 1)) {
	ss_code[current_base] = 1;
	ss_code[connected[current_base]] = 4;
      } else {/* single base pairs are assigned 12 */
	ss_code[current_base] = 12;
	ss_code[connected[current_base]] = 12;
      }
      current_base = connected[start_base];
    } else {
      current_base = start_base + 1;
    }
  } else
    current_base = start_base + 1;
  while (current_base < end_base) {
    if (connected[current_base] > current_base) {
      if (connected[current_base + 1] == (connected[current_base] - 1)) {
	ss_code[current_base] = 1;
	ss_code[connected[current_base]] = 4;
      } else {/* single base pairs are assigned 12 */
	ss_code[current_base] = 12;
	ss_code[connected[current_base]] = 12;
      }
      current_base = connected[current_base];
      ss_code[current_base] = 4;
    } else {
      current_base++;
    }
  }
}

void set_ss_code_starting_single(int start_base, int *connected,
				 int *ss_code, int length) {
  int             current_base;
  if (connected[start_base] > start_base)
    return;
  ss_code[start_base] = 8;
  current_base = start_base;
  while (current_base <= length) {
    current_base++;
    if (connected[current_base] > current_base) {
      ss_code[current_base - 1] = 9;
      return;
    }
  }
}

void set_ss_code_ending_single(int start_base, int end_base, int *connected,
			       int *ss_code) {
  int             current_base;
  if (connected[end_base] > 0)
    return;
  ss_code[end_base] = 9;
  current_base = end_base;
  while (current_base >= start_base) {
    current_base--;
    if (connected[current_base] > 0) {
      ss_code[current_base + 1] = 8;
      return;
    }
  }
}

void  set_ss_out_codes(int start_base, int end_base, int *ss_code,
		       int *connected, int length) {
  int             current_base;
  
  for (current_base = start_base; current_base <= end_base;
       current_base++) {
    ss_code[current_base] = 0;
    
  }
  set_ss_code_starting_single(start_base, connected, ss_code, length);
  traverse_ss_loop_set_code(start_base, end_base, start_base, connected,
			    ss_code);
  current_base = start_base + 1;
  connected[0] = 0;
  connected[length + 1] = 0;
  while (current_base <= end_base) {	/* If find last helix of base
					 * pair, traverse loop */
    if ((connected[current_base] > current_base) &&
	(connected[current_base + 1] != (connected[current_base] - 1))) {
      ss_code[current_base] = 2;
      ss_code[connected[current_base]] = 3;
      traverse_ss_loop_set_code(current_base, connected[current_base], 
				start_base, connected, ss_code);
    }
    current_base++;
  }
  set_ss_code_ending_single(start_base, end_base, connected, ss_code);
}

void general_read_ct_file(char *ct_filename, int *length, char *first_line) {
  int         *ss_code;
  int         *connected;
  int          tit_len; /* M. Zuker, June 27, 2006. Use figure title from the
			   command line if it exists */
  open_ct(ct_filename);
  read_ct_file(length, first_line);
  printf("Sequence length: %d\n", *length);
  /* Added by M. Zuker, June 27, 2006. */ 
  tit_len = strlen(g_fig_title);
  if (tit_len > 0) {
    strcpy(first_line, g_fig_title);
    printf("Sequence title replaced from command line\n");
  }
  /* --- */
  fclose(ctfp);
  ss_code = g_ss_code;
  connected = g_connected;
  set_ss_out_codes(1, *length, ss_code, connected, *length);
}
