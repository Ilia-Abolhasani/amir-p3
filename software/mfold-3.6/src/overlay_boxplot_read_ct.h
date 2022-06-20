
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int *g_ct_basepair;

char sequence[MAXIMUM_SIZE];

int zuker_distance(int row1,int col1,int row2,int col2) {
  int hor_dis;
  int ver_dis;
  hor_dis=(int)(abs(row1-row2));
  ver_dis=(int)(abs(col1-col2));
  if(hor_dis<ver_dis)
    return ver_dis;
  else
    return hor_dis;
}

int distance_plot_ct(int ct_len,int plot_row,int plot_col,int max_look) {
  /* return true if a ct dot is within max_look distance of plot_row,col
   * max look is 0 to 5
   * check 0 */
  int check_row;
  int distance;
  /* check rows at or below the plot dot */
  for(check_row=plot_row;check_row<=(plot_row+max_look);check_row++) {
    if(check_row<=ct_len) {
      distance = 
	zuker_distance( plot_row,plot_col,check_row,g_ct_basepair[check_row]);
      if(distance<=max_look) {
	return TRUE;
      }
    }
  }
  /* check rows above the plot dot */
  for(check_row=plot_row-max_look;check_row<plot_row;check_row++) {
    if(check_row>0) {
      distance = 
	zuker_distance(plot_row,plot_col,check_row,g_ct_basepair[check_row]);
      if(distance<=max_look) {
	return TRUE;
      }
    }
  }
  return FALSE;    
}

void ct_print_helix(int number_of_helices,struct ct_helix *diag) {
  int i;
  for(i=0;i<number_of_helices;i++) {
    printf("diag=%d, row=%d, energy=%d, column=%d, length=%d color=%d\n",
	   diag[i].diagonal, diag[i].row, diag[i].energy, diag[i].column,
	   diag[i].length, diag[i].color);
  }
}

void fix_name(char *sequence_name) {
  /* remove new line and other annoying characters 
   * Jan 24, 2007. M. Zuker removes tabs as well.
   */
  int i;
  for(i=0;i<strlen(sequence_name);i++) {
    if(sequence_name[i]=='\n' || sequence_name[i]=='\t')
      sequence_name[i]=' ';
  }
}

int ct_compare_helices(const void *h1,const void *h2) {
  int dif;
  struct ct_helix *helix1,*helix2;
  helix1=(struct ct_helix *)h1;
  helix2=(struct ct_helix *)h2; 
  dif=(*helix1).diagonal-(*helix2).diagonal; /* return <0 if diag1 <
					      * diag2 on same diagonal
					      */ 
  if(dif==0) 
    dif=(*helix1).row-(*helix2).row; /* return <0 if row1<row2 or >0 if > */
  return dif;
}   

void ct_sort_helices(struct ct_helix *diag,int helices_in_ct_file) {
  qsort(diag,helices_in_ct_file, sizeof(struct ct_helix),ct_compare_helices);
}

void append_structure_name(char *data_name, char *file_name) {
  /* place name from file at start of description, data_name */
  char temp_name[120];
  int i, k, length;
  char number[12];
  strcpy(temp_name,file_name);
  chop_suffix(temp_name, ".ct"); 
  length = strlen(temp_name) ;
  /* look for _#, _##, _### etc. at end of name */
  for (k=length-1; isdigit(temp_name[k]) && k > 0; k--) { } 
  if ( (k < length - 1) && (temp_name[k]=='_') ) {
    number[0] = ' ';
    number[1] = '[';
    for (i=k+1; i < length && i <= 8+k ; i++)
      number[i+1-k] = temp_name[i] ;
    number[length-k+1] = ']';
    number[length-k+2] = '\0';
    strcat(data_name,number);
  }
}

void read_ct_file(FILE *fp, int file, int *len, struct ct_helix *diag,
		  char file_data_ct[MAXIMUM_CT_FILES][80],
		  int *helices_in_ct_file, int *first_sequence_length,
		  char *g_sequence_name, int sequence_name_set) {
  char record[80];
  int file_length;
  int row,previous,next,column,junk;
  int helix_len;
  int helices_in_this_file; 
  char base;
  char energy_string[80];
  int first_char;
  int helix_started;
  int j,i;
  int sequence_length;
  int previous_column;
  int name_start;
  int energy_start;
  int record_length;
  char junk_length[80];
  char sequence_name[120];
  int previous_energy,energy;
  int color = COLOR_TEXT ;
  int plotlen = 0;
  sequence_length=0;
  helices_in_this_file=*helices_in_ct_file;
  /* diag is the diagonal that each row,col corresponds to */
  if(fgets(record,90,fp)==NULL) {
    printf("Error!\tThe ct file is completely empty\n");
    try_exit(13);
  } else {
    /* Record should look like " 741 dG =  -238.3    mok-lacZ + gtgtaa
     * ", for example, where "741" is the sequence length and "mok-lacZ
     * + gtgtaa+" (i.e. all the rest) is the name (or title).
     * The name should be be the same in all files on a given run
     * The energy varies and is stored as file_data_ct
     */
    sscanf(record,"%d",&file_length); /* assign length */
    /* allocate ct_basepairs */
    g_ct_basepair=(int *)malloc((file_length+1)*sizeof(int));
    if(g_ct_basepair==NULL) {
      printf("Error!\tCannot allocate space for ct_basepair array.\n");
      try_exit(14);
    }           
    sscanf(record,"%s",junk_length); /* assign length string */
    first_char = 0;
    name_start = 0;
    energy_start = 0;
    record_length = strlen(record);
    while((isspace(record[first_char]))&&(first_char<record_length))
      first_char++;/* Advance through spaces of record */
    while(isdigit(record[first_char])) /*advance through numbers of record */
      first_char++;
    while((isspace(record[first_char]))&&(first_char<record_length))
      first_char++; /* advance through spaces */
    energy_start=first_char;
    while((!isspace(record[first_char]))&&(first_char<record_length))
      first_char++; /* advance through energy */
    while((isspace(record[first_char]))&&(first_char<record_length))
      first_char++; /* advance through spaces */
    while((!isspace(record[first_char]))&&(first_char<record_length))
      first_char++; /* advance through = */
    while((isspace(record[first_char]))&&(first_char<record_length))
      first_char++; /* advance through spaces */
    while((!isspace(record[first_char]))&&(first_char<record_length))
      first_char++; /* advance through = number*/
    name_start = first_char;
    /* this should copy energy into file_data_ct[file]*/
    strcpy(file_data_ct[file]," ");
    /* if the filename ends with _# append that number to file_data_ct */
     if(name_start>30) {
       name_start=6;
       strcpy(energy_string,"Energy is Undef.");
     } else  {
       for(i=energy_start,j=0;i<name_start;i++,j++) 
	 energy_string[j]=record[i];
       energy_string[j]='\0';
     }
     strcat(file_data_ct[file]," ");
     strcat(file_data_ct[file],energy_string);
     if((file==0)&&(!sequence_name_set)) {
       for(i=name_start,j=0;i<record_length;i++,j++) 
	 sequence_name[j]=record[i];
       sequence_name[j]='\0';
       fix_name(sequence_name);
       strcpy(g_sequence_name,sequence_name);
     }
  }
  /* Read the ct file */
  helix_started = FALSE;
  helix_len = 0;
  previous_column = -10;
  sequence_length = 0;
  previous_energy = INT_MAX; 
  while(fgets(record,90,fp)!=NULL) {
    sscanf(record,"%d %c %d %d %d %d", &row, &base, &previous, &next, &column,
	   &junk);
    if (row!=(sequence_length+1)) {
      printf("Error in ct file. Record %d of the file contains\n%s,\n",
	     sequence_length+1,record);
      printf("but it should start with %d\n",sequence_length+1);
      try_exit(15);
    }
    if(file==0)
      sequence[sequence_length]=base;
    else {
      if(sequence_length>*first_sequence_length) {
	printf("Error!\tThe sequences are not the same length.");
	printf("Lengths are %d and at least %d", *first_sequence_length,
	       sequence_length);
        try_exit(16);
      }
      if(base!=sequence[sequence_length]) {
	printf("Warning!\tAt position %d, sequence %s has a %c,\n",
	       sequence_length+1, file_data_ct[0], sequence[sequence_length]);
	printf("but sequence %s has a %c.\n", file_data_ct[file],base);
	printf("Program continues anyway.\n");
      }
    }
    g_ct_basepair[row]=column;
    /* printf("\n row is %d column is %d",row,column);*/
    sequence_length++;
    if(helix_started==TRUE) {
      if(column==0) { /* finish off helix */  
	diag[*helices_in_ct_file].length=helix_len;
	helix_started = FALSE;
	(*helices_in_ct_file)++;
      } else
	if(column<row) { /* col<row has already been added */
	  /* finish off helix. all helices appear twice */
	  diag[*helices_in_ct_file].length = helix_len;
	  helix_started=FALSE; /* just finish current helix */
	  (*helices_in_ct_file)++;
	} else {
	  if(column!=0) {
	    energy = energy_read(row,column);
	    plotlen = chain_length_read(row,column);
	    if(energy==INT_MIN) {
	      color=COLOR_NOT_OPTIMAL;
	    } else 
	      color=find_color(energy);
	  }  else
	    energy=INT_MAX;
	  if((column!=(previous_column-1))||(energy!=previous_energy)) {
	    /* finish off helix */
	    diag[*helices_in_ct_file].length=helix_len;
	    (*helices_in_ct_file)++;
	    helix_started=FALSE;
	    /* start new helix */
	    previous_column=column;
	    helix_len=1;
	    if(*len<column)
	      *len = column;
	    diag[*helices_in_ct_file].row=row;
	    diag[*helices_in_ct_file].column=column;
	    diag[*helices_in_ct_file].diagonal=row+column-1;
	    previous_energy=energy;
	    diag[*helices_in_ct_file].energy=energy;
	    diag[*helices_in_ct_file].plotlen=plotlen;
	    diag[*helices_in_ct_file].color=color;
	    helix_started=TRUE;
	  } else {
	    helix_len++; /* add one to length */
	    previous_column=column;
	  }
	}
    } else { /* helix has not been started */
      if(column>row) { /* start a helix */
        previous_column = column;
	helix_len = 1;
	diag[*helices_in_ct_file].row=row;
	diag[*helices_in_ct_file].column=column;
	diag[*helices_in_ct_file].diagonal=row+column-1;
	energy=energy_read(row,column);
	plotlen=chain_length_read(row,column);
	if(energy==INT_MIN) {
	  color=COLOR_NOT_OPTIMAL;
	} else
	  color=find_color(energy);
	diag[*helices_in_ct_file].energy=energy;
	diag[*helices_in_ct_file].plotlen=plotlen;
	diag[*helices_in_ct_file].color=color;
	previous_energy=energy;
	helix_started=TRUE;
	if(*len<column)
	  *len=column;
      }
    }
  }
  /* if a helix was started, record its length */
  if(helix_started==TRUE) {
    diag[*helices_in_ct_file].length=helix_len;
    diag[*helices_in_ct_file].diagonal=row+column-1;
    (*helices_in_ct_file)++;
  }
  if(*helices_in_ct_file>=(MAXIMUM_HELICES-1)) {
    printf("Error!\n1.\tThere are more than %d helices in the plot file.\n",
	   MAXIMUM_HELICES);
    printf("Increase the constant MAXIMUM_HELICES in %s.c", PROGRAM_NAME);
    printf(" and recompile.\n");
    try_exit(17);
  }
  helices_in_this_file=*helices_in_ct_file-helices_in_this_file;
  printf("There are %d helices in this ct file with length %d.\n",
	 helices_in_this_file,sequence_length);
  if(file==0)
    *first_sequence_length=sequence_length;
  else {
    if(sequence_length!=*first_sequence_length) {
      printf("Error!\tThe sequences have different lengths.\n");
      printf("The first sequence length is %d, the current one is %d.\n",
	     *first_sequence_length,sequence_length);
      try_exit(18);
    }
  }
  if(file_length!=sequence_length) {
    printf("Error!\tThere are %d records, but the first line ",
	   sequence_length);
    printf("indicates that there should be %d records.\n",
	   file_length);
    try_exit(19);
  }
  /* set the base pairs to -10 to avoid distance problems with 0 */
  for(i=1;i<=file_length;i++) {
    if(g_ct_basepair[i]==0)
      g_ct_basepair[i]=-10;
  }
}

void ct_fix_counts(int helices_in_ct_file, int *diag_count, 
		   struct ct_helix *diag) {
  int i;
  for(i=0;i<helices_in_ct_file;i++) {
    diag_count[diag[i].diagonal]++;
    /* update the count for this diagonal */
  }
} 

void sort_all_helices_set_start(int *diag_count, int *diag_start,
				int helices_in_ct_file,
				struct ct_helix *diag,int len) {
  int i;
 if(helices_in_ct_file==0) {
   printf("Error!\tNo helices have bee found.\n");
   try_exit(20);
 }
 /* diag_count counts the total number of entries in each diagonal.
  * Sort the plot file */
  ct_sort_helices(diag,helices_in_ct_file);
  for(i=0;i<=(2*len-1);i++)
    diag_count[i]=0;
  /* fix the count of diagonal for each entry in plot file */
  diag_start[1]=0;
  ct_fix_counts(helices_in_ct_file,diag_count,diag);
  for(i=2;i<=(2*len-1);i++) {
    diag_start[i]=diag_start[i-1]+diag_count[i-1];
  } 
  /* diag_start counts forward into each array for the start of each
     diagonal */ 
}

/* add overlaps for an individual diagonal named current_diagonal */
void add_overlap_diag(int *diagonal_fill, int current_diagonal, 
		      int *diag_start, int *diag_count, struct ct_helix *diag,
		      int number_of_ct_files) {
  /* diagonal_fill is a diagonal line that is filled to indicate how
     many dots at each position diagonal_count will be updated */
  int count,i,starting_position,farthest_position;
  int row,length,current_row;
  int closest_position;
  count=diag_count[current_diagonal];
  if(number_of_ct_files>1) {
    if(count<2)
      return;
  } else {
    if(count<1)
      return;
  }
  /* add each entry along the diagonal */
  starting_position=diag_start[current_diagonal];
  /* process each helix on the diagonal */
  farthest_position=0;
  closest_position=MAXIMUM_HELICES;
  for(i=0;i<count;i++) {
    row = diag[starting_position+i].row;
    length = diag[starting_position+i].length;
    /* process a single helix */
    if(farthest_position<(row+length-1))
      farthest_position=row+length-1;
    if(closest_position>row)
      closest_position=row;
    for(current_row=row;current_row<=(row+length-1);current_row++) {
      diagonal_fill[current_row]++;
    }
  }
}

/* add overlaps for all diagonals */
void  add_overlaps(int *diag_count, int *diag_start, int *helices_in_ct_file,
                   struct ct_helix *diag, int len, int number_of_ct_files) {
  int new_helices_in_ct_file;
  int i;
  int *diagonal_fill;
  diagonal_fill = (int *)malloc(len*sizeof(int));
  /* The length of longest diagonal, this was too big. Should probably
   * use square_root of 2 * len/2+1. Then clear it.
   */ 
  if(diagonal_fill==NULL) {
    printf("Error!\tInsufficient memory to allocate diagonal_fill.\n");
    try_exit(21);
  }
  for(i=0;i<len;i++)
    diagonal_fill[i]=0;
  /* Save old information */
  new_helices_in_ct_file = *helices_in_ct_file;
  /* add overlap helices */ 
  for(i=1;i<=(2*len-1);i++) {
    add_overlap_diag(diagonal_fill, i, diag_start, diag_count, diag,
		     number_of_ct_files);
  }
  /* update for new starts */
  diag_start[1]=0;
  for(i=2;i<=(2*len-1);i++)
    diag_start[i]=diag_start[i-1]+diag_count[i-1];
  *helices_in_ct_file=new_helices_in_ct_file;
  free(diagonal_fill);
}

void read_of_all_ct_files(FILE *fp, int *diag_count, int *diag_start,
			  int number_of_ct_files, int *helices_in_ct_file,
			  int *len, char filename_ct[MAXIMUM_CT_FILES][80],
			  struct ct_helix *diag, 
			  char file_data_ct[MAXIMUM_CT_FILES][80],
			  char *sequence_name, int sequence_name_set) {
  char filename[80];
  int first_sequence_len;
  int i,successful_opens;
  successful_opens=0;  
  for(i=0;i<number_of_ct_files;i++) {
    strcpy(filename,filename_ct[i]); 
    printf("Trying to open  %s.\n",filename);   
    if((fp=fopen(filename,"r"))==NULL) {
      printf("Error!\tCould not open file %s.\n", filename);
      strcpy(file_data_ct[i],"File not found: ");
      strcat(file_data_ct[i],filename);
    } else {
      successful_opens++;
      read_ct_file(fp, i, len, diag, file_data_ct, helices_in_ct_file,
		   &first_sequence_len, sequence_name, sequence_name_set);
      fclose(fp);
    }
  }
  if((successful_opens==0) || (*len==0)) {
    printf("Error!\tUnable to open the ct file successfully.\n");
    try_exit(22);
  }  
 printf("The farthest base pair in the ct file is %d.\n",*len);
 sort_all_helices_set_start(diag_count, diag_start, *helices_in_ct_file,
			    diag, *len); 
 /* The rest of this function could be removed ? */
 add_overlaps(diag_count, diag_start, helices_in_ct_file, diag, *len,
	      number_of_ct_files);
 ct_sort_helices(diag, *helices_in_ct_file);
}
