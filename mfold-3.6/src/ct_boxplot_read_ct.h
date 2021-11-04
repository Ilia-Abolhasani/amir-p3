/* Insert file for ct_boxplot */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char sequence[MAXIMUM_SIZE];

void print_helix(int number_of_helices,struct helix *diag) {
  int i;
  for(i=0;i<number_of_helices;i++) {
    printf("diag=%d, row=%d, color=%d, column=%d, length=%d\n",
	   diag[i].diagonal,
	   diag[i].row,
	   diag[i].color,
	   diag[i].column,
	   diag[i].length);
  }
}

void fix_name(char *sequence_name) {
 /* remove new line and other annoying characters 
  * Jan 24, 2007. M. Zuker removes tabs as well
  */
  int i;
  int len;
  int non_space;
  for(i=0;i<strlen(sequence_name);i++) 	{
    if(sequence_name[i]=='\n' || sequence_name[i]=='\t') 
      sequence_name[i]=' ';
  }
  len=strlen(sequence_name);
  i=len-1;
  /* remove extra spaces at end of name */
  while( (i>4)&&(sequence_name[i]==' ')) {
    sequence_name[i]='\0';
    i--;
  }
  /* remove leading spaces at start of name
   * find first non-space */
 len=strlen(sequence_name);
 len--;
 i=0;
 non_space=0;
 while((sequence_name[i]==' ')&&(i<len)) {
   non_space++;
   i++;
 }
 if(non_space>0) {
   for(i=non_space;i<=strlen(sequence_name);i++)
     sequence_name[i-non_space]=sequence_name[i];
 }
}

int compare_helices(const void *h1,const void *h2) {
  int dif;
  struct helix *helix1,*helix2;
  helix1=(struct helix *)h1;
  helix2=(struct helix *)h2; 
  dif=(*helix1).diagonal-(*helix2).diagonal; /* return <0 if diag1 < diag2 */
  if(dif==0)                            /* on same diagonal */
    dif=(*helix1).row-(*helix2).row; /* return <0 if row1<row2 or >0 if > */
  return dif;
}   

void sort_helices(struct helix *diag,int helices_in_plot_file) {
  qsort(diag, helices_in_plot_file, sizeof(struct helix), compare_helices);
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

void read_ct_file(FILE *fp,int file,int *len,struct helix *diag,
		  char file_data_ct[MAXIMUM_CT_FILES][255],
		  int *helices_in_plot_file, int *first_sequence_length,
		  char *g_sequence_name,int sequence_name_set,char *filename) {
  char record[160];
  int file_length;
  int row,previous,next,column,junk;
  int helix_len;
  int helices_in_this_file; 
  char base;
  char energy_string[140];
  int first_char;
  int helix_started;
  int j,i;
  int sequence_length;
  int previous_column;
  int name_start;
  int energy_start;
  int record_length;
  char c1,c2;
  char junk_length[150];
  int last_char;
  char sequence_name[150];
 sequence_length=0;
 helices_in_this_file=*helices_in_plot_file;
 /* diag is the diagonal that each row,col corresponds to */
 if(fgets(record,90,fp)==NULL) {
   printf("Error!\tct file is empty\n");
   try_exit(10);
 } else {
/* The top record should look like
 * 741 dG =  -238.3    mok-lacZ + gtgtaa 
 * where '741' is the sequence length,
 * dG = -238.3 is the free energy and the rest of the record
 * is the sequence name. Only the free energies may vary from one
 * ct file to another. The free energies are stored in file_data_ct 
 */
   sscanf(record,"%d",&file_length); /* assign length */
   sscanf(record,"%s",junk_length); /* assign length string */
   first_char=0;
   name_start=0;
   energy_start=0;
   record_length=strlen(record);
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
   name_start=first_char;
   /* this should copy energy into file_data_ct[file]*/
   sprintf(file_data_ct[file],"%d",file+1);
   /* if the filename ends with _# append that number to */
   /* file_data_ct */
   append_structure_name(file_data_ct[file],filename);
   if(name_start>30) { /* failed to find energy */
     name_start=6;
     strcpy(energy_string,"Energy is Undef.");
   } else {
     for(i=energy_start,j=0;i<name_start;i++,j++)
       energy_string[j]=record[i];
     energy_string[j]='\0';
   }
   last_char=strlen(energy_string)-1;
   if(energy_string[last_char]=='\n')
     energy_string[last_char]='\0';
   strcat(file_data_ct[file]," ");
   strcat(file_data_ct[file],energy_string);
   if((file==0)&&(!sequence_name_set)) { /* copy the name in */
     for(i=name_start,j=0;(i<record_length)&&(j<79);i++,j++)
       sequence_name[j]=record[i];
     sequence_name[j]='\0';
     fix_name(sequence_name);
     strcpy(g_sequence_name,sequence_name);
   }
 }
 /* Read the ct file */
 helix_started=FALSE;
 helix_len=0;
 previous_column=-10;
 sequence_length=0; 
 while(fgets(record,90,fp)!=NULL) {
   sscanf(record,"%d %c %d %d %d %d",&row,&base,&previous,&next,&column,&junk);
   if(row!=(sequence_length+1)) {
     printf("Error in ct file! Format is wrong\n");
       printf("Record %d of the file contains %s,\n",
	      sequence_length+1,record);
       printf("but it should start with %d\n",sequence_length+1);
       try_exit(11);
   }
   if(file==0)
     sequence[sequence_length]=base;
   else {
     if(sequence_length>*first_sequence_length) {
       printf("The sequences have different lengths!\n");
       printf("Lengths are %d and at least %d\n",
	      *first_sequence_length, sequence_length);
       printf("Error!\tProgram abnormally terminated\n");
       try_exit(12);
     }
     if(base!=sequence[sequence_length]) {
       c1=toupper(base);
       c2=toupper(sequence[sequence_length]);
       if( !(((c1=='U')&&(c2=='T'))||((c1=='T')&&(c2=='U')))) {
	 printf("Warning!\tAt position %d, sequence %s has a %c,\n",
		sequence_length+1,file_data_ct[0],sequence[sequence_length]);
	 printf("but sequence %s has a %c\n",
                file_data_ct[file],base);
	 printf("Program continues anyway\n");
       }
     }
   }
   sequence_length++;
   if(helix_started==TRUE) {
     if(column==0) { /* finish off helix */  
       diag[*helices_in_plot_file].length=helix_len;
       helix_started=FALSE;
       (*helices_in_plot_file)++;
     } else
       if(column<row) { /* col<row has already been added 
			 * finish off helix; all helices appear twice */
	 diag[*helices_in_plot_file].length=helix_len;
	 helix_started=FALSE; /* just finish current helix */
	 (*helices_in_plot_file)++;
       } else
	 if(column!=(previous_column-1)) { /* finish off helix */
	   diag[*helices_in_plot_file].length=helix_len;
	   (*helices_in_plot_file)++;
	   helix_started=FALSE;
	   /* start new helix */
	   previous_column=column;
	   helix_len=1;
	   if(*len<column)
	     *len=column;
	   diag[*helices_in_plot_file].row=row;
	   diag[*helices_in_plot_file].column=column;
	   diag[*helices_in_plot_file].diagonal=row+column-1;
	   diag[*helices_in_plot_file].color=file+2;
	   helix_started=TRUE;
	 } else {
	   helix_len++; /* add one to length */
	   previous_column=column;
	 }
      } else {  /* helix has not been started */
     if(column>row) { /* start a helix */
        previous_column=column;
	helix_len=1;
	diag[*helices_in_plot_file].row=row;
	diag[*helices_in_plot_file].column=column;
	diag[*helices_in_plot_file].diagonal=row+column-1;
	diag[*helices_in_plot_file].color=file+2;
	helix_started=TRUE;
	if(*len<column)
	  *len=column;
     }
   }
 }
 /* if a helix was started, record its length */
 if(helix_started==TRUE) {
   diag[*helices_in_plot_file].length=helix_len;
   diag[*helices_in_plot_file].diagonal=row+column-1;
   (*helices_in_plot_file)++;
 }
 if(*helices_in_plot_file==MAXIMUM_HELICES) {
   printf("Error!\tThere are more than %d helices in the plot file.\n",
	  MAXIMUM_HELICES);
   printf("Increase the constant MAXIMUM_HELICES in %s.c", PROGRAM_NAME);
   printf(" and recompile.\n");
   try_exit(13);
 }
 helices_in_this_file=*helices_in_plot_file-helices_in_this_file;
 printf("There are %d helices in this ct file with length %d.\n",
	helices_in_this_file,sequence_length);
 if(file==0)
   *first_sequence_length=sequence_length;
 else {
   if(sequence_length!=*first_sequence_length) {
     printf("Error!\tThe sequences have different lengths.\n");
     printf("The first sequence length is %d. Another is %d.\n",
	    *first_sequence_length,sequence_length);
     try_exit(14);
   }
 }
 if(file_length!=sequence_length) {
   printf("Error!\tThere are %d records, but the first line ", 
	  sequence_length);
   printf("indicates that there should be %d records.\n", file_length);
   try_exit(15);
 }
}

void fix_counts(int helices_in_plot_file, int *diag_count,struct helix *diag) {
  int i;
  for(i=0;i<helices_in_plot_file;i++) {
    diag_count[diag[i].diagonal]++;
  }
} 

void sort_all_helices_set_start(int *diag_count, int *diag_start,
				int helices_in_plot_file, struct helix *diag,
				int len) {
  int i;
  if(helices_in_plot_file==0) {
    printf("Error!\tNo helices_have_been_read.\n");
    try_exit(16);
  } else {
    printf("There are %d  helices in the ct files.\n",helices_in_plot_file);
  }

 /* diag_count counts the total number of entries in each diagonal
  *
  * Sort the plot file
  */
  sort_helices(diag,helices_in_plot_file);
  for(i=0;i<=(2*len-1);i++) 
    diag_count[i]=0;
  /* fix the count of diagonal for each entry in plot file */
  diag_start[1]=0;
  fix_counts(helices_in_plot_file,diag_count,diag);
  for(i=2;i<=(2*len-1);i++) {
    diag_start[i] = diag_start[i-1]+diag_count[i-1];
  } 
  /* diag_start counts forward into each array for the start of each
     diagonal */ 
}

/* add overlaps for an individual diagonal named current_diagonal */
void add_overlap_diag(int *diagonal_fill,int current_diagonal, int *diag_start,
		      int *diag_count, struct helix *diag,
		      int number_of_ct_files,int *new_helices_in_plot_file) {
  /* diagonal_fill is a diagonal line that is filled to indicate how
   * many dots at each position diagonal_count will be updated */
  int count,i,starting_position,farthest_position;
  int row,length,current_row,helix_started;
  int start_row_of_helix = 0;
  int closest_position;
  int length_of_current_helix;
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
    row=diag[starting_position+i].row;
    length=diag[starting_position+i].length;
    /* process a single helix */
    if(farthest_position<(row+length-1))
      farthest_position=row+length-1;
    if(closest_position>row)
      closest_position=row;
    for(current_row=row;current_row<=(row+length-1);current_row++) {
      diagonal_fill[current_row]++;
    }
  }
  /* The number at each entry along the diagonal indicates some
   * overlap ,  0 or 1 is no overlap
   * Value equal to number_of_ct_files indicates full overlap
   * Other indicates partial overlap
   * Make overlaping helices for partial overlap*/
  if(number_of_ct_files>2) {
    helix_started=FALSE;
    length_of_current_helix=0;
    for(i=closest_position;i<=(farthest_position+1);i++) {
      if((diagonal_fill[i]>1)&&(diagonal_fill[i]<number_of_ct_files)) {
	/* continue helix or start new helix */
	if(helix_started==TRUE)
	  length_of_current_helix++;
	else {
	  length_of_current_helix=1;
	  start_row_of_helix=i;
	  helix_started=TRUE;
	}
      } else { /* end helix if one was started or do nothing*/
	if(helix_started==TRUE) {
	  diag_count[current_diagonal]++;
	  diag[*new_helices_in_plot_file].row=start_row_of_helix;
	  diag[*new_helices_in_plot_file].diagonal=current_diagonal;
	  diag[*new_helices_in_plot_file].color=1;
	  diag[*new_helices_in_plot_file].column=current_diagonal-
	    start_row_of_helix+1;
	  diag[*new_helices_in_plot_file].length=length_of_current_helix;
	  (*new_helices_in_plot_file)++;
	  helix_started=FALSE;
	  length_of_current_helix=0;
	}
      }
    }
  }
  /* do full overlap last */
  helix_started=FALSE;
  length_of_current_helix=0;
  for(i=closest_position;i<=(farthest_position+1);i++) {
    if(diagonal_fill[i]==number_of_ct_files) {
      /* continue helix or start new helix */
      if(helix_started==TRUE) {
	length_of_current_helix++;
      } else {
	length_of_current_helix=1;
	start_row_of_helix=i;
	helix_started=TRUE;
      }
    } else {
      if(helix_started==TRUE) {
	diag_count[current_diagonal]++;
	diag[*new_helices_in_plot_file].row=start_row_of_helix;
	diag[*new_helices_in_plot_file].diagonal=current_diagonal;
	diag[*new_helices_in_plot_file].color=0;
	diag[*new_helices_in_plot_file].column = 
	  current_diagonal-start_row_of_helix+1;
	diag[*new_helices_in_plot_file].length=length_of_current_helix;
	(*new_helices_in_plot_file)++;
	helix_started=FALSE;
	length_of_current_helix=0;
      }
    }
    diagonal_fill[i]=0; /* clear the array for the next pass */
  }
}

/* add overlaps for all diagonals */
void  add_overlaps(int *diag_count, int *diag_start, int *helices_in_plot_file,
                   struct helix *diag, int len, int number_of_ct_files) {
  int new_helices_in_plot_file;
  int i;
  int *diagonal_fill;
  diagonal_fill=(int *)malloc(len*sizeof(int));
  /* The length of longest diagonal, this was too big
   * Should probably use square_root of 2 * len/2+1
   * Then clear it */
  if(diagonal_fill==NULL) {
    printf("Error!\tInsufficient memory to allocate diagonal_fill\n");
    try_exit(17);
  }
  for(i=0;i<len;i++)
    diagonal_fill[i]=0;
  /* Save old information */
  new_helices_in_plot_file = *helices_in_plot_file;
  /* add overlap helices */ 
  for(i=1;i<=(2*len-1);i++) {
    add_overlap_diag(diagonal_fill, i, diag_start, diag_count, diag,
		     number_of_ct_files, &new_helices_in_plot_file);
  }
  /* update for new starts */
  diag_start[1] = 0;
  for(i=2;i<=(2*len-1);i++)
    diag_start[i]=diag_start[i-1]+diag_count[i-1];
  *helices_in_plot_file=new_helices_in_plot_file;
  free(diagonal_fill);
}

void read_of_all_ct_files(FILE *fp, int *diag_count, int *diag_start,
			  int number_of_ct_files, int *helices_in_plot_file,
			  int *len, char filename_ct[MAXIMUM_CT_FILES][255],
			  struct helix *diag, 
			  char file_data_ct[MAXIMUM_CT_FILES][255],
			  char *sequence_name, int sequence_name_set) {
  char filename[255];
  int first_sequence_len;
  int i,successful_opens;
  successful_opens=0;  
  printf("number_of_ct_files = %d\n", number_of_ct_files);
  for(i=0;i<number_of_ct_files;i++) {
    strcpy(filename,filename_ct[i]); 
    printf("Trying to open %s.\n",filename);   
    if((fp=fopen(filename,"r"))==NULL) {
      printf("Could not open file %s.\n", filename);
      strcpy(file_data_ct[i],"File not found: ");
      strcat(file_data_ct[i],filename);
    } else {
      successful_opens++;
      read_ct_file(fp, i, len, diag, file_data_ct, helices_in_plot_file,
		   &first_sequence_len, sequence_name, sequence_name_set,
		   filename);
      fclose(fp);
    }
  }
  printf("successful_opens = %d and *len = %d\n", successful_opens, *len);
  if((successful_opens==0) || (*len==0)) {
    printf("Error!\tUnable to open any input files successfully.\n");
    try_exit(18);
  }  
  printf("Farthest base pair is %d\n",*len);
  sort_all_helices_set_start(diag_count, diag_start, *helices_in_plot_file,
			     diag, *len); 
  add_overlaps(diag_count, diag_start, helices_in_plot_file, diag, *len, 
	       number_of_ct_files);
  sort_helices(diag,*helices_in_plot_file);
  printf("After overlap there are %d helices.\n", *helices_in_plot_file);
}
