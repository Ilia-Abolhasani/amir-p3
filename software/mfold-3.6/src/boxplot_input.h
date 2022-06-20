/* boxplot_input.h for boxplot(_ng) */

int count_helices(void){
  int helices;
  char rec[100];
  helices=0;
  while(fgets(rec,90,fp)!=NULL) {
    helices++;
  }
  if(helices<1) {
    printf("Error!\tNo data in plot file.\n");
    try_exit(6);
  }
  printf("The number of helices is %d\n",helices-1);
  return helices;
}

void check_parameters(int,char **);
void initialize_data2(void);

void open_file(int argc, char **argv) {
  int helix_total; 
  helix_total=0;
  check_parameters(argc,argv);
  printf("Trying to open %s\n",g_plot_filename);   
 if((fp=fopen(g_plot_filename,"r"))==NULL) {
   printf("Error!\tCould not open plot file: %s.\n",
	  g_plot_filename);
   try_exit(7);
 } 
 helix_total=count_helices();
 fclose(fp);
 g_helix_array_size=helix_total;
 g_diag=malloc(sizeof(struct helix)*(helix_total+1));
 if(g_diag==NULL) {
   printf("Error in boxplot.\tInsufficient memory to ");
   printf("execute with %d helices in plot file: %s\n",
	  helix_total, g_plot_filename);
   try_exit(8);
 }
 printf("Helix total is %d\n",helix_total);
 if((fp=fopen(g_plot_filename,"r"))==NULL) {
   printf("Error!\tCould not open %s.",g_plot_filename);
   try_exit(9);
 }  
}

int compare_helices(const void *h1,const void *h2) {
  int dif;
  struct helix *helix1,*helix2;
  helix1 = (struct helix *)h1;
  helix2 = (struct helix *)h2; 
  dif=(*helix1).diagonal-(*helix2).diagonal; /* return <0 if diag1 < diag2 */
  if(dif==0) /* on same diagonal */
    dif = (*helix1).row-(*helix2).row; /* return <0 if row1<row2 or >0 if > */
  return dif;
} 

void sort_helices(void) {
  qsort(g_diag,g_helices_in_plot_file, sizeof(struct helix),compare_helices);
}

/* Not used. For debug only. 
   void print_helix(int number_of_helices) {
   int i;
   for(i=0;i<number_of_helices;i++) {
   printf("\n diag=%d, row=%d, color=%d, column=%d, energy=%d, length=%d",
   g_diag[i].diagonal,
   g_diag[i].row,
   g_diag[i].color,
   g_diag[i].column,
   g_diag[i].energy,
   g_diag[i].length);
   }
   }
*/

void fix_counts(void) {
  int i;
  for(i=0;i<g_helices_in_plot_file;i++) {
    g_diag_count[g_diag[i].diagonal]++; /* Update the diagonal count. */
  }
}

int find_color(int, float, int);

void fix_colors(void) { /* set the colors for each line of input */
  int i;
  for (i=0;i<g_helices_in_plot_file;i++) {
    g_diag[i].color = find_color(g_diag[i].energy, g_color_increment,
				 g_optimal_energy);
  }
}

void initialize_len(int adjust_col_flag, int adjust_amount) {
  /* g_worst_energy and g_optimal_energy are also set */
  int level,length,row,col,energy,diag,i;
  /* diag is the diagonal corresponding to each row and col
   * zuker_convert_factor is not used, but is useful if units change
   * in the energy dot plot. The current assumption is that units are
   * 1/10 kcal/mol.
   */
  float rec_energy,zuker_convert_factor=1.0;
  g_opt_prob_flag=FALSE;
  g_helices_in_plot_file=0;
  
  if(fgets(rec,90,fp)==NULL) {
    printf("Error!\tPlot file is empty\n");
    try_exit(10);
  }
  g_length=0;
  /* Read the plot file */
  g_worst_energy=INT_MIN;
  g_optimal_energy=INT_MAX;
  while(fgets(rec,90,fp)!=NULL) {
    sscanf(rec,"%d%d%d%d%f",&level,&length,&row,&col,&rec_energy);
    if(adjust_col_flag) {
      printf("(%d,%d) converted to ",row,col);
      col = col - adjust_amount;
      if ((col-length+1)<1) {
	printf("Error!\tColumn %d rounded up to %d\n", 
	       col + adjust_amount, length);
	col = length;
      }
      printf("(%d,%d)\n",row,col);
      if(g_length<(row+length-1))
	g_length=row+length-1;
    }
    if(g_length<col)
      g_length=col;
    g_diag[g_helices_in_plot_file].length=length;
    g_diag[g_helices_in_plot_file].row=row;
    g_diag[g_helices_in_plot_file].column=col;
    if(g_prob)
      energy = - (int)(1000000.*rec_energy + .50); /* M. Zuker, May 30, 2007 */
    else
      /* July 27, 2005. Zuker adds scaling for energy dot plot It was
       * always assumed that units were tenths of a kcal/mol so
       * truncation of last digit occurred with old software. 
       * This allows for higher precision. No truncation is needed.
       */ 
      energy=(int)(rec_energy/zuker_convert_factor);
    g_diag[g_helices_in_plot_file].energy=energy;
    if( (g_prob) && (energy<(-1*2*1000000+1)) ) {
      g_opt_prob_flag=TRUE;
    } else {
      if(energy<g_optimal_energy)
	g_optimal_energy=energy;
    } 
    if(energy>g_worst_energy)
      g_worst_energy=energy;
    diag = row + col - 1;
    g_diag[g_helices_in_plot_file].diagonal=diag;
    g_helices_in_plot_file++;
    if(g_helices_in_plot_file==g_helix_array_size) {
      printf("Error!\tThere was an error reading the plot file.\n");
      printf("%d helices were expected, but",g_helix_array_size);
      printf("the number of lines in the plot file was less.\n");
      printf("Strange error\tLook in boxplot.c for clues.\n");
      try_exit(11);
    }
  }
  if(g_helices_in_plot_file==0) {
    printf("Error!\tThe plot file contains no helices\n");
    try_exit(12);
  }
  printf("There are %d helices in the plot file.\n",
	 g_helices_in_plot_file);
  /* glob_diag_count counts the total number of entries in each diagonal */
  fclose(fp);
  /* Sort the plot file */
  sort_helices();
  g_diag_start=(int *)malloc((2*g_length)*sizeof(int));
  if(g_diag_start==NULL)
    printf("Error!\tInsufficient memory for g_diag_start\n");
  g_diag_count=(int *)malloc((2*g_length)*sizeof(int));
  if(g_diag_count==NULL)
    printf("Error!\tInsufficient memory for g_diag_count\n"); 
  for(i=0;i<=(2*g_length-1);i++)
    g_diag_count[i]=0;
  /* fix the count of diagonal for each entry in plot file */
  initialize_data2();
  if(g_opt_prob_flag==TRUE)
    printf("Optimal structure found in data\n");
  fix_colors();
  g_diag_start[1]=0;
  fix_counts();
  for(i=2;i<=(2*g_length-1);i++) {
    g_diag_start[i]=g_diag_start[i-1]+g_diag_count[i-1];
    /* g_diag_start counts forward into each array for the start of
     *  each diagonal 
     */ 
  }
}
 
/* Conversion Routines to display integers and floats */
char g_str[30];

char *num_string_int(int x) { /* converts integer to string */
  sprintf(g_str,"%d",x);
  return g_str;
}

char *num_string_float(float x) { /* convert float to string */
  sprintf(g_str,"%0.6f",x);
  return g_str;
}

char *num_prob_string(int x) { 
  float x_float;
  x_float = (float)x;
  if (x_float < 1.0e-6) {
    x_float = 0.0;
  } else {
    x_float = -x_float/1000000. ;   /* round to million position */
  }
  sprintf(g_str,"%.6f",x_float);
  return g_str;
} 

char *num_energy_string(int x) { /* Divide by 10 and convert to string */ 
  float x_float;
  x_float = (float)x;
  x_float = x_float/10.;
  sprintf(g_str,"%7.1f",x_float);   /* round to tenth position */
  return g_str;
}

char *num_prob_string_float(float x) { 
  x = -x/1000000.;
  if ( x < 1.0e-6 )
    x = 0.0;
  sprintf(g_str,"%.6f",x);     /* round to millionth position */
  return g_str;
} 

char *num_energy_string_float(float x) { /* Divide by 10 & convert to string */
  x = x/10.;      /* round to tenth position */
  sprintf(g_str,"%7.1f",x);
  return g_str;
}

char *num_string_fancy_int(int x) { /* Integer to string for prob or energies*/
  float float_x;
  float_x = (float)x;
  if(g_prob==TRUE)
    return num_prob_string_float(float_x);
  else
    return num_energy_string_float(float_x);
}

char *num_string_fancy_float(float x) { /* float to string */
  if(g_prob==TRUE)
    return num_prob_string_float(x);
  else
    return num_energy_string_float(x);
}

int step_fun(int dis_l,int dis_r) {
  int step;
  step = (dis_r-dis_l + 1)/8; /* make 8 steps  (9 tick marks total) */
  if (step<40) {
    if (step<3) {
      step=2;
    } else if (step>5) {
      step = (step + 5)/10; /* make a large step a multiple of ten */
      step = step*10;
    } else if (step==4) /* make steps of 4 look like 5 */
      step = 5;
  } else if(step>90) {
    step = (step + 50)/100; /* round to nearest 100 */
    step = step*100;
  } else {
    step = (step + 25)/50; /* round to nearest 50 */
    step = step*50;
  }
  return step;
}

/* make the tick marks start at a nice place 
 * ends are treated separately and always drawn 
 */
int start_fun(int dis_l,int step) {
  int start;
  start=dis_l+step;
  if(step>90) {
    start = (start + 60)/100; /* start at a multiple of 100 */
   start = start*100;
  } else if(step>=40) {
    start = (start + 32)/50; /* start at a multiple of 50 */
    start = start*50;
  } else if (step>5) {
    start = (start + 5)/10; /* start at a multiple of ten */
    start = start*10;
  }
  return start;
}


/* other functions _______________________________*/

void finish_setting_energy_cutoff(float temp_filter) {
  g_energy_cutoff=(int)(temp_filter*10);
  if(g_energy_cutoff>=0)
    g_energy_cutoff=g_optimal_energy+g_energy_cutoff;
  if(g_energy_cutoff>g_worst_energy)
    g_energy_cutoff=g_worst_energy;
  if(g_energy_cutoff<g_optimal_energy)
    g_energy_cutoff=g_optimal_energy;
  g_color_increment = (float) (g_energy_cutoff - g_optimal_energy)
    /(float)(g_number_of_colors-1);
}

/* given a row and column, return the energy of that position */
int energy_read(int row,int col) {
  int diag, diag_start, diag_end, i, temp_row;
  diag = col + row - 1;
  diag_start = g_diag_start[diag];
  diag_end = g_diag_start[diag+1];
 for(i=diag_start;i<diag_end;i++) {
   temp_row=g_diag[i].row;
   if (temp_row > row) 
     return 0;   /*  entries in diagonal are sorted */
   if( (!g_opt_prob_flag) || (g_diag[i].energy >= -20000) ) {
     if((temp_row<=row) && (g_diag[i].row+g_diag[i].length-1>=row))  
       return g_diag[i].energy;
   }
 }
 return 0;
} 

int energy_read_opt_prob(int row,int col) {
  int diag,diag_start,diag_end,i,temp_row;
  diag = col + row - 1;
  diag_start = g_diag_start[diag];
  diag_end = g_diag_start[diag+1];
  for(i=diag_start;i<diag_end;i++) {
    temp_row=g_diag[i].row;
    if(temp_row>row) 
      return 0;   /*  entries in diagonal are sorted */
    if(g_diag[i].energy<-20000) {
      if((temp_row<=row) && (g_diag[i].row + g_diag[i].length-1>=row))  
	return g_diag[i].energy;
    }
  }
  return 0;
} 

FILE *engfp;

int open_energy(char *filename) {    /* open specified output file */
  if ((engfp = fopen(filename, "w")) == NULL) {
    printf ("Could not open file: %s\n", filename);
    return(1);
  }
  return (0);
}

void make_energy_file(char *filename, int label_i, int label_j) {
  /* make energy file */
  float energy, probability;
  printf("gifeng is %s at location (%d,%d).\n", filename, label_i, label_j);
  if(open_energy(filename)) {
    printf("Error creating .gifeng file %s.", filename);
    return;
  }
  if((label_i<1) || (label_i>g_length) || (label_j<1) || (label_j>g_length))
    energy = 9999.9;
  else {
    if (g_prob) {
      probability = -(float) energy_read(label_i, label_j)/1000000.;
      fprintf(engfp, "%.6f\n", probability);
    } else {
      energy = ((float) energy_read(label_i, label_j) + .0001)/10.;
      fprintf(engfp,"%.1f\n", energy);
    }
    fclose(engfp);
  }
}
