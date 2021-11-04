FILE *postfp;
#include "header_ps.h"

int open_post(char *post_filename) {
  /* open specified output file */
  if ((postfp = fopen(post_filename, "w")) == NULL) { /* open a file */
    printf ("Error!\tCould not open file %s,\n",post_filename);
    return(1);
  }
  return (0);
}
void fixcolor_plot_postscript(int color) { /*  postscript colors */
  if((color>=0)&&(color<TOTAL_COLORS))
    fprintf(postfp,"\n color_%d\n",color);
}

void fixcolor_plot_postscript_define(float color_set[TOTAL_COLORS][3]) {
  /*  postscript colors */
  int color;
  fprintf(postfp,"\n%s Set colors","%%");
  fprintf(postfp,"\n%s 0 is for text, 1 is optimal dots in ct mode","%%");
  fprintf(postfp,"\n%s 2,3,4 are higher energy dots from plot","%%");
  fprintf(postfp,"\n%s 5 is optimal,ct overlap ","%%");
  fprintf(postfp,"\n%s 6 is ct overlap energies near optimal in ct mode","%%");
  fprintf(postfp,"\n%s 6 is ct within a distance of optimal in opt mode","%%");
  fprintf(postfp,"\n%s 7 is ct, not optimal in ct mode","%%");
  fprintf(postfp,"\n%s 7 is optimal without ct in opt mode","%%");
  fprintf(postfp,"\n%s 9 is Labels, 8 is background, 10 is grid lines","%%");
  fprintf(postfp,"\n%s 11 is other ct dots in optimal mode","%%"); 
  fprintf(postfp,"\n%s original colors from  boxplot.col\n\n","%%");
  for(color=0;color<TOTAL_COLORS;color++) {
    fprintf(postfp,"/color_%d {  %.3f %.3f %.3f setrgbcolor} def \n",color,
	    color_set[color][0],color_set[color][1],color_set[color][2]);
  }  
}

void single_hor_tic_zoom_post(int j_col,float scale, int left,int grid_flag) {
  fprintf(postfp,"\n %f 1. scale",scale);
  fprintf(postfp,"\n\n %f setlinewidth",0.40/scale);
  fprintf(postfp,"\n %d 540 moveto",j_col-left+1);
  fprintf(postfp,"\n 0 8 rlineto");
  if(grid_flag==TRUE) {
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_GRID);
    fprintf(postfp,"\n %d 539 moveto",j_col-left+1);
    fprintf(postfp,"\n 0 -538 rlineto");
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_TEXT);
  }
  fprintf(postfp,"\n %d 550 moveto",j_col-left+1);
  fprintf(postfp,"\n %f 1. scale",1./scale);
  fprintf(postfp,"\n -4 0 rmoveto");
  fprintf(postfp,"\n (%d) show",j_col);
}

void plot_hor_tics_zoom_post(float scale,int left, int right,int grid_flag) { 
  int j_col,step;
  int start;
  step=step_fun(left,right);
  start=start_fun(left,step);
  single_hor_tic_zoom_post(right,scale,left,FALSE);
  single_hor_tic_zoom_post(left,scale,left,FALSE);
  for(j_col=start;j_col<(right-2*step/3); j_col=j_col+step) {
    single_hor_tic_zoom_post(j_col,scale,left,grid_flag);
  }
}

void single_ver_tic_zoom_post(int i_row,float scale,int bottom,int grid_flag) {
  fprintf(postfp,"\n 1. %f scale",scale);
  fprintf(postfp,"\n\n %f setlinewidth",0.40/scale);
  fprintf(postfp,"\n 540.5 %d moveto",bottom-i_row+1);
  fprintf(postfp,"\n 8 0 rlineto");
  if(grid_flag==TRUE) {
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_GRID);
    fprintf(postfp,"\n 539.5 %d moveto",bottom-i_row+1);
    fprintf(postfp,"\n -538 0 rlineto");
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_TEXT);
  }
  fprintf(postfp,"\n 550.5 %d moveto",bottom-i_row+1);
  fprintf(postfp,"\n 1. %f scale",1./scale);
  fprintf(postfp,"\n (%d) show",i_row);
}

void plot_ver_tics_zoom_post(float scale, int top,int bottom,int grid_flag) {
  int i_row,step,start;
  step=step_fun(top,bottom);
  start=start_fun(top,step);
  single_ver_tic_zoom_post(top,scale,bottom,FALSE);
  single_ver_tic_zoom_post(bottom,scale,bottom,FALSE);
  for(i_row=start;i_row<(bottom-step/2);i_row=i_row+step) {
    single_ver_tic_zoom_post(i_row,scale,bottom,grid_flag);
  }
}

void draw_diagonal_line_post(int display_l,int display_r,int display_t,
			     int display_b) {
  int left,top,bottom,right;
  if((display_r>=display_t)&&(display_b>=display_l)) {
    top=display_t;
    bottom=display_b;
    right=display_r;
    left=display_l;
    if(display_r>=display_b) {
      bottom=display_b;
      right=display_b;
    } else {
      right=display_r;
      bottom=display_r;
    }
    if(display_l>=display_t) {
      left=display_l;
      top=display_l;
    } else {
      top=display_t;
      left=top;
    }
    fixcolor_plot_postscript(COLOR_TEXT);
    fprintf(postfp,"\n %f  %f moveto", left-display_l+0.0, display_b-top+2.0);
    fprintf(postfp,"\n %f %f  lineto stroke", right-display_l+2.0,
	    display_b-bottom+0.0);
  } 
}

void display_time(void) {
  time_t now;
  now=time(NULL);
  fprintf(postfp,"\n%s %s","%%",ctime(&now));
}

void draw_postscript_dots_zoom(int left, int right, int top, int bottom,
			       int number_of_colors, int prob_flag,
			       int *diag_start, int *diag_count, 
			       struct helix *diag, int energy_cutoff,
			       int chain_len, int points_plotted[8],
			       int plot_len, int optimal_energy,
			       int ct_len, float dot_size[4], 
			       int *ct_diag_start, int *ct_diag_count, 
			       struct ct_helix *ct_diag, int worst_energy,
			       float close_optimal_energy) {
  int row,col,new_row,count;
  int current_diag,position,location,new_color;
  int color;
  float size; 
  int start_diag,end_diag;
  int plot_display_b,plot_display_r;
  int ct_display_b,ct_display_r;
  int type;
  float size_multiplier;
  float big_dot_size;
  float width,height,average,extra_multiplier;
  int current_type;
  int ct_type,energy;
  int i;
  start_diag = top + left - 1;
  /* 0 black, gray dots
   * 1 green dots
   * 2 yellow dots
   * 3 red dots
   * 4 red , gray overlap dots
   * 5 yellow, gray overlap
   * 6 Optimal dots
   * 7 blue ct missed dots, none here
   * take care of black,gray dots 
   */
  for(i=0;i<8;i++)
    points_plotted[i]=0;
  if (bottom>plot_len)
    plot_display_b=plot_len;
  else
    plot_display_b=bottom;
  if (right>plot_len)
    plot_display_r=plot_len;
  else
    plot_display_r=right;
  end_diag=plot_display_b+plot_display_r-1;
  for(color=number_of_colors;color>=1;color--) {
    fixcolor_plot_postscript(color);
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=diag_start[current_diag];
      for(count=0;( (count<diag_count[current_diag]) 
		    &&(diag[position+count].row<=bottom) ) ;count++) {
	location=position+count;
	if( ( diag[location].energy<=energy_cutoff &&
	      diag[location].length>=chain_len ) ||
	    diag[location].energy==optimal_energy) {
	  new_color=diag[location].color;
	  if(new_color==color) {
	    new_row = diag[location].row;
	    col = diag[location].column;
	    for(row=new_row; row<(new_row+diag[location].length);row++) {
	      if((row>=top) && (row<=bottom) && (col>=left) && (col<= right)) {
		points_plotted[0]++;
		if(prob_flag!=TRUE) {
		  fprintf(postfp,"\n%d %d b", row-top+1,col-left+1);
		} else {
		  size=sqrt((-1.0*(float)diag[location].energy)/1000.);
		  fprintf(postfp,"\n%d %d %.5f pb", row-top+1,col-left+1,size);
		}
		if(color==1)
		  points_plotted[6]++;  
	      }
	      col--;
	    }
	  }
	}
      }
    }
  }
  /* plot black in lower triangle */
  fixcolor_plot_postscript(1);
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=diag_start[current_diag];
    for(count=0;( (count<diag_count[current_diag]) 
		  &&(diag[position+count].row<=right) ) ;count++) {
      location=position+count;
      if(diag[location].energy==optimal_energy) {
	new_row = diag[location].row;
	col = diag[location].column;
	if(col>=top) {
	  for(row=new_row;  row<(new_row+diag[location].length);row++) {
	    if((col>=top)&&(col<=bottom) && (row>=left) && (row<= right)) {
	      if(prob_flag!=TRUE) {
		fprintf(postfp,"\n%d %d b", col-top+1,row-left+1);
	      }	else { 
		size=sqrt((-1.0*(float)diag[location].energy)/1000.);
		fprintf(postfp,"\n%d %d %.5f pb", col-top+1,row-left+1,size);
	      }
	    }
	    col--;
	  }                        
	}
      }
    }
  }
  /* Make green, red, yellow for ct */
  if(bottom>ct_len)
    ct_display_b=ct_len;
  else
    ct_display_b=bottom;
  if(right>ct_len)
    ct_display_r=ct_len;
  else
    ct_display_r=right;
  end_diag=ct_display_b+ct_display_r-1;
  /* keep dots from getting too small */
  width = right-left;
  height = bottom-top;
  average = (width+height)/2.0;
  if (average>420.0)
    extra_multiplier=average/420.0;
  else
    extra_multiplier=1.0;
  for(current_type=2;current_type>=0;current_type--) {
    size_multiplier=dot_size[current_type+1];
    fixcolor_plot_postscript(current_type+COLOR_OPTIMAL);
    big_dot_size=size_multiplier*extra_multiplier;
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=ct_diag_start[current_diag];
      for(count=0; (count<ct_diag_count[current_diag]) &&
	    ct_diag[position+count].row<=bottom; count++) {
	location=position+count;
	new_row=ct_diag[location].row;
	col=ct_diag[location].column;
	type=set_ct_color(ct_diag[location].energy,
			  optimal_energy,worst_energy,close_optimal_energy);
	if(type==current_type) {
	  if(col>=left) {
	    for(row=new_row;row<new_row+ct_diag[location].length; row++) { 
	      if((row>=top)&&(row<=bottom) &&(col>=left)&&(col<=right)) {
		fprintf(postfp,"\n%d %d %.5f pb",
			row-top+1,col-left+1,big_dot_size);
		if(type==0)
		  points_plotted[1]++;
		else {
		  if(type==1)
		    points_plotted[2]++;
		  else
		    points_plotted[3]++;
		} 
	      }
	      col--;
	    }
	  }
	}
      }
    }
  }
  /* put black,gray dots on top of green,yellow,red */
  for(color=4;color>=1;color--) {
    fixcolor_plot_postscript(color);
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=ct_diag_start[current_diag];
      for(count=0;( (count<ct_diag_count[current_diag])
		    &&(ct_diag[position+count].row<=bottom) ); count++) {
	location=position+count;
	new_row=ct_diag[location].row;
	col=ct_diag[location].column;
	new_color=ct_diag[location].color;
	if(new_color==color) {
	  if(col>=left) {
	    energy=ct_diag[location].energy; 
	    if(energy<=energy_cutoff) {
	      for(row=new_row;row<new_row+ct_diag[location].length; row++) { 
		if((row>=top)&&(row<=bottom) &&(col>=left)&&(col<=right)) {
		  if((energy==optimal_energy)||
		     (ct_diag[location].plotlen>=chain_len)) {
		    ct_type=set_ct_color(energy, optimal_energy, worst_energy,
					 close_optimal_energy);
		    fprintf(postfp,"\n%d %d b",  row-top+1, col-left+1);
		    /* plot black on top of green,yellow,red */
		    if(ct_type==2)
		      points_plotted[4]++;/* red overlap */
		    else if(ct_type==1)
		      points_plotted[5]++;/* yellow overlap */
		  }
		}
		col--;
	      }
	    }
	  }
	}
      }
    }
  }
}


void opt_draw_postscript_dots_zoom(int left, int right, int top, int bottom,
				   int number_of_colors, int prob_flag,
				   int *diag_start, int *diag_count,
				   struct helix *diag, int energy_cutoff,
				   int chain_len, int points_plotted[8],
				   int plot_len, int optimal_energy, 
				   int ct_len, float dot_size[4], 
				   int *ct_diag_start, int *ct_diag_count,
				   struct ct_helix *ct_diag, int opt_distance,
				   int *ct_basepair) {
  int row,col,new_row,count;
  int current_diag,position,location,new_color;
  int color;
  float size; 
  int start_diag,end_diag;
  int plot_display_b,plot_display_r;
  int ct_display_b,ct_display_r;
  float size_multiplier;
  float big_dot_size;
  float width,height,average,extra_multiplier;
  int current_color;
  int energy;
  int i;
  start_diag=top+left-1;
  for(i=0;i<8;i++)
    points_plotted[i]=0;
  /* black, gray dots
   * green dots
   * yellow dots
   * red dots
   * red , gray overlap dots
   * yellow, gray overlap
   * Optimal dots
   * ct missed optimal 
   */
  if(bottom>plot_len)
    plot_display_b=plot_len;
  else
    plot_display_b=bottom;
  if(right>plot_len)
    plot_display_r=plot_len;
  else
    plot_display_r=right;
  end_diag=plot_display_b+plot_display_r-1;
  /* plot gray in upper triangle */
  for(color=number_of_colors;color>=2;color--) {
    fixcolor_plot_postscript(color);
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=diag_start[current_diag];
      for(count=0;( (count<diag_count[current_diag]) 
		    &&(diag[position+count].row<=bottom) ) ;count++) {
	location=position+count;
	if((diag[location].energy<=energy_cutoff)&&
	   (diag[location].length>=chain_len)) {
	  new_color=diag[location].color;
	  if(new_color==color) {
	    new_row=diag[location].row;
	    col=diag[location].column;
	    for(row=new_row; row<(new_row+diag[location].length);row++) {
	      if((row>=top)&&(row<=bottom) &&(col>=left)&&(col<= right)) {
		points_plotted[0]++;
		if(prob_flag!=TRUE) {
		  fprintf(postfp,"\n%d %d b", row-top+1,col-left+1);
		} else {
		  size = sqrt((-1.0*(float)diag[location].energy)/1000.);
		  fprintf(postfp,"\n%d %d %.5f pb", row-top+1,col-left+1,size);
		}
	      }
	      col--;
	    }                        
	  }
	}
      }
    }
  }
  /* plot black in lower triangle */
  fixcolor_plot_postscript(1);
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=diag_start[current_diag];
    for(count=0;( (count<diag_count[current_diag]) 
		  &&(diag[position+count].row<=right) ) ;count++) {
      location=position+count;
      if(diag[location].energy==optimal_energy) {
	new_row=diag[location].row;
	col=diag[location].column;
	if(col>=top) {
	  for(row=new_row; row<(new_row+diag[location].length);row++) {
	    if((col>=top)&&(col<=bottom) && (row>=left)&&(row<= right)) {
	      if(prob_flag!=TRUE)
		fprintf(postfp,"\n%d %d b",  col-top+1,row-left+1);
	      else {
		size=sqrt((-1.0*(float)diag[location].energy)/1000.);
		fprintf(postfp,"\n%d %d %.5f pb", col-top+1,row-left+1,size);
	      }
	    }
	    col--;
	  }                        
	}
      }
    }
  }
  /* Make cyan for missed  ct */
  if(bottom>ct_len)
    ct_display_b=ct_len;
  else
    ct_display_b=bottom;
  if(right>ct_len)
    ct_display_r=ct_len;
  else
    ct_display_r=right;
  end_diag=ct_display_b+ct_display_r-1;
  /* keep dots from getting too small */
  width=right-left;
  height=bottom-top;
  average=(width+height)/2.0;
  if(average>420.0)
    extra_multiplier=average/420.0;
  else
    extra_multiplier=1.0;
  size_multiplier=dot_size[0];
  fixcolor_plot_postscript(COLOR_CT_MISSED_OPTIMAL);
  big_dot_size=size_multiplier*extra_multiplier*.5;
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=ct_diag_start[current_diag];
    for(count=0;( (count<ct_diag_count[current_diag])
		  &&(ct_diag[position+count].row<=bottom) ); count++) {
      location=position+count;
      new_row=ct_diag[location].row;
      col=ct_diag[location].column;
      energy=ct_diag[location].energy;
      if(energy!=optimal_energy) {
	if(col>=left) {
	  for(row=new_row;row<new_row+ct_diag[location].length; row++) { 
	    if((row>=top)&&(row<=bottom) && (col>=left)&&(col<=right)) {
	      fprintf(postfp,"\n%d %d %.5f pb",
		      row-top+1,col-left+1,big_dot_size);
	      points_plotted[7]++;
	    } 
	    col--;
	  }
	}
      }
    }
  }
  /* Draw red, yellow, green optimal dots */
  end_diag=plot_display_b+plot_display_r-1;

for(current_color=2;current_color>=0;current_color--) {
  big_dot_size=dot_size[current_color+1]*extra_multiplier;
  fixcolor_plot_postscript(COLOR_OPTIMAL+current_color);
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=diag_start[current_diag];
    for(count=0;( (count<diag_count[current_diag]) 
       &&(diag[position+count].row<=bottom) ) ;count++) {
      location=position+count;
      if(diag[location].energy==optimal_energy) {
	new_row=diag[location].row;
	col=diag[location].column;
	for(row=new_row; row<(new_row+diag[location].length);row++) {
	  if((row>=top)&&(row<=bottom) && (col>=left)&&(col<= right)) {
	    if(row<=ct_len) {
	      if(ct_basepair[row]==col) {
		new_color=0;
	      } else {
		if(distance_plot_ct(ct_len,row,col, opt_distance)) {
		  new_color=1;
		} else {
		  new_color=2;
		}
	      }
	    } else {
	      new_color=2;
	    }
	    if(new_color==current_color) {
	      points_plotted[1+current_color]++;
	      if(prob_flag!=TRUE) {
		fprintf(postfp,"\n%d %d %.5f pb",
			row-top+1,col-left+1,big_dot_size);
	      } else {
		size=sqrt((-1.0*(float)diag[location].energy)/1000.);
		fprintf(postfp,"\n%d %d %.5f pb", row-top+1,col-left+1,size);
	      }
	    }
	  }
	  col--;
	}                        
      }
    }
  }
 }
 points_plotted[6]=points_plotted[1]+points_plotted[2]+points_plotted[3];
 points_plotted[0]=points_plotted[0]+points_plotted[6];
}

void general_post(int opt_flag, char *post_filename, int post_file, int width,
		  int height, char *plot_filename, char *ct_filename,
		  int optimal_energy, int energy_cutoff, int zoom_labels_count,
		  int prob_flag, int chain_len, int number_of_colors,
		  float color_increment, float prob[4][6],
		  char *zoom_labels_string[50], int zoom_labels_row[50],
		  int zoom_labels_col[50], int left, int right, int top,
		  int bottom, int *diag_start, int *diag_count,
		  struct helix *diag, int grid_flag, int plot_len, int ct_len,
		  float dot_size[4], int *ct_diag_start, int *ct_diag_count, 
		  struct ct_helix *ct_diag, char *sequence_name,
		  int sequence_name_set, float close_optimal_energy,
		  int worst_energy, float color_set[TOTAL_COLORS][3],
		  int opt_distance, int *ct_basepair) {
  float hscale,vscale,x_float;
  int error;
  char str[80];
  int color,offset,i;
  float close_energy_formatted;
  char close_to_optimal_string[30];
  int points_plotted[8]; /* 0: all, 1: green, 2: yellow, 3: red, 4:
			  * red/gray overlap, 5: yellow/gray overlap,
			  * points_plotted[6]: optimal dots, points_plotted[7]:
			  * ct_missed */  
  float center;         
  error = open_post(post_filename);
  if(error) {
    printf("Error creating PostScript file %s.\n",post_filename);
    return;
  }
  vscale = 540.0/(height+1);
  hscale = 540.0/(width+1);
  fprintf(postfp,"%c!PS-Adobe-3.0 EPSF-2.0\n",'%');
  fprintf(postfp,"%sBoundingBox: 0 0 612 792\n","%%");
  fprintf(postfp,"\n%s Created by %s, %s","%%",PROGRAM_NAME, PACKAGE_STRING);
  display_time();
  fprintf(postfp,"%s\n%s\n","%%","%%");
  fprintf(postfp," %f %f scale\n",1.0,1.0);
  fprintf(postfp,"%s\n","%%");
  fprintf(postfp,"%s To change image size ....\n","%%");
  fprintf(postfp,"%s For full size image Use Bound ngBox:0 0 612 792\n","%%");
  fprintf(postfp,"%s                     and 1.0 1.0 scale\n","%%");
  fprintf(postfp,"%s For image .25 of full size, replace 612 with .25*612\n","%%");
  fprintf(postfp,"%s        replace 792 with .25*792 and replace\n","%%");
  fprintf(postfp,"%s        1.0 with .25: Other scales are similar.\n","%%");
  fprintf(postfp,"%s                Valid scales are .05 to 1.00\n","%%");
  fprintf(postfp,"%s ONLY scale and Bound ngBox lines above should be changed.\n",
       "%%");
  fprintf(postfp,"%s\n\n\n\n","%%");
  fixcolor_plot_postscript_define(color_set);
  /* Make smalll hexagon */
  fprintf(postfp,"\n/sbox {"); /* for labels */
  fprintf(postfp,"\n exch len exch sub 1 add ");
  fprintf(postfp,"\n moveto");
  fprintf(postfp,"\n   0.25  .25 rlineto");
  fprintf(postfp,"\n   -.25   0  rlineto");
  fprintf(postfp,"\n  -.25    -.25 rlineto");
  fprintf(postfp,"\n   0   -.25  rlineto");
  fprintf(postfp,"\n   .25    0 rlineto");
  fprintf(postfp,"\n   .25  0.25  rlineto");
  fprintf(postfp,"\n   0    0.25  rlineto");
  fprintf(postfp,"\n   closepath");
  fprintf(postfp,"\n   fill");
  fprintf(postfp,"\n   stroke");
  fprintf(postfp,"\n}def");
  /* Make big hexagon */
  fprintf(postfp,"\n/b {"); /* for non-probability or main window */
  fprintf(postfp,"\n  exch len exch sub 1 add ");
  fprintf(postfp,"\n moveto");
  fprintf(postfp,"\n   .5  .5  rlineto");
  fprintf(postfp,"\n   -.5   0 rlineto");
  fprintf(postfp,"\n  -.5    -.5 rlineto");
  fprintf(postfp,"\n   0   -.5  rlineto");
  fprintf(postfp,"\n   .5    0 rlineto");
  fprintf(postfp,"\n   .5  0.5  rlineto");
  fprintf(postfp,"\n   0    .5  rlineto");
  fprintf(postfp,"\n   closepath");
  fprintf(postfp,"\n   fill");
  fprintf(postfp,"\n   stroke");
  fprintf(postfp,"\n}def");
  /* Make scalable squares */
  fprintf(postfp,"\n/pb {"); /* for boxes whose size depends on probability */
  fprintf(postfp,"\n /s exch def");
  fprintf(postfp,"\n exch len exch sub 1 add bx");
  fprintf(postfp,"\n} def");
  fprintf(postfp,"\n/bx {");
  fprintf(postfp,"\n newpath");
  fprintf(postfp,"\n moveto");
  fprintf(postfp,"\n 0.5 s mul  0 rlineto");
  fprintf(postfp,"\n 0 0.5  s mul  rlineto");
  fprintf(postfp,"\n -1 s mul   0 rlineto");
  fprintf(postfp,"\n 0     -1 s mul rlineto");
  fprintf(postfp,"\n 1 s mul     0 rlineto");
  fprintf(postfp,"\n 0   0.5 s mul rlineto");
  fprintf(postfp,"\n closepath");
  fprintf(postfp,"\n fill");
  fprintf(postfp,"\n stroke");
  fprintf(postfp,"\n} def");
  fprintf(postfp,"\n/len { %d } def",height);
  /* set background color */
  fixcolor_plot_postscript(COLOR_BACK);
  fprintf(postfp,"\n0 0 moveto");
  fprintf(postfp,"\n611 0 rlineto");
  fprintf(postfp,"\n0 791 rlineto");
  fprintf(postfp,"\n-611 0 rlineto");
  fprintf(postfp,"\n0 -791 rlineto");
  fprintf(postfp,"\nclosepath");
  fprintf(postfp,"\nfill");
  fixcolor_plot_postscript(COLOR_TEXT);
  display_header_message_ps(postfp,TRUE);
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 20 scalefont");
  fprintf(postfp,"\n setfont");
  /* center the title */
  if(sequence_name_set) {
    center=72*8.5/2-(strlen(sequence_name))/2.*20.*7./16.;
    fprintf(postfp,"\n %f 720 moveto",center);
    fprintf(postfp,"\n ( %s) show", sequence_name);
  } else {
    center = 72*8.5/2. - (20.+strlen(plot_filename) + strlen(ct_filename))/
      2.*20.*7./16.;
    fprintf(postfp,"\n %f 720 moveto",center);
    if (prob_flag) 
    fprintf(postfp,"\n (Probability dotplot for %s and %s) show", 
	    plot_filename, ct_filename);
    else
      fprintf(postfp,"\n (Energy dotplot for %s and %s) show", plot_filename,
	      ct_filename);
  }
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 10 scalefont");
  fprintf(postfp,"\n setfont");
  if(prob_flag) {
    fprintf(postfp,"\n 70 90 moveto");
    strcpy(str,num_string_float(-1.*((float)energy_cutoff)/1000.));
    strcat(str,"    <= Probability <= 1");
  } else {
    fprintf(postfp,"\n 130 690 moveto");
    strcpy(str,"      Energy Increment: ");
    strcat(str,num_string_fancy_int(energy_cutoff-optimal_energy,prob_flag));
    strcat(str," kcal/mol");
  }
  if(chain_len>1) {
    strcat(str,"                 Filter: ");
    strcat(str,num_string_int(chain_len));
  }
  fprintf(postfp,"\n (%s) show",str); 
  /* display colors */
  fixcolor_plot_postscript(1);
  fprintf(postfp,"\n 35 90  moveto ");
  if(!prob_flag) {
    strcpy(str,"Optimal energy = ");
    strcat(str,num_string_fancy_int(optimal_energy,prob_flag));
    strcat(str,"kcal/mol");
    fprintf(postfp,"\n ( %s ) show",str);
    x_float=(float)optimal_energy;
    for(color=2,offset=0;color<=number_of_colors;color++,offset=offset+15) {
      fixcolor_plot_postscript(color);
      strcpy(str,num_string_fancy_float(x_float+(color-2)*color_increment,
					prob_flag));
      strcat(str," < energy <=");
      strcat(str,num_string_fancy_float(x_float+(color-1)*color_increment,
					prob_flag));
      strcat(str,"  kcal/mol");
      fprintf(postfp,"\n 30 %d  moveto",75-offset);
      fprintf(postfp,"\n ( %s ) show",str);
    }
  } else {
    fixcolor_plot_postscript(1);
    strcpy(str, num_string_float(prob[number_of_colors-4][0]));
    strcat(str," < prob.");
    fprintf(postfp,"\n 30 75 moveto");
    fprintf(postfp,"\n (%s) show",str);
    for(color=2,offset=12;color<=number_of_colors-1;color++,offset=offset+13) {
      fixcolor_plot_postscript(color);
      strcpy(str, num_string_float(prob[number_of_colors-4][color-1]));
      strcat(str," < prob. <=");
      strcat(str, num_string_float(prob[number_of_colors-4][color-2])); 
      fprintf(postfp,"\n 30 %d moveto",75-offset);
      fprintf(postfp,"\n (%s) show",str);
    }
    fixcolor_plot_postscript(color);
    strcpy(str,"             prob. <=");
    strcat(str,num_string_float(prob[number_of_colors-4][number_of_colors-2]));
    fprintf(postfp,"\n 30 %d moveto",75-offset);
    fprintf(postfp,"\n (%s) show",str);
  }
  fixcolor_plot_postscript(COLOR_TEXT);
  fprintf(postfp,"\n 25 100 translate");
  /* label top ticks */
  fprintf(postfp,"\n\n 1.0 setlinewidth");
  plot_hor_tics_zoom_post(hscale,left,right,grid_flag);
  /* label side ticks */
  plot_ver_tics_zoom_post(vscale,top,bottom,grid_flag);
  /* draw box or triangle */ 
  fprintf(postfp,"\n\n 1.0 setlinewidth");
  fprintf(postfp,"\n 0 540.5 moveto 0 0 lineto stroke");
  fprintf(postfp,"\n 0 0.0 moveto 540.5 0 lineto stroke");
  /* draw right and top */
  fprintf(postfp,"\n 0 540.5 moveto 540.5 540.5 lineto stroke");
  fprintf(postfp,"\n 540.5 0  moveto 540.5 540.5 lineto stroke");
  /* set scale */
  fprintf(postfp,"\n %f %f scale",hscale,vscale);
  fprintf(postfp,"\n\n %f setlinewidth",1.0/(hscale+vscale)/2.0);
  draw_diagonal_line_post(left,right,top,bottom);
  if(opt_flag)
    opt_draw_postscript_dots_zoom(left, right, top, bottom, number_of_colors,
				  prob_flag, diag_start, diag_count, diag,
				  energy_cutoff, chain_len, points_plotted,
				  plot_len, optimal_energy, ct_len, dot_size,
				  ct_diag_start, ct_diag_count, ct_diag,
				  opt_distance, ct_basepair);
  else
    draw_postscript_dots_zoom(left, right, top, bottom, number_of_colors,
			      prob_flag, diag_start, diag_count, diag, 
			      energy_cutoff, chain_len, points_plotted,
			      plot_len, optimal_energy, ct_len, dot_size,
			      ct_diag_start, ct_diag_count, ct_diag,
			      worst_energy, close_optimal_energy);
  if((post_file!=MAIN_POST)&&(zoom_labels_count>0)) {
    fixcolor_plot_postscript(COLOR_LABEL);
    fprintf(postfp,"\n /Helvetica findfont");
    fprintf(postfp,"\n 8 scalefont");
    fprintf(postfp,"\n setfont");
    for(i=0;i<zoom_labels_count;i++) {
      fprintf(postfp,"\n%d %d sbox", zoom_labels_row[i]-top+1, 
	      zoom_labels_col[i]-left+1);
      fprintf(postfp,"\n %d %d moveto", zoom_labels_col[i]-left+2,
              bottom-zoom_labels_row[i]+1);
      fprintf(postfp,"\n %f %f scale",1/hscale,1/vscale);
      fprintf(postfp,"\n 0 -1.5 rmoveto");
      fprintf(postfp,"\n (%s) show ",zoom_labels_string[i]);
      fprintf(postfp,"\n %f %f scale",hscale,vscale);
    }  
  }
  fprintf(postfp,"\n %f %f scale",1/hscale,1/vscale);
  fprintf(postfp,"\n -25 -100 translate");
  /* indicate how many points were plotted */
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 10 scalefont");
  fprintf(postfp,"\n setfont");
  fprintf(postfp,"\n 300 90 moveto");
  strcpy(str,"Base Pairs for Plot file: ");
  fixcolor_plot_postscript(COLOR_TEXT);
  strcat(str,num_string_int(points_plotted[0]));
  fprintf(postfp,"\n ( %s ) show",str); 
  fprintf(postfp,"\n 300 75 moveto");
  strcpy(str,"Optimal CT Overlap:     ");
  fixcolor_plot_postscript(COLOR_OPTIMAL);
  strcat(str,num_string_int(points_plotted[1]));
  fprintf(postfp,"\n ( %s ) show",str);
  fixcolor_plot_postscript(1);
  fprintf(postfp,"\n 460 75 moveto");
  fprintf(postfp,"\n ( of %d ) show",points_plotted[6]);  
  if(!opt_flag) {
    if (close_optimal_energy!=0.0) {
      fprintf(postfp,"\n 300 60 moveto");
      strcpy(str,"Near Optimal:    ");
      fixcolor_plot_postscript(COLOR_NEAR_OPTIMAL);
      strcat(str,num_string_int(points_plotted[2]));
      strcat(str," ,Overlap=");
      strcat(str,num_string_int(points_plotted[5]));
      strcat(str," within ");
      close_energy_formatted=(float)(worst_energy-optimal_energy);
      close_energy_formatted=close_energy_formatted*close_optimal_energy;
      close_energy_formatted=close_energy_formatted*.1;
      sprintf(close_to_optimal_string,"%.1f",close_energy_formatted);
      strcat(str,close_to_optimal_string);
      strcat(str," kcal/mol");
      fprintf(postfp,"\n ( %s ) show",str);
    }
  } else if (opt_distance>0) {
      fprintf(postfp,"\n 300 60 moveto");
      strcpy(str,"Within    ");
      strcat(str,num_string_int(opt_distance));
      strcat(str,"  of Optimal: ");
      fixcolor_plot_postscript(COLOR_NEAR_OPTIMAL);
      strcat(str,num_string_int(points_plotted[2]));
      fprintf(postfp,"\n ( %s ) show",str);
  }
  fprintf(postfp,"\n 300 45 moveto");
  if(opt_flag)
    strcpy(str,"Optimal Only: ");
  else
    strcpy(str,"Not Optimal:     ");
  strcat(str,num_string_int(points_plotted[3]));
  fixcolor_plot_postscript(COLOR_NOT_OPTIMAL);
  if(!opt_flag) {
    strcat(str," , Overlap=");
    strcat(str,num_string_int(points_plotted[4]));
  }
  fprintf(postfp,"\n ( %s ) show",str);
  if(opt_flag) {
    fixcolor_plot_postscript(COLOR_CT_MISSED_OPTIMAL);
    strcpy(str,"Other CT Base Pairs: ");
    strcat(str,num_string_int(points_plotted[7]));
    fprintf(postfp,"\n 300 30 moveto");
    fprintf(postfp,"\n ( %s ) show",str);
  }
  fprintf(postfp,"\n showpage\n");
  fprintf(postfp,"%sEOF\n","%%");
  fclose(postfp);
}
