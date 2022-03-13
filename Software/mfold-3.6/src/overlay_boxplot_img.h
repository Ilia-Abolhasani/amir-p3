#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
#include "strings_img.h"
#endif

int img_color[TOTAL_COLORS];

int open_img(char *img_tempname) {    /* open specified output file */
  if ((giffp = fopen(img_tempname, "w")) == NULL) { /* open a file */
    printf ("Could not open file %s.", img_tempname);
    printf ("This is a temporary file used to create an image file.\n");
    return(1);
  }
  return (0);
}

/* Create an image
 * convert a row to a position in pixels down from the top */
int img_row_to_y(int dis_t,int y_start,float y_pic_per_dot,int row) { 
   return y_start+(int)(y_pic_per_dot*((float)(row-dis_t))+.500001);
}

/* convert a column to position in pixels right from the left edge */
int img_col_to_x(int dis_l,int x_start,float x_pic_per_dot,int column) {
  return x_start+(int)(x_pic_per_dot*((float)(column-dis_l))+.50001);
}

/* draw a single dot */
void img_plot_dot(int x_start,int y_start, float x_pic_per_dot,
		  float y_pic_per_dot,int dis_l,int dis_t, int color,int row,
		  int column,int energy, float ct_adjust,int prob_flag) {
  int x_pos; /* center of dot */
  int y_pos;
  int ix_dot_size,iy_dot_size;
  float scale_adjust_x,scale_adjust_y; /* make large dots smaller
					* make small dots larger */
  int lx,rx,ty,by; /* positions for corners of dot */
  x_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,column);
  y_pos=img_row_to_y(dis_t,y_start,y_pic_per_dot,row);
  if((ct_adjust<=1.0)&&((x_pic_per_dot<1.2)||(y_pic_per_dot<1.2))) {
    /* draw a dot here */
    gdImageSetPixel(image, x_pos, y_pos, color);
  } else { /* draw a rectangle here 94% of the proper width */
    if(prob_flag==TRUE) { /* adjust size for value of probability */
      scale_adjust_x=sqrt((-1.0*(float)energy)/1000.)*ct_adjust;
      scale_adjust_y=scale_adjust_x;
    } else {
      scale_adjust_x=ct_adjust;
      scale_adjust_y=ct_adjust;
      if(ct_adjust<=1.0) { /* make big dots smaller */
	if(x_pic_per_dot>3.)
	  scale_adjust_x=.85*ct_adjust;
	if(y_pic_per_dot>3.)
	  scale_adjust_y=.85*ct_adjust;
      } else {
	if(x_pic_per_dot<=3.0)
	  scale_adjust_x=3.0*1.0/x_pic_per_dot;
	if(y_pic_per_dot<=3.0)
	  scale_adjust_y=3.0*1.0/y_pic_per_dot;
      }
    }
    ix_dot_size=(int)(x_pic_per_dot*scale_adjust_x+.7);
    if(ix_dot_size<1)
      ix_dot_size=1;
    iy_dot_size=(int)(y_pic_per_dot*scale_adjust_y+.7);
    if(iy_dot_size<1)
      iy_dot_size=1;
    if(ix_dot_size>100)
      ix_dot_size=100;
    if(iy_dot_size>100)
      iy_dot_size=100;
    lx=x_pos-ix_dot_size/2;
    rx=lx+ix_dot_size-1;
    ty=y_pos-iy_dot_size/2;
    by=ty+iy_dot_size-1;
    gdImageFilledRectangle(image, lx, ty, rx, by, color);
  } 
}

/* create the tic marks */
void img_single_hor_tic(int x_start,int j_col,float x_pic_per_dot,int ty,
                        int by,int dis_l,int grid_flag) {
  int x_pos;
  x_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,j_col);
  if(grid_flag==TRUE)
    gdImageLine(image, x_pos, ty, x_pos, by - 2, img_color[COLOR_GRID]);
  gdImageLine(image, x_pos, ty, x_pos, ty - 15, img_color[COLOR_TEXT]);
  gdImagePrintf(image,gdFontLarge,x_pos-10,ty-36,img_color[COLOR_TEXT],
		"%d",j_col);
}

void img_plot_hor_tics(int x_start,float x_pic_per_dot,int dis_l,
                       int dis_r,int ty,int by,int grid_flag) { 
  int j_col,step;
  int start;
  step=step_fun(dis_l,dis_r);
  start=start_fun(dis_l,step);
  img_single_hor_tic(x_start,dis_l,x_pic_per_dot,ty,by,dis_l,FALSE);
  img_single_hor_tic(x_start,dis_r,x_pic_per_dot,ty,by,dis_l,FALSE);
  for(j_col=start;j_col<(dis_r-2*step/3); j_col=j_col+step) {
    img_single_hor_tic(x_start,j_col,x_pic_per_dot,ty,by,dis_l, grid_flag);
  }
}

void img_single_ver_tic(int y_start,int i_row,float y_pic_per_dot,int lx,
			int rx,int dis_t,int grid_flag) {
  int y_pos; 
  y_pos=img_row_to_y(dis_t,y_start,y_pic_per_dot,i_row);
  if(grid_flag==TRUE)
    gdImageLine(image, rx - 2, y_pos, lx + 2, y_pos, img_color[COLOR_GRID]);
  gdImageLine(image, rx, y_pos, rx + 15, y_pos, img_color[COLOR_TEXT]);
  gdImagePrintf(image, gdFontLarge, rx + 18, y_pos - 12, img_color[COLOR_TEXT],
		"%d", i_row);
}

void img_plot_ver_tics(int y_start,float y_pic_per_dot,int dis_t,
                       int dis_b,int lx,int rx,int grid_flag) {
  int i_row,step,start;
  step=step_fun(dis_t,dis_b);
  start=start_fun(dis_t,step);
  img_single_ver_tic(y_start,dis_t,y_pic_per_dot,lx,rx,dis_t,FALSE);
  img_single_ver_tic(y_start,dis_b,y_pic_per_dot,lx,rx,dis_t,FALSE);
  for(i_row=start;i_row<(dis_b-step/2);i_row=i_row+step) {
    img_single_ver_tic(y_start,i_row,y_pic_per_dot,lx,rx,dis_t,
		       grid_flag);  
  }
}

void  img_tics(int x_start,int y_start,float x_pic_per_dot,float y_pic_per_dot,
               int dis_l,int dis_r,int dis_t,int dis_b,int lx,
               int rx,int ty,int by,int grid_flag) {
  img_plot_ver_tics(y_start,y_pic_per_dot,dis_t,dis_b,lx,rx,grid_flag);
  img_plot_hor_tics(x_start,x_pic_per_dot,dis_l,dis_r,ty,by,grid_flag);
}

  /* make the base pairs */
void  img_points(float dot_size[4],int x_start,int y_start,float x_pic_per_dot,
		 float y_pic_per_dot,int left,int right,int top,int bottom,
		 int points_plotted[8],int optimal_energy,int *diag_start,
		 int *diag_count,struct helix *diag,int energy_cutoff,
		 int chain_len,int number_of_colors,int plot_len,
		 int *ct_diag_start,int *ct_diag_count, 
		 struct ct_helix *ct_diag,int ct_len,float close_optimal_energy,
		 int worst_energy,int prob_flag) {
  int row,column,new_row,location,new_color,color,count,start_diag,end_diag;
  int plot_display_b,plot_display_r,ct_display_b,ct_display_r,type;
  int ct_type,energy,position,current_type,current_diag;
  float size_multiplier;
  start_diag=top+left-1;
  if(bottom>plot_len)
    plot_display_b=plot_len;
  else
    plot_display_b=bottom;
  if(right>plot_len)
    plot_display_r=plot_len;
  else
    plot_display_r=right;
  end_diag=plot_display_b+plot_display_r-1; 
  /* black,gray dots in upper triangle */
  for(color=number_of_colors;color>=1;color--) {
    for (current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=diag_start[current_diag];
      for (count=0;((count<diag_count[current_diag]) &&
		   (diag[position+count].row<=bottom));count++) {
	location=position+count;
	if ( (diag[location].energy<=energy_cutoff) &&
	   ( (diag[location].length>=chain_len) ||
	    (diag[location].energy==optimal_energy) ) ) {
	  new_color=diag[location].color;
	  if(new_color==color) {
	    new_row=diag[location].row;
	    column=diag[location].column;
	    for(row=new_row;
		row<(new_row+diag[location].length);row++) {
	      if((row>=top)&&(row<=bottom)&&(column>=left)&&(column<= right)) {
		img_plot_dot(x_start,y_start,x_pic_per_dot,y_pic_per_dot,left,
			     top,img_color[color],row,column,
			     diag[location].energy,1.0, prob_flag);
		points_plotted[0]++;
		if(color==1)
		  points_plotted[6]++;
	      }
	      column--;
	    }
	  }                        
	}
      }
    }
  }
  /* draw optimal dots in lower triangle */
  color=1;
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=diag_start[current_diag];
    for(count=0;((count<diag_count[current_diag]) &&
		 (diag[position+count].row<=right));count++) {
      location=position+count;
      if((diag[location].energy==optimal_energy)) {
	new_row=diag[location].row;
	column=diag[location].column;
	if(column>=top) {
	  for(row=new_row;row<(new_row+diag[location].length);row++) {
	    if((column>=top)&&(column<=bottom)&&(row>=left)&&(row<=right)) {
	      img_plot_dot(x_start,y_start,x_pic_per_dot,y_pic_per_dot,left,top,
			   img_color[1],column,row,diag[location].energy,1.0,
			   prob_flag);
	    }
	    column--;
	  }			     		    
	}
      }
    }
  }
  /* draw red,yellow,green dots */
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
  for(current_type=2;current_type>=0;current_type--) {
    size_multiplier=dot_size[current_type+1];
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=ct_diag_start[current_diag];
      for(count=0;( (count<ct_diag_count[current_diag])
		    &&(ct_diag[position+count].row<=bottom) ); count++) {
	location=position+count;
	new_row=ct_diag[location].row;
	column=ct_diag[location].column;
	type=set_ct_color(ct_diag[location].energy,
			  optimal_energy,worst_energy,close_optimal_energy);
	if(type==current_type) {
	  if(column>=left) {
	    for(row=new_row;row<new_row+ct_diag[location].length;row++) { 
	      if((row>=top)&&(row<=bottom)&&(column>=left)&&
		 (column<=right)){
		img_plot_dot(x_start,y_start,x_pic_per_dot,y_pic_per_dot,
			     left,top,img_color[current_type+COLOR_OPTIMAL],
			     row,column,-1001,size_multiplier,prob_flag);
		if(type==0)
		  points_plotted[1]++;
		else {
		  if(type==1)
		    points_plotted[2]++;
		  else
		    points_plotted[3]++;
		} 
	      }
	      column--;
	    }
	  }
	}
      }
    }
  }
  /* draw gray on top of red,green,yellow */ 
  for(color=4;color>=1;color--) {
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=ct_diag_start[current_diag];
      for(count=0;( (count<ct_diag_count[current_diag])
		    &&(ct_diag[position+count].row<=bottom));count++) {
	location=position+count;
	new_row=ct_diag[location].row;
	column=ct_diag[location].column;
	new_color=ct_diag[location].color;
	if(new_color==color) {
	  if(column>=left) {
	    energy=ct_diag[location].energy;
	    if(energy<=energy_cutoff) {
	      for(row=new_row;row<new_row+ct_diag[location].length;
		  row++) { 
		if((row>=top)&&(row<=bottom)&&(column>=left)
		   &&(column<=right)) {
		  if((energy==optimal_energy)||
		     (ct_diag[location].plotlen>=chain_len)) {
		    ct_type=
		      set_ct_color(energy,optimal_energy,
				   worst_energy,close_optimal_energy);
		    if((x_pic_per_dot>3)&&(y_pic_per_dot>3))
		      img_plot_dot(x_start,y_start,x_pic_per_dot,
				   y_pic_per_dot,left,top,
				   img_color[color],row,column,
				   optimal_energy,1.0,
				   prob_flag);
		    /* plot black on top of green,yellow,red */
		    if(ct_type==2)
		      points_plotted[4]++;/* red overlap */
		    else if(ct_type==1)
		      points_plotted[5]++;/* yellow overlap */
		  }
		}
		column--;
	      }
	    }
	  }
	}
      }
    }
  }
}

void opt_img_points(float dot_size[4],int x_start,int y_start,
		    float x_pic_per_dot,float y_pic_per_dot,int left,int right,
		    int top,int bottom,int points_plotted[8],
		    int optimal_energy,int *diag_start,int *diag_count,
		    struct helix *diag,int energy_cutoff,int chain_len,
		    int number_of_colors,int plot_len,int *ct_diag_start,
		    int *ct_diag_count,struct ct_helix *ct_diag,int ct_len,
		    int prob_flag,int opt_distance,int *ct_basepair) {
  int row,column,new_row,location,new_color,color,count,start_diag,end_diag;
  int plot_display_b,plot_display_r,ct_display_b,ct_display_r, energy;
  int position,current_diag;
  float size_multiplier;
  start_diag=top+left-1;
  if(bottom>plot_len)
    plot_display_b=plot_len;
  else
    plot_display_b=bottom;
  if(right>plot_len)
    plot_display_r=plot_len;
  else
    plot_display_r=right;
  end_diag=plot_display_b+plot_display_r-1; 
  /* black,gray dots in upper triangle */
  size_multiplier=dot_size[0]*1.11;
  for(color=number_of_colors;color>=2;color--) {
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=diag_start[current_diag];
      for(count=0;((count<diag_count[current_diag])&&
		   (diag[position+count].row<=bottom));count++) {
	location=position+count;
	if((diag[location].energy<=energy_cutoff)&&
	   (diag[location].length>=chain_len)) {
	  new_color=diag[location].color;
	  if(new_color==color) {
	    new_row=diag[location].row;
	    column=diag[location].column;
	    for(row=new_row;row<(new_row+diag[location].length);row++) {
	      if((row>=top)&&(row<=bottom)&&(column>=left)&&(column<= right)) {
		img_plot_dot(x_start,y_start,x_pic_per_dot,y_pic_per_dot,left,
			     top,img_color[color],row,column,
			     diag[location].energy,size_multiplier,prob_flag);
		points_plotted[0]++;
	      }
	      column--;
	    }
	  }                        
	}
      }
    }
  }
  /* draw optimal black dots in lower triangle */
  color=1;
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=diag_start[current_diag];
    for(count=0;((count<diag_count[current_diag])&&
		 (diag[position+count].row<=right));count++) {
      location=position+count;
      if((diag[location].energy==optimal_energy)) {
	new_row=diag[location].row;
	column=diag[location].column;
	if(column>=top) {
	  for(row=new_row;row<(new_row+diag[location].length);row++) {
	    if((column>=top)&&(column<=bottom)&&(row>=left)&&(row<=right)) {
	      img_plot_dot(x_start,y_start,x_pic_per_dot,y_pic_per_dot,left,
			   top,img_color[1],column,row,
			   diag[location].energy,1.0,prob_flag);
	    }
	    column--;
	  }			     		    
	}
      }
    }
  }
  /* Make cyan for missed ct */
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
  size_multiplier=dot_size[0]*.6;
  for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
    position=ct_diag_start[current_diag];
    for(count=0;( (count<ct_diag_count[current_diag])
		  &&(ct_diag[position+count].row<=bottom) ); count++) {
      location=position+count;
      new_row=ct_diag[location].row;
      column=ct_diag[location].column;
      energy=ct_diag[location].energy;
      if(energy!=optimal_energy) {
	if(column>=left) {
	  for(row=new_row;row<new_row+ct_diag[location].length; row++) { 
	    if((row>=top)&&(row<=bottom)&&(column>=left)&&(column<=right)) {
	      img_plot_dot(x_start,y_start,x_pic_per_dot,y_pic_per_dot,left,top,
			   img_color[COLOR_CT_MISSED_OPTIMAL],row,column,-1001,
			   size_multiplier,prob_flag);
	      points_plotted[7]++;
	    }
	    column--;
	  }
	}
      }
    }
  }
  /* Draw red, yellow, green optimal dots */
  end_diag=plot_display_b+plot_display_r-1;
  for(color=2;color>=0;color--) {
    size_multiplier=dot_size[color+1];
    for(current_diag=start_diag;current_diag<=end_diag;current_diag++) {
      position=diag_start[current_diag];
      for(count=0;((count<diag_count[current_diag])&&
		   (diag[position+count].row<=bottom));count++) {
	location=position+count;
	if(diag[location].energy==optimal_energy) {
	  new_row=diag[location].row;
	  column=diag[location].column;
	  for(row=new_row;row<(new_row+diag[location].length);row++) {
	    if((row>=top)&&(row<=bottom)&&(column>=left)&&(column<= right)) {
	      if(row<=ct_len) {
		if(ct_basepair[row]==column) {
		  new_color=0;
		}
		else {
		  if(distance_plot_ct(ct_len,row,column,opt_distance)) {
		    new_color=1;
		  } else {
		    new_color=2;
		  }
		}
	      } else {
		new_color=2;
	      }
	      if(new_color==color) {
		points_plotted[1+color]++; 
		img_plot_dot(x_start,y_start,x_pic_per_dot,
			     y_pic_per_dot,left,top,
			     img_color[color + COLOR_OPTIMAL],row,column,
			     diag[location].energy,size_multiplier,
			     prob_flag);
	      }
	    }
	    column--;
	  }
	}                        
      }
    }
  }
  points_plotted[6]=points_plotted[1]+points_plotted[2]+points_plotted[3];
  points_plotted[0]=points_plotted[0]+points_plotted[6];
}

/* make diagonal line */
void img_make_diag(int display_l,int display_r,int display_t,int display_b,
		   int  x_start,int start_y,float x_pic_per_dot,
		   float y_pic_per_dot) {
  int ty,rx,lx,by,left,top,bottom,right;
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
    ty=img_row_to_y(display_t,start_y,y_pic_per_dot,top) -
      (int) (y_pic_per_dot/2.+.51);
    lx=img_col_to_x(display_l,x_start,x_pic_per_dot,left) -
      (int) (x_pic_per_dot/2.+.51);
    rx=img_col_to_x(display_l,x_start,x_pic_per_dot,right) +
      (int) (x_pic_per_dot/2.+.51);
    by=img_row_to_y(display_t,start_y,y_pic_per_dot,bottom) +
      (int) (y_pic_per_dot/2.+.51);    
    gdImageLine(image, rx, by, lx, ty, gdBrushed);
    /* The above draws a wide diagonal line */
  }
} 

void img_label_points(int x_start,int y_start,float x_pic_per_dot,
		      float y_pic_per_dot,int dis_l,int dis_t,
		      int zoom_labels_row[50],int zoom_labels_col[50],
		      int zoom_labels_count) {
  int i,hor_pos,ver_pos,label_hor_pos,lx,ty,rx,by,hor_shift,ver_shift;
  char string[80];
  hor_shift=(int)(x_pic_per_dot/3);
  ver_shift=(int)(y_pic_per_dot/3);
  for(i=0;i<zoom_labels_count;i++) {
    hor_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,
			 zoom_labels_col[i]);
    ver_pos=img_row_to_y(dis_t,y_start,y_pic_per_dot,
			 zoom_labels_row[i]);
    label_hor_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,
			       zoom_labels_col[i]+1);
    lx=hor_pos-hor_shift;
    rx=hor_pos+hor_shift;
    ty=ver_pos-ver_shift;
    by=ver_pos+ver_shift;
    strcpy(string,num_string_int(zoom_labels_row[i]));
    strcat(string,",");
    strcat(string,num_string_int(zoom_labels_col[i]));
    gdImageFilledRectangle(image, lx, ty, rx, by, img_color[COLOR_LABEL]);
    gdImagePrintf(image, gdFontSmall, label_hor_pos + 1, ver_pos - 5, 
		  img_color[COLOR_LABEL], "<%s>", string);
  }
}

void img_create_colors(float color_set[TOTAL_COLORS][3]) {
  int color[TOTAL_COLORS][3],i,j;
  gdImagePtr bkim;
  /* convert colors to integers */
  for(j=0;j<TOTAL_COLORS;j++) {
    for(i=0;i<=2;i++)
      color[j][i]=(int)(255.*color_set[j][i]);
  }
  bkim = gdImageCreate(2, 2);
  img_color[8] = gdImageColorAllocate(bkim, color[8][0], color[8][1], 
				      color[8][2]);
  img_color[0] = gdImageColorAllocate(bkim, color[0][0], color[0][1], 
				      color[0][2]);
  gdImageLine(bkim, 0, 0, 1, 0, img_color[0]);
  gdImageLine(bkim, 0, 1, 1, 1, img_color[8]);
  /* end of brush definition */
  img_color[8] = gdImageColorAllocate(image, color[8][0], color[8][1], 
				      color[8][2]);
  /* first is background color */
  img_color[0] = gdImageColorAllocate(image, color[0][0], color[0][1], 
				      color[0][2]);
  img_color[1] = gdImageColorAllocate(image, color[1][0], color[1][1], 
				      color[1][2]);
  /* declare colors */
  for(i=2;i<TOTAL_COLORS;i++)
    img_color[i] = gdImageColorAllocate(image, color[i][0], color[i][1], 
					color[i][2]);
  gdImageSetBrush(image, bkim); /* Sets the brush  for wide lines*/
  gdImageInterlace(image, 1); /* $$$$$ turn on interlace */
}

void finish_making_img(int opt_flag, int clear_background_flag,
		       char *img_filename, char *plot_filename,
		       char *ct_filename, int optimal_energy, 
		       int energy_cutoff, int zoom_labels_count, int prob_flag,
		       int chain_len, int number_of_colors, 
		       float color_increment, float prob[4][6],
		       int zoom_labels_row[50], int zoom_labels_col[50],
		       int left, int right, int top, int bottom,
		       int *diag_start, int *diag_count, struct helix *diag,
		       int grid_flag, int plot_len, int ct_len, 
		       float dot_size[4], int *ct_diag_start, 
		       int *ct_diag_count, struct ct_helix *ct_diag,
		       char *sequence_name, int sequence_name_set,
		       int resolution, float close_optimal_energy, 
		       int worst_energy, float color_set[TOTAL_COLORS][3], 
		       int opt_distance, int *ct_basepair) {
  /* 0 for dot_plot, 1 for green, 2 for yellow 3 for red, 4 for
   * overlap red, 5 for yellow overlap 6 for total optimal dots 
   */
  int points_plotted[8]; 
  int error; /* resolution is 72, 110, or 200 ? */
  int dis_l,dis_r,dis_t,dis_b; /* local versions of display_l,display_r */
  char string[80];
  char opt_string[15];
  int x_float;
  int color,offset,offset_adjust; 
  int color_index;
  /* use 8.5 wide by 11 inches wide gif ? */
  int img_width; /* width of image in pixels */
  int img_height; /* height of image in pixels */
  int lx,ty,rx,by; /* top,bottom,right,left corn of rec */
  float x_pic_per_dot,y_pic_per_dot;
  float close_energy_formatted;
  int x_start,y_start;
  char img_tempname[120];
  int i;
  for(i=0;i<8;i++)
    points_plotted[i]=0;
  dis_l=left;
  dis_r=right;
  dis_t=top;
  dis_b=bottom;
  /* correct message for text window */
  img_width = (int)(8.5*resolution) - 1;
  img_height = 11*resolution-1; 
  strcpy(img_tempname,img_filename);
  error=open_img(img_filename);
  if(error) {
    printf("Error creating image file.\n");
    return;
  }                         /* hor by ver or x by y */
  image = gdImageCreate(img_width + 1, img_height + 1);
  /* declare image size in pixels */ 
  printf("Creating image %d by %d\n",img_width+1,img_height+1);
  img_create_colors(color_set); /* declare gif colors */
  if(clear_background_flag==TRUE)
    gdImageColorTransparent(image, img_color[COLOR_BACK]);/* make box 7x7 in */
  ty = (int) (1.97*resolution); /* come down 3 inches */
  by = (int) (img_height-1.97*resolution);
  lx = (int) (.47*resolution); /* left side of box */ 
  rx = img_width - (.95*resolution);
  /* 6.5 inches wide drawing region
   * an extra .1 inch was used to make the line outside the region
   * Using a wide line too 
   */
  display_header_message_img(image, img_width, 0);
  gdImageLine(image, lx, ty, rx, ty, gdBrushed);
  gdImageLine(image, lx, by, rx, by, gdBrushed);
  gdImageLine(image, lx, ty, lx, by, gdBrushed);
  gdImageLine(image, rx, ty, rx, by, gdBrushed);
  if(sequence_name_set) {
    gdImageStringCenter(image, gdFontGiant, img_width/2, (int) resolution/2,
			sequence_name, img_color[COLOR_TEXT]);
  } else {
    if (prob_flag==TRUE) {
      gdImagePrintfCenter(image, gdFontGiant, img_width/2, (int)resolution/2,
			  img_color[COLOR_TEXT], 
			  "Probability dotplot for %s and %s",
			  plot_filename, ct_filename);
    } else {  
      gdImagePrintfCenter(image, gdFontGiant, img_width/2, (int) resolution/2,
		    img_color[COLOR_TEXT], "Energy dotplot for %s and %s", 
		    plot_filename, ct_filename);
    }
  }
  x_pic_per_dot=resolution*7.0/(dis_r-dis_l+1);
  /* width of base pair*/
  y_pic_per_dot=resolution*7.0/(dis_b-dis_t+1); /* height of base pair*/
  x_start=(int)(.5*resolution+x_pic_per_dot/2+.5);
  /* distance from left edge for first base pair */
  y_start=(int)(2.0*resolution+y_pic_per_dot/2+.5);
  /* distance from top for first base pair both above are in pixels
   * create the tic marks 
   */
  img_tics(x_start, y_start, x_pic_per_dot, y_pic_per_dot, dis_l,
	   dis_r,dis_t,dis_b,lx,rx,ty,by,grid_flag);
  /* make the diagonal line */
  img_make_diag(dis_l, dis_r, dis_t, dis_b, x_start, y_start, x_pic_per_dot,
		y_pic_per_dot);
  /* create the points */
  if(opt_flag)
    opt_img_points(dot_size,x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
		   dis_r,dis_t,dis_b,points_plotted,optimal_energy,diag_start,
		   diag_count,diag,energy_cutoff,chain_len,number_of_colors,
		   plot_len,ct_diag_start,ct_diag_count,ct_diag,ct_len,
		   prob_flag,opt_distance,ct_basepair);
  else
    img_points(dot_size,x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
	       dis_r,dis_t,dis_b,points_plotted,optimal_energy,diag_start,
	       diag_count,diag,energy_cutoff,chain_len,number_of_colors,
	       plot_len,ct_diag_start,ct_diag_count,ct_diag,ct_len,
	       close_optimal_energy, worst_energy,prob_flag);
  /* State number of points plotted */
  gdImagePrintf(image, gdFontLarge, lx, by + 45, img_color[COLOR_TEXT], 
		"Base pairs Plotted:  %d", points_plotted[0]);
  if(opt_flag)
    strcpy(opt_string,"Optimal & CT");
  else
    strcpy(opt_string,"Optimal");
  gdImagePrintf(image, gdFontLarge, lx, by + 64, img_color[COLOR_OPTIMAL], 
		"%s: %d of ", opt_string, points_plotted[1]);
  gdImagePrintf(image, gdFontLarge, lx + 200, by + 64, 
		img_color[COLOR_PLOT_OPTIMAL], "%d", points_plotted[6]);
  if(!opt_flag) {
    if(close_optimal_energy!=0.0) {
      close_energy_formatted=(float)(worst_energy-optimal_energy);
      close_energy_formatted=close_energy_formatted*close_optimal_energy;
      close_energy_formatted=close_energy_formatted*.1;
      gdImagePrintf(image,gdFontLarge,lx,by+83,img_color[COLOR_NEAR_OPTIMAL], 
		    "Near Optimal: %d, Overlap=%d within %.1f kcal/mol", 
		    points_plotted[2],points_plotted[5],
		    close_energy_formatted);
    }
    gdImagePrintf(image,gdFontLarge,lx,by+102,img_color[COLOR_NOT_OPTIMAL],
		  "Not Optimal: %d, Overlap=%d", points_plotted[3], 
		  points_plotted[4]);
  } else {
    if(opt_distance>0) {
      gdImagePrintf(image,gdFontLarge,lx,by+83,img_color[COLOR_NEAR_OPTIMAL],
		    "Near Optimal: %d  within %d ",points_plotted[2],
		    opt_distance);
    }
    gdImagePrintf(image,gdFontLarge,lx,by+102,img_color[COLOR_NOT_OPTIMAL],
		  "Optimal Only: %d ", points_plotted[3]);
    gdImagePrintf(image,gdFontLarge,lx+200,by+102,
		  img_color[COLOR_CT_MISSED_OPTIMAL], "CT Only: %d ", 
		  points_plotted[7]);
  }
  /* display energy cutoff and filter*/
  if(!prob_flag) {
    strcpy(string,
	   num_string_fancy_float(energy_cutoff-optimal_energy,prob_flag));
    strcat(string," kcal/mol");
    if(chain_len>1) {
      strcat(string,"         Filter: ");
      strcat(string,num_string_int(chain_len));
    }
    gdImagePrintf(image,gdFontLarge,resolution/2+5,resolution, 
		  img_color[COLOR_TEXT], " Energy Increment:%s", string);
  } else {
    strcpy(string,num_string_fancy_float(energy_cutoff,prob_flag));
    strcat(string," <= Probability <= ");
    strcat(string, num_string_fancy_float(optimal_energy,prob_flag));
    gdImagePrintf(image, gdFontLarge, resolution / 2 + 5, resolution, 
		  img_color[COLOR_TEXT], " %s ", string);
  }
  /* display optimal energy */
  if(!prob_flag) {
    strcpy(string,num_string_fancy_int(optimal_energy,prob_flag));
    gdImagePrintf(image, gdFontLarge, lx, by + 20, img_color[1], 
		  "Optimal Energy: %s ", string);
  }
  /* display color map */
  if(!prob_flag) {
    x_float=(float)optimal_energy;
    for(color=2,offset=12;color<=number_of_colors;color++,offset+=22) {
      strcpy(string,num_string_fancy_float(x_float+(float)(color-2)*
					   color_increment,prob_flag));
      strcat(string,"< energy <=");
      strcat(string,num_string_fancy_float(x_float+(float)(color-1)*
					   color_increment,prob_flag));
      gdImagePrintf(image, gdFontLarge, (int) (4.5 * resolution), by + offset,
		    img_color[color], "%s ", string);
    }
  } else {
    strcpy(string, num_string_float(prob[number_of_colors-4][0]));
    strcat(string," < Prob. ");
    gdImagePrintf(image, gdFontLarge,(int) (4.2*resolution), by+6,img_color[8],
		  "%s ", string);
    color_index = img_color[7];
    if((number_of_colors<6)||(resolution>90))
      offset_adjust=23;
    else
      offset_adjust=20;
    for(color=2,offset=6+offset_adjust;color<=number_of_colors-1;
	color++,offset+=offset_adjust) {
      strcpy(string, num_string_float(prob[number_of_colors-4][color-1]));
      strcat(string," < Prob. <=");
      strcat(string, num_string_float(prob[number_of_colors-4][color-2])); 
      gdImagePrintf(image, gdFontLarge, (int) (4.2 * resolution), 
		    by + offset, color_index, "%s ", string);
      color_index = img_color[color + 1];
    }
    strcpy(string,"        Prob. <=");
    strcat(string,
	   num_string_float(prob[number_of_colors-4][number_of_colors-2]));
    gdImagePrintf(image,gdFontLarge,(int)(4.2*resolution),
		  by+offset,color_index, "%s ", string);
  }
  /* make and label the labeled dots */
  if(zoom_labels_count>0)
    img_label_points(x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
		     dis_t,zoom_labels_row,zoom_labels_col,zoom_labels_count);
#if defined(HAVE_LIBGD) && defined(HAVE_GDIMAGEGIF)
  gdImageGif(image, giffp);
#endif
}
