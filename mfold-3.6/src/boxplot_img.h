/* boxplot_img.h */


#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
#include "strings_img.h"
#endif

void clean_title_img(char *title) {
  int i;
  int length;
  length = strlen(title);
  /* replace double with single quotes */
  for(i=0;i<length;i++) {
    if(title[i] == 34) 
      title[i] = 39;
  }
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
		  float y_pic_per_dot,int dis_l,int dis_t,
		  int color,int row,int column,int energy) {
  int x_pos;/* positions for center of dot */
  int y_pos;
  int ix_dot_size,iy_dot_size;
  float scale_adjust_x,scale_adjust_y;/* make large dots smaller */
  /* make small dots larger */
  int lx,rx,ty,by; /* positions for corners of dot */
  x_pos = img_col_to_x(dis_l,x_start,x_pic_per_dot,column);
  y_pos = img_row_to_y(dis_t,y_start,y_pic_per_dot,row);
  if((x_pic_per_dot<1.2)||(y_pic_per_dot<1.2))  { /* draw a dot here */ 
    gdImageSetPixel(image, x_pos, y_pos, color); 
  } else { /* draw a rectangle here 94% of the proper width */ 
    if(g_prob==TRUE) { /* adjust for value of prob */ 
      if (energy <= -1200000) 
	energy = -1200000; 
      scale_adjust_x = sqrt((-1.0*(float)energy)/1000000.); 
      scale_adjust_y = scale_adjust_x; 
    } else {
      if(x_pic_per_dot>3.) {
	scale_adjust_x=.85; 
      } else {
	scale_adjust_x=1; 
      }
      if(y_pic_per_dot>3.) {
	scale_adjust_y=.80; /* make big dots smaller */ 
      } else {
	scale_adjust_y=1.1; /* make small dots bigger */ 
      }
    }
    ix_dot_size=(int)(x_pic_per_dot*scale_adjust_x+.7); 
    if(ix_dot_size<1) 
      ix_dot_size=1; 
    iy_dot_size=(int)(y_pic_per_dot*scale_adjust_y+.7); 
    if(iy_dot_size<1) 
      iy_dot_size=1; 
    lx=x_pos-ix_dot_size/2; 
    rx=lx+ix_dot_size-1; 
    ty=y_pos-iy_dot_size/2; 
    by=ty+iy_dot_size-1; 
    gdImageFilledRectangle(image, lx, ty, rx, by, color); 
  }
} 

/* create the tic marks */
void img_single_hor_tic(int x_start,int j_col,float x_pic_per_dot,int ty,
 int by,int dis_l,int make_grid_lines) {
  int x_pos;
  x_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,j_col);
  if(make_grid_lines==TRUE)
    gdImageLine(image, x_pos, ty, x_pos, by - 2, color_gl);
  gdImageLine(image, x_pos, ty, x_pos, ty - 15, color_gl);
  gdImagePrintf(image, gdFontLarge, x_pos-10, ty-36, color_b0, "%d", j_col);
}

void img_plot_hor_tics(int x_start,float x_pic_per_dot,int dis_l, 
		       int dis_r,int ty,int by,int make_grid_lines) { 
  int j_col,step;
  int start;
  step=step_fun(dis_l,dis_r);
  start=start_fun(dis_l,step);
  img_single_hor_tic(x_start,dis_l,x_pic_per_dot,ty,by,dis_l,FALSE);
  img_single_hor_tic(x_start,dis_r,x_pic_per_dot,ty,by,dis_l,FALSE);
  for(j_col=start;j_col<(dis_r-2*step/3); j_col=j_col+step) {
    img_single_hor_tic(x_start,j_col,x_pic_per_dot,ty,by,dis_l,
		       make_grid_lines);
  }
}

void img_single_ver_tic(int y_start, int i_row, float y_pic_per_dot, int lx,
			int rx, int dis_t, int make_grid_lines) {
  int y_pos; 
  y_pos = img_row_to_y(dis_t,y_start,y_pic_per_dot,i_row);
  if(make_grid_lines==TRUE)
    gdImageLine(image, rx - 2, y_pos, lx + 2, y_pos, color_gl);
  gdImageLine(image, rx, y_pos, rx + 15, y_pos, color_b0);
  gdImagePrintf(image, gdFontLarge, rx+18, y_pos-12, color_b0, "%d", i_row);
}

void img_plot_ver_tics(int y_start,float y_pic_per_dot,int dis_t,
                       int dis_b,int lx,int rx,int make_grid_lines) {
  int i_row,step,start;
  step=step_fun(dis_t,dis_b);
  start=start_fun(dis_t,step);
  img_single_ver_tic(y_start,dis_t,y_pic_per_dot,lx,rx,dis_t,FALSE);
  img_single_ver_tic(y_start,dis_b,y_pic_per_dot,lx,rx,dis_t,FALSE);
  for(i_row=start;i_row<(dis_b-step/2);i_row=i_row+step) {
    img_single_ver_tic(y_start, i_row, y_pic_per_dot, lx, rx, dis_t,
		       make_grid_lines);  
  }
}

void  img_tics(int x_start,int y_start,float x_pic_per_dot,float y_pic_per_dot,
               int dis_l,int dis_r,int dis_t,int dis_b,int lx,
               int rx,int ty,int by,int make_grid_lines) {
  img_plot_ver_tics(y_start,y_pic_per_dot,dis_t,dis_b,lx,rx,make_grid_lines);
  img_plot_hor_tics(x_start,x_pic_per_dot,dis_l,dis_r,ty,by,make_grid_lines);
}

/* make the base pairs */
int img_points(int x_start, int y_start, float x_pic_per_dot, 
	       float y_pic_per_dot, int dis_l,int dis_r, int dis_t,int dis_b) {
  int row,column,new_row;
  int diag,position,location,new_color;
  int color;
  int count;
  int energy;
  int start_diag,end_diag;
  int points_plotted;
  start_diag = dis_t+dis_l-1;
  end_diag = dis_b+dis_r-1;
  points_plotted = 0; 
  for(color = g_number_of_colors;color>=1;color--) {
    for(diag = start_diag;diag<end_diag;diag++) {
      position = g_diag_start[diag];
      for(count = 0;count<g_diag_count[diag];count++) {
	location = position + count;
	if((!g_opt_prob_flag) || (g_diag[location].energy>=-20000)) {
	  if(g_diag[location].row<=dis_b) {
	    if((g_diag[location].energy==g_optimal_energy)
	       ||((g_diag[location].energy<=g_energy_cutoff)&&
		  (g_diag[location].length>=g_chain_len))) {
	      new_color=g_diag[location].color;
	      if(new_color==color) {
		new_row = g_diag[location].row;
		column = g_diag[location].column;
		for(row=new_row; row<(new_row+g_diag[location].length);
		    row++) {
		  if((row>=dis_t)&&(row<=dis_b)
		     &&(column>=dis_l)&&(column<= dis_r)) {
		    img_plot_dot(x_start,y_start,x_pic_per_dot,
				 y_pic_per_dot,dis_l,dis_t,*img_color[color],
				 row,column,g_diag[location].energy);
		    points_plotted++;
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
  /* plot dots below main diagonal */
  color = 1;
  for(diag=start_diag;diag<end_diag;diag++) {
    position = g_diag_start[diag];
    for(count=0;count<g_diag_count[diag];count++) {
      location = position + count;
      if(g_diag[position+count].row<=dis_r && (g_diag[location].color==1)) {
	location = position + count;
	new_row = g_diag[location].row;
	column = g_diag[location].column;
	if(column>=dis_t)  {
	  for(row=new_row; row<(new_row+g_diag[location].length); row++) {
	    if((column>=dis_t) && (column<=dis_b) && (row>=dis_l) 
	       && (row<= dis_r)) {
	      energy = g_diag[location].energy;
	      if(g_prob) {
		if(g_opt_prob_flag) {
		  if(energy<=-20001)
		    img_plot_dot(x_start, y_start, x_pic_per_dot,
				 y_pic_per_dot, dis_l, dis_t, color_OP,
				 column, row, -1000000);
		} else {
		  img_plot_dot(x_start, y_start, x_pic_per_dot, y_pic_per_dot,
			       dis_l, dis_t, color_OP, column, row, energy);
		}
	      } else {
		img_plot_dot(x_start, y_start, x_pic_per_dot, y_pic_per_dot,
			     dis_l, dis_t, color_OP, column, row, energy);
	      }
	    }
	    column--;
	  }
	}
      }
    }
  }
  return points_plotted;
}

/* make diagonal line */
void img_make_diag(int dis_l, int dis_r, int dis_t, int dis_b, int x_start,
		   int start_y, float x_pic_per_dot, float y_pic_per_dot,
		   int resolution) {
  int ty,rx;
  int lx,by;
  int right,left,top,bottom;
  if((!((dis_r>=dis_t)&&(dis_b>=dis_l)))) /* cases all below, all above */
    return; /* There is no diagonal line */
  /* find right end of line */
  if(dis_r>dis_b) { /* diagonal crosses bottom border */
    bottom = dis_b;
    right = bottom;
  } else { /* diagonal crosses right border */
    right = dis_r;
    bottom = right;
  }
  /* find left end of line */
  if(dis_l>dis_t) { /* diagonal crosses left side */
    left = dis_l;
    top = dis_l;
  } else { /* diagonal crosses top side */
    top = dis_t;
    left = top;
  }
  ty = img_row_to_y(dis_t, start_y, y_pic_per_dot, top) - 
    (int) (y_pic_per_dot/2. + .5);
  lx = img_row_to_y(dis_l, x_start, x_pic_per_dot, left) -
    (int) (x_pic_per_dot/2. + .51);
  rx = img_col_to_x(dis_l, x_start, x_pic_per_dot, right) +
    (int) (x_pic_per_dot/2.+.51);
  by = img_row_to_y(dis_t, start_y, y_pic_per_dot, bottom) +
    (int) (y_pic_per_dot/2. + .5);
  if(resolution<70) /* draw a narrow or wide line */
    gdImageLine(image, lx - 1, ty - 1, rx, by, color_b0);
  else
    gdImageLine(image, lx - 1, ty - 1, rx, by, gdBrushed);
} 

void img_create_colors(char *file_type,int grid_flag) {
  int color[TOTAL_COLORS][3];
  int i , j ;
  gdImagePtr bkim;
  /* convert colors to integers */
  for(j=0;j<TOTAL_COLORS;j++) {
    for(i=0;i<=2;i++)
      color[j][i]=(int)(255.*g_color[j][i]);
  }
  if(g_number_of_colors>=7) {
    for(j=1;j<=g_number_of_colors;j++) {
      for(i=0;i<=2;i++)
        color[j][i]=(int)(255.*g_prob_color[j][i]);
    }
  }
  bkim = gdImageCreate(2, 2);
  color_bak = gdImageColorAllocate(bkim, color[COLOR_BACKGROUND][0],
				   color[COLOR_BACKGROUND][1],
				   color[COLOR_BACKGROUND][2]);
  color_b0 = gdImageColorAllocate(bkim, color[0][0], color[0][1], color[0][2]);
  gdImageLine(bkim, 0, 0, 1, 0, color_b0);
  gdImageLine(bkim, 0, 1, 1, 1, color_b0);
  /* end of brush definition */
  color_bak = gdImageColorAllocate(image, color[COLOR_BACKGROUND][0],
				   color[COLOR_BACKGROUND][1],
				   color[COLOR_BACKGROUND][2]);
  /* first is background color */
  color_b0 = gdImageColorAllocate(image,color[0][0], color[0][1], color[0][2]);
  color_bk = gdImageColorAllocate(image,color[1][0], color[1][1], color[1][2]);
  /* special optimal color */
  color_OP = gdImageColorAllocate(image,
				  (int) (255. * g_color[1][0]), 
				  (int) (255. * g_color[1][0]), 
				  (int) (255. * g_color[1][0]));
  /* declare colors */
  color_re = gdImageColorAllocate(image, color[2][0],color[2][1],color[2][2]);
  color_gr = gdImageColorAllocate(image, color[3][0],color[3][1],color[3][2]);
  color_ye = gdImageColorAllocate(image, color[4][0],color[4][1],color[4][2]);
  color_pu = gdImageColorAllocate(image, color[5][0],color[5][1],color[5][2]);
  color_br = gdImageColorAllocate(image, color[6][0],color[6][1],color[6][2]);
  color_bl = gdImageColorAllocate(image, color[7][0],color[7][1],color[7][2]);
  color_fl = gdImageColorAllocate(image, color[8][0],color[8][1],color[8][2]);
  color_la = gdImageColorAllocate(image, color[COLOR_LABEL][0],
				  color[COLOR_LABEL][1],
				  color[COLOR_LABEL][2]);
  color_co = gdImageColorAllocate(image, color[COLOR_COPY][0],
				  color[COLOR_COPY][1],
				  color[COLOR_COPY][2]);
  if(grid_flag==TRUE)
    color_gl = gdImageColorAllocate(image,
				    color[COLOR_GRID][0],
				    color[COLOR_GRID][1],
				    color[COLOR_GRID][2]);
               /* Gray for grid lines */
  gdImageSetBrush(image, bkim);
  if(strcmp(file_type,"jpg")) {
    gdImageInterlace(image, 1);  /* $$$$$ turn on interlace */
  }
}

void img_label_points(int label_i, int label_j, int x_start, int y_start,
		      float x_pic_per_dot, float y_pic_per_dot, int dis_l,
		      int dis_t) {
  int hor_pos,ver_pos,label_hor_pos,lx,ty,rx,by;
  int hor_shift,ver_shift;
  char string[100];
  hor_shift = (int)(x_pic_per_dot/3);
  ver_shift = (int)(y_pic_per_dot/3);
  hor_pos = img_col_to_x(dis_l,x_start,x_pic_per_dot, label_j);
  ver_pos = img_row_to_y(dis_t,y_start,y_pic_per_dot, label_i);
  label_hor_pos = img_col_to_x(dis_l,x_start,x_pic_per_dot, label_j+1);
  lx = hor_pos-hor_shift;
  rx = hor_pos+hor_shift;
  ty = ver_pos-ver_shift;
  by = ver_pos+ver_shift;
  strcpy(string,num_string_int(label_i));
  strcat(string,",");
  strcat(string,num_string_int(label_j)); 
  gdImageFilledRectangle(image, lx, ty, rx, by, color_la);
  gdImagePrintf(image, gdFontSmall, label_hor_pos + 1, ver_pos - 5, color_la, 
		"(%s) ", string);
}

int open_img(char *filename) {
   if ((giffp = fopen(filename, "wb")) == NULL) {
     printf ("Error!\tCould not open file: %s\n", filename);
     return(1);
   }
   return (0);
}   

void make_imgdat(char *filename, int resolution, int x_start, int y_start,
		 float x_pic_per_dot, float y_pic_per_dot, int dis_l,
		 int dis_t) {
  strcat(filename, ".gifdat");
  printf("gifdat file is %s\n",filename);
  if( open_img(filename) ) {
    printf("Error creating .gifdat file %s.\n",filename);
    return;
  }
  fprintf(giffp, "%d %d %d %f %f %d %d %d \n", resolution, x_start, y_start,
	  x_pic_per_dot, y_pic_per_dot, g_length, dis_l, dis_t);
  fclose(giffp); 
}

void finish_making_img(int clear_flag, int resolution, char *name_for_img, 
		       char *title, int make_imgdat_flag, int make_label_flag,
		       int label_i, int label_j, int make_grid_flag, 
		       int png_mode, int jpg_mode, int mi_flag, 
		       int opt_prob_flag, int label_array_count,
		       int *label_row_array, int *label_col_array) {
  int hor_pos; /* resolution is 72, 110, or 200 ? */
  int color;
  int i;
  char file_type[6];
  char temp_filename[120];
  char gifeng_filename[120];
  int offset,offset_adjust;
  int dis_l,dis_r,dis_t,dis_b;/* local versions of display_l,display_r */
  char string[120],img_filename[120];
  int x_float; 
  /* use 8.5 x 11 in for image */
  int img_width; /* width of image in pixels */
  int img_height; /* height of image in pixels */
  int lx, ty, rx, by; /* top, bottom, right, left corner of rectangle  */
  float x_pic_per_dot, y_pic_per_dot;
  int x_start,y_start;
  int points_plotted;
  dis_l = display_l;
  dis_r = display_r;
  dis_t = display_t;
  dis_b = display_b;
  /* correct message for text window */
  strcpy(img_filename,name_for_img);
  if(strlen(img_filename)==0) {
    printf("Error!\nImage file name is empty.\n");
    return; /* abort making file */
  }
  img_width = (int)(8.5*resolution)-1;
  img_height = 11*resolution-1;
  strcpy(temp_filename, img_filename);
  chop_suffix(temp_filename, ".img"); 
  if(make_label_flag==TRUE) {
    strcpy(gifeng_filename, temp_filename);
    strcat(gifeng_filename, ".gifeng");
  }
  if(png_mode) {
    strcpy(file_type, "png");
  } else if(jpg_mode) {
    strcpy(file_type, "jpg");
  } else {
    strcpy(file_type, "gif");
  }
  if(open_img(img_filename)) printf("Error for %s file\n", file_type);
  /* hor by ver or x by y */
  image = gdImageCreate(img_width + 1, img_height + 1);
  /* declare image size in pixels */ 
  printf("Creating %s %d by %d\n",file_type,img_width+1,img_height+1);
  img_create_colors(file_type,make_grid_flag);
  if(clear_flag==TRUE)
    gdImageColorTransparent(image, color_bak);  /* make box  7 by 7 inches */
  ty = (int)(1.97*resolution); /* come down 3 inches */
  by = (int)(img_height-1.97*resolution);
  lx = (int)(.47*resolution); /* left side of box */ 
  rx = img_width - .95*resolution; /* 6.5 inches wide drawing region
				    * an extra .1 inch was used to
				    * make the line outside the region */ 
  if(resolution<70)
    gdImageRectangle(image, lx, ty, rx, by, color_b0);
  else { /* Make a wide line */
      gdImageLine(image, lx, ty, rx, ty, gdBrushed);
      gdImageLine(image, lx, by, rx, by, gdBrushed);
      gdImageLine(image, lx, ty, lx, by, gdBrushed);
      gdImageLine(image, rx, ty, rx, by, gdBrushed);
  }
  hor_pos = (int)(4.25*resolution - (strlen(title) + 2)/2*10);
  if(hor_pos<1)
     hor_pos=1;
  clean_title_img(title);
  gdImageString(image, gdFontLarge, hor_pos, resolution / 2,
		(unsigned char*) title, color_b0); 
  x_pic_per_dot = resolution*7.0/(dis_r-dis_l+1); /* width of base pair */
  y_pic_per_dot = resolution*7.0/(dis_b-dis_t+1); /* height of base pair */
  x_start=(int)(.5*resolution+x_pic_per_dot/2+.5); /* distance from
						      left edge for
						      first base pair*/  
  y_start=(int)(2.0*resolution+y_pic_per_dot/2+.5); /* distance from
						       top for first
			 			       base pair both
						       above are in
						       pixels create
						       the tic marks */  
  img_tics(x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
	   dis_r,dis_t,dis_b,lx,rx,ty,by,make_grid_flag);
 /* make the diagonal line */
  img_make_diag(dis_l, dis_r, dis_t, dis_b, x_start, y_start, x_pic_per_dot,
		y_pic_per_dot, resolution);
  /* both above are in pixels 
   * create the points */
  points_plotted = img_points(x_start,y_start,x_pic_per_dot,y_pic_per_dot,
			      dis_l,dis_r,dis_t,dis_b);
  gdImagePrintf(image,gdFontLarge,lx,by+15, color_OP,"Lower Triangle Shows"); 
  display_header_message_img(image, img_width, 0 );
  if(g_prob) {
    if(opt_prob_flag)
      gdImagePrintf(image, gdFontLarge, lx, by + 35, color_OP, 
		    "Optimal Structure");
    else
      gdImagePrintf(image, gdFontLarge, lx, by + 35, color_OP, 
		    "Highest Color Range");
  } else
    gdImagePrintf(image, gdFontLarge, lx, by + 35, color_OP, "Optimal Energy");
  gdImagePrintf(image, gdFontLarge, lx, by + 68, color_b0, "Upper Triangle");
  gdImagePrintf(image, gdFontLarge, lx, by + 91, color_b0, 
		"Base pairs Plotted: %d", points_plotted);
  if(g_prob!=TRUE) {
    strcpy(string, num_string_fancy_float(g_energy_cutoff-g_optimal_energy));
    strcat(string," kcal/mol");
    if(g_chain_len>1) {
      strcat(string,"     Filter: ");
      strcat(string,num_string_int(g_chain_len));
    }
    gdImagePrintf(image, gdFontLarge, resolution / 2 + 5, resolution, color_b0,
		  "deltaG in Plot File =%s ", string);
  } else {
    strcpy(string,num_string_fancy_float(g_energy_cutoff));
    if(mi_flag)
      strcat(string," <= M.I. <= 2");
    else
      strcat(string," <= Probability <= 1");
    gdImagePrintf(image, gdFontLarge, resolution / 2 + 5, resolution, color_b0,
		  " %s ", string);
  }
  /* display optimal energy */
  if(g_prob!=TRUE) {
    strcpy(string,num_string_fancy_int(g_optimal_energy));
    gdImagePrintf(image, gdFontLarge, (int) (4.2 * resolution) + 14, by + 3,
		  color_bk, "Optimal energy: %s ", string);
  }
  if(g_prob!=TRUE)
    strcpy(string, num_string_fancy_float(g_energy_cutoff-g_optimal_energy));
  else
    strcpy(string, num_string_fancy_float(g_energy_cutoff));
  if(g_prob!=TRUE)
    strcat(string," kcal/mol");
  if(g_chain_len>1) {
    strcat(string,"         Filter: ");
    strcat(string,num_string_int(g_chain_len));
  }
  /* print the color scheme */
 if(g_prob!=TRUE) {
   x_float=(float)g_optimal_energy;
   if((g_number_of_colors<6)||(resolution>90))
     offset_adjust=23;
   else
     offset_adjust=17;
   for(color=2,offset=22; color<=g_number_of_colors; 
       color++,offset+=offset_adjust) {
     strcpy(string,num_string_fancy_float(x_float+(float)(color-2)*
					  g_color_increment));
     strcat(string,"< energy <= ");
     strcat(string,num_string_fancy_float(x_float+(float)(color-1)*
					  g_color_increment));
     gdImagePrintf(image, gdFontLarge, (int) (4.2 * resolution), by + offset, 
		   *img_color[color], "%s ", string);
   }
 } else {
   strcpy(string,num_string_float(g__prob[g_number_of_colors-4][0]));
   strcat(string," < Prob. < 1");
   gdImagePrintf(image, gdFontLarge, (int) (4.2 * resolution), by + 6, 
		 color_bk, "%s ", string);
   if((g_number_of_colors<6)||(resolution>90))
     offset_adjust=23;
   else
     offset_adjust=17;
   for(color=2,offset=6+offset_adjust;color<=g_number_of_colors-1;
       color++,offset+=offset_adjust) {
     strcpy(string, num_string_float(g__prob[g_number_of_colors-4][color-1]));
     strcat(string," < Prob. <= ");
     strcat(string, num_string_float(g__prob[g_number_of_colors-4][color-2])); 
     gdImagePrintf(image, gdFontLarge, (int) (4.2 * resolution), by + offset, 
		   *img_color[color], "%s ", string);
   }
   strcpy(string,"           Prob. <= ");
   strcat(string,num_string_float(g__prob[g_number_of_colors-4]
				  [g_number_of_colors-2]));
   gdImagePrintf(image, gdFontLarge, (int)(4.2 * resolution), by + offset, 
		 *img_color[g_number_of_colors], "%s ", string);
 }
 if(make_label_flag==TRUE) {
   img_label_points(label_i, label_j, x_start, y_start, x_pic_per_dot,
		    y_pic_per_dot, dis_l, dis_t);
 }
 if(label_array_count>0) {
   for(i=0;i<label_array_count;i++)
     img_label_points(label_row_array[i], label_col_array[i], x_start, y_start,
		      x_pic_per_dot, y_pic_per_dot, dis_l, dis_t);
 }
 if (png_mode)
#if HAVE_LIBPNG
   gdImagePng(image, giffp)
#endif
     ;
 else if (jpg_mode)
#if HAVE_LIBJPEG
   gdImageJpeg(image, giffp, -1)
#endif
     ;
 else
#if defined(HAVE_LIBGD) && defined(HAVE_GDIMAGEGIF)
   gdImageGif(image, giffp)
#endif
     ;
 fclose(giffp);
 printf("%s file is %s.\n", file_type, img_filename);
 if(make_imgdat_flag==TRUE) {
   make_imgdat(temp_filename, resolution, x_start, y_start, x_pic_per_dot,
	       y_pic_per_dot, dis_l,dis_t);
 }
 if(make_label_flag==TRUE)
   make_energy_file(gifeng_filename, label_i, label_j);
} 
