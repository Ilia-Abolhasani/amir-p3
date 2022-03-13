#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#if defined(HAVE_LIBGD) && defined(HAVE_GD_H)
#include "strings_img.h"
#endif

FILE *gifdatfp;

int img_color[COLOR_TABLE_SIZE];

/* create image output for ct_boxplot */

int open_img(char *image_file) {    /* open specified output file */
  if ((giffp = fopen(image_file, "w")) == NULL) {
    printf ("Error!\tCould not open image file %s.\n", image_file);
    return(1);
  }
  return (0);
}

/* create an image */

/* convert a row to a position in pixels down from the top */
int img_row_to_y(int dis_t,int y_start,float y_pic_per_dot,int row) { 
  return y_start+(int)(y_pic_per_dot*((float)(row-dis_t))+.500001);
}

/* convert a column to position in pixels right from the left edge */
int img_col_to_x(int dis_l,int x_start,float x_pic_per_dot,int column) {
  return x_start+(int)(x_pic_per_dot*((float)(column-dis_l))+.50001);
} 

/* draw a single dot */
void img_plot_dot(int x_start, int y_start, float x_pic_per_dot, 
		  float y_pic_per_dot, int dis_l, int dis_t, int color_index,
		  int color, int row, int column, float dot_magnifier,
                  int stripe_flag, int last_structure) {
  int x_pos; /* positions for center of dot */
  int y_pos;
  int ix_dot_size,iy_dot_size;
  int upper_left_x,upper_left_y,lower_right_x,lower_right_y;
  float adjuster;
  float scale_adjust_x,scale_adjust_y;/* make large dots smaller */
  /* make small dots larger */
  int lx,rx,ty,by; /* positions for corners of dot */
  gdPoint points[4];
  x_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,column);
  y_pos=img_row_to_y(dis_t,y_start,y_pic_per_dot,row);
  /* draw a rectangle here */
  scale_adjust_x=dot_magnifier;
  scale_adjust_y=dot_magnifier; /* make small dots bigger */
  ix_dot_size=(int)(x_pic_per_dot*scale_adjust_x+.7);
  iy_dot_size=(int)(y_pic_per_dot*scale_adjust_y+.7);
  if(ix_dot_size<2)
    ix_dot_size=2;
  if(iy_dot_size<2)
    iy_dot_size=2;
  lx=(int)(x_pos-ix_dot_size/2+.5);
  rx=(int)(lx+ix_dot_size-1+.5);
  ty=(int)(y_pos-iy_dot_size/2+.5);
  by=(int)(ty+iy_dot_size-1+.5);
  if(!stripe_flag)
    gdImageFilledRectangle(image, lx, ty, rx, by, color_index);
  else {
    adjuster=(color-1)*1./(float)(last_structure-1);
    upper_left_x=(int)(x_pos-adjuster*ix_dot_size/2+.5);
    upper_left_y=(int)(y_pos-adjuster*iy_dot_size/2+.5);
    lower_right_x=(int)(x_pos+adjuster*ix_dot_size/2+.5-1.);
    lower_right_y=(int)(y_pos+adjuster*iy_dot_size/2+.5-1.);
    points[0].x = upper_left_x; points[0].y = upper_left_y;
    points[1].x = rx; points[1].y = ty;
    points[2].x = lower_right_x; points[2].y = lower_right_y;
    points[3].x = lx; points[3].y = by;
    gdImageFilledPolygon(image, points, 4, color_index);
  }
}

/* create the tic marks */
void img_single_hor_tic(int x_start,int j_col,float x_pic_per_dot,int ty,
                        int by,int dis_l,int grid_flag) {
  int x_pos;
  x_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot,j_col);
  if(grid_flag==TRUE)
    gdImageLine(image, x_pos, ty, x_pos, by - 2, img_color[COL_GRID]);
  gdImageLine(image, x_pos, ty, x_pos, ty - 15, img_color[COL_TEXT]);
  gdImagePrintf(image, gdFontLarge, x_pos - 10, ty - 36, img_color[COL_TEXT], 
		"%d", j_col);
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
    img_single_hor_tic(x_start,j_col,x_pic_per_dot,ty,by,dis_l,
		       grid_flag);
  }
}

void img_single_ver_tic(int y_start,int i_row,float y_pic_per_dot,int lx,
			int rx, int dis_t,int grid_flag) {
  int y_pos; 
  y_pos=img_row_to_y(dis_t,y_start,y_pic_per_dot,i_row);
  if(grid_flag==TRUE)
    gdImageLine(image, rx - 2, y_pos, lx + 2, y_pos, img_color[COL_GRID]);
  gdImageLine(image, rx, y_pos, rx + 15, y_pos, img_color[COL_TEXT]);
  gdImagePrintf(image, gdFontLarge, rx + 18, y_pos - 12, img_color[COL_TEXT], 
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
    img_single_ver_tic(y_start,i_row,y_pic_per_dot,lx,rx,dis_t, grid_flag);  
  }
}

void  img_tics(int x_start,int y_start,float x_pic_per_dot,float y_pic_per_dot,
               int dis_l,int dis_r,int dis_t,int dis_b,int lx,
               int rx,int ty,int by,int grid_flag) {
  img_plot_ver_tics(y_start,y_pic_per_dot,dis_t,dis_b,lx,rx,grid_flag);
  img_plot_hor_tics(x_start,x_pic_per_dot,dis_l,dis_r,ty,by,grid_flag);
}

void  img_points(int counter[MAXIMUM_CT_FILES+1], int stripe_flag,
		 int x_start, int y_start, float x_pic_per_dot,
		 float y_pic_per_dot, int dis_l, int dis_r, int dis_t,
		 int dis_b, int last_structure, struct helix *diag,
		 int lower_structure, int *diag_start, int *diag_count,
		 float dot_magnifier, int nongray_for_3_flag,
		 int *gray3_content, int gray3[3]) {
  int row,column,n_column,new_row;
  int location,new_color;
  int color,starting_color,ending_color;
  int diagonal,start_diag,end_diag,position,count;
  int old_color;
  int color_index;
  int gray3_column;
  start_diag=dis_t+dis_l-1;
  end_diag=dis_b+dis_r-1;
  /* single structures */
  gray3[0]=0;
  gray3[1]=0;
  gray3[2]=0;
  for(color=last_structure;color>=2;color--) {
    color_index = img_color[color];
    for(diagonal=start_diag;diagonal<=end_diag;diagonal++) {
      position=diag_start[diagonal];
      for(count=0;((count<diag_count[diagonal])&&
		   (diag[position+count].row<=dis_b) ) ;count++) {
	location=position+count;
	new_color=diag[location].color;
	if(new_color==color) {
	  new_row=diag[location].row;
	  column=diag[location].column;
	  for(row=new_row; row<(new_row+diag[location].length);row++) {
	    if((row>=dis_t) && (row<=dis_b) && (column>=dis_l) && 
	       (column<= dis_r)) {
	      img_plot_dot(x_start, y_start, x_pic_per_dot, y_pic_per_dot,
			   dis_l, dis_t, color_index, color, row, column,
			   dot_magnifier, stripe_flag, last_structure);
	      counter[color]++;  
	    }
	    column--;
	  }
	}                        
      }    
    }
  }
  /* black or gray dots above main diagonal */
  if((last_structure>3)&&(!stripe_flag))
    starting_color=1;
  else
    starting_color=0;
  for(color=starting_color;color>=0;color--) {
    color_index = img_color[color];
    if((nongray_for_3_flag)&&(color==1))
      color_index = img_color[5];
    for(diagonal=start_diag;diagonal<=end_diag;diagonal++) {
      position=diag_start[diagonal];
      for(count=0;((count<diag_count[diagonal])&&
		   (diag[position+count].row<=dis_b) ) ;count++) {
	location=position+count;
	new_color=diag[location].color;
	if(new_color==color) {
	  new_row=diag[location].row;
	  column=diag[location].column;
	  for(row=new_row; row<(new_row+diag[location].length);row++) {
	    if((row>=dis_t)&&(row<=dis_b)&&(column>=dis_l)&&(column<= dis_r)) {
	      img_plot_dot(x_start, y_start, x_pic_per_dot, y_pic_per_dot,
			   dis_l, dis_t, color_index, color, row, column,
			   dot_magnifier, FALSE, last_structure); 
	      counter[color]++;
	    }
	    column--;
	  }
	}                      
      }    
    }
  }
  /* plot partial overlap for 3 sequences for gray interpreted */
  if(nongray_for_3_flag) {
    for(row=dis_t;row<=dis_b;row++) {
      column = gray3_content[row];
      if(column!=0) {
	if(column<0)
          n_column=-1*column;
        else
          n_column=column;
        if((n_column>=dis_l)&&(n_column<=dis_r)) { 
	  if(column>0) { /* for seq 2 nd 3 */
	    color_index = img_color[7];
	    gray3[2]++;
	  } else { /* for seq 1 & 3 */
	    color_index = img_color[6];
	    gray3[1]++;           
	  }
	  img_plot_dot(x_start, y_start, x_pic_per_dot, y_pic_per_dot, dis_l,
		       dis_t, color_index, 1, row, n_column, dot_magnifier,
		       FALSE, last_structure);
	}
      }
    }
    gray3[0]=counter[1]-gray3[1]-gray3[2];
  }
  /* subdiagonal
   * points below the main diagonal
   * draw dots on lower diagonal */
  if (lower_structure != -1) { /* -1 is no dots below */
    if (lower_structure != -10) {
      if (lower_structure == -2) { /* -2 is black and gray dots below */
	starting_color = 1;
	ending_color = 0;
      } else if (lower_structure == -12) {
	starting_color=0;
	ending_color=0;
      } else {
	starting_color=lower_structure; /* single structure color dot below */
	ending_color=lower_structure;
      }
      for(color=starting_color;color>=ending_color;color--) {
	color_index = img_color[color];
	for(diagonal=start_diag;diagonal<=end_diag;diagonal++) {
	  position = diag_start[diagonal];
	  for(count=0;((count<diag_count[diagonal])&&
		       (diag[position+count].row<=dis_r) ) ;count++) {
	    location=position+count;
	    new_color=diag[location].color;
	    if(new_color==color) {
	      new_row=diag[location].row;
	      column=diag[location].column;
	      if(column>=dis_t) {
		for(row=new_row; row<(new_row+diag[location].length);row++) {
		  if((column>=dis_t) && (column<=dis_b) && (row>=dis_l) && 
		     (row<= dis_r)) {
		    img_plot_dot(x_start, y_start, x_pic_per_dot,y_pic_per_dot,
				 dis_l, dis_t, color_index, color, column, row,
				 dot_magnifier, FALSE, last_structure);
		  }
		  column--;
		}
	      }
	    }                        
	  }    
	}
      }
    }
    if(lower_structure<=-10) {
      color=1; /* take care of gray interpreted in lower triangle */
      for(diagonal=start_diag;diagonal<=end_diag;diagonal++) {
	position=diag_start[diagonal];
	for(count=0;((count<diag_count[diagonal]) &&
		     (diag[position+count].row<=dis_r) ) ;count++) {
	  location = position + count;
	  new_color = diag[location].color;
	  if(new_color==color) {
	    new_row = diag[location].row;
	    column = diag[location].column;
	    if(column >= dis_t) {
	      for (row=new_row;row<(new_row+diag[location].length);row++) {
		if ((column>=dis_t) && (column<=dis_b) && (row>=dis_l) &&
		   (row<= dis_r)) {
		  gray3_column = gray3_content[row]; 
		  if (gray3_column == 0) {
		    color_index = img_color[5];
		    old_color=5;
		  } else if (gray3_column<0) {
		    color_index = img_color[6];
		    old_color=6;
		  } else {
		    color_index = img_color[7];
		    old_color=7;
		  }
		  if(!nongray_for_3_flag)
		    gray3[old_color-5]++;
		  img_plot_dot(x_start, y_start, x_pic_per_dot, y_pic_per_dot,
			       dis_l, dis_t, color_index, color, column,row,
			       dot_magnifier, FALSE, last_structure);
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

/* make diagonal line */
void img_make_diag(int display_l, int display_r, int display_t, int display_b,
		   int x_start, int start_y, float x_pic_per_dot, 
		   float y_pic_per_dot) {
  int ty,rx;
  int lx,by;
  int right,top,bottom,left;
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
    
    ty=img_row_to_y(display_t,start_y,y_pic_per_dot,top)-
              (int)(y_pic_per_dot/2.+.51);
    lx=img_row_to_y(display_l,x_start,x_pic_per_dot,left)-
              (int)(x_pic_per_dot/2.+.51);
    rx=img_col_to_x(display_l,x_start,x_pic_per_dot,right)+
               (int)(x_pic_per_dot/2.+.51);
    by=img_row_to_y(display_t,start_y,y_pic_per_dot,bottom)+
               (int)(y_pic_per_dot/2.+.51);
    gdImageLine(image, rx - 1, by, lx - 1, ty, gdBrushed);
      /* The above draws a wide diagonal line */
   }
} 

void img_label_points(int x_start, int y_start, float x_pic_per_dot,
		      float y_pic_per_dot, int dis_l, int dis_t,
		      int zoom_labels_count, int *zoom_labels_row,
		      int *zoom_labels_col) {
  int i, hor_pos, ver_pos, label_hor_pos, lx, ty, rx, by, hor_shift, ver_shift;
  char string[255];
  char temp_string[50];
  hor_shift=(int)(x_pic_per_dot/3);
  ver_shift=(int)(y_pic_per_dot/3);
  for(i=0;i<zoom_labels_count;i++) {
    hor_pos=img_col_to_x(dis_l,x_start,x_pic_per_dot, zoom_labels_col[i]);
    ver_pos=img_row_to_y(dis_t,y_start,y_pic_per_dot, zoom_labels_row[i]);
    label_hor_pos=img_col_to_x(dis_l, x_start, x_pic_per_dot, 
			       zoom_labels_col[i]+1);
    lx=hor_pos-hor_shift;
    rx=hor_pos+hor_shift;
    ty=ver_pos-ver_shift;
    by=ver_pos+ver_shift;
    strcpy(string, num_string_int(zoom_labels_row[i], temp_string));
    strcat(string,",");
    strcat(string,num_string_int(zoom_labels_col[i], temp_string));
    gdImageFilledRectangle(image, lx, ty, rx, by, img_color[COL_DOTS]);
    gdImagePrintf(image, gdFontSmall, label_hor_pos + 1, ver_pos - 5, 
		  img_color[COL_DOTS], "(%s) ", string);
  }
}

void img_create_colors(float color_table[][3], int last_structure,
		       int nongray_for_3_flag, int lower_structure,
		       int jpg_mode) {
  int color[COLOR_TABLE_SIZE][3];
  int i,j;
  gdImagePtr bkim;
  /* convert colors to integers */
  if((nongray_for_3_flag)||(lower_structure<=-10))
    last_structure=last_structure+3;
  for(j=0;j<COLOR_TABLE_SIZE;j++) {
    for(i=0;i<=2;i++)
      color[j][i]=(int)(255.*color_table[j][i]);
  }
  bkim = gdImageCreate(2, 2);
    img_color[COL_BACKGROUND] = 
      gdImageColorAllocate(bkim, color[COL_BACKGROUND][0],
			   color[COL_BACKGROUND][1], color[COL_BACKGROUND][2]);
    img_color[COL_TEXT] = 
      gdImageColorAllocate(bkim, color[0][0], color[0][1], color[0][2]);
    gdImageLine(bkim, 0, 0, 1, 0, img_color[COL_TEXT]);
    gdImageLine(bkim, 0, 1, 1, 1, img_color[COL_TEXT]);
  /* end of brush definition */
    img_color[COL_BACKGROUND] = 
      gdImageColorAllocate(image, color[COL_BACKGROUND][0],
			   color[COL_BACKGROUND][1], color[COL_BACKGROUND][2]);
    /* first is background color */
    for(j=0;j<=last_structure;j++)
      img_color[j] = gdImageColorAllocate(image, color[j][0], color[j][1], 
					  color[j][2]);
    for(j=202;j<COLOR_TABLE_SIZE;j++)
      img_color[j] = gdImageColorAllocate(image, color[j][0], color[j][1], 
					  color[j][2]);
    /* Gray for grid lines */ 
    gdImageSetBrush(image, bkim); /* Sets the brush  for wide lines*/
    if(!jpg_mode)
      gdImageInterlace(image, 1); /* $$$$$ turn on interlace */
}
/* Make the image */

int open_imgdat(char *filename) {    /* open specified output file */
  if ((gifdatfp = fopen(filename, "w")) == NULL) {
    printf ("Error!\tCould not open output image file %s.", filename);
    return(1);
  }
  return (0);
}
void make_imgdat(char *filename,int resolution,int x_start,int y_start,
      float x_pic_per_dot,float y_pic_per_dot,int dis_l,int dis_t,
      int max_len) {
  char gifdat_name[120];
  strcpy(gifdat_name, filename);
  chop_suffix(gifdat_name, ".img");
  strcat(gifdat_name, ".gifdat");
  printf("The name for the .gifdat file is %s.\n", gifdat_name);
  if(open_imgdat(gifdat_name)) {
    printf("Error creating .gifdat file %s.",gifdat_name);
    return;
  }
  fprintf(gifdatfp,"%d %d %d %f %f %d %d %d \n", resolution, x_start, y_start,
	  x_pic_per_dot, y_pic_per_dot, max_len, dis_l,dis_t); 
  fclose(gifdatfp); 
}

void finish_making_img(int stripe_flag, char *img_filename, int resolution,
		       int zoom_labels_count, int *zoom_labels_row,
		       int *zoom_labels_col, int grid_flag, int last_structure,
		       struct helix *diag, float color_table[][3], 
		       int display_l, int display_r, int display_t,
		       int display_b, int lower_structure,
		       char *sequence_name, float dot_magnifier,
		       int *diag_start, int *diag_count, 
		       char file_data_ct[][255], int clear_flag, 
		       int gifdat_flag, int max_len, int forced_stripe_flag,
		       int nongray_for_3_flag,int *gray3_content,int png_mode,
		       int jpg_mode) {
  /* resolution is 72, 110, or 200 ? */
  int dis_l, dis_r, dis_t, dis_b;/* local versions of display_l,display_r */
  int hor_pos,ver_pos;
  int color,offset;
  int other_stop_point;
  int counter[MAXIMUM_CT_FILES+1];
  int stop_point;
  /* use 8.5 in x 11 in for image */
  int img_width; /* width of image in pixels */
  int img_height; /* height of image in pixels */
  int lx,ty,rx,by; /* top,bottom,right,left corn of rec */
  int i;
  int gray3[3];
  float x_pic_per_dot,y_pic_per_dot;
  int x_start,y_start;
  for(i=0; i<=last_structure; i++)
     counter[i]=0;
  dis_l = display_l;
  dis_r = display_r;
  dis_t = display_t;
  dis_b = display_b;
  img_width = (int) (8.5*resolution) - 1;
  img_height = 11*resolution - 1; 
  if(open_img(img_filename)) {
    if(png_mode)
      printf("Error creating png file.\n");
    else if(jpg_mode)
      printf("Error creating jpg file.\n"); 
    else
      printf("Error creating gif file.\n");
    return;
  }    /* hor by ver or x by y */
  image = gdImageCreate(img_width + 1, img_height + 1);
  /* declare image size in pixels */
  if(png_mode)
    printf("Creating png\n");
  else if(jpg_mode)
    printf("Creating jpg\n");
  else
    printf("Creating gif\n");
  printf("%d by %d \n", img_width+1, img_height+1);
  img_create_colors(color_table, last_structure, nongray_for_3_flag,
		    lower_structure, jpg_mode);
  /* declare gif colors */
  if(clear_flag==TRUE)
    gdImageColorTransparent(image, img_color[COL_BACKGROUND]);
  /* make box 7 x 7 in^2 */
  ty = (int)(1.57*resolution); /* come down 3 inches */
  by = (int)(img_height-2.37*resolution);
  lx = (int)(.47*resolution); /* left side of box */ 
  rx = img_width - (.95*resolution); /* 6.5 inches wide drawing region an
				      * extra .1 inch was used to make the
				      * line outside the region. Using a
				      * wide line too */ 
  display_header_message_img(image, img_width, 0 );
  gdImageLine(image, lx, ty, rx, ty, gdBrushed);
  gdImageLine(image, lx, by, rx, by, gdBrushed);
  gdImageLine(image, lx, ty, lx, by, gdBrushed);
  gdImageLine(image, rx, ty, rx, by, gdBrushed);
  gdImagePrintfCenter(image, gdFontLarge, img_width/2, (int) (resolution / 2),
		      img_color[COL_TEXT],
		      "Structure dot plot for %s", sequence_name);
  x_pic_per_dot = resolution*7.0/(dis_r-dis_l+1);
  /* width of base pair*/
  y_pic_per_dot = resolution*7.0/(dis_b-dis_t+1); /* height of base pair*/
  x_start = (int) (0.5*resolution + x_pic_per_dot/2 + 0.5);
  /* distance from left edge for first base pair*/
  y_start = (int) (1.59*resolution+y_pic_per_dot/2+.5);
  if(gifdat_flag == TRUE) {
    make_imgdat(img_filename, resolution, x_start, y_start, x_pic_per_dot,
		y_pic_per_dot, dis_l, dis_t, max_len);
  }
  /* distance from top for first base pair both above are in pixels
   * create the tic marks  
   */
  img_tics(x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
	   dis_r,dis_t,dis_b,lx,rx,ty,by,grid_flag);
  /* make the diagonal line */
  img_make_diag(dis_l,dis_r,dis_t,dis_b,x_start,
		y_start,x_pic_per_dot,y_pic_per_dot);
  /* create the points */
  if(forced_stripe_flag)
    stripe_flag=TRUE;
  else {
    if((nongray_for_3_flag)||(last_structure<4))
      stripe_flag=FALSE;
    if(stripe_flag==TRUE) {
      if( !( (x_pic_per_dot>(STRIPE_GIF*(last_structure-1)))&&
	     (y_pic_per_dot>(STRIPE_GIF*(last_structure-1))) ) )
	stripe_flag=FALSE;
    }
  }
  img_points(counter,stripe_flag,x_start,y_start,x_pic_per_dot,y_pic_per_dot,
	     dis_l,dis_r,dis_t,dis_b,last_structure,diag,
	     lower_structure,diag_start,diag_count,dot_magnifier,
	     nongray_for_3_flag,gray3_content,gray3);
  /* Display colors */
  hor_pos = resolution/2;
  ver_pos = img_height - 2.2*resolution;
  gdImageString(image, gdFontLarge, hor_pos, ver_pos,
		(unsigned char*) "Overlap", img_color[COL_TEXT]);
  gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + 18, 
		img_color[COL_COMPLETE_OVERLAP], "Full:   %4d", counter[0]);
  if(!stripe_flag) {
    if(last_structure>3) {
      gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + 36, 
		    img_color[COL_PART_OVERLAP], "Partial: %3d ", counter[1]);
       } 
   } 
   stop_point = last_structure/2 + 1;
   if(last_structure>16) {
     stop_point = 9;
     other_stop_point = 15; 
     hor_pos = resolution*6 - 10;
     gdImagePrintf(image, gdFontLarge, hor_pos, (int) 
		   (img_height - .47*resolution + 2), img_color[16], 
		   "%d others", last_structure - 15);
   } else {
     other_stop_point=last_structure;
   }
   hor_pos=resolution*3-50; 
   for(color=2,offset=0;color<=stop_point; color++,offset=offset+18)
     gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + offset, 
		   img_color[color], "%s %3d", file_data_ct[color-2],
		   counter[color]);
   hor_pos=resolution*6-30; 
   for(color=(stop_point+1),offset=0;color<=other_stop_point;
       color++,offset=offset+18)
     gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + offset, 
		   img_color[color], "%s %3d", file_data_ct[color-2],
		   counter[color]);
   if((nongray_for_3_flag)||(lower_structure<=-10)) {
     gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + 38, img_color[5], 
		   "1 & 2     %4d", gray3[0]);
     gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + 56, img_color[6],
		   "1 & 2     %4d", gray3[1]);
     gdImagePrintf(image, gdFontLarge, hor_pos, ver_pos + 74, img_color[7],
		   "2 & 3     %4d", gray3[2]);
   }
   /* Display labels */
   if(zoom_labels_count>0)
     img_label_points(x_start, y_start, x_pic_per_dot, y_pic_per_dot, dis_l,
		      dis_t, zoom_labels_count, zoom_labels_row, 
		      zoom_labels_col);
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
   if(png_mode)
     printf("The name for the png file is %s\n.", img_filename);
   else if(jpg_mode)
     printf("The name for the jpg file is %s.\n", img_filename);
   else
     printf("The name for the gif file is %s.\n", img_filename); 
}
