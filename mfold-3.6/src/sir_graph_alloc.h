void report_mem_error(char *string) {
	printf("Error!\tInsufficient memory to allocate %s.", string);
	try_exit(22);
}

void memory_allocator(int length) {
  int             size, size2;
  size = length + 4;
  size2 = length + 3;
  
  g_base_x = (float *) malloc((size) * sizeof(float));
  if (g_base_x == NULL)
    report_mem_error("g_base_x");
  g_base_y = (float *) malloc((size) * sizeof(float));
  if (g_base_y == NULL)
    report_mem_error("g_base_y");
  g_structure_angle = (float *) malloc((size) * sizeof(float));
  if (g_structure_angle == NULL)
    report_mem_error("g_structure_angle");
  g_structure_label_x = (float *) malloc((size) * sizeof(float));
  if (g_structure_label_x == NULL)
    report_mem_error("g_structure_label_x");
  g_structure_label_y = (float *) malloc((size) * sizeof(float));
  if (g_structure_label_y == NULL)
    report_mem_error("g_structure_label_y");
  
  g_structure_label_value = (int *) malloc((size) * sizeof(int));
  if (g_structure_label_value == NULL)
    report_mem_error("g_structure_label_value");
  g_history = (int *) malloc((size) * sizeof(int));
  if (g_history == NULL)
    report_mem_error("g_history");
  g_next = (int *) malloc((size) * sizeof(int)); /* 2003-06-21 Nick Markham */
  if (g_next == NULL)                            /* 2003-06-21 Nick Markham */
    report_mem_error("g_next");            /* 2003-06-21 Nick Markham */
  g_prev = (int *) malloc((size) * sizeof(int)); /* 2003-06-21 Nick Markham */
  if (g_prev == NULL)                            /* 2003-06-21 Nick Markham */
    report_mem_error("g_prev");            /* 2003-06-21 Nick Markham */
  
  g_undo1_base_x = (float *) malloc((size) * sizeof(float));
  if (g_undo1_base_x == NULL)
    printf("Insufficient memory to allocate g_undo1_base_x\n");
  g_undo1_base_y = (float *) malloc((size) * sizeof(float));
  if (g_undo1_base_y == NULL)
    report_mem_error("g_undo1_base_y");
  g_undo1_structure_angle = (float *) malloc((size) * sizeof(float));
  if (g_undo1_structure_angle == NULL)
    report_mem_error("g_undo1_structure_angle");
  g_undo2_base_x = (float *) malloc((size) * sizeof(float));
  if (g_undo2_base_x == NULL)
    printf("Insufficient memory to allocate g_undo2_base_x\n");
  g_undo2_base_y = (float *) malloc((size) * sizeof(float));
  if (g_undo2_base_y == NULL)
    report_mem_error("g_undo2_base_y");
  g_undo2_structure_angle = (float *) malloc((size) * sizeof(float));
  if (g_undo2_structure_angle == NULL)
    report_mem_error("g_undo2_structure_angle");
  
  g_previous_loop = (int *) malloc((size) * sizeof(int));

  if (g_previous_loop == NULL)
    report_mem_error("g_previous_loop");
  
  g_ss_code = (int *) malloc((size) * sizeof(int));
  if (g_ss_code == NULL)
    report_mem_error("g_ss_code");
  g_connected = (int *) malloc((size) * sizeof(int));
  if (g_connected == NULL)
    report_mem_error("g_connected");
  g_oldconnect = (int *) malloc((size) * sizeof(int));
  if (g_oldconnect == NULL)
    report_mem_error("g_oldconnect");
  g_bases = (char *) malloc((size) * sizeof(char));
  if (g_bases == NULL)
    report_mem_error("g_bases");
  g_base_list = (int *) malloc((size) * sizeof(int));
  if (g_base_list == NULL)
    report_mem_error("g_base_list");
  g_base_advance = (float *) malloc((size) * sizeof(float));
  if (g_base_advance == NULL)
    report_mem_error("g_base_advance");
  
  g_loop_center_x = (float *) malloc((size2) * sizeof(float));
  if (g_loop_center_x == NULL)
    report_mem_error("g_loop_center_x");
  
  g_loop_center_y = (float *) malloc((size2) * sizeof(float));
  if (g_loop_center_y == NULL)
    report_mem_error("g_loop_center_y");
  g_undo1_loop_center_x = (float *) malloc((size2) * sizeof(float));
  if (g_undo1_loop_center_x == NULL)
    report_mem_error("g_undo1_loop_center_x");
  g_undo1_loop_center_y = (float *) malloc((size2) * sizeof(float));
  if (g_undo1_loop_center_y == NULL)
    report_mem_error("g_undo1_loop_center_y");
  g_undo2_loop_center_x = (float *) malloc((size2) * sizeof(float));
  if (g_undo2_loop_center_x == NULL)
    report_mem_error("g_undo2_loop_center_x");
  g_undo2_loop_center_y = (float *) malloc((size2) * sizeof(float));
  if (g_undo2_loop_center_y == NULL)
    report_mem_error("g_undo2_loop_center_y");
  
  g_loop_radius = (float *) malloc((size2) * sizeof(float));
  if (g_loop_radius == NULL)
    report_mem_error("g_loop_radius");
  
  
  g_stack = (int *) malloc((size) * sizeof(int));
  if (g_stack == NULL)
    report_mem_error("g_stack");
  
  
  g_midpoint_x = (float *) malloc((size2) * sizeof(float));
  if (g_midpoint_x == NULL)
    report_mem_error("g_midpoint_x");
  g_midpoint_y = (float *) malloc((size2) * sizeof(float));
  if (g_midpoint_y == NULL)
    report_mem_error("g_midpoint_y");
  
  
  g_ann_to_color = (int *) malloc((size) * sizeof(int));
  if (g_ann_to_color == NULL)
    report_mem_error("g_ann_to_color");
}
