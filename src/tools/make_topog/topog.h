void create_rectangular_topog(int nx, int ny, double basin_depth, double *depth);
void create_bowl_topog(int nx, int ny, const double *x, const double *y, double bottom_depth,
		       double min_depth, double bowl_east, double bowl_south, double bowl_west,
		       double bowl_north, double *depth);
void create_gaussian_topog(int nx, int ny, const double *x, const double *y, double bottom_depth,
			   double min_depth, double gauss_amp, double gauss_scale, double slope_x,
			   double slope_y, double *depth);
void create_idealized_topog( int nx, int ny, const double *x, const double *y,
			     double bottom_depth, double min_depth, double *depth);
void create_realistic_topog(int nx_dst, int ny_dst, const double *x_dst, const double *y_dst, const char *vgrid_file, 
			    const char* topog_file, const char* topog_field, double scale_factor,
			    int tripolar_grid, int cyclic_x, int cyclic_y, 
			    int fill_first_row, int filter_topog, int num_filter_pass,
			    int smooth_topo_allow_deepening, int round_shallow, int fill_shallow,
			    int deepen_shallow, int full_cell, int flat_bottom, int adjust_topo,
			    int fill_isolated_cells, int dont_change_landmask, int kmt_min, double min_thickness,
			    int open_very_this_cell, double fraction_full_cell, double *depth, 
                            int *num_levels, domain2D domain, int debug, int great_circle_algorithm );

void create_box_channel_topog(int nx, int ny, double basin_depth,
			      double jwest_south, double jwest_north, double jeast_south,
			      double jeast_north, double *depth);
void create_dome_topog(int nx, int ny, const double *x, const double *y, double dome_slope,
		       double dome_bottom, double dome_embayment_west, double dome_embayment_east,
		       double dome_embayment_south, double dome_embayment_depth, double *depth);
