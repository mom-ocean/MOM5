#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include "mpp.h"
#include "mpp_domain.h"
#include "mpp_io.h"
#include "tool_util.h"
#include "topog.h"

char *usage[] = {
  "",
  "                                                                                         ",
  "                                                                                         ",
  "                    Usage of make_topog                                                  ",
  "                                                                                         ",
  "   make_topog --mosaic mosaic_file [--topog_type topog_type]  [x_refine #] [y_refine #]  ",
  "              [--bottom_depth #] [--min_depth #]                                         ",
  "              [--topog_file topog_file]  [--topog_field topog_field]                     ",
  "              [--scale_factor #] [--num_filter_pass #] [--gauss_amp #]                   ",
  "              [--gauss_scale #] [--slope_x #] [--slope_y #] [--bowl_south #]             ",
  "              [--bowl_north #] [--bowl_west #] [--bowl_east #]                           ",
  "              [--fill_first_row] [--filter_topog] [--round_shallow]                      ",
  "              [--fill_shallow] [--deepen_shallow] [--smooth_topo_allow_deepening]        ",
  "              [--dome_slope #] [--dome_bottom #], [--dome_embayment_west #]              ",
  "              [--dome_embayment_east # [--dome_embayment_south #]                        ",
  "              [--dome_embayment_depth #] [--vgrid_file vgrid_file]                       ",
  "              [--flat_bottom] [--full_cell] [--fill_isolated_cells #]                    ",
  "              [--dont_change_landmask] [--kmt_min #] [--dont_adjust_topo ]               ",
  "              [--fraction_full_cell #] [--dont_open_very_this_cell ] [--min_thickness #] ",
  "              [--output output_file [--help] [--verbose]                                 ",
  "                                                                                         ",
  "   make_topog can generate topography for any Mosaic. The output file                    ",
  "   will contains the topography for each tile in the Mosaic. The field name in           ",
  "   the output topography file will be depth_tile# and it is positive down.               ",
  "   The topography data will be defined on model grid, the model grid size will be        ",
  "   supergrid grid size divided by refinement (x_refine, y_refine, default is 2).         ",
  "   --mosaic is a required option and all other options are optional, but                 ",
  "   some options are required depending on the choice of topog_type.                      ",
  "   Below specify the option (required or non-required) needed for every kind             ",
  "   of topog_type. when topog_type (--topog_type) is                                      ",
  "                                                                                         ",
  "     1. 'realistic':          Remap the topography onto the current grid from some       ",
  "                              source data file.                                          ",
  "                              --topog_file and --topog_field must be specified,          ",
  "                              --min_depth, --scale_factor,                               ",
  "                              --num_filter_pass, --flat_bottom                           ",
  "                              --fill_first_row  --filter_topog                           ",
  "                              --round_shallow] --fill_shallow --deepen_shallow]          ",
  "                              --smooth_topo_allow_deepening --vgrid_file                 ",
  "                              --full_cell --dont_fill_isolated_cells                     ",
  "                              --dont_change_landmask --kmt_min  --dont_adjust_topo       ",
  "                              --fraction_full_cell --dont_open_very_this_cell            ",
  "                              --min_thickness  are optional arguments.                   ",
  "                                                                                         ",
  "     2.  'rectangular_basin': Constructing a rectangular basin with a flat bottom.       ",
  "                              --bottome_depth are optional arguments. Set bottom_depth   ",
  "                              to 0 to get all land topography.                           ",
  "                                                                                         ",
  "     3.  'gaussian':          Construct gaussian bump on a sloping bottom.               ",
  "                              --bottom_depth, --min_depth --gauss_amp, --gauss_scale,    ",
  "                              --slope_x, --slope_y are optional arguments.               ",
  "                                                                                         ",  
  "     4.. 'bowl':              --bottom_depth, --min_depth, --bowl_south, --bowl_north,   ",
  "                              --bowl_west, --bowl_east are optional arguments.           ",
  "                                                                                         ",
  "     5. 'idealized':          Generates an 'idealized' not very realistic topography.    ",
  "                              --bottom_depth, --min_depth are optional arguments.        ",
  "                                                                                         ",
  "     6. 'box_channel':        Generate a box_channel topography. The interiol of the     ",
  "                              grid box is a flot bottom. The boundary of the grid        ",
  "                              box is land except points [jwest_south:jwest_north]        ",
  "                              and [ieast_south:ieast_north]. --jwest_south, jwest_north, ",
  "                              --ieast_south and ieast_north need to be specified.        ",
  "                              --bottom_depth are optional arguments.                     ",
  "                                                                                         ",
  "     7. 'dome':               similar (not identical) to DOME configuration of           ",
  "                              Legg etal Ocean Modelling (2005).  --dome_slope,           ",
  "                              --dome_bottom, --dome_embayment_west,                      ",
  "                              --dome_embayment_east, --dome_embayment_south and          ",
  "                              --dome_embayment_depth are optional arguments.             ",
  "                                                                                         ",
  "   generate_mosaic_topo take the following flags                                         ",
  "                                                                                         ",
  "   --mosaic mosaic_file          Specify the mosaic file where topography data located.  ",
  "                                                                                         ",
  "   --topog_type topog_type       Specify type of topography. Its value can be            ",
  "                                 'realistic', 'rectangular_basin', 'gaussian', 'bowl'    ",
  "                                 or 'idealized'. The default value is 'realistic'.       ",
  "                                                                                         ",
  "   --x_refine #                  specify the refinement ratio of model grid vs supergrid ",
  "                                 ins x-directin. default value 2.                        ",
  "                                                                                         ",
  "   --y_refine #                  specify the refinement ratio of model grid vs supergrid ",
  "                                 ins y-directin. default value 2.                        ",
  "                                                                                         ",
  "   --topog_file topog_file       Specify name of topograhy file (e.g. scripps,           ",
  "                                 navy_topo, ...)                                         ",
  "                                                                                         ",
  "   --topog_field topog_field     Specify name of topography field name in topog_file.    ",
  "                                                                                         ",
  "   --bottom_depth #              Specify maximum depth (or bottom depth) of ocean.       ",
  "                                 default value is 5000 meter.                            ",
  "                                                                                         ",
  "   --min_depth #                 Specify minimum depth of ocean.                         ",
  "                                 default value is 10 meter.                              ",
  "                                                                                         ",  
  "   --scale_factor #              Specify scaling factor for topography data (e.g. -1 to  ",
  "                                 flip sign or 0.01 to convert from centimeters).         ",
  "                                 default value is 1.                                     ",
  "                                                                                         ",
  "   --num_filter_pass #           Specify number of passes of spatial filter              ",
  "                                 default value is 1.                                     ",  
  "                                                                                         ",
  "   --gauss_amp #                 specify height of gaussian bump as percentage of ocean  ",
  "                                 depth. default value is 0.5.                            ",
  "                                                                                         ",
  "   --gauss_scale #               Specify width of gaussian bump as percentag e of        ",
  "                                 basin width. Default value is 0.25.                     ",
  "                                                                                         ",
  "   --slope_x #                   Specify rise of the ocean floor to the east for         ",
  "                                 the gaussian bump. Default value is 0.                  ",
  "                                                                                         ",
  "   --slope_y #                   Specify rise of the ocean floor to the north for        ",
  "                                 the gaussian bump. Default value is 0.                  ",
  "                                                                                         ",
  "   --bowl_south #                Specify southern boundary of Winton bowl.               ",
  "                                 Default value is 60.                                    ",
  "                                                                                         ",
  "   --bowl_north #                Specify northern boundary of Winton bowl.               ",
  "                                 Default value is 70.                                    ",  
  "                                                                                         ",
  "   --bowl_west #                 Specify western boundary of Winton bowl.                ",
  "                                 Default value is 0.                                     ",
  "                                                                                         ",  
  "   --bowl_east #                 Specify eastern boundary of Winton bowl.                ",
  "                                 Default value is 20.                                    ",
  "                                                                                         ",
  "   --fill_first_row              when specified, make first row of ocean model all       ",
  "                                 land points for ice model.                              ",
  "                                                                                         ",
  "   --filter_topog                When specified, apply filter to topography.             ",
  "                                                                                         ",
  "   --round_shallow               When specified, Make cells land if depth is less        ",
  "                                 than 1/2 mimumim depth, otherwise make ocean.           ",
  "                                                                                         ",
  "   --fill_shallow                When specified, Make cells less than minimum            ",
  "                                 depth land.                                             ",
  "                                                                                         ",
  "   --deepen_shallow              When specified, Make cells less than minimum            ",
  "                                 depth equal to minimum depth.                           ",
  "                                                                                         ",
  "   --smooth_topo_allow_deepening when specified, allow filter to deepen cells.           ",
  "                                                                                         ",
  "   --output output_file          The created netcdf file that contains mosaic            ",
  "                                 topography. Default value is 'topog.nc'                 ",
  "                                                                                         ",
  "   --jwest_south #               Specify the starting j-index on the west boundary of    ",
  "                                 the box channel.                                        ",
  "                                                                                         ",
  "   --jwest_north #               Specify the ending j-index on the west boundary of      ",
  "                                 the box channel.                                        ",  
  "                                                                                         ",
  "   --jeast_south #               Specify the starting j-index on the east boundary of    ",
  "                                 the box channel.                                        ",
  "                                                                                         ",
  "   --jeast_north #               Specify the ending j-index on the east boundary of      ",
  "                                 the box channel.                                        ",
  "                                                                                         ",
  "   --dome_slope #                Slope for the dome configuration.  Default = 0.01       ",
  "                                                                                         ",
  "   --dome_bottom #               Bottom of the dome configuration.  Default=3600.0       ",
  "                                                                                         ",
  "   --dome_embayment_west #       western edge of dome embayment. Default=19.0            ",
  "                                                                                         ",
  "   --dome_embayment_east #       eastern edge of dome embayment. Default=21.0            ",
  "                                                                                         ",
  "   --dome_embayment_south #      southern edge of dome embayment. Default=69.0           ",
  "                                                                                         ",
  "   --dome_embayment_depth #      Depth of the embayment. Default=600.0                   ",
  "                                                                                         ",
  "   --vgrid_file vgrid_file       Specify the vertical grid file. When vgrid_file is      ",
  "                                 specified, number of levels will be written out and     ",
  "                                 topography could be ajusted when set --dont_adjust_topo(",
  "                                 remove_isolated_cells, restrict_partial_cells and       ",
  "                                 enforce_min_depth ) is not set.                         ",
  "                                                                                         ",
  "   --flat_bottom                 generate flat bottom over ocean points when specified.  ",
  "                                                                                         ",
  "   --full_cell                   do not generate partial bottom cells when specified     ",
  "                                                                                         ",
  "   --dont_fill_isolated_cells    Do not allow non-advective tracer cells (strongly       ",
  "                                 recommended) when it is not set.                        ",
  "                                                                                         ",
  "   --dont_change_landmask        Do not change land/sea mask when filling isolated cells ",
  "                                 when specified                                          ",
  "                                                                                         ",
  "   --kmt_min #                   minimum number of vertical levels. default = 2          ",
  "                                                                                         ",
  "   --dont_adjust_topo            Topography will not be adjusted when it is set. When    ",
  "                                 it is not set, adjust topography                        ",
  "                                 (enforce_min_depth, remove_isolated_cells,              ",
  "                                 restrict_partial_cells). Strongly recommend not setting ",
  "                                 this option                                             ",
  "                                                                                         ",
  "   --fraction_full_cell #        Fraction of the associated full cell that a             ",
  "                                 corresponding partial cell thickness is no smaller than.",
  "                                 That is, we maintain partial_cell_min_dht(i,j,k) =      ",
  "                                 fraction_full_cell*full_cell_dzt(k). If                 ",
  "                                 fraction_full_cell=0.0, then partial_cell_min_dht       ",
  "                                 = min(zw(1), 50.0). default is 0.2.                     ",
  "                                                                                         ",
  "   --dont_open_very_this_cell    When it is not set, check which change is larger,       ",
  "                                 opening or closing the cell, and to do that with smaller",
  "                                 effect in depth_t.                                      ",
  "                                                                                         ",
  "   --min_thickness #             inimum vertical thickness allowed. with default value   ",
  "                                 0.1. Increase or decrease this number as needed.        ",
  "                                                                                         ",
  "   --help                        Print out this message and then exit.                   ",
  "                                                                                         ",
  "   --verbose                     Will print out running time message when this           ",
  "                                 option is set. Otherwise the run will be silent         ",
  "                                 when there is no error.                                 ",
  
  "   Example                                                                               ", 
  "                                                                                         ",
  "   1. Generate 'realistic' topography                                                    ",
  "   > make_topog --mosaic mosaic.nc --topog_type realistic                                ",
  "                           --topog_file /archive/fms/mom4/input_data/OCCAM_p5degree.nc   ",
  "                           --topog_field TOPO --scale_factor -1                          ",
  "                                                                                         ",
  "   2. Generate 'rectangular_basin' topography (with uniform topography 200 meters).      ",
  "   > make_topog --mosaic mosaic.nc --topog_type  rectangular_basin                       ",
  "                           --bottome_depth 200                                           ",
  "                                                                                         ",
  "   3. Generate 'box_channel' topography (with flat bottom topography 5000 meters).       ",
  "   > make_topog --mosaic mosaic.nc --topog_type  box_channel --jwest_south 20            ",
  "                --jwest_north 40 --jeast_south 30 --jeast_north 60 --bottom_depth 5000   ",
  "                                                                                         ",
  "   4. Generate 'bowl' topography                                                         ",
  "     make_topog --mosaic ocean_mosaic.nc --topog_type bowl --bottom_depth 5500           ",
  "                                                                                         ",
  "   5. Generate 'dome' topography                                                         ",
  "     make_topog --mosaic ocean_mosaic.nc --topog_type dome --dome_slope 0.01             ",
  "                --dome_embayment_south 69 --dome_embayment_west 19.25                    ",
  "                --dome_embayment_east 21.25 --dome_bottom 5500                           ",
  "",
  NULL };

const int REALISTIC         = 1;
const int RECTANGULAR_BASIN = 2;
const int GAUSSIAN          = 3;
const int BOWL              = 4;
const int IDEALIZED         = 5;
const int BOX_CHANNEL       = 6;
const int DOME              = 7;

char grid_version[] = "0.2";
char tagname[] = "$Name: tikal $";

int main(int argc, char* argv[])
{
  char   *mosaic_file = NULL, *topog_file = NULL, *topog_field = NULL;
  char   topog_type[32] = "realistic", output_file[32] = "topog.nc";
  int    num_filter_pass = 1;
  int    x_refine = 2, y_refine = 2;
  double bottom_depth = 5000, min_depth = 10, scale_factor = 1;
  double gauss_amp = 0.5, gauss_scale = 0.25, slope_x = 0, slope_y = 0;
  double bowl_south = 60, bowl_north = 70, bowl_west = 0, bowl_east = 20;
  int    flat_bottom = 0, fill_first_row = 0;
  int    filter_topog = 0, round_shallow = 0, fill_shallow = 0;
  int    deepen_shallow = 0, smooth_topo_allow_deepening = 0;
  int    jwest_south=0, jwest_north=0, jeast_south=0, jeast_north=0;
  double dome_slope=0.01;
  double dome_bottom=3600.0;
  double dome_embayment_west=19.0; 
  double dome_embayment_east=21.0;
  double dome_embayment_south=69.0;
  double dome_embayment_depth=600.0;
  char   *vgrid_file = NULL;
  int    full_cell = 0;
  int    fill_isolated_cells = 1;
  int    dont_change_landmask = 0;
  int    kmt_min = 2;
  int    adjust_topo = 1;
  double fraction_full_cell = 0.2;
  int    open_very_this_cell = 1;
  double min_thickness = 0.1;
  int    my_topog_type;
  int    use_great_circle_algorithm=0;
  int    cyclic_x, cyclic_y, tripolar_grid;
  int    errflg = (argc == 1);
  int    option_index, i, c;
  unsigned int verbose = 0;
  
  /*
   * process command line
   */

  static struct option long_options[] = {
    {"mosaic",                      required_argument, NULL, 'a'},
    {"topog_type",                  required_argument, NULL, 'b'},
    {"x_refine",                    required_argument, NULL, 'X'},
    {"y_refine",                    required_argument, NULL, 'Y'},
    {"topog_file",                  required_argument, NULL, 'd'},
    {"topog_field",                 required_argument, NULL, 'e'},
    {"bottom_depth",                required_argument, NULL, 'f'},
    {"min_depth",                   required_argument, NULL, 'g'},
    {"scale_factor",                required_argument, NULL, 'i'},
    {"num_filter_pass",             required_argument, NULL, 'j'},
    {"gauss_amp",                   required_argument, NULL, 'k'},
    {"gauss_scale",                 required_argument, NULL, 'l'},
    {"slope_x",                     required_argument, NULL, 'm'},    
    {"slope_y",                     required_argument, NULL, 'n'},
    {"bowl_south",                  required_argument, NULL, 'p'},
    {"bowl_north",                  required_argument, NULL, 'q'},
    {"bowl_west",                   required_argument, NULL, 'r'},
    {"bowl_east",                   required_argument, NULL, 's'},
    {"fill_first_row",              no_argument,       NULL, 't'},
    {"filter_topog",                no_argument,       NULL, 'u'},    
    {"round_shallow",               no_argument,       NULL, 'v'},
    {"fill_shallow",                no_argument,       NULL, 'w'},
    {"deepen_shallow",              no_argument,       NULL, 'x'},
    {"smooth_topo_allow_deepening", no_argument,       NULL, 'y'},
    {"output",                      required_argument, NULL, 'o'},
    {"jwest_south",                 required_argument, NULL, 'A'},
    {"jwest_north",                 required_argument, NULL, 'B'},
    {"jeast_south",                 required_argument, NULL, 'C'},
    {"jeast_north",                 required_argument, NULL, 'D'},
    {"dome_slope",                  required_argument, NULL, 'E'},
    {"dome_bottom",                 required_argument, NULL, 'F'},
    {"dome_embayment_west",         required_argument, NULL, 'G'},
    {"dome_embayment_east",         required_argument, NULL, 'H'},
    {"dome_embayment_south",        required_argument, NULL, 'I'},
    {"dome_embayment_depth",        required_argument, NULL, 'J'},
    {"vgrid_file",                  required_argument, NULL, 'K'},
    {"flat_bottom",                 no_argument,       NULL, 'L'},
    {"full_cell",                   no_argument,       NULL, 'M'},
    {"dont_fill_isolated_cells",    no_argument,       NULL, 'N'},
    {"dont_change_landmask",        no_argument,       NULL, 'O'},
    {"kmt_min",                     required_argument, NULL, 'P'},
    {"dont_adjust_topo",            no_argument,       NULL, 'Q'},
    {"fraction_full_cell",          required_argument, NULL, 'R'},
    {"dont_open_very_this_cell",    no_argument,       NULL, 'S'},
    {"min_thickness",               required_argument, NULL, 'T'},
    {"help",                        no_argument,       NULL, 'h'},
    {"verbose",                     no_argument,       NULL, 'V'},
    {0, 0, 0, 0},
  };

  /* start parallel */

  mpp_init(&argc, &argv);

  mpp_domain_init();  
   
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
    switch (c) {
    case 'a':
      mosaic_file = optarg;
      break;      
    case 'b':
      strcpy(topog_type, optarg);
      break;
    case 'X':
      x_refine = atoi(optarg);
      break;
    case 'Y':
      y_refine = atoi(optarg);
      break;      
    case 'd':
      topog_file = optarg;
      break;
    case 'e':
      topog_field = optarg;
      break;
    case 'f':
      bottom_depth = atof(optarg);
      break; 
    case 'g':
      min_depth = atof(optarg);
      break;
    case 'i':
      scale_factor = atof(optarg);
      break;
    case 'j':
      num_filter_pass = atoi(optarg);
      break; 
    case 'k':
      gauss_amp = atof(optarg);
      break;
    case 'l':
      gauss_scale = atof(optarg);
      break;
    case 'm':
      slope_x = atof(optarg);
      break;
    case 'n':
      slope_y = atof(optarg);
      break;
    case 'p':
      bowl_south = atof(optarg);
      break;
    case 'q':
      bowl_north = atof(optarg);
      break;
    case 'r':
      bowl_west  = atof(optarg);
      break;
    case 's':
      bowl_east  = atof(optarg);
      break;
    case 't':
      fill_first_row = 1;
      break;
    case 'u':
      filter_topog = 1;
      break;
    case 'v':
      round_shallow = 1;
      break;
    case 'w':
      fill_shallow = 1;
      break;
    case 'x':
      deepen_shallow = 1;
      break;
    case 'y':
      smooth_topo_allow_deepening = 1;
      break;
    case 'o':
      strcpy(output_file,optarg);
      break;
    case 'A':
      jwest_south = atoi(optarg);
      break;
    case 'B':
      jwest_north = atoi(optarg);
      break;
    case 'C':
      jeast_south = atoi(optarg);
      break;
    case 'D':
      jeast_north = atoi(optarg);
      break;
    case 'E':
      dome_slope = atof(optarg);
      break;
    case 'F':
      dome_bottom = atof(optarg);
      break;
    case 'G':
      dome_embayment_west = atof(optarg);
      break;
    case 'H':
      dome_embayment_east = atof(optarg);
      break;
    case 'I':
      dome_embayment_south = atof(optarg);
      break;
    case 'J':
      dome_embayment_depth = atof(optarg);
      break;
    case 'K':
      vgrid_file = optarg;
      break;
    case 'L':
      flat_bottom = 1;
      break;
    case 'M':
      full_cell = 1;
      break;
    case 'N':
      fill_isolated_cells = 0;
      break;
    case 'O':
      dont_change_landmask = 1;
      break;
    case 'P':
      kmt_min = atoi(optarg);
      break;
    case 'Q':
      adjust_topo = 0;
      break;
    case 'R':
      fraction_full_cell = atof(optarg);
      break;
    case 'S':
      open_very_this_cell = 0;
      break;
    case 'T':
      min_thickness = atof(optarg);
      break;
    case 'V':
      verbose = 1;
      break;      
    case '?':
      errflg++;
      break;
    }

  if (errflg || !mosaic_file ) {
    if( mpp_pe() == mpp_root_pe() ) {
      char **u = usage;
      while (*u) { fprintf(stderr, "%s\n", *u); u++; }
      mpp_error("make_topog: Check your arguments");
    }
  }


  /* Write out arguments value  */
  if(mpp_pe() == mpp_root_pe() && verbose) printf("NOTE from make_topog ==> the topog_type is: %s\n",topog_type);
  if(x_refine != 2 || y_refine != 2 ) mpp_error("Error from make_topog: x_refine and y_refine should be 2, contact developer");
  if(mpp_pe() == mpp_root_pe() && verbose) printf("NOTE from make_topog ==> x_refine is %d, y_refine is %d\n",
				       x_refine, y_refine);

  /* vgrid_file can only be passed in when topog_type is realistic */ 
  if(vgrid_file && strcmp(topog_type,"realistic"))
      mpp_error("make_topog: --vgrid_file should not be specified when topog_type = realistic");

  if (strcmp(topog_type,"rectangular_basin") == 0) {
    my_topog_type = RECTANGULAR_BASIN;
    if(mpp_pe() == mpp_root_pe() && verbose) printf("NOTE from make_topog ==> the basin depth is %f\n",bottom_depth);
  }
  else if (strcmp(topog_type,"gaussian") == 0) {
    my_topog_type = GAUSSIAN;
    if(mpp_pe() == mpp_root_pe() && verbose){
      printf("NOTE from make_topog ==> bottom_depth is: %f\n", bottom_depth );
      printf("NOTE from make_topog ==> min_depth is: %f\n", min_depth );
      printf("NOTE from make_topog ==> gauss_amp is: %f\n", gauss_amp );
      printf("NOTE from make_topog ==> gauss_scale is: %f\n", gauss_scale );
      printf("NOTE from make_topog ==> slope_x is: %f\n", slope_x );
      printf("NOTE from make_topog ==> slope_y is: %f\n", slope_y );      
    }
  }
  else if(strcmp(topog_type,"bowl") == 0) {
    my_topog_type = BOWL;
    if(mpp_pe() == mpp_root_pe() && verbose){
      printf("NOTE from make_topog ==> bottom_depth is: %f\n",bottom_depth);
      printf("NOTE from make_topog ==> min_depth is: %f\n",min_depth);
      printf("NOTE from make_topog ==> bowl_south is: %f\n",bowl_south);
      printf("NOTE from make_topog ==> bowl_north is: %f\n",bowl_north);
      printf("NOTE from make_topog ==> bowl_west is: %f\n",bowl_west);
      printf("NOTE from make_topog ==> bowl_east is: %f\n",bowl_east);
    }
  }
  else if(strcmp(topog_type,"idealized") == 0) {
    my_topog_type = IDEALIZED;
    if(mpp_pe() == mpp_root_pe() && verbose){
      printf("NOTE from make_topog ==> bottom_depth is: %f\n",bottom_depth);
      printf("NOTE from make_topog ==> min_depth is: %f\n",min_depth);
    }
  }
  else if(strcmp(topog_type,"realistic") == 0) {
    my_topog_type = REALISTIC;
    if(!topog_file || !topog_field)
      mpp_error("Error from make_topog: when topog_type is realistic, topog_file and topog_field must be specified.");
    if(mpp_pe() == mpp_root_pe() && verbose){
      printf("\n\n ************************************************************\n\n");
      printf("NOTE from make_topog ==> input arguments\n\n");
      printf("NOTE from make_topog ==> min_depth is: %f\n",min_depth);
      printf("NOTE from make_topog ==> topog_file is: %s\n", topog_file);
      printf("NOTE from make_topog ==> topog_field is: %s\n", topog_field);
      printf("NOTE from make_topog ==> scale_factor is: %f\n", scale_factor);
      printf("NOTE from make_topog ==> num_filter_pass is: %d\n", num_filter_pass);
      printf("NOTE from make_topog ==> kmt_min is %d\n", kmt_min);
      printf("NOTE from make_topog ==> fraction_full_cell is %f\n", fraction_full_cell);
      printf("NOTE from make_topog ==> min_thickness is %f\n", min_thickness);
      if(vgrid_file)
	printf("NOTE from make_topog ==> vgrid_file is %s\n", vgrid_file);
      else
	printf("NOTE from make_topog ==> no vgrid_file is specified\n");
      
      if(fill_first_row) printf("NOTE from make_topog ==>make first row of ocean model all land points.\n");
      if(filter_topog) printf("NOTE from make_topog ==>will apply filter to topography.\n");
      if(round_shallow) printf("NOTE from make_topog ==>Make cells land if depth is less than 1/2 "
			       "mimumim depth, otherwise make ocean.\n");
      if(fill_shallow) printf("NOTE from make_topog ==>Make cells less than minimum depth land.\n");
      if(deepen_shallow) printf("NOTE from make_topog ==>Make cells less than minimum depth equal to minimum depth.\n");
      if(smooth_topo_allow_deepening) printf("NOTE from make_topog ==>allow filter to deepen cells.\n");
      if(flat_bottom) printf("NOTE from make_topog ==> generate flat bottom over ocean points.\n");
      if(full_cell) printf("NOTE from make_topog ==> not generate partial bottom cells.\n");
      if(fill_isolated_cells) printf("NOTE from make_topog ==> not allow non-advective tracer cells\n");
      if(dont_change_landmask) printf("NOTE from make_topog ==> not change land/sea mask when filling isolated cells\n");
      if(open_very_this_cell) printf("NOTE from make_topog ==> open this cell\n");
      if(adjust_topo) printf("NOTE from make_topog ==> adjust topography\n");
      printf("\n\n ************************************************************\n\n");
    }
  }
  else if(strcmp(topog_type,"box_channel") == 0) {
    my_topog_type = BOX_CHANNEL;
    if( jwest_south <= 0) mpp_error("make_topog: jwest_south must a positive integer when topog_type = box_channel");
    if( jwest_north <= 0) mpp_error("make_topog: jwest_north must a positive integer when topog_type = box_channel");
    if( jeast_south <= 0) mpp_error("make_topog: jeast_south must a positive integer when topog_type = box_channel");
    if( jeast_north <= 0) mpp_error("make_topog: jeast_north must a positive integer when topog_type = box_channel");
    if( jwest_south > jwest_north ) mpp_error("make_topog: jwest_south > jwest_north when topog_type = box_channel");
    if( jeast_south > jeast_north ) mpp_error("make_topog: jeast_south > jeast_north when topog_type = box_channel");
  }
  else if(strcmp(topog_type,"dome") == 0) {
    my_topog_type = DOME;
  }
  else {
    mpp_error("make_topog: topog_type should be rectangular_basin, gaussian, bowl, idealized, realistic, box_channel or dome");
  }
  
  if(mpp_pe() == mpp_root_pe() && verbose) {
    printf("**************************************************\n");
    printf("Begin to generate topography \n");
  }

  {
    const int STRING = 255;
    int m_fid, g_fid, vid;
    int ntiles, fid, dim_ntiles, n, dims[2];
    size_t start[4], nread[4], nwrite[4];
    int *nx=NULL, *ny=NULL, *nxp=NULL, *nyp=NULL;
    int *id_depth=NULL, *id_level=NULL;
    double *depth=NULL, *x=NULL, *y=NULL;
    int *num_levels=NULL;
    char **tile_files=NULL;
    char history[512], dimx_name[128], dimy_name[128], depth_name[128], level_name[128];
    char gridfile[256], griddir[256];
    
    /* history will be write out as global attribute
       in output file to specify the command line arguments
    */

    strcpy(history,argv[0]);

    for(i=1;i<argc;i++) {
      strcat(history, " ");
      strcat(history, argv[i]);
    }

    /* grid should be located in the same directory of mosaic file */
    get_file_path(mosaic_file, griddir);
    
    /* get mosaic dimension */
    m_fid = mpp_open(mosaic_file, MPP_READ);
    ntiles = mpp_get_dimlen( m_fid, "ntiles");
    
    tile_files = (char **)malloc(ntiles*sizeof(double *));
    id_depth = (int *)malloc(ntiles*sizeof(int));
    id_level = (int *)malloc(ntiles*sizeof(int));
    /* loop through each tile to get tile information and set up meta data for output file */
    fid = mpp_open(output_file, MPP_WRITE);

    dim_ntiles = mpp_def_dim(fid, "ntiles", ntiles);
    nx = (int *)malloc(ntiles*sizeof(int));
    ny = (int *)malloc(ntiles*sizeof(int));
    nxp = (int *)malloc(ntiles*sizeof(int));
    nyp = (int *)malloc(ntiles*sizeof(int));   
    for( n = 0; n < ntiles; n++ ) {
      int use_great_circle_algorithm_prev=0;
      
      tile_files[n] = (char *)malloc(STRING*sizeof(double));
      start[0] = n;
      start[1] = 0;
      nread[0] = 1;
      nread[1] = STRING;
      vid = mpp_get_varid(m_fid, "gridfiles");
      mpp_get_var_value_block(m_fid, vid, start, nread, gridfile);
      sprintf(tile_files[n], "%s/%s", griddir, gridfile);
      g_fid = mpp_open(tile_files[n], MPP_READ);
      nx[n] = mpp_get_dimlen(g_fid, "nx");
      ny[n] = mpp_get_dimlen(g_fid, "ny");
      if( nx[n]%x_refine != 0 ) mpp_error("make_topog: supergrid x-size can not be divided by x_refine");
      if( ny[n]%y_refine != 0 ) mpp_error("make_topog: supergrid y-size can not be divided by y_refine");
      nx[n] /= x_refine;
      ny[n] /= y_refine;
      nxp[n] = nx[n] + 1;
      nyp[n] = ny[n] + 1;
      if(ntiles == 1) {
	strcpy(dimx_name, "nx");
	strcpy(dimy_name, "ny");
	strcpy(depth_name, "depth");
	if(vgrid_file)strcpy(level_name, "num_levels");
      }
      else {
	sprintf(dimx_name, "nx_tile%d", n+1);
	sprintf(dimy_name, "ny_tile%d", n+1);
	sprintf(depth_name, "depth_tile%d", n+1);
	if(vgrid_file)sprintf(level_name, "num_levels_tile%d", n+1);
      }

      dims[1] = mpp_def_dim(fid, dimx_name, nx[n]); 
      dims[0] = mpp_def_dim(fid, dimy_name, ny[n]);
      id_depth[n] = mpp_def_var(fid, depth_name, NC_DOUBLE, 2, dims,  2, "standard_name",
				"topographic depth at T-cell centers", "units", "meters");
      if(vgrid_file) id_level[n] = mpp_def_var(fid, level_name, NC_INT, 2, dims, 2, "standard_name",
					       "number of vertical T-cells", "units", "none");
      /* when topog_type is realistics, check if use great_circle_algorithm */
      use_great_circle_algorithm = get_great_circle_algorithm(g_fid);
      if(n>0) {
	if( use_great_circle_algorithm != use_great_circle_algorithm_prev)
	  mpp_error("make_topog: atribute 'great_circle_algorithm' of field 'tile' have different value for different tile");
      }
      use_great_circle_algorithm_prev = use_great_circle_algorithm;

      mpp_close(g_fid);
    }
    mpp_close(m_fid);
    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    if(use_great_circle_algorithm) mpp_def_global_att(fid, "great_circle_algorithm", "TRUE");
    mpp_def_global_att(fid, "history", history);
    
    mpp_end_def(fid);

    if(mpp_pe()==mpp_root_pe() && use_great_circle_algorithm)
      printf("\n NOTE from make_topog: use great circle algorithm\n");
    
    /* get the boundary condition for realistics topogrpahy, currently only support tripolar_grid,
     cyclic_x and cyclic_y*/
    if(my_topog_type == REALISTIC) get_boundary_type(mosaic_file, VERSION_2, &cyclic_x, &cyclic_y, &tripolar_grid);
    
    /* Generate topography and write out to the output_file */
    for(n=0; n<ntiles; n++) {
      int layout[2], isc, iec, jsc, jec, nxc, nyc, ni, i, j;
      double *gdata=NULL, *tmp=NULL;
      domain2D domain;
      
      /* define the domain, each tile will be run on all the processors. */
      mpp_define_layout( nx[n], ny[n], mpp_npes(), layout);
      mpp_define_domain2d( nx[n], ny[n], layout, 0, 0, &domain);
      mpp_get_compute_domain2d( domain, &isc, &iec, &jsc, &jec);
      nxc = iec - isc + 1;
      nyc = jec - jsc + 1;

      if(my_topog_type == DOME ) {
	x   = (double *)malloc(nxc*nyc*sizeof(double));
	y   = (double *)malloc(nxc*nyc*sizeof(double));
      }
      else {
	x   = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));
	y   = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));
      }
      depth = (double *)malloc(nxc*nyc*sizeof(double));
      if(vgrid_file) num_levels = (int* )malloc(nxc*nyc*sizeof(int)); 
      tmp   = (double *)malloc((nxc*x_refine+1)*(nyc*y_refine+1)*sizeof(double));
      start[0] = jsc*y_refine; start[1] = isc*x_refine;
      nread[0] = nyc*y_refine+1; nread[1] = nxc*x_refine+1;
      ni       = nxc*x_refine+1;
      g_fid = mpp_open(tile_files[n], MPP_READ);
      vid = mpp_get_varid(g_fid, "x");
      mpp_get_var_value_block(g_fid, vid, start, nread, tmp);
      if(my_topog_type == DOME ) {
	for(j = 0; j < nyc; j++) for(i = 0; i < nxc; i++)
	  x[j*nxc+i] = tmp[(j*y_refine+1)*ni+i*x_refine+1];
      }
      else {
	for(j = 0; j < nyc+1; j++) for(i = 0; i < nxc+1; i++)
	  x[j*(nxc+1)+i] = tmp[(j*y_refine)*ni+i*x_refine];
      }
      vid = mpp_get_varid(g_fid, "y");
      mpp_get_var_value_block( g_fid, vid, start, nread, tmp);
      mpp_close(g_fid);
      if(my_topog_type == DOME ) {
	for(j = 0; j < nyc; j++) for(i = 0; i < nxc; i++)
	  y[j*nxc+i] = tmp[(j*y_refine+1)*ni+i*x_refine+1];
      }
      else {
	for(j = 0; j < nyc+1; j++) for(i = 0; i < nxc+1; i++)
	  y[j*(nxc+1)+i] = tmp[(j*y_refine)*ni+i*x_refine];
      }
      switch (my_topog_type) {
      case RECTANGULAR_BASIN:
	create_rectangular_topog(nx[n], ny[n], bottom_depth, depth);
	break;
      case GAUSSIAN:
	create_gaussian_topog(nx[n], ny[n], x, y, bottom_depth, min_depth,
			      gauss_amp, gauss_scale, slope_x, slope_y, depth);
	break;
      case BOWL:
	create_bowl_topog(nx[n], ny[n], x, y, bottom_depth, min_depth, bowl_east,
			  bowl_south, bowl_west, bowl_north, depth);
	break;
      case IDEALIZED:
	create_idealized_topog( nx[n], ny[n], x, y, bottom_depth, min_depth, depth);
	break;
      case REALISTIC:	
	create_realistic_topog(nxc, nyc, x, y, vgrid_file, topog_file, topog_field, scale_factor,
			       tripolar_grid, cyclic_x, cyclic_y, fill_first_row, filter_topog, num_filter_pass,
			       smooth_topo_allow_deepening, round_shallow, fill_shallow,
			       deepen_shallow, full_cell, flat_bottom, adjust_topo,
			       fill_isolated_cells, dont_change_landmask, kmt_min, min_thickness, open_very_this_cell,
			       fraction_full_cell, depth, num_levels, domain, verbose, use_great_circle_algorithm );
	break;
      case BOX_CHANNEL:
	create_box_channel_topog(nx[n], ny[n], bottom_depth,
				 jwest_south, jwest_north, jeast_south, jeast_north, depth);
	break;
      case DOME:
	create_dome_topog(nx[n], ny[n], x, y, dome_slope, dome_bottom, dome_embayment_west,
			  dome_embayment_east, dome_embayment_south, dome_embayment_depth, depth);
	break;
      }
      gdata = (double *)malloc(nx[n]*ny[n]*sizeof(double));
      mpp_global_field_double(domain, nxc, nyc, depth, gdata);
      mpp_put_var_value(fid, id_depth[n], gdata);
      if(vgrid_file) {
	double *tmp_double=NULL;
	int    *gdata_int=NULL;
	tmp_double = (double *)malloc(nxc*nyc*sizeof(double));
	gdata_int = (int *)malloc(nx[n]*ny[n]*sizeof(int));
	for(i=0; i<nxc*nyc; i++) tmp_double[i] = num_levels[i];
	mpp_global_field_double(domain, nxc, nyc, tmp_double, gdata);
	for(i=0; i<nx[n]*ny[n]; i++) gdata_int[i] = gdata[i];
	mpp_put_var_value(fid, id_level[n], gdata_int);
	free(gdata_int);
	free(tmp_double);
	free(num_levels);
      }
      free(x);
      free(y);
      free(tmp);
      free(depth);
      free(gdata);
      mpp_delete_domain2d(&domain);
    }
    mpp_close(fid);
  
    /*release memory */
    free(id_depth);
    free(id_level);
    for(n=0; n<ntiles; n++) free(tile_files[n]);
    free(tile_files);
    free(nx);
    free(ny);
    free(nxp);
    free(nyp);
  }
    
  if(mpp_pe() == mpp_root_pe() && verbose ) printf("Successfully generate %s\n",output_file);

  mpp_end();

  return 0;
}; //main
