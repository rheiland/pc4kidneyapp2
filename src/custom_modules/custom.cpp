/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <cstdlib>

// double pbm_grad_x[75][88];
// double pbm_grad_y[75][88];

// static constexpr int NXG=88;
// static constexpr int NYG=75;
static constexpr int NXG=875;
static constexpr int NYG=750;
double pbm_grad_x[NYG][NXG];
double pbm_grad_y[NYG][NXG];


void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = basic_volume_model; // standard_volume_update_function;
//	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_velocity = heterotypic_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	// cell_defaults.functions.custom_cell_rule = plasto_elastic_mechanics; 
	// cell_defaults.functions.custom_cell_rule = epithelial_special_mechanics; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

	// let's add custom data not in the XML 
	std::vector<double> zero = {0,0,0};
	cell_defaults.custom_data.add_vector_variable( "rest_position" , "microns" , zero ); 
	cell_defaults.custom_data.add_vector_variable( "BM_attach_point" , "microns", zero ); 

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

    Cell_Definition* pCD; 
    pCD = find_cell_definition("parietal_epithelial"); 
	pCD->functions.custom_cell_rule = parietal_epithelial_mechanics; 
    // pCD = find_cell_definition("mesangial_matrix"); 
	// pCD->functions.custom_cell_rule = mesangial_matrix_mechanics; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void read_pbm_membrane_gradient_data( void )
{
    std::ifstream grad_file;
    std::string line;

    double dval;
    int n = 0;
    float vmin = 1.e6;
    float vmax = -1.e6;
    float v;

    int idy = 0;
    int idx = 0;

    // std::string fname1 = "./config/grad_x_pbm_875x750.dat";
    // std::string fname2 = "./config/grad_y_pbm_875x750.dat";

    // std::string fname1 = "../data/grad_x_pbm_875x750.dat";   // on nanoHUB, use ".."
    // std::string fname2 = "../data/grad_y_pbm_875x750.dat";

    const char* env_p = std::getenv("KIDNEY_DATA_PATH");
    std::cout << "       read_pbm_membrane_gradient_data(): KIDNEY_DATA_PATH: " << env_p << std::endl;
    // std::string full_filename = cwd_string + "/../" + filename;   // bloody nanoHUB
    // std::string full_filename = env_p + filename;
    std::string tmpname = "/grad_x_pbm_875x750.dat";   
    std::string fname1 = env_p + tmpname;

    tmpname = "/grad_y_pbm_875x750.dat";   
    std::string fname2 = env_p + tmpname;

    try {
        std::cout << " ---- Reading " << fname1 << std::endl;
        // grad_file.open("./config/grad_x_pbm.dat");
        grad_file.open(fname1);

        while (std::getline(grad_file,line))
        {
            // std::cout << "-- " << count1 << std::endl;
            // count1++;
            std::istringstream s(line);
            std::string field;
            // if (count1 > 5) break;
            idx = 0;
            while (std::getline(s,field,' ')) // beware, not comma separated!
            {
                // std::cout << count2 << ": " << field << std::endl;
                v = std::stof(field);
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
                pbm_grad_x[idy][idx] = v;
                n++;
                // count2++;
                // if (count2 > 5) break;
                idx += 1;
            }
            idy += 1;
        }
        std::cout << "-------- grad_x vmin,vmax= " << vmin << ", "  << vmax << std::endl;
    }
    catch (const std::ifstream::failure& e) {
      std::cout << ">>>>>>> Exception opening/reading " << fname1 << std::endl;
    }

    grad_file.close();

    try {
        std::cout << " ---- Reading " << fname2 << std::endl;
        // grad_file.open("./config/grad_x_pbm.dat");
        grad_file.open(fname2);
        idx = 0;
        idy = 0;
        while (std::getline(grad_file,line))
        {
            // std::cout << "-- " << count1 << std::endl;
            // count1++;
            std::istringstream s(line);
            std::string field;
            // if (count1 > 5) break;
            idx = 0;
            while (std::getline(s,field,' '))   // beware, not comma separated!
            {
                // std::cout << count2 << ": " << field << std::endl;
                v = std::stof(field);
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
                pbm_grad_y[idy][idx] = v;
                n++;
                // count2++;
                // if (count2 > 5) break;
                idx += 1;
            }
            idy += 1;
        }
        std::cout << "-------- grad_y vmin,vmax= " << vmin << ", "  << vmax << std::endl;
    }
    catch (const std::ifstream::failure& e) {
      std::cout << ">>>>>>> Exception opening/reading " << fname2 << std::endl;
    }

    // std::exit(1);

}

void read_membrane_distance_data( void )
{
    //---------------------------------
    // read glom BM distances (signed?)
	int gbm_index = microenvironment.find_density_index( "pbm_gbm_distance" ); 
    std::cout << "\n---------- gbm_index = " << gbm_index << std::endl;

    const char* env_p = std::getenv("KIDNEY_DATA_PATH");
    std::cout << "       read_membrane_distance_data(): KIDNEY_DATA_PATH: " << env_p << std::endl;
    std::string tmpname = "/pbm_gbm_dist.dat";   
    std::string fname1 = env_p + tmpname;

    tmpname = "/vessels4_dist.dat";   
    std::string fname2 = env_p + tmpname;


    // for( int n=0; n < microenvironment.number_of_voxels(); n++ )
    std::ifstream gbm_file;
    std::string line;
    try {
        // gbm_file.open("./data/pbm_gbm_dist.dat");
        gbm_file.open(fname1);

        double dval;
        int n = 0;
        float vmin = 1.e6;
        float vmax = -1.e6;
        float v;

        while (std::getline(gbm_file,line))
        {
            // std::cout << "-- " << count1 << std::endl;
            // count1++;
            std::istringstream s(line);
            std::string field;
            // if (count1 > 5) break;
            while (std::getline(s,field,','))
            {
                // std::cout << count2 << ": " << field << std::endl;
                v = std::stof(field);
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
                microenvironment(n)[gbm_index] = v;
                n++;
                // count2++;
                // if (count2 > 5) break;
            }
        }
        std::cout << "-------- vmin,vmax= " << vmin << ", "  << vmax << std::endl;
    }
    catch (const std::ifstream::failure& e) {
      std::cout << "Exception opening/reading file";
    }

    //---------------------------------
    // read glom capillary (4 blood vessel regions) distances (signed)
	int gcap_index = microenvironment.find_density_index( "blood_vessel_distance" ); 
    std::cout << "\n---------- gcap_index = " << gcap_index << std::endl;

    std::ifstream gcap_file;
    try {
        // gcap_file.open("./data/vessels4_dist.dat");
        gcap_file.open(fname2);

        double dval;
        int n = 0;
        float vmin = 1.e6;
        float vmax = -1.e6;
        float v;

        while (std::getline(gcap_file,line))
        {
            // std::cout << "-- " << count1 << std::endl;
            // count1++;
            std::istringstream s(line);
            std::string field;
            // if (count1 > 5) break;
            while (std::getline(s,field,','))
            {
                // std::cout << count2 << ": " << field << std::endl;
                v = std::stof(field);
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
                microenvironment(n)[gcap_index] = v;
                n++;
                // count2++;
                // if (count2 > 5) break;
            }
        }
        std::cout << "-------- vmin,vmax= " << vmin << ", "  << vmax << std::endl;
    }
    catch (const std::ifstream::failure& e) {
      std::cout << "Exception opening/reading file";
    }

}

void setup_microenvironment( void )
{
    // static std::string pbm_gbm_dist_file = parameters.strings( "pbm_gbm_distance_file" ); 
    // static std::string vessels_dist_file = parameters.strings( "vessels_distance_file" ); 

    // std::cout << "\n---------- setup_microenvironment(): dist1 file " << pbm_gbm_dist_file << std::endl;
    // std::cout << "\n---------- setup_microenvironment(): dist2 file" << vessels_dist_file << std::endl;

    // set domain parameters 

    // put any custom code to set non-homogeneous initial conditions or 
    // extra Dirichlet nodes here. 

    // initialize BioFVM 

    initialize_microenvironment();

    //---------------------------------
    // read glom BM distances (signed)
	// int pbm_gbm_index = microenvironment.find_density_index( "pbm_gbm_distance" ); 
    // std::cout << "\n---------- pbm_gbm_index = " << pbm_gbm_index << std::endl;

    // // for( int n=0; n < microenvironment.number_of_voxels(); n++ )
    // std::ifstream pbm_gbm_fp;
    // std::string line;
    // try {
    //     // gbm_file.open("../data/pbm_gbm_dist.dat");
    //     // gbm_file.open(parameters.strings("pbm_gbm_distance_file"));
    //     pbm_gbm_fp.open(pbm_gbm_dist_file);

    //     double dval;
    //     int n = 0;
    //     float vmin = 1.e6;
    //     float vmax = -1.e6;
    //     float v;

    //     while (std::getline(pbm_gbm_fp,line))
    //     {
    //         // std::cout << "-- " << count1 << std::endl;
    //         // count1++;
    //         std::istringstream s(line);
    //         std::string field;
    //         // if (count1 > 5) break;
    //         while (std::getline(s,field,','))
    //         {
    //             // std::cout << count2 << ": " << field << std::endl;
    //             v = std::stof(field);
    //             if (v < vmin) vmin = v;
    //             if (v > vmax) vmax = v;
    //             microenvironment(n)[pbm_gbm_index] = v;
    //             n++;
    //             // count2++;
    //             // if (count2 > 5) break;
    //         }
    //     }
    //     std::cout << "-------- vmin,vmax= " << vmin << ", "  << vmax << std::endl;
    // }
    // catch (const std::ifstream::failure& e) {
    //   std::cout << "Exception opening/reading file";
    // }

    // //---------------------------------
    // // read glom capillary (4 blood vessel regions) distances (signed)
	// int vessel_index = microenvironment.find_density_index( "blood_vessel_distance" ); 
    // std::cout << "\n---------- vessel_index = " << vessel_index << std::endl;

    // std::ifstream vessel_dist_fp;
    // try {
    //     // gcap_file.open("../data/vessels4_dist.dat");
    //     vessel_dist_fp.open(vessels_dist_file);

    //     double dval;
    //     int n = 0;
    //     float vmin = 1.e6;
    //     float vmax = -1.e6;
    //     float v;

    //     while (std::getline(vessel_dist_fp,line))
    //     {
    //         // std::cout << "-- " << count1 << std::endl;
    //         // count1++;
    //         std::istringstream s(line);
    //         std::string field;
    //         // if (count1 > 5) break;
    //         while (std::getline(s,field,','))
    //         {
    //             // std::cout << count2 << ": " << field << std::endl;
    //             v = std::stof(field);
    //             if (v < vmin) vmin = v;
    //             if (v > vmax) vmax = v;
    //             microenvironment(n)[vessel_index] = v;
    //             n++;
    //             // count2++;
    //             // if (count2 > 5) break;
    //         }
    //     }
    //     std::cout << "-------- vmin,vmax= " << vmin << ", "  << vmax << std::endl;
    // }
    // catch (const std::ifstream::failure& e) {
    //   std::cout << "Exception opening/reading file";
    // }

    return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	// for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	// {
	// 	Cell_Definition* pCD = cell_definitions_by_index[k]; 
	// 	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	// 	for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
	// 	{
	// 		std::vector<double> position = {0,0,0}; 
	// 		position[0] = Xmin + UniformRandom()*Xrange; 
	// 		position[1] = Ymin + UniformRandom()*Yrange; 
	// 		position[2] = Zmin + UniformRandom()*Zrange; 
			
	// 		pC = create_cell( *pCD ); 
	// 		pC->assign_position( position );
	// 	}
	// }
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
    if (load_subcells_from_pugixml() == false)
    {
        std::cout << "NOTE: ------- <cell_positions> in .xml are NOT enabled" << std::endl;
        return;
    }
	

	// more custom setup 
	// // for each cell, set its relaxed position to current position 
	int nRP = 0; // cell_defaults.custom_data.find_vector_variable( "rest_position");
	for( int n=0 ; n < (*all_cells).size(); n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		pC->custom_data.vector_variables[nRP].value = pC->position; 
	}
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring_wrapped(pCell); } 
// paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


bool load_subcells_from_pugixml( pugi::xml_node root )
{
	pugi::xml_node node = root.child( "initial_conditions" ); 
	if( !node )
	{ 
		std::cout << "Warning: XML-based cell positions has wrong formating. Ignoring!" << std::endl; 
		return false;
	}

	node = node.child( "cell_positions" ); 
	if( !node )
	{
		std::cout << "Warning: XML-based cell positions has wrong formating. Ignoring!" << std::endl; 
		 return false;
	}

	// enabled? 
	if( node.attribute("enabled").as_bool() == false )
	{ return false; }

	// get filename 

	std::string folder = xml_get_string_value( node, "folder" ); 
	std::string filename = xml_get_string_value( node, "filename" ); 
	std::string input_filename = folder + "/" + filename; 

	std::string filetype = node.attribute("type").value() ; 

	// what kind? 
	if( filetype == "csv" || filetype == "CSV" )
	{
		std::cout << "Loading cells from CSV file " << input_filename << " ... " << std::endl; 
		load_subcells_csv( input_filename );
		system("sleep 1");
		return true; 
	}
	if( filetype == "matlab" || filetype == "mat" || filetype == "MAT" )
	{
		std::cout << "Error: Load cell positions from matlab not yet supported. Try CSV." << std::endl; 
		exit(-1); 
		std::cout << "Loading cells from matlab file " << input_filename << " ... " << std::endl; 
		return false; 
	}
	if( filetype == "scene" )
	{
		std::cout << "Error: load cell positions from scene not yet supported. Try CSV." << std::endl; 
		exit(-1); 
		std::cout << "Loading cells from scene file " << input_filename << " ... " << std::endl; 
		return false; 
	}
	if( filetype == "physicell" || filetype == "PhysiCell" )
	{
		std::cout << "Error: load cell positions from PhysiCell snapshot not yet supported. Try CSV." << std::endl; 
		exit(-1); 
		std::cout << "Loading cells from PhysiCell file " << input_filename << " ... " << std::endl; 
		return false; 
	}

	return false; 
}

void load_subcells_csv( std::string filename )
{
    char tmp[256];
    getcwd(tmp, 256);
    std::cout << "custom.cpp: load_subcells_csv(): Current working directory: " << tmp << std::endl;
    std::cout << "custom.cpp: load_subcells_csv(): (param) filename= " << filename << std::endl;
    std::string cwd_string = tmp;   
    // if(const char* env_p = std::getenv("KIDNEY_DATA_PATH"))
    const char* env_p = std::getenv("KIDNEY_DATA_PATH");
    std::cout << "        KIDNEY_DATA_PATH: " << env_p << std::endl;
    // std::string full_filename = cwd_string + "/../" + filename;   // bloody nanoHUB
    // std::string full_filename = env_p + filename;
    // std::string tmpname = "/cells.csv";  //rwh: NO! 
    // std::string full_filename = env_p + tmpname;
    // std::string sep_str = "/";
    std::string full_filename = env_p + filename;
    std::cout << "custom.cpp: load_subcells_csv(): full_filename= " << full_filename << std::endl;


	// std::ifstream file( filename, std::ios::in );
	std::ifstream file( full_filename, std::ios::in );
	if( !file )
	{ 
		// std::cout << "Error: " << filename << " not found during cell loading. Quitting." << std::endl; 
		std::cout << "Error: " << full_filename << " not found during cell loading. Quitting." << std::endl; 
		exit(-1);
	}
	
	static int nCellID = cell_defaults.custom_data.find_variable_index( "cell_ID" ); 

	std::string line;
    int prev_type = -1;
	while (std::getline(file, line))
	{
		std::vector<double> data;
		csv_to_vector( line.c_str() , data ); 

		if( data.size() != 5 )
		{
			std::cout << "Error! Importing subcells from a CSV file expects each row to be x,y,z,typeID,cellID." << std::endl;
			exit(-1);
		}

		std::vector<double> position = { data[0] , data[1] , data[2] };

		int my_type = (int) data[3]; 
		int my_ID = (int) data[4]; 
		Cell_Definition* pCD = find_cell_definition( my_type );
		if( pCD != NULL )
		{
            if (pCD->type != prev_type)
            {
                std::cout << "Creating " << pCD->name << " (type=" << pCD->type << ") at " 
                    << position << std::endl; 
                prev_type = pCD->type;
            }
			Cell* pCell = create_cell( *pCD ); 
			pCell->assign_position( position ); 
			
			pCell->custom_data[nCellID] = my_ID; 
		}
		else
		{
			std::cout << "Warning! No cell definition found for index " << my_type << "!" << std::endl
			<< "\tIgnoring cell in " << full_filename << " at position " << position << std::endl; 
//			<< "\tIgnoring cell in " << filename << " at position " << position << std::endl; 
		}

	}

	file.close(); 	
}

bool load_subcells_from_pugixml( void )
{ return load_subcells_from_pugixml( physicell_config_root ); }

std::vector<std::string> paint_by_number_cell_coloring_wrapped( Cell* pCell )
{
	static std::vector< std::string > colors(0); 
	static bool setup_done = false; 
	if( setup_done == false )
	{
		colors.push_back( "grey" ); // default color will be grey 

		colors.push_back( "red" );
		colors.push_back( "yellow" ); 	
		colors.push_back( "green" ); 	
		colors.push_back( "blue" ); // 4
		
		colors.push_back( "magenta" ); 	
		colors.push_back( "orange" ); 	
		colors.push_back( "lime" ); 	
		colors.push_back( "cyan" ); // 8
		
		colors.push_back( "hotpink" ); 	
		colors.push_back( "peachpuff" ); 	
		colors.push_back( "darkseagreen" ); 	
		colors.push_back( "lightskyblue" ); // 12

		colors.push_back( "darkred" );
		colors.push_back( "goldenrod" ); 	
		colors.push_back( "darkgreen" ); 	
		colors.push_back( "darkblue" );  // 16
		
		setup_done = true; 
	}
	
	int number_of_colors = colors.size(); 
	
	// start all black 
	
	std::vector<std::string> output = { "black", "black", "black", "black" }; 
	
	// paint by number -- by cell type 
	
	// modular arithmetic  
	// int type = (int) ( ((unsigned int) pCell->type) % 13 ); 
	static int nCellID = cell_defaults.custom_data.find_variable_index( "cell_ID" ); 
	// int agent_type = pCell->type % number_of_colors; 
	int agent_type = (int) (pCell->custom_data[nCellID]) % number_of_colors; 
	
//	std::cout << type << std::endl; 
	
	std::string interior_color = "white"; 
	if( agent_type < number_of_colors )
	{ interior_color = colors[ agent_type ]; }
	
	output[0] = interior_color; // set cytoplasm color 
	
	if( pCell->phenotype.death.dead == false ) // if live, color nucleus same color 
	{
		output[2] = interior_color; 
		output[3] = interior_color; 
	}
	else
	{
		// apoptotic cells will retain a black nucleus 
		// if necrotic, color the nucleus brown 
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "rgb(139,69,19)";
			output[3] = "rgb(139,69,19)";
		}
	}
	
	return output; 
}
