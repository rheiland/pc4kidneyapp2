#include "./mechanics.h"


void epithelial_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
    static bool first_time = true;

    // these need to match that in custom.cpp
    static constexpr int NXG=875;  
    static constexpr int NYG=750;
    extern double pbm_grad_x[NYG][NXG];
    extern double pbm_grad_y[NYG][NXG];

    static double x_min = microenvironment.mesh.bounding_box[0]; 
	static double x_max = microenvironment.mesh.bounding_box[3]; 
	static double x_diff = microenvironment.mesh.bounding_box[3] - microenvironment.mesh.bounding_box[0]; 

	static double y_min = microenvironment.mesh.bounding_box[1]; 
	static double y_max = microenvironment.mesh.bounding_box[4]; 
	static double y_diff = microenvironment.mesh.bounding_box[4] - microenvironment.mesh.bounding_box[1]; 

    // std::cout << "------ epithelial_special_mechanics" << std::endl;

    double xpos = pCell->position[0];
    double ypos = pCell->position[1];

    double vmin = 1.e30;
    double vmax = -vmin;
    for (int idy=0; idy<NYG; idy++)
    for (int idx=0; idx<NXG; idx++)
    {
        if (pbm_grad_x[idy][idx] > vmax)
        {
            vmax = pbm_grad_x[idy][idx];
            // std::cout << "updating vmax="<<vmax<<" at idx,idy=" << idx <<"," <<idy<<std::endl;
        }
        if (pbm_grad_x[idy][idx] < vmin)
        {
            vmin = pbm_grad_x[idy][idx];
            // std::cout << "updating vmin="<<vmin<<" at idx,idy=" << idx <<"," <<idy<<std::endl;
        }
    }
    // std::cout << "------ epithelial_special_mechanics: dx min,max = " << vmin << ","<<vmax << std::endl;


    if (first_time) 
    {
        std::cout << "------ epithelial_special_mechanics: ID = " << pCell->ID << ": " << xpos << ", " << ypos << "; " << "x_diff = " << x_diff << std::endl;
    }

	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	// static int nRP = 0; // "rest_position"
	// std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value ; 

    int ix = (int)(((xpos - x_min)/x_diff) * NXG);
    // int iy = 75 - (int)(((ypos - y_min)/y_diff) * 75);
    int iy = (int)(((ypos - y_min)/y_diff) * NYG);

    if (first_time) 
    {
        std::cout << " ------ ID="<<pCell->ID <<": ix,iy= " << ix << "," << iy << "; dx,dy = " << pbm_grad_x[iy][ix] <<"," << pbm_grad_y[iy][ix] << std::endl;
        first_time = false;
    }

    // static double dscale = 0.1;
    // remember to negate the direction (gradient points toward higher values)
    // double vnorm = 0.1;
    // pCell->position[0] -= pbm_grad_x[iy][ix] * dscale;   
    // pCell->position[1] -= pbm_grad_y[iy][ix] * dscale;

    // adopted from update_motility_vector (in core/PhysiCell_cell.cpp)

    // choose a uniformly random unit vector 
    double temp_angle = 6.28318530717959 * UniformRandom();
    double temp_phi = 3.1415926535897932384626433832795 * UniformRandom();
    
    double sin_phi = sin(temp_phi);
    double cos_phi = cos(temp_phi);
    
    if( phenotype.motility.restrict_to_2D == true )
    { 
        sin_phi = 1.0; 
        cos_phi = 0.0;
    }
    
    std::vector<double> dvec; 
    dvec.resize(3,0.0); 
    
    // invert to point *to* PBM
    dvec[0] = -pbm_grad_x[iy][ix]; 
    dvec[1] = -pbm_grad_y[iy][ix]; 
    dvec[2] = 0.0; //  assuming 2D model
    
    // if the update_bias_vector function is set, use it  
    if( pCell->functions.update_migration_bias )
    {
        pCell->functions.update_migration_bias( pCell, phenotype, dt ); 
    }
    
    phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
    phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector 

    double one_minus_bias = 1.0 - phenotype.motility.migration_bias; 
		
    axpy( &(phenotype.motility.motility_vector), one_minus_bias, dvec ); // motility = (1-bias)*randvec + bias*bias_vector
		
    normalize( &(phenotype.motility.motility_vector) ); 
		
    phenotype.motility.motility_vector *= phenotype.motility.migration_speed;
    pCell->update_motility_vector( dt );
	return; 
}

void plasto_elastic_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	static int nRP = 0; // "rest_position"
	
	// displacement 
	std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value - pCell->position; 
	
	static int nEConst = pCell->custom_data.find_variable_index( "cell_elasticity" );
	static int nPConst = pCell->custom_data.find_variable_index( "cell_plasticity" );

	// first, update the agent's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[nEConst] , displacement );

	// now, plastic mechanical relaxation

	double plastic_temp_constant = -dt * pCell->custom_data[nPConst];
	axpy( &(pCell->custom_data.vector_variables[nRP].value) , plastic_temp_constant , displacement );
	
	return; 
}

void plasto_elastic_mechanics_v0( Cell* pCell, Phenotype& phenotype, double dt )
{
    std::cout << "------ plasto_elastic_mechanics\n";
	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	static int nRP = 0; // "rest_position"
	
	// displacement 
	std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value - pCell->position; 
	
	static int nEConst = pCell->custom_data.find_variable_index( "cell_elasticity" );
	static int nPConst = pCell->custom_data.find_variable_index( "cell_plasticity" );

	// first, update the agent's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[nEConst] , displacement );

	// now, plastic mechanical relaxation

	double plastic_temp_constant = -dt * pCell->custom_data[nPConst];
	axpy( &(pCell->custom_data.vector_variables[nRP].value) , plastic_temp_constant , displacement );
	
	return; 
}

// specialized potential function 
void add_heterotypic_potentials(Cell* my_cell , Cell* other_agent)
{
	static int nCellID = my_cell->custom_data.find_variable_index( "cell_ID" ); 
	
	// if( this->ID == other_agent->ID )
	if( my_cell == other_agent )
	{ return; }
	
	static int nRAOC = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cells" ); 
	static int nRAOCT = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cell_types" ); 
	
	double rel_heterotypic_adhesion = my_cell->custom_data[nRAOCT];  
	double rel_other_cells_adhesion = my_cell->custom_data[nRAOC];  
	
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		my_cell->displacement[i] = my_cell->position[i] - (*other_agent).position[i]; 
		distance += (my_cell->displacement[i]) * (my_cell->displacement[i]); 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = my_cell->phenotype.geometry.radius + (*other_agent).phenotype.geometry.radius; 
	
	double RN = my_cell->phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	else
	{
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		my_cell->state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( my_cell->phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion; 
	
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	double max_interactive_distance = my_cell->phenotype.mechanics.relative_maximum_adhesion_distance * my_cell->phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( my_cell->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		
		int my_id = (int) my_cell->custom_data[nCellID] ; 
		int other_id = (int) other_agent->custom_data[nCellID] ; 
	
		if( my_id != other_id )
		{ effective_adhesion *= rel_other_cells_adhesion; }
		
		if( my_cell->type != other_agent->type )
		{ effective_adhesion *= rel_heterotypic_adhesion; }
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	axpy( &(my_cell->velocity) , temp_r , my_cell->displacement ); 
	
	return;
}

void heterotypic_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// pCell->add_potentials(*neighbor);
		add_heterotypic_potentials( pCell, *neighbor ); 
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			// pCell->add_potentials(*neighbor);
			add_heterotypic_potentials( pCell, *neighbor ); 
		}
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}

void BM_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	return; 
}

