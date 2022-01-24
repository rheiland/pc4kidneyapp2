#ifndef __mechanics_h__
#define __mechanics_h__

#include "../core/PhysiCell.h" 

using namespace PhysiCell;
using namespace BioFVM; 

// viscoplastic model for BM elements -- they "anchor" to a preferred location 

void plasto_elastic_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); // done 

void add_heterotypic_potentials(Cell* my_cell , Cell* other_agent); // done 
void heterotypic_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt); // done 



void parietal_epithelial_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
// void mesangial_matrix_mechanics( Cell* pCell, Phenotype& phenotype, double dt );

void BM_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt );


#endif 