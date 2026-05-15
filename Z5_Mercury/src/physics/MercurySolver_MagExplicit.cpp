
#include "MercurySolver.h"

void MercurySolver::Build_E_explicit_edge_()
{
    AddIdealEdgeEMF_();

    AddPoleResistiveEdgeEMF_FromJcell_();

    // AddAmbipolarEdgeEMF_();
}