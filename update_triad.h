#ifndef include_update_triad
#define included_update_triad

#include "CartesianPatchGeometry.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchLevel.h"
#include "ibamr/IBRodForceSpec.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/app_namespaces.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
void update_triad(Pointer<PatchHierarchy<NDIM> > hierarchy,
                  LDataManager* lag_manager,
                  const double current_time,
                  const double dt,
                  std::vector<double>& params,
                  const int loop_num);
void read_from_input(Pointer<Database> input_db, std::vector<double>& params);

#endif
