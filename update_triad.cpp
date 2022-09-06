#include "ibamr/IBAnchorPointSpec.h"
#include "ibamr/IBSpringForceSpec.h"

#include "update_triad.h"

void
update_triad(Pointer<PatchHierarchy<NDIM> > hierarchy,
             LDataManager* lag_manager,
             const double current_time,
             const double dt,
             std::vector<double>& params,
             const int loop_num)
{
    const double& angw = params[0];
    double rot = params[1];
    const int finest_ln = hierarchy->getFinestLevelNumber();
    Pointer<LData> D_data = lag_manager->getLData("D", finest_ln);
    Pointer<LData> X_data = lag_manager->getLData("X", finest_ln);
    Vec D_vec = D_data->getVec();
    Vec X_vec = X_data->getVec();
    double* D_vals;
    int ierr = VecGetArray(D_vec, &D_vals);
    IBTK_CHKERRQ(ierr);
    double* X_vals;
    ierr = VecGetArray(X_vec, &X_vals);
    IBTK_CHKERRQ(ierr);
    // Get mesh
    Pointer<LMesh> l_mesh = lag_manager->getLMesh(finest_ln);
    const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
    // Loop over mesh nodes and update triads
    for (const auto& node : local_nodes)
    {
        const int& curr_idx = node->getLagrangianIndex();
        if (curr_idx == 0)
        {
            const int petsc_idx = node->getLocalPETScIndex();
            const int D1_offset = 0;
            Eigen::Map<Vector3d> D1(&D_vals[petsc_idx * 3 * 3 + D1_offset]);
            const int D2_offset = 3;
            Eigen::Map<Vector3d> D2(&D_vals[petsc_idx * 3 * 3 + D2_offset]);
            const int D3_offset = 6;
            Eigen::Map<Vector3d> D3(&D_vals[petsc_idx * 3 * 3 + D3_offset]);
            if (current_time >= 0.2) rot = -rot;
            D1(0) = std::cos(angw * dt * loop_num);
            D1(1) = rot * std::sin(angw * dt * loop_num);
            D1(2) = 0.0;
            D2(0) = -rot * std::sin(angw * dt * loop_num);
            D2(1) = std::cos(angw * dt * loop_num);
            D2(2) = 0.0;
            D3(0) = 0.0;
            D3(1) = 0.0;
            D3(2) = 1.0;
        }
    }
    return;
}

void
read_from_input(Pointer<Database> input_db, std::vector<double>& params)
{
    double angw = input_db->getDouble("Ang_Freq");
    double rot = input_db->getDouble("Rot_Dir");
    params.resize(2);
    params[0] = angw;
    params[1] = rot;
    return;
}
