// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBRodForceSpec.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#include "petscmat.h"
#include "petscvec.h"
#include <petsclog.h>

#include "ibamr/namespaces.h" // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Geometry"
#include "unsupported/Eigen/MatrixFunctions"
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <string>
#include <vector>

// Local includes
#include "IBBistableRodForceGen.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_compute_lagrangian_force_and_torque;
static Timer* t_initialize_level_data;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBistableRodForceGen::IBBistableRodForceGen(Pointer<Database> input_db)
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    IBAMR_DO_ONCE(t_compute_lagrangian_force_and_torque =
                      TimerManager::getManager()->getTimer("IBAMR::IBBistableRodForceGen::"
                                                           "computeLagrangianForceAndTorque()");
                  t_initialize_level_data =
                      TimerManager::getManager()->getTimer("IBAMR::IBBistableRodForceGen::initializeLevelData()"););
    return;
} // IBBistableRodForceGen

IBBistableRodForceGen::~IBBistableRodForceGen()
{
    int ierr;
    for (auto& D_next_mat : d_D_next_mats)
    {
        if (D_next_mat)
        {
            ierr = MatDestroy(&D_next_mat);
            IBTK_CHKERRQ(ierr);
        }
    }
    for (auto& X_next_mat : d_X_next_mats)
    {
        if (X_next_mat)
        {
            ierr = MatDestroy(&X_next_mat);
            IBTK_CHKERRQ(ierr);
        }
    }

    for (auto& D_next_2 : d_D_next_2_mats)
    {
        if (D_next_2)
        {
            ierr = MatDestroy(&D_next_2);
            IBTK_CHKERRQ(ierr);
        }
    }
    for (auto& D_prev_mat : d_D_prev_mats)
    {
        if (D_prev_mat)
        {
            ierr = MatDestroy(&D_prev_mat);
            IBTK_CHKERRQ(ierr);
        }
    }
    return;
} // ~IBBistableRodForceGen

void
IBBistableRodForceGen::setUniformBodyForce(IBTK::Vector F, int structure_id, int level_number)
{
    d_uniform_body_force_data[level_number][structure_id] = std::move(F);
    return;
} // setUniformBodyForce

void
IBBistableRodForceGen::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const int level_number,
                                           const double /*init_data_time*/,
                                           const bool /*initial_time*/,
                                           LDataManager* const l_data_manager)
{
    // The primary use of this function is to prepare for the transfer of data.
    // To compute forces and torques, we need to compute derivatives involving director vectors and point motion. This
    // data is stored on each node, and is distributed across processors. To ensure we have the data locally to compute
    // with, we need to communicate the data to processors. Here, we create PETSc matrices that are permutations of the
    // identity matrix to transfer the data for us. When needed in other functions, we make appropriate calls to matrix
    // multiplication routines with the matrices set up in this function.
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    IBAMR_TIMER_START(t_initialize_level_data);

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    int ierr;

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int level_num = level->getLevelNumber();
    const int new_size = std::max(level_num + 1, static_cast<int>(d_is_initialized.size()));

    d_D_next_mats.resize(new_size);
    d_D_next_2_mats.resize(new_size);
    d_D_prev_mats.resize(new_size);
    d_X_next_mats.resize(new_size);
    d_petsc_curr_node_idxs.resize(new_size);
    d_petsc_next_2_node_idxs.resize(new_size);
    d_petsc_prev_node_idxs.resize(new_size);
    d_petsc_next_node_idxs.resize(new_size);
    d_material_params.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_next_mat = d_D_next_mats[level_num];
    Mat& X_next_mat = d_X_next_mats[level_num];
    Mat& D_next_2_mat = d_D_next_2_mats[level_num];
    Mat& D_prev_mat = d_D_prev_mats[level_num];
    std::vector<int>& petsc_curr_node_idxs = d_petsc_curr_node_idxs[level_num];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_num];
    std::vector<int>& petsc_next_2_node_idxs = d_petsc_next_2_node_idxs[level_num];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_num];
    std::vector<std::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& material_params =
        d_material_params[level_num];

    if (D_next_mat)
    {
        ierr = MatDestroy(&D_next_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (X_next_mat)
    {
        ierr = MatDestroy(&X_next_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (D_next_2_mat)
    {
        ierr = MatDestroy(&D_next_2_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (D_prev_mat)
    {
        ierr = MatDestroy(&D_prev_mat);
        IBTK_CHKERRQ(ierr);
    }
    petsc_curr_node_idxs.clear();
    petsc_next_node_idxs.clear();
    petsc_next_2_node_idxs.clear();
    petsc_prev_node_idxs.clear();
    material_params.clear();

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Determine the "next", "previous", and "next_2" node indices for all rods associated with the
    // present MPI process.
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBRodForceSpec* const force_spec = node_idx->getNodeDataItem<IBRodForceSpec>();
        if (force_spec)
        {
            const int& curr_idx = node_idx->getLagrangianIndex();
            const unsigned int num_rods = force_spec->getNumberOfRods();
#if !defined(NDEBUG)
            TBOX_ASSERT(curr_idx == force_spec->getMasterNodeIndex());
#endif
            const std::vector<int>& next_idxs = force_spec->getNextNodeIndices();
            const std::vector<std::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& params =
                force_spec->getMaterialParams();
#if !defined(NDEBUG)
            TBOX_ASSERT(num_rods == next_idxs.size());
#endif
            for (unsigned int k = 0; k < num_rods; ++k)
            {
                petsc_curr_node_idxs.push_back(curr_idx);
                petsc_next_node_idxs.push_back(next_idxs[k]);
                material_params.push_back(params[k]);
                if (curr_idx > 0) petsc_prev_node_idxs.push_back(curr_idx - 1);
                if (curr_idx + 2 < static_cast<int>(l_data_manager->getNumberOfNodes(level_num)))
                    petsc_next_2_node_idxs.push_back(curr_idx + 2);
            }
        }
    }

    // Map the Lagrangian node indices to the PETSc indices corresponding to the
    // present data distribution.
    l_data_manager->mapLagrangianToPETSc(petsc_curr_node_idxs, level_num);
    l_data_manager->mapLagrangianToPETSc(petsc_next_node_idxs, level_num);
    l_data_manager->mapLagrangianToPETSc(petsc_prev_node_idxs, level_num);
    l_data_manager->mapLagrangianToPETSc(petsc_next_2_node_idxs, level_num);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_num);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_num);

    // Determine the non-zero structure for the matrices.
    const int local_sz = static_cast<int>(petsc_curr_node_idxs.size());
    const int prev_local_sz = static_cast<int>(petsc_prev_node_idxs.size());
    const int next_local_sz = static_cast<int>(petsc_next_2_node_idxs.size());
    std::vector<int> next_d_nz(local_sz, 1), next_o_nz(local_sz, 0);
    std::vector<int> prev_d_nz(prev_local_sz, 1), prev_o_nz(prev_local_sz, 0);
    std::vector<int> next_2_d_nz(next_local_sz, 1), next_2_o_nz(next_local_sz, 0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& next_idx = petsc_next_node_idxs[k];
        if (next_idx >= global_node_offset && next_idx < global_node_offset + num_local_nodes)
        {
            ++next_d_nz[k]; // a "local"    next index
        }
        else
        {
            ++next_o_nz[k]; // a "nonlocal" next index
        }
    }
    for (int k = 0; k < prev_local_sz; ++k)
    {
        const int& prev_idx = petsc_prev_node_idxs[k];
        if (prev_idx >= global_node_offset && prev_idx < global_node_offset + num_local_nodes)
        {
            ++prev_d_nz[k];
        }
        else
        {
            ++prev_o_nz[k];
        }
    }
    for (int k = 0; k < next_local_sz; ++k)
    {
        const int& next_2_idx = petsc_next_2_node_idxs[k];
        if (next_2_idx >= global_node_offset && next_2_idx < global_node_offset + num_local_nodes)
        {
            ++next_2_d_nz[k];
        }
        else
        {
            ++next_2_o_nz[k];
        }
    }

    // Create new MPI block AIJ matrices and set the values of the non-zero
    // entries.
    {
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                             3 * 3,
                             3 * 3 * local_sz,
                             3 * 3 * num_local_nodes,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             0,
                             local_sz ? &next_d_nz[0] : NULL,
                             0,
                             local_sz ? &next_o_nz[0] : NULL,
                             &D_next_mat);
        IBTK_CHKERRQ(ierr);
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                             3 * 3,
                             3 * 3 * prev_local_sz,
                             3 * 3 * num_local_nodes,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             0,
                             prev_local_sz ? &prev_d_nz[0] : NULL,
                             0,
                             prev_local_sz ? &prev_o_nz[0] : NULL,
                             &D_prev_mat);
        IBTK_CHKERRQ(ierr);
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                             3 * 3,
                             3 * 3 * next_local_sz,
                             3 * 3 * num_local_nodes,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             0,
                             next_local_sz ? &next_2_d_nz[0] : NULL,
                             0,
                             next_local_sz ? &next_2_o_nz[0] : NULL,
                             &D_next_2_mat);
        IBTK_CHKERRQ(ierr);

        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> curr_vals(
            Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> next_vals(
            Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> prev_vals(
            Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> next_2_vals(
            Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
        for (unsigned int d = 0; d < 3 * 3; ++d)
        {
            curr_vals(d, d) = 0.0;
            next_vals(d, d) = +1.0;
            prev_vals(d, d) = +1.0;
            next_2_vals(d, d) = +1.0;
        }

        int i_offset;

        ierr = MatGetOwnershipRange(D_next_mat, &i_offset, NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= 3 * 3;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_curr = petsc_curr_node_idxs[k];
            int j_next = petsc_next_node_idxs[k];
            ierr = MatSetValuesBlocked(D_next_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(D_next_mat, 1, &i, 1, &j_next, next_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }

        ierr = MatGetOwnershipRange(D_next_2_mat, &i_offset, NULL);
        i_offset /= 3 * 3;
        IBTK_CHKERRQ(ierr);
        for (int k = 0; k < next_local_sz; ++k)
        {
            int i = i_offset + k;
            int j_curr = petsc_curr_node_idxs[k];
            int j_next_2 = petsc_next_2_node_idxs[k];
            ierr = MatSetValuesBlocked(D_next_2_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(D_next_2_mat, 1, &i, 1, &j_next_2, next_2_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }

        ierr = MatGetOwnershipRange(D_prev_mat, &i_offset, NULL);
        i_offset /= 3 * 3;
        IBTK_CHKERRQ(ierr);
        for (int k = 0; k < prev_local_sz; ++k)
        {
            int i = i_offset + k;
            int j_prev = petsc_prev_node_idxs[k];
            int j_curr = petsc_curr_node_idxs[k];
            ierr = MatSetValuesBlocked(D_prev_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(D_prev_mat, 1, &i, 1, &j_prev, prev_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    {
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                             NDIM,
                             NDIM * local_sz,
                             NDIM * num_local_nodes,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             0,
                             local_sz ? &next_d_nz[0] : NULL,
                             0,
                             local_sz ? &next_o_nz[0] : NULL,
                             &X_next_mat);
        IBTK_CHKERRQ(ierr);

        Matrix curr_vals(Matrix::Zero());
        Matrix next_vals(Matrix::Zero());
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            curr_vals(d, d) = 0.0;
            next_vals(d, d) = +1.0;
        }

        int i_offset;

        ierr = MatGetOwnershipRange(X_next_mat, &i_offset, NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= NDIM;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_curr = petsc_curr_node_idxs[k];
            int j_next = petsc_next_node_idxs[k];
            ierr = MatSetValuesBlocked(X_next_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(X_next_mat, 1, &i, 1, &j_next, next_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(D_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(D_next_2_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(D_prev_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(X_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_2_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_prev_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(X_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_num] = true;

    IBAMR_TIMER_STOP(t_initialize_level_data);
    return;
} // initializeLevelData

void
IBBistableRodForceGen::computeLagrangianForceAndTorque(Pointer<LData> F_data,
                                                       Pointer<LData> N_data,
                                                       Pointer<LData> X_data,
                                                       Pointer<LData> D_data,
                                                       const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                       const int level_number,
                                                       const double data_time,
                                                       LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;
    d_torque.setZero();

    IBAMR_TIMER_START(t_compute_lagrangian_force_and_torque);

#if !defined(NDEBUG)
    TBOX_ASSERT(level_number < static_cast<int>(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    const std::pair<int, int>& helix_idxs = l_data_manager->getLagrangianStructureIndexRange(0, level_number);

    const int global_offset = l_data_manager->getGlobalNodeOffset(level_number);

    int ierr;

    // We need to gather all the data needed to do computations. That data might be on other processors. We do a matrix
    // multiplication to gather the data needed on this processor.

    // Create appropriately sized temporary vectors.
    int i_start, i_stop;

    ierr = MatGetOwnershipRange(d_D_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec = D_data->getVec();
    Vec D_next_vec;

    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &D_next_vec);
    IBTK_CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(d_D_next_2_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_next_2_vec;

    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &D_next_2_vec);
    IBTK_CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(d_D_prev_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_prev_vec;

    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &D_prev_vec);
    IBTK_CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(d_X_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec X_vec = X_data->getVec();
    Vec X_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &X_next_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_next_mats[level_number], D_vec, D_next_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_D_next_2_mats[level_number], D_vec, D_next_2_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_D_prev_mats[level_number], D_vec, D_prev_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_X_next_mats[level_number], X_vec, X_next_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the rod forces acting on the nodes of the Lagrangian mesh.
    double* D_vals;
    ierr = VecGetArray(D_vec, &D_vals);
    IBTK_CHKERRQ(ierr);

    double* D_next_vals;
    ierr = VecGetArray(D_next_vec, &D_next_vals);
    IBTK_CHKERRQ(ierr);

    double* D_prev_vals;
    ierr = VecGetArray(D_prev_vec, &D_prev_vals);
    IBTK_CHKERRQ(ierr);

    double* D_next_2_vals;
    ierr = VecGetArray(D_next_2_vec, &D_next_2_vals);
    IBTK_CHKERRQ(ierr);

    double* X_vals;
    ierr = VecGetArray(X_vec, &X_vals);
    IBTK_CHKERRQ(ierr);

    double* X_next_vals;
    ierr = VecGetArray(X_next_vec, &X_next_vals);
    IBTK_CHKERRQ(ierr);

    std::vector<int>& petsc_curr_node_idxs = d_petsc_curr_node_idxs[level_number];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_number];
    const std::vector<std::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& material_params =
        d_material_params[level_number];

    const size_t local_sz = petsc_curr_node_idxs.size();
    std::vector<double> F_curr_node_vals(NDIM * local_sz, 0.0);
    std::vector<double> N_curr_node_vals(NDIM * local_sz, 0.0);
    std::vector<double> F_next_node_vals(NDIM * local_sz, 0.0);
    std::vector<double> N_next_node_vals(NDIM * local_sz, 0.0);

    std::vector<int> lag_node_idxs(petsc_curr_node_idxs);
    l_data_manager->mapPETScToLagrangian(lag_node_idxs, level_number);
    bool contains_pt0 = false;
    unsigned int kk = -1;
    // Forces are computed on nodes. We compute quantities at the center of a rod, then discretize the derivative so
    // that we compute forces on the nodes.
    for (unsigned int k = 0; k < local_sz; ++k)
    {
        // Note that on the last rod, we compute the regularization term in a one-sided manner. In the below notation,
        // this means that D2_next_2 is not needed (further, it does not exist!). This causes issues in our data
        // structures. The array the keeps track of the D2_next_2 data is now one vector shorter than the rest of the
        // arrays. To complicate matters, because of the PETSc ordering of indices, we aren't guaranteed that the last
        // rod will be at the end of the loop. Therefore, we use a separate index to keep track of which vector we need
        // to pull out for D2_next_2.
        const bool last_pt = lag_node_idxs[k] == (helix_idxs.second - 2);
        if (!last_pt) kk++;
        const bool on_hook = lag_node_idxs[k] < d_hook_length;
        const double ds = material_params[k][0];
        const double a1 = material_params[k][1];
        const double a2 = material_params[k][2];
        const double a3 = material_params[k][3];
        const double b1 = material_params[k][4];
        const double b2 = material_params[k][5];
        const double b3 = material_params[k][6];
        const double kappa1 = material_params[k][7];
        const double kappa2 = material_params[k][8];
        const double tau1 = material_params[k][9];
        const double tau2 = on_hook ? d_tau_2_hook : d_tau_2_main;
        const double gamma = d_gamma;
        // I'm not sure what this is for? Is this a bug??
        // I'm inclined to believe we need another index to keep track of the previous index, similar to that of the
        // next_2 index. There may be something related to the fact that regularization is not used on the hook.
        contains_pt0 = contains_pt0 || (lag_node_idxs[k] == 0);
        // Compute the forces applied by the rod to the "current" and "next"
        // nodes.
        // The following notation is used:
        // Di : The ith director vector at the beginning of the rod
        // Di_next: The ith director vector at the end of the rod (equivalently: beginning of next rod)
        // Di_next_2: The ith director vector at the end of the next rod
        // Di_prev: The ith director vector at the beginning of the previous rod.
        const int D1_offset = 0;
        Eigen::Map<const Vector3d> D1(&D_vals[(petsc_curr_node_idxs[k] - global_offset) * 3 * 3 + D1_offset]);
        Eigen::Map<const Vector3d> D1_next(&D_next_vals[k * 3 * 3 + D1_offset]);
        Vector3d D1_prev(Vector3d::Zero());
        Vector3d D1_next_2(Vector3d::Zero());
        if (!on_hook)
        {
            if (!last_pt) D1_next_2 = Eigen::Map<Vector3d>(&D_next_2_vals[kk * 3 * 3 + D1_offset]);
            D1_prev = contains_pt0 ? Eigen::Map<Vector3d>(&D_prev_vals[(k - 1) * 3 * 3 + D1_offset]) :
                                     Eigen::Map<Vector3d>(&D_prev_vals[k * 3 * 3 + D1_offset]);
        }

        const int D2_offset = 3;
        Eigen::Map<const Vector3d> D2(&D_vals[(petsc_curr_node_idxs[k] - global_offset) * 3 * 3 + D2_offset]);
        Eigen::Map<const Vector3d> D2_next(&D_next_vals[k * 3 * 3 + D2_offset]);
        Vector3d D2_prev(Vector3d::Zero());
        Vector3d D2_next_2(Vector3d::Zero());
        if (!on_hook)
        {
            if (!last_pt) D2_next_2 = Eigen::Map<Vector3d>(&D_next_2_vals[kk * 3 * 3 + D2_offset]);
            D2_prev = contains_pt0 ? Eigen::Map<Vector3d>(&D_prev_vals[(k - 1) * 3 * 3 + D2_offset]) :
                                     Eigen::Map<Vector3d>(&D_prev_vals[k * 3 * 3 + D2_offset]);
        }

        const int D3_offset = 6;
        Eigen::Map<const Vector3d> D3(&D_vals[(petsc_curr_node_idxs[k] - global_offset) * 3 * 3 + D3_offset]);
        Eigen::Map<const Vector3d> D3_next(&D_next_vals[k * 3 * 3 + D3_offset]);
        Vector3d D3_prev(Vector3d::Zero());
        Vector3d D3_next_2(Vector3d::Zero());
        if (!on_hook)
        {
            if (!last_pt) D3_next_2 = Eigen::Map<Vector3d>(&D_next_2_vals[kk * 3 * 3 + D3_offset]);
            D3_prev = contains_pt0 ? Eigen::Map<Vector3d>(&D_prev_vals[(k - 1) * 3 * 3 + D3_offset]) :
                                     Eigen::Map<Vector3d>(&D_prev_vals[k * 3 * 3 + D3_offset]);
        }

        Eigen::Map<const Vector3d> X(&X_vals[(petsc_curr_node_idxs[k] - global_offset) * NDIM]);
        Eigen::Map<const Vector3d> X_next(&X_next_vals[k * NDIM]);
        std::array<Eigen::Map<const Vector3d>*, 3> D = { { &D1, &D2, &D3 } };
        std::array<Eigen::Map<const Vector3d>*, 3> D_next = { { &D1_next, &D2_next, &D3_next } };
        std::array<Vector3d*, 3> D_next_2 = { { &D1_next_2, &D2_next_2, &D3_next_2 } };
        std::array<Vector3d*, 3> D_prev = { { &D1_prev, &D2_prev, &D3_prev } };

        // We need director vectors at rod midpoints. "Interpolate" between the first and second using the square root
        // of the rotation matrix.
        Matrix3d A(Matrix3d::Zero());
        Matrix3d A_next(Matrix3d::Zero());
        Matrix3d A_prev(Matrix3d::Zero());
        for (int i = 0; i < 3; ++i)
        {
            A += (*D_next[i]) * (*D[i]).transpose();
            if (!on_hook)
            {
                if (!last_pt) A_next += (*D_next_2[i]) * (*D_next[i]).transpose();
                A_prev += (*D[i]) * (*D_prev[i]).transpose();
            }
        }
        Matrix3d sqrt_A = A.sqrt();

        Matrix3d sqrt_A_next(Matrix3d::Zero());
        Matrix3d sqrt_A_prev(Matrix3d::Zero());
        if (!on_hook)
        {
            if (!last_pt) sqrt_A_next = A_next.sqrt();
            sqrt_A_prev = A_prev.sqrt();
        }
        Vector3d D1_half, D2_half, D3_half;
        Vector3d D2_half_next;
        Vector3d D2_half_prev;
        std::array<Vector3d*, 3> D_half = { { &D1_half, &D2_half, &D3_half } };
        // Now compute the director vectors at rod mid-points.
        for (int i = 0; i < 3; ++i) *D_half[i] = sqrt_A * (*D[i]);
        // Note we only need D2 on the previous and next rods.
        if (!on_hook)
        {
            if (!last_pt) D2_half_next = sqrt_A_next * D2_next;
            D2_half_prev = sqrt_A_prev * D2_prev;
        }

        const Vector3d dX_ds((X_next - X) / ds);
        const double F1 = b1 * D1_half.dot(dX_ds);
        const double F2 = b2 * D2_half.dot(dX_ds);
        const double F3 = b3 * (D3_half.dot(dX_ds) - 1.0);
        const Vector3d F_half = F1 * D1_half + F2 * D2_half + F3 * D3_half;

        const Vector3d dD1_ds((D1_next - D1) / ds);
        Vector3d dD1_ds_next, dD1_ds_prev;
        if (!on_hook)
        {
            if (!last_pt) dD1_ds_next = Vector3d((D1_next_2 - D1_next) / ds);
            dD1_ds_prev = Vector3d((D1 - D1_prev) / ds);
        }
        const Vector3d dD2_ds((D2_next - D2) / ds);
        const Vector3d dD3_ds((D3_next - D3) / ds);

        const double N1 = a1 * (dD2_ds.dot(D3_half) - kappa1);
        const double N2 = a2 * (dD3_ds.dot(D1_half) - kappa2);
        double omega_3 = dD1_ds.dot(D2_half);
        double N3 = a3 * (omega_3 - tau1) * (omega_3 - tau2) * (omega_3 - 0.5 * (tau1 + tau2));
        if (!on_hook)
        {
            // We only include the regularization term on the body of the flagellum.
            if (last_pt)
            {
                // We are at the boundary. Use a one-sided difference
                const double omega_3_prev = dD1_ds_prev.dot(D2_half_prev);
                N3 -= gamma * gamma / (ds * ds) * (omega_3_prev - omega_3);
            }
            else
            {
                // Use a centered difference.
                const double omega_3_prev = dD1_ds_prev.dot(D2_half_prev);
                const double omega_3_next = dD1_ds_next.dot(D2_half_next);
                N3 -= gamma * gamma / (ds * ds) * (omega_3_next - 2.0 * omega_3 + omega_3_prev);
            }
        }
        const Vector3d N_half = N1 * D1_half + N2 * D2_half + N3 * D3_half;

        // Here we compute FORCES (torques), not force (torque) densities. Essentially, we discretise the derivative,
        // then multiply by ds to convert force densities. The ds cancels.
        Eigen::Map<Vector3d> F_curr(&F_curr_node_vals[k * NDIM]);
        Eigen::Map<Vector3d> F_next(&F_next_node_vals[k * NDIM]);
        F_curr = F_half;
        F_next = -F_half;

        Eigen::Map<Vector3d> N_curr(&N_curr_node_vals[k * NDIM]);
        Eigen::Map<Vector3d> N_next(&N_next_node_vals[k * NDIM]);
        N_curr = N_half + 0.5 * ((X_next - X)).cross(F_half);
        N_next = -N_half + 0.5 * ((X_next - X)).cross(F_half);

        // If we are on lag index 0, output the torque
        if (lag_node_idxs[k] == 0)
            for (int d = 0; d < NDIM; ++d) d_torque[d] = N_curr[d];
    }

    SAMRAI_MPI::sumReduction(d_torque.data(), NDIM);

    ierr = VecRestoreArray(D_vec, &D_vals);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_next_vec, &D_next_vals);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&D_next_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_next_2_vec, &D_next_2_vals);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&D_next_2_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_prev_vec, &D_prev_vals);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&D_prev_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_vec, &X_vals);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_next_vec, &X_next_vals);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X_next_vec);
    IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getVec();
    Vec N_vec = N_data->getVec();
    if (local_sz > 0)
    {
        ierr = VecSetValuesBlocked(F_vec,
                                   static_cast<int>(petsc_curr_node_idxs.size()),
                                   &petsc_curr_node_idxs[0],
                                   &F_curr_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(F_vec,
                                   static_cast<int>(petsc_next_node_idxs.size()),
                                   &petsc_next_node_idxs[0],
                                   &F_next_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(N_vec,
                                   static_cast<int>(petsc_curr_node_idxs.size()),
                                   &petsc_curr_node_idxs[0],
                                   &N_curr_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(N_vec,
                                   static_cast<int>(petsc_next_node_idxs.size()),
                                   &petsc_next_node_idxs[0],
                                   &N_next_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(F_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(N_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(N_vec);
    IBTK_CHKERRQ(ierr);

    computeLagrangianBodyForce(F_data, hierarchy, level_number, data_time, l_data_manager);

    IBAMR_TIMER_STOP(t_compute_lagrangian_force_and_torque);
    return;
} // computeLagrangianForceAndTorque

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBBistableRodForceGen::getFromInput(Pointer<Database> db)
{
    if (db)
    {
        d_tau_2_main = db->getDouble("tau_2_main");
        d_tau_2_hook = db->getDouble("tau_2_hook");
        d_gamma = db->getDouble("gamma");
        d_hook_length = db->getInteger("hook_length");
    }
    return;
} // getFromInput

void
IBBistableRodForceGen::computeLagrangianBodyForce(Pointer<LData> F_data,
                                                  Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                  const int level_number,
                                                  const double /*data_time*/,
                                                  LDataManager* const l_data_manager)
{
    if (d_uniform_body_force_data[level_number].empty()) return;

    double* const F_node = F_data->getLocalFormVecArray()->data();
    SAMRAI::tbox::Pointer<LMesh> l_mesh = l_data_manager->getLMesh(level_number);
    for (const auto& n : l_mesh->getLocalNodes())
    {
        const auto lag_idx = n->getLagrangianIndex();
        const auto structure_id = l_data_manager->getLagrangianStructureID(lag_idx, level_number);
        if (d_uniform_body_force_data[level_number].count(structure_id))
        {
            const IBTK::Vector& F = d_uniform_body_force_data[level_number][structure_id];
            const int local_petsc_idx = n->getLocalPETScIndex();
            for (int d = 0; d < NDIM; ++d)
            {
                F_node[NDIM * local_petsc_idx + d] += F(d);
            }
        }
    }
    F_data->restoreArrays();
    return;
} // computeLagrangianBodyForce

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
