// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include "ibamr/CFINSForcing.h"
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/GeneralizedIBMethod.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBKirchhoffRodForceGen.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include "update_triad.h"

#include <stdio.h>
// Set up application namespace declarations
#include <unsupported/Eigen/MatrixFunctions>

#include <ibamr/app_namespaces.h>

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

void calculatePitchAndRadius(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                             LDataManager* l_data_manager,
                             const double loop_time,
                             const int interation_num);

void calculateTorque(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                     LDataManager* l_data_manager,
                     const double loop_time,
                     const int iteration_num);

void Build_mpi_type_PR(MPI_Datatype* type);

static std::ofstream pitch_file;

struct struct_PR
{
    int idx;
    double R;
    double P;
};
bool
comparePR(const struct_PR& lhs, const struct_PR& rhs)
{
    return lhs.idx < rhs.idx;
}

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    SAMRAIManager::setMaxNumberPatchDataEntries(512);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool is_from_restart = app_initializer->isFromRestart();
        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<GeneralizedIBMethod> ib_method_ops = new GeneralizedIBMethod(
            "GeneralizedIBMethod", app_initializer->getComponentDatabase("GeneralizedIBMethod"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBKirchhoffRodForceGen> ib_force_and_torque_fcn = new IBKirchhoffRodForceGen();
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);
        ib_method_ops->registerIBKirchhoffRodForceGen(ib_force_and_torque_fcn);

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when
        // necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }
        if (input_db->keyExists("ComplexFluid"))
        {
            Pointer<CFINSForcing> complex_fluid =
                new CFINSForcing("ComplexFluidForcing",
                                 app_initializer->getComponentDatabase("ComplexFluid"),
                                 navier_stokes_integrator,
                                 grid_geometry,
                                 adv_diff_integrator,
                                 visit_data_writer);
            time_integrator->registerBodyForceFunction(complex_fluid);
        }

        std::vector<double> params;
        read_from_input(app_initializer->getComponentDatabase("Motor"), params);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write restart data before starting main time integration loop.
        if (dump_restart_data && !is_from_restart)
        {
            pout << "\nWriting restart files...\n\n";
            RestartManager::getManager()->writeRestartFile(restart_dump_dirname, 0);
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            update_triad(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, dt, params, iteration_num);
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            /*
                        if (last_step)
                        {
                            output_new_files(patch_hierarchy,
               ib_method_ops->getLDataManager(), loop_time, dt, params,
               iteration_num);
                        }
            */
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
                calculatePitchAndRadius(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, iteration_num);
                calculateTorque(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, iteration_num);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy,
                            navier_stokes_integrator,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }

        pout << "\nFinished simulation. Wrapping up...\n\n";

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            LDataManager* l_data_manager,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
    Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    VecView(X_lag_vec, viewer);
    PetscViewerDestroy(&viewer);
    VecDestroy(&X_lag_vec);
    return;
} // output_data

void
calculatePitchAndRadius(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                        LDataManager* l_data_manager,
                        const double loop_time,
                        const int iteration_num)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    const int global_offset = l_data_manager->getGlobalNodeOffset(ln);
    const int local_node_num = l_data_manager->getNumberOfLocalNodes(ln);
    Pointer<LData> D_data = l_data_manager->getLData("D", ln);
    Vec D_vec = D_data->getVec();
    double* D_vals;
    int ierr = VecGetArray(D_vec, &D_vals);
    IBTK_CHKERRQ(ierr);
    Pointer<LMesh> l_mesh = l_data_manager->getLMesh(ln);
    std::pair<int, int> helix_idxs = l_data_manager->getLagrangianStructureIndexRange(0, ln);
    const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
    const std::vector<LNode*>& ghost_nodes = l_mesh->getGhostNodes();

    std::vector<LNode*> nodes = local_nodes;
    std::vector<int> petsc_curr_node_idxs, petsc_next_idxs, lag_idxs;
    std::vector<double> ds_vec;
    petsc_curr_node_idxs.clear();
    petsc_next_idxs.clear();
    lag_idxs.clear();
    ds_vec.clear();
    for (std::vector<LNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        const LNode* const node_idx = *it;
        const int& curr_idx = node_idx->getLagrangianIndex();
        if ((curr_idx >= helix_idxs.first) && (curr_idx <= helix_idxs.second) &&
            node_idx->getNodeDataItem<IBRodForceSpec>())
        {
            const IBRodForceSpec* const force_spec = node_idx->getNodeDataItem<IBRodForceSpec>();
            ds_vec.push_back(force_spec->getMaterialParams()[0][0]);
            lag_idxs.push_back(curr_idx);
            petsc_curr_node_idxs.push_back(curr_idx);
            const std::vector<int>& next_idxs = force_spec->getNextNodeIndices();
            if (next_idxs.size() != 0)
            {
                TBOX_ASSERT(next_idxs.size() == 1);
                petsc_next_idxs.push_back(next_idxs[0]);
            }
            else
            {
                TBOX_ERROR("ROD has no next term");
            }
        }
    }
    l_data_manager->mapLagrangianToPETSc(petsc_curr_node_idxs, ln);
    l_data_manager->mapLagrangianToPETSc(petsc_next_idxs, ln);
    const int local_sz = static_cast<int>(petsc_curr_node_idxs.size());
    std::vector<int> next_d_nz(local_sz, 1), next_o_nz(local_sz, 0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& next_idx = petsc_next_idxs[k];
        if (next_idx >= global_offset && next_idx < global_offset + local_node_num)
        {
            ++next_d_nz[k];
        }
        else
        {
            ++next_o_nz[k];
        }
    }
    Mat D_next_mat;
    ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                         3 * 3,
                         3 * 3 * local_sz,
                         3 * 3 * local_node_num,
                         PETSC_DETERMINE,
                         PETSC_DETERMINE,
                         0,
                         local_sz ? &next_d_nz[0] : NULL,
                         0,
                         local_sz ? &next_o_nz[0] : NULL,
                         &D_next_mat);
    IBTK_CHKERRQ(ierr);
    Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> next_vals(
        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
    Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> curr_vals(
        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
    for (unsigned int d = 0; d < 3 * 3; ++d)
    {
        next_vals(d, d) = 1.0;
        curr_vals(d, d) = 0.0;
    }
    int i_offset;
    ierr = MatGetOwnershipRange(D_next_mat, &i_offset, NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= 3 * 3;
    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_curr = petsc_curr_node_idxs[k];
        int j_next = petsc_next_idxs[k];
        ierr = MatSetValuesBlocked(D_next_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(D_next_mat, 1, &i, 1, &j_next, next_vals.data(), INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(D_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    int i_start, i_stop;
    ierr = MatGetOwnershipRange(D_next_mat, &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &D_next_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(D_next_mat, D_vec, D_next_vec);
    IBTK_CHKERRQ(ierr);
    double* D_next_vals;
    ierr = VecGetArray(D_next_vec, &D_next_vals);
    IBTK_CHKERRQ(ierr);
    const int D1_offset = 0, D2_offset = 3, D3_offset = 6;
    std::vector<struct_PR> pitch_rad_vals(petsc_curr_node_idxs.size());
    for (int k = 0; k < petsc_curr_node_idxs.size(); ++k)
    {
        const int idx1 = petsc_curr_node_idxs[k];
        const double ds = ds_vec[k];

        Eigen::Map<const Vector3d> D1(&D_vals[(idx1 - global_offset) * 3 * 3 + D1_offset]);
        Eigen::Map<const Vector3d> D2(&D_vals[(idx1 - global_offset) * 3 * 3 + D2_offset]);
        Eigen::Map<const Vector3d> D3(&D_vals[(idx1 - global_offset) * 3 * 3 + D3_offset]);

        Eigen::Map<const Vector3d> D1_next(&D_next_vals[k * 3 * 3 + D1_offset]);
        Eigen::Map<const Vector3d> D2_next(&D_next_vals[k * 3 * 3 + D2_offset]);
        Eigen::Map<const Vector3d> D3_next(&D_next_vals[k * 3 * 3 + D3_offset]);

        boost::array<Eigen::Map<const Vector3d>*, 3> D = { { &D1, &D2, &D3 } };
        boost::array<Eigen::Map<const Vector3d>*, 3> D_next = { { &D1_next, &D2_next, &D3_next } };
        Matrix3d A(Matrix3d::Zero());
        for (int i = 0; i < 3; ++i)
        {
            A += (*D_next[i]) * (*D[i]).transpose();
        }
        Matrix3d sqrt_A = A.sqrt();

        Vector3d D1_half, D2_half, D3_half;
        boost::array<Vector3d*, 3> D_half = { { &D1_half, &D2_half, &D3_half } };
        for (int i = 0; i < 3; ++i)
        {
            *D_half[i] = sqrt_A * (*D[i]);
        }

        const Vector3d dD1_ds((D1_next - D1) / ds);
        const Vector3d dD2_ds((D2_next - D2) / ds);
        const Vector3d dD3_ds((D3_next - D3) / ds);

        double om1 = dD2_ds.dot(D3_half), om2 = dD3_ds.dot(D1_half), om3 = dD1_ds.dot(D2_half);

        struct_PR val;
        val.idx = lag_idxs[k];
        double tau = om3;
        double kappa = std::sqrt(om1 * om1 + om2 * om2);
        val.R = kappa / (kappa * kappa + tau * tau);
        val.P = 2 * M_PI * tau / (kappa * kappa + tau * tau);
        pitch_rad_vals[k] = val;
    }
    // Sort vals, Reduce to one processor, print to file.
    MPI_Datatype tuple;
    Build_mpi_type_PR(&tuple);
    int tuple_size;
    MPI_Type_size(tuple, &tuple_size);
    int rank = SAMRAI_MPI::getRank();
    int num_nodes = SAMRAI_MPI::getNodes();
    std::vector<int> vals_on_procs(num_nodes);
    vals_on_procs[rank] = pitch_rad_vals.size();
    SAMRAI_MPI::sumReduction(vals_on_procs.data(), num_nodes);
    int tot = 0;
    for (int i = 0; i < num_nodes; ++i) tot += vals_on_procs[i];

    std::vector<struct_PR> pitch_rad_vals_red;
    pitch_rad_vals_red.resize(tot);

    std::vector<int> displacements(num_nodes);
    displacements[0] = 0;
    for (int i = 1; i < num_nodes; ++i) displacements[i] = displacements[i - 1] + /*tuple_size**/ vals_on_procs[i - 1];

    MPI_Gatherv(pitch_rad_vals.data(),
                pitch_rad_vals.size(),
                tuple,
                pitch_rad_vals_red.data(),
                vals_on_procs.data(),
                displacements.data(),
                tuple,
                0,
                SAMRAI_MPI::getCommunicator());

    if (SAMRAI_MPI::getRank() == 0)
    {
        std::sort(pitch_rad_vals_red.begin(), pitch_rad_vals_red.end(), comparePR);
        ostringstream file_name;
        system("mkdir -p pitch");
        file_name << "pitch/Pitch_" << iteration_num << ".out";
        pitch_file.open(file_name.str().c_str());
        int i = 0;
        for (std::vector<struct_PR>::const_iterator it = pitch_rad_vals_red.begin(); it != pitch_rad_vals_red.end();
             ++it)
        {
            const struct_PR val = *it;
            pitch_file << val.idx << " " << val.R << " " << val.P << "\n";
        }
        pitch_file.close();
    }
    MPI_Type_free(&tuple);
    return;
}

void
Build_mpi_type_PR(MPI_Datatype* type)
{
    int rank = SAMRAI_MPI::getRank();
    struct_PR object;
    int struct_length = 3;
    int blocklengths[struct_length];
    MPI_Datatype types[struct_length];
    MPI_Aint displacements[struct_length];
    blocklengths[0] = blocklengths[1] = blocklengths[2] = blocklengths[3] = 1;
    types[0] = MPI_INT;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_UB;
    displacements[0] = (size_t) & (object.idx) - (size_t)&object;
    displacements[1] = (size_t) & (object.R) - (size_t)&object;
    displacements[2] = (size_t) & (object.P) - (size_t)&object;
    displacements[3] = sizeof(object);
    MPI_Type_create_struct(struct_length, blocklengths, displacements, types, type);
    MPI_Type_commit(type);
    return;
}

void
calculateTorque(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                LDataManager* l_data_manager,
                const double loop_time,
                const int iteration_num)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    const int global_offset = l_data_manager->getGlobalNodeOffset(ln);
    const int local_node_num = l_data_manager->getNumberOfLocalNodes(ln);
    Pointer<LData> N_data = l_data_manager->getLData("N", ln);
    Vec N_vec = N_data->getVec();
    double* N_vals;
    int ierr = VecGetArray(N_vec, &N_vals);
    IBTK_CHKERRQ(ierr);
    Pointer<LMesh> l_mesh = l_data_manager->getLMesh(ln);
    std::pair<int, int> helix_idxs = l_data_manager->getLagrangianStructureIndexRange(0, ln);
    const std::vector<LNode*>& nodes = l_mesh->getLocalNodes();
    std::vector<double> torque(NDIM);
    torque[0] = torque[1] = torque[2] = 0.0;
    int petsc_idx;

    for (std::vector<LNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        const LNode* const node_idx = *it;
        const int curr_idx = node_idx->getLagrangianIndex();
        if (curr_idx == 1)
        {
            petsc_idx = node_idx->getGlobalPETScIndex();
            Eigen::Map<const Vector3d> N(&N_vals[(petsc_idx - global_offset) * 3]);
            for (int i = 0; i < NDIM; ++i)
            {
                torque[i] = N(i);
            }
        }
    }
    SAMRAI_MPI::sumReduction(&torque[0], NDIM);

    if (SAMRAI_MPI::getRank() == 0)
    {
        std::ofstream torque_file;
        ostringstream file_name;
        system("mkdir -p torque");
        file_name << "torque/Torque_" << iteration_num << ".out";
        torque_file.open(file_name.str().c_str());
        torque_file << torque[0] << " " << torque[1] << " " << torque[2] << "\n";
        torque_file.close();
    }
    return;
}
