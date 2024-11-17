
/**
*    \author Mohit Mehindratta / Hakim Amer
*    \date   1/9/2022
Acado ccode for mpc implementation for controlling the bluerov2 in 2-D. The model used can be found in the paper "Super-twisting integral sliding mode control for trajectory tracking of an
Unmanned Underwater Vehicle"
*/


#include <acado_code_generation.hpp>
#include <ros/ros.h>
#include <ros/package.h>
#include <boost/algorithm/string.hpp>

USING_NAMESPACE_ACADO


int main()
{

    //State Variables:
    DifferentialState x;  // the body position w.r.t X_I
    DifferentialState y;  // the body position w.r.t Y_I
    DifferentialState z;  // the body position w.r.t Z_I

    DifferentialState u;  // the translation velocity along X_B
    DifferentialState v;  // the translation velocity along Y_B
    DifferentialState w;  // the translation velocity along Z_B
    
    DifferentialState psi;  // yaw angle 
    DifferentialState r;   // yaw rate

    
    DifferentialState aux_state_px;
    DifferentialState aux_state_py;
    DifferentialState aux_state_pz;


    DifferentialState aux_state_ox;
    DifferentialState aux_state_oy;
    DifferentialState aux_state_oz;

    DifferentialState aux_slack;

    OnlineData pe_x;  // entagelemnent point_x
    OnlineData pe_y;  // entagelemnent point_y
    OnlineData pe_z;  // entagelemnent point_z

    OnlineData ob_x;  // obstacle centre_x
    OnlineData ob_y;  // obstacle centre_y
    OnlineData ob_z;  // obstacle centre_z



    //OnlineData Fx_dist;  // the external disturbance force along X_B
    //OnlineData Fy_dist;  // the external disturbance force along Y_B
    //OnlineData Fz_dist;  // the external disturbance force along Z_B
    // MPC control input 
    Control X;  // Force along X_B
    Control Y;  // Force along Y_B
    Control Z;  // Force along Z_B
    Control M_z;  // Torque about Z_B (Yawing moment)

    Control slack;

     

    // BlueROV2 Model Parameters 
    
    //const double F_bouy = 114.8; // Bouyancy force (N)
    const double m = 11.4;    // BlueROV2 mass (kg)  
    const double g = 9.82;  // gravitational field strength (m/s^2)

    const double F_bouy = 1026 * 0.0115 * g; // Bouyancy force (N)
    const double eps = 0.00001;
    //const double F_bouy = 114.8; // Buoyancy force (N)

    const double X_ud = -2.6 ; // Added mass in x direction (kg)
    const double Y_vd = -18.5 ; // Added mass in y direction (kg)
    const double Z_wd = -13.3 ; // Added mass in z direction (kg)
    const double N_rd = -0.28 ; // Added mass for rotation about z direction (kg)

    const double I_xx = 0.21 ; // Moment of inertia (kg.m^2)
    const double I_yy = 0.245 ; // Moment of inertia (kg.m^2)
    const double I_zz = 0.245 ; // Moment of inertia (kg.m^2)

    const double X_u = -0.09 ; // Linear damping coefficient in x direction (N.s/m)
    const double Y_v  = -0.26 ; // Linear damping coefficient  in y direction (N.s/m)
    const double Z_w = -0.19; // Linear damping coefficient  in z direction (N.s/m)
    const double N_r = -4.64 ;  // Linear damping coefficient for rotation about z direction (N.s/rad)

    const double X_uc = -34.96 ; // quadratic damping coefficient in x direction (N.s^2/m^2)
    const double Y_vc = -103.25 ; // quadratic damping coefficient  in y direction (N.s^2/m^2)
    const double Z_wc = -74.23 ; // quadratic damping coefficient  in z direction (N.s^2/m^2)
    const double N_rc = - 0.43 ; // quadratic damping coefficient for rotation about z direction (N.s^2/rad^2)
    // Model equations: 2-D model, assuming no roll or pitch
    DifferentialEquation f;

    f << dot(x) == cos(psi) * u - sin(psi) * v;
    f << dot(y) == sin(psi) * u +  cos(psi) * v;
    f << dot(z) ==  w;
 
    f << dot(u) == (X + (m * v - Y_vd * v) * r + (X_u + X_uc *sqrt( u * u + eps) ) * u)/(m - X_ud);
    f << dot(v) == (Y - (m * u - X_ud * u) * r + (Y_v + Y_vc *sqrt( v * v + eps) ) * v)/(m - Y_vd)  ;
    f << dot(w) == (Z + (Z_w + Z_wc * sqrt(w * w + eps)) * w + (m * g - F_bouy))/(m - Z_wd)  ;

    f << dot(psi) ==  r;
    f << dot(r) == (M_z - (m * v - Y_vd * v) * u - (X_ud * u - m * u) * v + (N_r + N_rc * sqrt(r * r + eps)) * r)/(I_zz - N_rd);

    f << dot(aux_state_px) == pe_x ;
    f << dot(aux_state_py) == pe_y;
    f << dot(aux_state_pz) == pe_z;

    f << dot(aux_state_ox) == ob_x;
    f << dot(aux_state_oy) == ob_y;
    f << dot(aux_state_oz) == ob_z;

    f << dot(aux_slack) == slack;


    // calculate entagelement cost

    //IntermediateState p_ref_x = (pe_x - ob_x);  //
    //IntermediateState p_ref_y = (pe_y - ob_y);  //

    //IntermediateState p_traj_x = (x - ob_x);  //
    //IntermediateState p_traj_y = (x - ob_y);  //
    
    //IntermediateState s =  p_ref_x * p_traj_x + p_ref_y * p_traj_y;


    //IntermediateState s =  (x - pe_x) * (x - pe_x) + (y - pe_y) * (y - pe_y);

    //IntermediateState p_ref_norm = sqrt(p_ref_x * p_ref_x + p_ref_y  * p_ref_y);  // Constant added for numerical stability
    //IntermediateState p_traj_norm = sqrt(p_traj_x * p_traj_x + p_traj_y  * p_traj_y + 0.00001);  // Constant added for numerical stability



    // Reference functions and weighting matrices:
    Function h, hN;

    //IntermediateState distance_to_obstacle;
    //distance_to_obstacle = sqrt((x - ob_x) * (x - ob_x) + (y - ob_y) * (y - ob_y));

    //h << x << y << z << u << v << w << psi << s - 1 << r << X << Y<< Z << M_z;
    
    h << x << y << z << u << v << w << psi << r << x << X << Y<< Z << M_z << slack;
    hN << x << y << z << u << v << w << psi << r;


    BMatrix W = eye<bool>(h.getDim());
    BMatrix WN = eye<bool>(hN.getDim());

    //
    // Optimal Control Problem
    //
    double N = 30;
    double Ts = 0.01;
    OCP ocp(0.0, N * Ts, N);

    ocp.subjectTo(f);

    ocp.minimizeLSQ(W, h);
    ocp.minimizeLSQEndTerm(WN, hN);


   // Constraints on Inputs 
    ocp.subjectTo(-80 <= X <= 80);
    ocp.subjectTo(-80  <= Y <= 80);
    ocp.subjectTo(-160 <= Z <= 160);
    ocp.subjectTo(-160 <= M_z <= 160);    //in Nm
    ocp.subjectTo( slack >= 0 );

    //ocp.subjectTo(-160 <= x <= 160);    //in Nm

    double d_min = 2.0;  // Minimum distance from each side of the square

// Constraints to keep (x, y) outside a square centered at (ob_x, ob_y)
//ocp.subjectTo((y - (ob_y + d_min)) >= 0);  // Above the top side
//ocp.subjectTo((y - (ob_y - d_min)) <= 0);  // Below the bottom side
//ocp.subjectTo((x - (ob_x + d_min)) >= 0);  // Right of the right side
//ocp.subjectTo((x - (ob_x - d_min)) <= 0);  // Left of the left side

    //ocp.subjectTo((x - ob_x) * (x - ob_x) + (y - ob_y) * (y - ob_y) + eps - (d_min - slack) * (d_min - slack) >= 0);
    //ocp.subjectTo( slack >= 0 );

   // double eps = 0.0001;  // Small tolerance value
    // Export the code:
    OCPexport mpc(ocp);

    mpc.set(HESSIAN_APPROXIMATION, GAUSS_NEWTON);
    mpc.set(DISCRETIZATION_TYPE, MULTIPLE_SHOOTING);
    mpc.set(INTEGRATOR_TYPE, INT_RK4);     // RK4 for best comprimise between accuracy and computational effort 
    mpc.set(NUM_INTEGRATOR_STEPS, 2 * N);

    mpc.set(SPARSE_QP_SOLUTION, CONDENSING);
    mpc.set(QP_SOLVER, QP_QPOASES);
    mpc.set(MAX_NUM_QP_ITERATIONS, 1000);

    // 	mpc.set( SPARSE_QP_SOLUTION, SPARSE_SOLVER );
    // 	mpc.set( QP_SOLVER, QP_QPDUNES );

    //mpc.set(HOTSTART_QP, YES);

    mpc.set(CG_HARDCODE_CONSTRAINT_VALUES, NO);  // Possible to Change Constraints Afterwards (only with qpOASES)

    mpc.set(GENERATE_TEST_FILE, NO);
    mpc.set(GENERATE_MAKE_FILE, NO);
    mpc.set(GENERATE_MATLAB_INTERFACE, NO);
    mpc.set(GENERATE_SIMULINK_INTERFACE, NO);

    // Optionally set custom module name:
    mpc.set(CG_MODULE_NAME, "nmpc");
    mpc.set(CG_MODULE_PREFIX, "NMPC");


     
    std::string path = ros::package::getPath("acado_ccode_generation");  //Set the name of folder of the acado c-code generator ros node  
    std::string path_dir = path + "/solver/NMPC_BlueRov2_2-D";           // Set the name of the folder that contains the acado c-code solver files
    ROS_INFO("%s", path_dir.c_str());

    try
    {
        ROS_WARN("TRYING TO EXPORT");
        if (mpc.exportCode(path_dir) != SUCCESSFUL_RETURN)
            ROS_ERROR("FAIL EXPORT CODE");
    }
    catch (...)
    {
        ROS_ERROR("FAIL TO EXPORT");
    }

    mpc.printDimensionsQP();

    ROS_WARN("DONE CCODE");

    return EXIT_SUCCESS;
}
