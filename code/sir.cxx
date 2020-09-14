//////////////////////////////////////////////////////////////////////////
////////////////        SIR.cxx         /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         SIR control problem				  ////////////////
//////// Last modified: 14 September 2020                     ////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

////////////////////////////////////////////////////////////////////////////
///////////////////  Define data values ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

typedef struct{

	double R0; // the reproduction number

} Mydata;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   adouble i = final_states[1];

   //return tf; // We want to minimize the time tf needed to bring S (#susceptibles) below a critical value S*.
   return 0.0;

}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
	adouble u = controls[0]; // the control parameter u reduces the infection rate via social distancing

	adouble t = time;

	//if (u >= 0)
	//	return sqrt(u);
	//else
	//	return sqrt(-u);
	
	//return 0.0; 
    //return exp(u);

	return u*u; // quadratic cost function
	
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

   // The SIR model
   // \f[
   // \dot S = -(1-\alpha)SI
   // \dot I = (1-\alpha)SI-I/R_0
   // \f]
 
   adouble sdot, idot;

   adouble s = states[ CINDEX(1) ];
   adouble i = states[ CINDEX(2) ];
   adouble u = controls[ CINDEX(1) ];

   // ODE constraints (SIR model)

   Mydata* md = (Mydata*) workspace->problem->user_data;
   double R0 = md->R0; // the basic reproduction number

   sdot = -R0*(1-u)*s*i;
   idot = R0*(1-u)*s*i - i;

   /////////////////////////////

   derivatives[ CINDEX(1) ] = sdot;
   derivatives[ CINDEX(2) ] = idot;

   // path constraints
   //path[ CINDEX(1) ] = i; // I should stay below a threshold I_thresh (hospital capacity limit)
   //path[ CINDEX(2) ] = u; // \f$ 0 \leq u \leq 1 \f$
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

   adouble s0 = initial_states[ CINDEX(1) ];
   adouble i0 = initial_states[ CINDEX(2) ];
   adouble sf = final_states[ CINDEX(1) ];

   e[ CINDEX(1) ] = s0;
   e[ CINDEX(2) ] = i0;
   e[ CINDEX(3) ] = sf; // we want to bring S (#susceptibles) below a critical value S* at t=tf.

}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "SIR epidemic control";
    problem.outfilename                 = "sir.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2; // S,I
    problem.phases(1).ncontrols 		= 1; // u (alpha in the manuscript)
    problem.phases(1).nevents   		= 3; // S0, I0, Sf
    problem.phases(1).npath     		= 0;
    //problem.phases(1).nodes                     = "[200]"; // TODO: What is a reasonable number here? Maybe check Bachelor's thesis (KMH)
    problem.phases(1).nodes                     = "[500]"; // TODO: What is a reasonable number here? Maybe check Bachelor's thesis (KMH)

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem parameters //////////////////////////////
////////////////////////////////////////////////////////////////////////////


	Mydata* md = (Mydata*) malloc(sizeof(Mydata)); // allocate user data

	double R0 = 3.0; // estimate of the basic reproduction number for COVID-19
	md->R0 = R0; // set the basic reproduction number

	problem.user_data = (void*) md;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double sL = 0.0; 
    double sU = 1.0; 
    double iL = 0.0; 
    double iU = 0.01; // want to keep I below the hospital capacity limit i_thresh

    //double uL = 0.001; // 0 = no control, 1 = zero spreading of virus
    double uL = 0.0; // 0 = no control, 1 = zero spreading of virus
    double uU = 1.0;

    //double i0 = iU; // initial fraction of infected; here we start at Imax
    double i0 = 0.0025; // initial fraction of infected,

    double s0 = 1-i0;
    //double s0 = 0.4;
    double sf = 1/R0; // goal: herd immunity sf = 1/R0

	double TGuess = 50.0; // guess for the time horizon
	
    //double s0 = s_initial;
    //double i0 = i_initial;
    //double sf = s_final;

    problem.phases(1).bounds.lower.states(1) = sL;
    problem.phases(1).bounds.lower.states(2) = iL;

    problem.phases(1).bounds.upper.states(1) = sU;
    problem.phases(1).bounds.upper.states(2) = iU;


    problem.phases(1).bounds.lower.controls(1) = uL;
    problem.phases(1).bounds.upper.controls(1) = uU;


    problem.phases(1).bounds.lower.events(1) = s0;
    problem.phases(1).bounds.lower.events(2) = i0;
    problem.phases(1).bounds.lower.events(3) = sf;

    problem.phases(1).bounds.upper.events(1) = s0;
    problem.phases(1).bounds.upper.events(2) = i0;
    problem.phases(1).bounds.upper.events(3) = sf;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.0; 
    problem.phases(1).bounds.upper.EndTime      = 500.0; //TODO: is this upper bound necessary



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= problem.phases(1).nodes(1);
    int ncontrols           = problem.phases(1).ncontrols;
    int nstates             = problem.phases(1).nstates;

    DMatrix x_guess    =  zeros(nstates,nnodes);

    x_guess(1,colon()) = s0*ones(1,nnodes);
    x_guess(2,colon()) = i0*ones(1,nnodes);
	// TODO: propagate ode as initial guess, or choose heuristic solution from manuscript

    //problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.controls       = 0.5*zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0,TGuess,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 5000;
    //algorithm.nlp_tolerance               = 1.e-2;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.mesh_refinement             = "automatic";
    algorithm.collocation_method = "Legendre";
//  algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-6;
    //algorithm.ode_tolerance               = 1.e-4;
	algorithm.mr_max_iterations = 1;




////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    //lambda = solution.get_dual_costates_in_phase(1);
    //H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

	//string dataDir(argv[1]); // pass outputdir via command line

    //x.Save((dataDir+"x.dat").c_str());
    //u.Save((dataDir+"u.dat").c_str());
    //t.Save((dataDir+"t.dat").c_str());
    //lambda.Save("lambda.dat");
    //H.Save("H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////



	//cerr << dataDir << endl;

    plot(t,x,t,u, problem.name+": states, control", "time (t / tau)", "states, control","S I u");
    plot(t,x(2,colon()), problem.name+": states, control", "time (t / tau)", "states, control","I");

    plot(x(1,colon()),x(2,colon()),problem.name+": I(S)","S/N", "I/N");

    //plot(t,x,t,u,problem.name+": states, control", "time (t*beta)", "states, control","S I u",
                             //"pdf", "sir_control.pdf");

    //plot(t,u,problem.name+": control","time (t*beta)", "control", "u",
                             //"pdf", "sir_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
