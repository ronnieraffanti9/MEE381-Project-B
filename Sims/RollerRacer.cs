//============================================================================
// RollerRacer.cs : Defines derived class for simulating a Roller Racer.
//       Equations of motion are derived in class notes.
//============================================================================
using System;

public class RollerRacer : Simulator
{
    LinAlgEq matrix;    //so we can solve the matrix
    // physical parameters, names are the same as that in the notes
    double m;   // mass of vehicle
    double Ig;  // moment of inertia (vertical axis) about center of mass
    double b;   // distance of com ahead of rear axle
    double c;   // distance of rear contact patch from symmetry axis
    double d;   // caster length
    double h;   // longitudinal distance between center of mass and steer axis

    double rW;  // radius of rear wheel, used for calculating rotation rate
    double rWs; // radius of steered wheel, used for calculating rotation rate

    double kPDelta;  // proportional gain for steer filter
    double kDDelta;  // derivative gain for steer filter
    double deltaDes; // desired steer angle

    double muS;      // static frict coeff, lower bound

    bool simBegun;   // indicates whether simulation has begun

    public RollerRacer() : base(11)
    {
        g = 9.81;
        muS = 0.9;
        SetInertia(25.0 /*mass*/, 0.3 /*radius of gyration*/);
        SetGeometry(1.3 /*wheel base*/, 0.6 /* cg dist from axle*/,
            0.3 /*caster dist*/, 1.0 /*wheel sep*/, 0.5*0.75 /*Rwheel radius*/,
            0.15 /*steered wheel radius*/);
        kPDelta = 10.0;
        kDDelta = 4.0;

        x[0] = 0.0;   // x coordinate of center of mass
        x[1] = 0.0;   // xDot, time derivative of x
        x[2] = 0.0;   // z coordinate of center of mass
        x[3] = 0.0;   // zDot, time derivative of z
        x[4] = 0.0;   // psi, heading angle
        x[5] = 0.0;   // psiDot, time derivative of heading, yaw rate
        x[6] = 0.0;   // rotation angle of left rear wheel
        x[7] = 0.0;   // rotation angle of right rear wheel
        x[8] = 0.0;   // rotation angle of front steered wheel
        x[9] = 0.0;   // delta, steer angle
        x[10] = 0.0;  // deltaDot, steer rate

        SetRHSFunc(RHSFuncRRacer);
        simBegun = false;
    }

    private void RHSFuncRRacer(double[] xx, double t, double[] ff)
    {
        matrix = new LinAlgEq(5);//5 equations for 5 unknowns
        // give names to some state variable so code is easier to read & write
        double xDot = xx[1];
        double zDot = xx[3];
        double psi  = xx[4];
        double psiDot = xx[5];
        double delta = xx[9];
        double deltaDot = xx[10];

        // calculate some trig functions here, so you only have to do it once
        double cosPsi = Math.Cos(psi);
        double sinPsi = Math.Sin(psi);
        double cosDelta = Math.Cos(delta);
        double sinDelta = Math.Sin(delta);
        double cosPsiPlusDelta = Math.Cos(psi + delta);
        double sinPsiPlusDelta = Math.Sin(psi + delta);

        // #### You will do some hefty calculations here

        //set first row of matrix:eq1
        matrix.A[0][0] = - m;//x
        matrix.A[0][1] = 0.0;//z
        matrix.A[0][2] = 0.0;//psi
        matrix.A[0][3] = sinPsi;//Fb
        matrix.A[0][4] = sinPsiPlusDelta;//Ff
        matrix.b[0] = 0;//remaining

        //set second row of matrix:eq2
        matrix.A[1][0] = 0;//x
        matrix.A[1][1] = - m;//z
        matrix.A[1][2] = 0;//psi
        matrix.A[1][3] = cosPsi;//Fb
        matrix.A[1][4] = cosPsiPlusDelta;//Ff
        matrix.b[1] = 0;//remaining

        //set third row of matrix:eq3
        matrix.A[2][0] = 0;//x
        matrix.A[2][1] = 0;//z
        matrix.A[2][2] = - Ig;//psi
        matrix.A[2][3] = b;//Fb
        matrix.A[2][4] = - h * cosDelta - d;//Ff
        matrix.b[2] = 0;//remaining

        //set fourth row of matrix:eq7
        matrix.A[3][0] = sinPsi;//x
        matrix.A[3][1] = cosPsi;//z
        matrix.A[3][2] = b;//psi
        matrix.A[3][3] = 0;//Fb
        matrix.A[3][4] = 0;//Ff
        matrix.b[3] = - xDot*psiDot*cosPsi + zDot*psiDot*sinPsi;//remaining

        //set fifth row of matrix:eq10
        matrix.A[4][0] = sinPsiPlusDelta;//x
        matrix.A[4][1] = cosPsiPlusDelta;//z
        matrix.A[4][2] = - h * cosPsi;//psi
        matrix.A[4][3] = 0;//Fb
        matrix.A[4][4] = 0;//Ff
        double deltaDDot = -kDDelta*deltaDot -kPDelta*(delta - deltaDes);
        matrix.b[4] = - (psiDot + deltaDDot)*d - xDot*(psiDot + deltaDot)*cosPsiPlusDelta + zDot*(psiDot + deltaDot)*sinPsiPlusDelta;//remaining

        matrix.SolveGauss();//solve the matrix

        // apply solutions to the simulation
        ff[0] = xDot;// x coordinate of center of mass
        ff[1] = matrix.sol[0];// time derivative of x
        ff[2] = zDot;// z coordinate of center of mass
        ff[3] = matrix.sol[1];// time derivative of z
        ff[4] = psiDot;// heading angle
        ff[5] = matrix.sol[2];// yaw rate
        ff[6] = -(xDot*cosPsi - zDot*sinPsi - c*psiDot) / rW;// rotation angle of left rear wheel
        ff[7] = -(xDot*cosPsi - zDot*sinPsi + c*psiDot) / rW;// rotation angle of right rear wheel
        ff[8] = -(xDot*cosPsiPlusDelta - zDot*sinPsiPlusDelta + h*psiDot*sinDelta) / rWs;// rotation angle of front steered wheel
        ff[9] = deltaDot;// steer angle
        ff[10] = deltaDDot;// steer rate
        simBegun = true;
    }

    //------------------------------------------------------------------------
    // SetInitialSpeed: Sets the initial speed of the vehicle. Must be set
    //          before simulation has begun.
    //------------------------------------------------------------------------
    public void SetInitalSpeed(double val)
    {
        if(simBegun) return;

        x[1] = val;
    }

    //------------------------------------------------------------------------
    // SetInertia: sets the two inertia properties of the vehicle. 
    //     mm: total mass in kilograms
    //     rgyr: radius of gyration in meters
    //------------------------------------------------------------------------
    public void SetInertia(double mm, double rgyr)
    {
        if(mm <= 0.1)   // check lower bound for mass
            return;     // return and not update parameters.

        if(rgyr < 0.03) // check lower bound for radius of gyration
            return;     // return and not update parameters.

        m = mm;
        Ig = m*rgyr*rgyr;
    }

    //------------------------------------------------------------------------
    // SetGeometry: Sets the geometry of the vehicle.
    //    wsb: distance between rear axle and steer axis
    //    dcg: distance from wheel axle to center of mass
    //    dcst: length of the caster
    //    wid: distance between rear wheels
    //    wRad: radius of rear wheel
    //    wRadS: radius of steered wheel
    //------------------------------------------------------------------------
    public void SetGeometry(double wsb, double dcg, double dcst, double wid, 
        double wRad, double wRadS)
    {
        // check lower bounds
        if(wsb < 0.01) return;
        if(dcg <= 0.0) return;
        if(dcst < 0.0) return;
        if(wid < 0.05) return;
        if(wRad < 0.05) return;
        if(wRadS < 0.05) return;

        if(wsb-dcst < dcg) return; //cg must be btw rear axle and steer contact

        b = dcg;
        c = 0.5*wid;
        d = dcst;
        h = wsb-dcg;

        rW = wRad;
        rWs = wRadS;
    }

    //------------------------------------------------------------------------
    // Getters/Setters
    //------------------------------------------------------------------------

    public double SteerAngleSignal
    {
        set{
            deltaDes = value;
        }
    }

    public double SteerAngle
    {
        get{
            return x[9];
        }
    }

    public double xG
    {
        get{
            return x[0];
        }
    }

    public double zG
    {
        get{
            return x[2];
        }
    }

    public double Heading
    {
        get{
            return x[4];
        }
    }

    public double WheelAngleL
    {
        get{
            return x[6];
        }
    }

    public double WheelAngleR
    {
        get{
            return x[7];
        }
    }

    public double WheelAngleF
    {
        get{
            return x[8];
        }
    }

    public double Speed
    {
        get{
            // ######## You have to write this part ################
            double Speed = Math.Sqrt(Math.Pow(x[1],2) + Math.Pow(x[3],2));
            return(Speed);//velocity x + velocity z
        }
    }

    public double KineticEnergy
    {
        get{
            // ######## You have to write this part ################
            double Speed = Math.Sqrt(Math.Pow(x[1],2) + Math.Pow(x[3],2));
            double kE = Math.Pow(Speed,2) * m * 0.5;
            return(kE);
        }
    }

    public double SlipRateRear
    {
        get{
            double xDot = x[1];
            double zDot = x[3];
            double psi  = x[4];
            double psiDot = x[5];
            double delta = x[9];
            double deltaDot = x[10];

            double cosDelta = Math.Cos(delta);
            double cosPsiPlusDelta = Math.Cos(psi + delta);
            double sinPsiPlusDelta = Math.Sin(psi + delta);
            
            double slipRate = xDot*sinPsiPlusDelta + zDot*cosPsiPlusDelta - h*psiDot*cosDelta + (psiDot+deltaDot)*d;

            return(slipRate);
        }
    }

    public double SlipRateFront
    {
        get{
            double xDot = x[1];
            double zDot = x[3];
            double psi  = x[4];
            double psiDot = x[5];
            double delta = x[9];

            double cosPsi = Math.Cos(psi);
            double sinPsi = Math.Sin(psi);

            double slipRate = (xDot*cosPsi - zDot*sinPsi + c*psiDot) +
            (xDot*sinPsi + zDot*cosPsi + b*psiDot);

            return(slipRate);
        }
    }

    public double FontFrictionFactor
    {
        get{
            // ######## You have to write this part ################

            return(-1.21212121);
        }
    }
}