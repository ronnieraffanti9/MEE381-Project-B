//============================================================================
// LinAlgEq.cs  Defines a class for a linear algebraic equations solver
//============================================================================
using System;

public class LinAlgEq
{
    private int n = 3;       // number of algebraic equations (number of unknowns too)
    private double[][] _A;   // coefficient matrix
    private double[] _b;     // right hand side
    private double[][] M;    // augmented matrix
    private double[] _x;     // solution

    //--------------------------------------------------------------------
    // Constructor for the class.
    //--------------------------------------------------------------------
    public LinAlgEq(int nn = 3)
    {
        _b = new double[1];   // these three lines get rid of the warning
        _A = new double[1][];
        M  = new double[1][];
        _x = new double[1];

        Resize(nn);
    }

    //--------------------------------------------------------------------
    // resize: resize the matrices to hold the right number of equations
    //         and unknowns.
    //--------------------------------------------------------------------
    public void Resize(int nn)
    {
        // ##### should check if nn is bigger than zero
        n = nn;

        _b = new double[n];
        _A = new double[n][];
        M  = new double[n][];
        _x = new double[n];

        int i,j;
        for(i=0;i<n;++i)
        {
            _A[i] = new double[n];
            M[i] = new double[n+1];
                
            for(j=0;j<n;++j)
            {
                _A[i][j] = 0.0;
            }
            _A[i][i] = 1.0;
            _b[i] = 0.0;
            _x[i] = 0.0;
        }
    }

    //--------------------------------------------------------------------
    // SolveGauss: Solve by Gauss Eliminition
    //--------------------------------------------------------------------
    public void SolveGauss()
    {
        // form augmented matrix
        int i, j, k;
        for(i=0;i<n;++i)
        {
            for(j=0;j<n;++j)
            {
                M[i][j] = _A[i][j];
            }
            M[i][n] = _b[i];
        }

        // perform Gauss elimination
        for (i = 0; i < (n - 1); ++i) {// loop through rows
            PivotRow(i); // move largest number into pivot element 
            for (j = i + 1; j < n; ++j) {// loop through rows below the pivot row
                double multiple = M[j][i] / M[i][i];// Calculate the multaple for row elimination
                for (k = i; k <= n; ++k) {//loop through values in row
                    M[j][k] = M[j][k] - multiple * M[i][k];//calculate new value
                }
            }
        }

        // perform back substitution
        for(i = n - 1; i >= 0; --i){// loop through the matrix starting from the bottom
            double multaple = M[i][n];// rightmost value in row (b matrix)
            for(j = n - 1; j > i; --j){//loops through values right of pivot
                multaple = multaple - M[i][j] * _x[j];//subtract known variables
            }
            _x[i] = multaple/M[i][i];//calculate the solution for the current row
        }
    }

    //--------------------------------------------------------------------
    // PivotRow
    //--------------------------------------------------------------------
    private void PivotRow(int j)
    {
        double[] holder;
        double maxElem = Math.Abs(M[j][j]);
        int rowIdx = j;
        int i;

        for(i = j+1; i<n; ++i)
        {
            // find largest element in jth column
            if(Math.Abs(M[i][j])>maxElem)
            {
                maxElem = Math.Abs(M[i][j]);
                rowIdx = i;
            }
        }

        // swap rows
        if(rowIdx != j)
        {
            holder = M[j];
            M[j] = M[rowIdx];
            M[rowIdx] = holder;
            //Console.WriteLine("Swap " + j.ToString());
        }
    }

    //--------------------------------------------------------------------
    // Check: checks the solution
    //--------------------------------------------------------------------
    public double Check()
    {
        double sum = 0.0;
        double sum2 = 0.0;

        int i,j;
        for(i=0;i<n;++i)
        {
            sum = 0.0;
            for(j=0;j<n;++j)
            {
                sum += _A[i][j] * _x[j];
            }
            double delta = sum - _b[i];
            sum2 += delta*delta; 
        }

        return(Math.Sqrt(sum2/(1.0*n)));
    }

    //--------------------------------------------------------------------
    // getters and setters 
    //--------------------------------------------------------------------
    public double[] b 
    {
        get
        {
            return _b;
        }            
        set
        {
            _b = value;
        }
    }

    public double[][] A
    {
        get
        {
            return _A;
        }
        set
        {
            _A = value;
        }
    }

    public double[] sol
    {
        get
        {
            return _x;
        }
    }
}