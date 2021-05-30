package utils;

public class Methods {

    // Vector one norm
    public static double normaNje( double[] x )
    {
        double sum = 0.0;
        for( double xi : x ) sum += Math.abs( xi );
        return sum;
    }

    // Vector second norm
    public static double normaDy( double[] x )
    {
        double sum = 0.0;
        for( double xi : x ) sum += Math.pow( xi, 2 );

        return Math.sqrt( sum );
    }

    // Vector infinite norm
    public static double normaInf( double[] x )
    {
        double max = 0.0;
        for( double xi : x ) if( Math.abs( xi ) > max ) max = Math.abs( xi );

        return max;
    }

    // Vector times Scalar
    public static double[] multVS(double[] a, double b) {
        double[] rez = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            rez[i] = a[i] / b;
        }
        return rez;
    }

    // Vector - Vector
    public static double[] differenceV(double[] a, double[] b) {
        double[] rez = new double[a.length];

        if (a.length == b.length) {
            for (int i = 0; i < a.length; i++) {
                rez[i] = a[i] - b[i];
            }
        }
        else {
            System.out.println("Vectors are not equal!");
        }
        return rez;
    }

    // Matrix - Matrix
    public static double[][] differenceM(double[][] a, double[][] b) {
        double[][] rez = new double[a.length][a[0].length];
        if (a.length == b[0].length) {
            for (int i = 0; i < a.length; i++) {
                for (int j = 0; j < a[0].length; j++) {
                    rez[i][j] = a[i][j] - b[i][j];
                }
            }
        }
        return rez;
    }

    // Matrix and Scalar multiplication
    public static double[][] multSM(double[][] a, double b) {
        double[][] rez = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                rez[i][j] = b * rez[i][j];
            }
        }
        return rez;
    }

    // Matrix and Vector multiplication
    public static double[] multMV( double[][] a, double[] b )
    {
        double[] rez = new double[b.length];

        for( int i = 0; i < a.length; i++ )
            for( int j = 0; j < b.length; j++ )
                rez[i] += a[i][j] * b[j];

        return rez;
    }

    // Matrix and Matrix multiplication
    public static double[][] multMM( double[][] A, double[][] B )
    {
        double[][] C = new double[A.length][A[0].length];

        for( int i = 0; i < A.length; i++ )
            for( int j = 0; j < A[0].length; j++ )
            {
                C[i][j] = 0;
                for( int k = 0; k < A[0].length; k++ )
                    C[i][j] += A[i][k] * B[k][j];
            }

        return C;
    }

    // Matrix one norm
    public static double normaMNje( double[][] A )
    {
        double max = 0.0;
        double[] colSum = new double[A.length];

        for( int i = 0; i < A.length; i++ )
            for( int j = 0; j < A[i].length; j++ )
                colSum[i] += Math.abs( A[j][i] );

        for( double k : colSum ) if( k > max ) max = k;

        return max;
    }

    // Matrix infinite norm
    public static double normaMInf( double[][] A )
    {
        double max = 0.0, sum = 0.0;

        for( double[] ai : A )
        {
            for( double aj : ai )
                sum += Math.abs( aj );

            if( sum > max ) max = sum;
            sum = 0.0;
        }

        return max;
    }

    // Solving system of equations (A - qI)y = x
    public static double[] solveSystem(double[][] A, double[] x, double q) {
        // Creating a unit matrix
        double[][] unitMatrix = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                if (i == j) {
                    unitMatrix[i][j] = 1;
                }
                else {
                    unitMatrix[i][j] = 0;
                }
            }
        }
        double[] y = new double[x.length];
        double[][] temp1 = differenceM(A, multSM(unitMatrix, q));
        double[][] temp = InverseMatrix.inverse(temp1);
        y = multMV(temp, x);

        return y;
    }








    //Gauss Elimination method for solving system of equations
    // public static double[] gaussianElimination(  int n, double[][] A, double[] b ) throws Exception
    // {
    //     double[][] AM = new double[n][n+1];    // the augmented matrix
    //     double[] x = new double[n];

    //     // Initializing the augmented matrix [A:b] from matrix A and b
    //     for( int i = 0; i < n; i++ )
    //         for( int j = 0; j < n + 1; j++ )
    //             AM[i][j] = ( j == n ) ? b[i] : A[i][j];

    //     for( int i = 0; i < n - 1; i++ )
    //     {
    //         // find p the smallest integer with i<=p<=n && A[p][i] != 0
    //         // if no integer p can be found, no unique solution exists
    //         int p = i;
    //         while( p < n && AM[p][i] == 0 ) p++;
    //         if( AM[p][i] == 0 ) throw new Exception( "The system has no unique solution" );

    //         // if p!=i then perform (Ep)<->(Ei)
    //         if( p != i )
    //         {
    //             for( int j = 0; j < n + 1; j++ )
    //             {
    //                 double temp = AM[i][j];
    //                 AM[i][j] = AM[p][j];
    //                 AM[p][j] = temp;
    //             }
    //         }

    //         for( int j = i + 1; j < n; j++ )
    //         {
    //             double mji = AM[j][i] / AM[i][i];
    //             for( int k = 0; k < n + 1; k++ )
    //                 AM[j][k] -= mji * AM[i][k];
    //         }
    //     }

    //     // if AM[n][n] == 0 no unique solution exists
    //     if( AM[n-1][n-1] == 0 ) throw new Exception( "The system has no unique solution" );

    //     x[n-1] = AM[n-1][n] / AM[n-1][n-1];
    //     for( int i = n - 2; i >= 0; i-- )
    //     {
    //         // Signum(j=1...n) AM[i][j]*x[j]
    //         double sum = 0.0;
    //         for( int j = i + 1; j < n; j++ )
    //             sum += AM[i][j] * x[j];

    //         x[i] = ( AM[i][n] - sum ) / AM[i][i];
    //     }

    //     return x;
    // }
}
