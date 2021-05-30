package main;

import utils.Methods;

import java.util.Arrays;

public class InversePowerMethod {
    public static void inversePowerMethod(int n, double[][] A, double[] x, double TOL, int maxIterations) {
        double[] y;

        //Step 1 - Set q = (x^t * Ax) / (x^t * x)
        double q = 0;
        double q1 = Methods.normaInf(x);
        double q2 = Methods.normaInf(Methods.multMV(A, x));
        double q3 = q1 * Methods.normaNje(x);
        q = (q1 * q2) / q3;
        System.out.println("q = " + q);

        //Step 2 - Set k = 1
        int k = 1;

        //Step 3 - Find the smallest integer p with 1 <= p <= n and |Xp| = ||x||inf
        int p = 0;
        for (int i = 1; i < n; i++) {
            if (x[i] == Methods.normaInf(x)) {
                p = i;
            }
        }

        //Step 4 - Set x =  x / Xp
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i] / x[p];
        }

//        int iterations = 0;
//        if (iterations > maxIterations) {
//            System.out.println("Maximim number of Iterations Exceeded!");
//        }
//        else {
            //Step 5 - While (k <= N) do Steps 6-12
            while (k <= maxIterations) {
                //Step 6 - Solve the linear system (A - qI)y = x
                y = Methods.solveSystem(A, x, q);

                //Step 7 - if the system doesnt have a unique solution then output: q is an eigenvalue
                if (y.length == 1) {
                    System.out.println("q = " + q + ", is an eigen value");
                }

                //Step 8 - Set u = y[p]
                double u = y[p];

                //Step 9 - Find the smallest intiger p with 1 <= p <= n and |Yp| = ||y||inf
                int p1 = 0;
                for (int i = 1; i < n; i++) {
                    if (y[i] == Methods.normaInf(y)) {
                        p1 = i;
                    }
                }

                //Step 10 - Set ERR = || x - (y/Yp) ||inf
                double[] ERR1 = Methods.differenceV(x, Methods.multVS(y, y[p]));
                for (int i = 0; i < x.length; i++) {
                    x[i] = y[i] / y[p1];
                }
                double ERR = Methods.normaInf(ERR1);

                //Step 11 - if ERR < TOL then set u = (1/u) + q
                if (ERR < TOL) {
                    u = (1 / u) + q;
                    System.out.println("u = " + u + ", x = " + Arrays.toString(x));
                }

                //Step 12 - Set k = k + 1
                k++;
            }
        }
//    }
    public static void main(String[] args) {
        double[][] A = {
                {-4, 14, 0},
                {-5, 13, 0},
                {-1, 0, 2}
        };
        double[] B = {1, 1, 1};
        double TOL = 10E-6;
        int maxIterations = 100;
        inversePowerMethod(3, A, B, TOL, maxIterations);
    }
}