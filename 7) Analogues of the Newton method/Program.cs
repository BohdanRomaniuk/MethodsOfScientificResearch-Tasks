using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _7__Analogues_of_the_Newton_method
{
    using Helpers;
    class Program
    {
        private static double left;
        private static double right;
        private static int n;
        private static double h;

        public static double F(double x)
        {
            return (-1) * ((9 + x * x) * (9 + x * x));
        }

        public static double[] GetDiscreteFunction()
        {
            var res = new double[n];
            for (int i = 0; i < n; ++i)
            {
                res[i] = F(left + i * h);
            }
            return res;
        }

        private static double[,] GetMatrixT()
        {
            var matrix = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                matrix[i, i] = -2;
                if (i + 1 < n)
                {
                    matrix[i, i + 1] = 1;
                    matrix[i + 1, i] = 1;
                }
            }
            return matrix;
        }

        private static double[,] GetMatrixA()
        {
            var f = GetDiscreteFunction();
            var diagonalMatrix = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                diagonalMatrix[i, i] = f[i];
            }
            return GetMatrixT().Multiply(diagonalMatrix);
        }

        private static double[,] GetMatrixD(double[,] A, double lambda)
        {
            //D=A-λI
            var D = new double[n, n];
            D = D.FillDiagonal(lambda);
            return A.Sub(D);
        }

        public static void LU(double[,] D, out double[,] L, out double[,] U)
        {
            L = new double[n, n];
            U = new double[n, n];

            for (int i = 0; i < n; ++i)
            {
                U[0, i] = D[0, i];
                for (int j = i; j < n; ++j)
                {
                    var sum = 0.0;
                    for (int k = 0; k < i; ++k)
                    {
                        sum += L[i, k] * U[k, j];
                    }
                    U[i, j] = D[i, j] - sum;
                    if (i <= j)
                    {
                        sum = 0.0;
                        for (int k = 0; k < i; ++k)
                        {
                            sum += L[j, k] * U[k, i];
                        }
                        L[j, i] = (D[j, i] - sum) / U[i, i];
                    }
                }
            }
        }

        public static void MV(double[,] B, out double[,] M, double[,] U, double[,] L, out double[,] V)
        {
            M = new double[n, n];
            V = new double[n, n];

            for (int i = 0; i < n; ++i)
            {
                V[0, i] = B[0, i];
                for (int j = i; j < n; ++j)
                {
                    var sum = 0.0;
                    for (int k = 0; k < i; ++k)
                    {
                        sum += M[i, k] * U[k, j] + L[i, k] * V[k, j];
                    }
                    V[i, j] = B[i, j] - sum;
                    if (i <= j)
                    {
                        sum = 0;
                        for (int k = 0; k < i; ++k)
                        {
                            sum += M[j, k] * U[k, i] + L[j, k] * V[k, i];
                        }
                        M[j, i] = (B[j, i] - sum - L[j, i] * V[i, i]) / U[i, i];
                    }
                }
            }
        }

        public static void WN(double[,] L, double[,] U, double[,] M, double[,] V, out double[,] W)
        {
            W = new double[n, n];
            W.FillDiagonal(0);

            var C = new double[n, n];
            var N = new double[n, n];
            C.FillDiagonal(1);
            N.FillDiagonal(1);

            for (int m = 0; m < n; ++m)
            {
                for (int i = m; i < n; ++i)
                {
                    var sum = 0.0;
                    for (int k = 0; k < m; ++k)
                    {
                        sum += N[m, k] * U[k, i] + 2 * M[m, k] * V[k, i] + L[m, k] * W[k, i];
                    }
                    W[m, i] = C[m, i] - sum;
                }

                for (int j = m + 1; j < n; ++j)
                {
                    var sum = 0.0;
                    for (int k = 0; k < m; ++k)
                    {
                        sum += N[j, k] * U[k, m] + 2 * M[j, k] * V[k, m] + L[j, k] * W[k, m];
                    }
                    N[j, m] = (C[j, m] - sum - (2 * M[j, m] * V[m, m]) - L[j, m] * W[m, m]) / U[m, m];
                }
            }
        }

        private static double FirstApproximation(double[,] D, int n)
        {
            var B = new double[n, n];
            var C = new double[n, n];
            B.FillDiagonal(1);

            var U = new double[n, n];
            var L = new double[n, n];
            var V = new double[n, n];
            var M = new double[n, n];
            var W = new double[n, n];

            LU(D, out L, out U);
            MV(B, out M, U, L, out V);
            WN(L, U, M, V, out W);


            ////////////Second approach
            double sum = 0;
            for (int i = 0; i < n; ++i)
            {
                sum += V[i, i] / U[i, i];
            }

            //Delta lamdba	
            return 1 / sum;
        }

        static void Main(string[] args)
        {
            left = 0;
            right = 1;
            n = 4;
            h = (right - left) / (n - 1);
            var deviation = 0.00000001;
            var A = GetMatrixA();

            double lambdaPrev = 32; //Initial approach
            double lambdaNext = 0;
            double deltaLambda = double.MaxValue;

            Console.WriteLine($"Lambda = {lambdaPrev}");
            do
            {
                var D = GetMatrixD(A, lambdaPrev);
                deltaLambda = FirstApproximation(D, n);

                lambdaNext = lambdaPrev + deltaLambda;
                lambdaPrev = lambdaNext;
                Console.WriteLine($"Lambda = {lambdaNext}");
            } while (Math.Abs(deltaLambda) > deviation);

            Console.WriteLine($"\nLambda final = {lambdaNext}");
            var decA = new Accord.Math.Decompositions.EigenvalueDecomposition(A, false, true);
            var eigValsA = decA.RealEigenvalues.OrderBy(c => c);
            Console.WriteLine($"\nReal values   = { string.Join("\t", eigValsA.Select(c => c))}");

            Console.ReadKey();
        }
    }
}
