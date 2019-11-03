using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;

namespace WeinsteinMethod
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

        private static double[,] GetMatrixB()
        {
            var matrix = GetMatrixT();
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    matrix[i, j] *= -1;
                }
            }
            return matrix;
        }

        private static double[,] GetMatrixC(double[,] A, double[,] B)
        {
            var c = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    c[i, j] = A[i, j] - B[i, j];
                }
            }
            return c;
        }

        private static double[] GetFk(int k)
        {
            var fk = new double[n];
            if (k == 0)
            {
                for (int i = 0; i < n; ++i)
                {
                    fk[i] = 1;
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    fk[i] = Math.Pow(left + i * h, k);
                }
            }
            return fk;
        }

        private static double[] GetGk(double[,] C, double[] fk)
        {
            return Matrix.Dot(C, fk);
        }

        private static double[,] GetCm(double[,] C, int m)
        {
            var fk_arr = new List<double[]>();
            for (int i = 0; i < m; ++i)
            {
                fk_arr.Add(GetFk(i));
            }

            var gk_arr = new List<double[]>();
            for (int i = 0; i < m; ++i)
            {
                gk_arr.Add(GetGk(C, fk_arr[i]));
            }

            var bij = new double[m, m];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    bij[i, j] += ScalarSum(fk_arr[j], gk_arr[i]);
                }
            }

            //bij = bij.Inverse();

            var Cm = new double[n, n];
            for (int i = 0; i < m; ++i)
            {
                // Copy
                var gi = new double[n];
                for (int k = 0; k < n; ++k)
                {
                    gi[k] = gk_arr[i][k];
                }

                for (int j = 0; j < m; ++j)
                {
                    // gi = gi * bij[i][j]
                    for (int z = 0; z < n; ++z)
                    {
                        gi[z] = gi[z] * bij[i, j];
                    }

                    Cm = Cm.Sum(MultiVecs(gi, gk_arr[j]));
                }
            }
            return Cm;
        }

        private static double ScalarSum(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; ++i)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }

        private static double[,] MultiVecs(double[] a, double[] b)
        {
            double[,] matr = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    matr[i, j] = a[i] * b[j];
                }
            }
            return matr;
        }

        private static double[] GetRealEigenValues(double[,] A)
        {
            var decA = new Accord.Math.Decompositions.EigenvalueDecomposition(A, false, true);
            return decA.RealEigenvalues.OrderBy(c => c).ToList().ToArray();
        }

        private static double[] GetCalculatedEigenValues(double[,] B, double[,] C, int m)
        {
            var Cm = GetCm(C, m);
            var Am = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    Am[i, j] = B[i, j] + Cm[i, j];
                }
            }
            var dec = new Accord.Math.Decompositions.EigenvalueDecomposition(Am, false, true);
            return dec.RealEigenvalues.OrderBy(c => c).ToList().ToArray();
        }

        private static double[] GetCalculatedEigenValues(double[,] A, int m)
        {
            var decAInverse = new Accord.Math.Decompositions.EigenvalueDecomposition(A.Inverse(), false, true);
            var eigValsAInverse = decAInverse.RealEigenvalues;
            for (int i = 0; i < eigValsAInverse.Length; ++i)
            {
                eigValsAInverse[i] = 1 / eigValsAInverse[i];
            }
            return eigValsAInverse.OrderBy(c => c).ToList().ToArray();
        }

        static void Main(string[] args)
        {
            left = 0;
            right = 2;
            n = 20;
            h = (right - left) / (n - 1);

            var m = 10;
            var A = GetMatrixA();
            var B = GetMatrixB();
            var C = GetMatrixC(A, B);

            Console.WriteLine("Real eigen values:");
            var lambdaReal = GetRealEigenValues(A);
            lambdaReal.Print();

            Console.WriteLine("\nCalculated eigen values:");
            var lambdaCalc = GetCalculatedEigenValues(A, m);
            lambdaCalc.Print();
            Console.ReadKey();
        }
    }
}
