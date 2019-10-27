using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2__Test
{
    class Program
    {
        static double left = 1;
        static double right = 2;

        public static double F(double x)
        {
            return (-1) * (9 + x * x) * (9 + x * x);
        }

        public static double[] get_discrete_function(int n, double h)
        {
            var res = new double[n - 1];
            for (int i = 1; i < n; ++i)
            {
                res[i - 1] = F(left + i * h);
            }
            return res;
        }

        private static double[,] get_T(int n, double h)
        {
            var matrix = new double[n - 1, n - 1];
            for (int i = 0; i < n - 1; ++i)
            {
                matrix[i, i] = -2 / (h * h);
                var val = left + (i * h);
                if (i + 1 < n - 1)
                {
                    matrix[i, i + 1] = 1 / (h * h);
                    matrix[i + 1, i] = 1 / (h * h);
                }
            }
            return matrix;
        }

        private static double[,] get_A(int n, double h)
        {
            var f = get_discrete_function(n, h);
            var matrixA = get_T(n, h);
            for (int i = 0; i < n - 1; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    matrixA[i, j] = matrixA[i, j] * f[j];
                }
            }
            return matrixA.Transpose();
        }

        private static double[,] get_B(int n, double h)
        {
            var matrix = get_T(n, h);
            for (int i = 0; i < n - 1; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    matrix[i, j] *= -1;
                }
            }
            return matrix;
        }

        private static double[,] get_C(double[,] A, double[,] B, int n)
        {
            var c = new double[n - 1, n - 1];
            for (int i = 0; i < n - 1; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    c[i, j] = A[i, j] - B[i, j];
                }
            }
            return c;
        }

        private static double[] get_fk(int k, int n, double h)
        {
            var fk = new double[n - 1];
            if (k == 0)
            {
                for (int i = 0; i < n - 1; ++i)
                {
                    fk[i] = 1;
                }
            }
            else
            {
                for (int i = 1; i < n; ++i)
                {
                    fk[i - 1] = Math.Pow(left + i * h, k);
                }
            }
            return fk;
        }

        private static double[] get_gk(double[,] C, double[] fk)
        {
            return Matrix.Dot(C, fk);
        }

        private static double[,] get_Cm(double[,] C, int n, double h, int m)
        {
            var fk_arr = new List<double[]>();
            for (int i = 0; i < m; ++i)
            {
                fk_arr.Add(get_fk(i, n, h));
            }

            var gk_arr = new List<double[]>();
            for (int i = 0; i < m; ++i)
            {
                gk_arr.Add(get_gk(C, fk_arr[i]));
            }

            var bij = new double[m, m];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    bij[i, j] += scalarSum(fk_arr[j], gk_arr[i]);
                }
            }

            bij = bij.Inverse();
            var Cm = new double[n - 1, n - 1];
            for (int i = 0; i < m; ++i)
            {
                var gi = new double[n - 1];
                for (int k = 0; k < n - 1; ++k)
                {
                    gi[k] = gk_arr[i][k];
                }

                for (int j = 0; j < m; ++j)
                {
                    for (int z = 0; z < n - 1; ++z)
                    {
                        gi[z] = gi[z] * bij[i, j];
                    }

                    for (int i1 = 0; i1 < n - 1; ++i1)
                    {
                        for (int i2 = 0; i2 < n - 1; ++i2)
                        {
                            Cm[i1, i2] += gi[i1] * gk_arr[j][i2];
                        }
                    }

                }
            }
            return Cm;
        }

        private static double scalarSum(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; ++i)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }

        private static double run(double[,] B, double[,] C, int n, double h, int m)
        {
            var Cm = get_Cm(C, n, h, m);
            var Am = new double[n - 1, n - 1];
            for (int i = 0; i < n - 1; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    Am[i, j] = B[i, j] + Cm[i, j];
                }
            }
            var dec = new Accord.Math.Decompositions.EigenvalueDecomposition(Am, false, true);
            var eigvals = dec.RealEigenvalues;
            return eigvals.Min();
        }


        static void Main(string[] args)
        {
            var n = 5;
            var h = (right - left) / (n - 1);

            var m = 10;
            var A = get_A(n, h);
            var B = get_B(n, h);
            var C = get_C(A, B, n);

            var lambda = run(B, C, n, h, m);
            Console.WriteLine($"lambda = {lambda}");
            Console.ReadKey();
        }
    }
}
