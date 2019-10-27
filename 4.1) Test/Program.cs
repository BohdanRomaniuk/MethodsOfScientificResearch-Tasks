using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _4._1__Test
{
    class Program
    {
        private static double left = 1;
        private static double right = 2;
        private static double deviation = 0.001;
        private static int n_max = 3000;

        public static double f(double x)
        {
            return -1 * (9 + x * x) * (9 + x * x);
        }

        public static double[] get_discrete_function(int n, double h)
        {
            var res = new double[n];
            for (int i = 0; i < n; ++i)
            {
                res[i] = f(left + i * h);
            }
            return res;
        }

        private static double[,] get_T(int n, double h)
        {
            var matrix = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                matrix[i, i] = -2 / (h * h);
                if (i + 1 < n)
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
            var tempM = get_T(n, h);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    tempM[i, j] = tempM[i, j] * f[j];
                }
            }

            return tempM.Transpose();
        }

        static void Main(string[] args)
        {
            int n = 10;
            double h = (right - left) / (n - 1);

            var lambda_previous = double.MaxValue;
            while (n <= n_max)
            {
                Console.WriteLine($"n: {n}");
                var A = get_A(n, h);
                var inverse = A.Inverse();

                var dec = new Accord.Math.Decompositions.EigenvalueDecomposition(inverse);
                var eigvals = dec.RealEigenvalues;
                var m = eigvals.Max();
                Console.WriteLine($"m: {m}");
                var lambda = 1 / m;
                Console.WriteLine($"lambda: {lambda}");
                if (Math.Abs(lambda - lambda_previous) < deviation)
                {
                    break;
                }
                lambda_previous = lambda;
                n += 100;
                Console.WriteLine();
            }
            Console.ReadKey();
        }
    }
}
