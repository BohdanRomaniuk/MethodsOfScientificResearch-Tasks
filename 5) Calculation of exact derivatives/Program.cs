using Accord.Math;
using System;
using System.Linq;

namespace Calculation_of_exact_derivatives
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

        private static void ddet(double[,] mat, int n, out double f, out double df, out double ddf)
        {
            f = 0;
            df = 0;
            ddf = 0;
        }

        private static double[,] sub_diag(double[,] mat, double value, int n)
        {
            var sub = new double[n, n];
            sub = sub.FillDiagonal(value);
            return mat.Sub(sub);
        }

        private static void n_eig(double l, int n, double[,] AA)
        {
            double l_n = 0;

            for (int i = 0; i < n; ++i)
            {
                var _mat = sub_diag(AA, l, n);
                var _mat_1 = sub_diag(AA.Inverse(), l, n);
                ddet(_mat, n, out var f, out var df, out var ddf);

                l_n = l + f / df;
                l = l_n;
            }


            Console.WriteLine($"First = {l_n}");
            var decA = new Accord.Math.Decompositions.EigenvalueDecomposition(AA, false, true);
            var eigValsA = decA.RealEigenvalues.OrderBy(c => c);
            Console.WriteLine($"Second = { string.Join("\t", eigValsA.Select(c => c))}");
        }

        static void Main(string[] args)
        {
            left = 0;
            right = 1;
            n = 4;
            h = (right - left) / (n - 1);
            var AA = GetMatrixA();
            n_eig(5, 30, AA);

            Console.ReadKey();
        }
    }
}
