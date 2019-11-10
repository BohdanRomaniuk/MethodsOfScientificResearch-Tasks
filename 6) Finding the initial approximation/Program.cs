using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Finding_the_initial_approximation
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

        static void Main(string[] args)
        {
        }
    }
}
