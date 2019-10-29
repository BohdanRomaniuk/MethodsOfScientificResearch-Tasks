using System;
using System.Linq;
using Accord.Math;

namespace FikerMethod
{
    class Program
    {
        private static double left;
        private static double right;
        private static int n;
        private static double h;

        public static double F(double x)
        {
            return 1 / ((9 + x * x) * (9 + x * x));
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
            left = 0;
            right = 2;
            n = 5;
            h = (right - left) / (n - 1);
            var A = GetMatrixA();
            var AInverse = A.Inverse();

            var decA = new Accord.Math.Decompositions.EigenvalueDecomposition(A, false, true);
            var eigValsA = decA.RealEigenvalues.OrderBy(c => c);
            Console.WriteLine($"Real eigen values: {string.Join("\t", eigValsA)}");

            var decAInverse = new Accord.Math.Decompositions.EigenvalueDecomposition(AInverse, false, true);
            var eigValsAInverse = decAInverse.RealEigenvalues;
            for (int i = 0; i < eigValsAInverse.Length; ++i)
            {
                eigValsAInverse[i] = 1 / eigValsAInverse[i];
            }
            Console.WriteLine($"Calculated values: {string.Join("\t", eigValsAInverse.OrderBy(c => c))}");

            Console.ReadKey();
        }
    }
}
