using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Kolats_and_KrylovBogoliubov_methods
{
    class Program
    {
        private static double left;
        private static double right;
        private static int n;
        private static double h;
        private static double l2;

        public static double F(double x)
        {
            return (-1) * ((9 + x * x) * (9 + x * x));
        }

        private static Vector GetInitialF()
        {
            var vector = new DenseVector(n);
            for (int i = 0; i < n; ++i)
            {
                vector[i] = F(left + (i * h));
            }
            return vector;
        }

        private static MathNet.Numerics.LinearAlgebra.Matrix<double> CreateMatrix(Vector f)
        {
            var rightMatrix = new DenseMatrix(n, n);
            var matrix = new DenseMatrix(n, n);
            for (int i = 0; i < n; ++i)
            {
                rightMatrix[i, i] = f[i];
                matrix[i, i] = -2;
                if (i + 1 < n)
                {
                    matrix[i, i + 1] = 1;
                    matrix[i + 1, i] = 1;
                }
            }
            return matrix.Multiply(rightMatrix);
        }

        private static double Integrate(Vector u, Vector v)
        {
            double sum = 0;
            for (int i = 0; i < n; ++i)
            {
                sum += u[i] * v[i] * h;
            }
            return sum;
        }

        static void Main(string[] args)
        {
            //Initialization
            left = 0;
            right = 1;
            n = 4;
            h = (right - left) / (n - 1);
            l2 = 34;

            //Iterations count
            int itersCount = 15;

            List<Vector> f = new List<Vector>();
            f.Add(GetInitialF());
            var T = CreateMatrix(f[0]);

            Console.WriteLine("T matrix:\t\t\t\t\t\t\t F vector:");
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    Console.Write($"{T[i, j],-5}{(j == n - 1 ? "\t\t\t" + f[0][i].ToString() : string.Empty)}");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            // Kolats method begin
            ++itersCount;
            for (int i = 0; i < itersCount - 1; ++i)
            {
                var fCurrent = f.Last(); //fi-1
                var fNext = T.Solve(fCurrent) as Vector; //fi
                f.Add(fNext);
            }

            // a0, a1, a2, ... (Shwartz constants)
            double[] a = new double[itersCount];
            for (int i = 0; i < itersCount; ++i)
            {
                a[i] = Integrate(f[0], f[i]);
            }

            // μ calculation
            double[] mu = new double[itersCount];
            for (int i = 1; i < itersCount; ++i)
            {
                mu[i] = a[i - 1] / a[i];
            }

            // ν calculation
            double[] nu = new double[itersCount];
            for (int i = 2; i < itersCount; ++i)
            {
                nu[i] = mu[i] - (mu[i] - mu[i - 1]) / (l2 / mu[i] - 1);
            }

            // λ calculation
            double[] lambdaList = new double[itersCount];
            for (int i = 2; i < itersCount; ++i)
            {
                lambdaList[i] = (mu[i] + nu[i]) / 2;
            }

            // Kolats method results
            Console.OutputEncoding = Encoding.UTF8;
            Console.WriteLine("Kolats method:");
            Console.WriteLine("{0,3} {1,10} {2,15} {3,15} ", "Iter.", "ν", "λ", "μ");
            for (int i = 2; i < itersCount; ++i)
            {
                Console.WriteLine("{0}:\t{1:F10}\t{2:F10}\t{3:F10}", i, nu[i], lambdaList[i], mu[i]);
            }
            Console.WriteLine();

            // Krylov-Bogoliubov method begin
            // μ calculation
            double[] lowerBound = new double[itersCount];
            double[] upperBound = new double[itersCount];
            for (int i = 2; i < itersCount; ++i)
            {
                lowerBound[i] = mu[i] - Math.Sqrt(Math.Abs((mu[i - 1] - mu[i]) * mu[i]));
                upperBound[i] = mu[i] + Math.Sqrt(Math.Abs((mu[i - 1] - mu[i]) * mu[i]));
                lambdaList[i] = (upperBound[i] + lowerBound[i]) / 2;
            }

            // Krylov-Bogoliubov method results
            Console.WriteLine("Krylov-Bogoliubov method:");
            Console.WriteLine("i:\t\tμ lower\t\tλ\t\tμ upper");
            for (int i = 2; i < itersCount; ++i)
            {
                Console.WriteLine("{0}:\t{1:F10}\t{2:F10}\t{3:F10}", i, lowerBound[i], lambdaList[i], upperBound[i]);
            }

            Console.ReadKey();
        }
    }
}
