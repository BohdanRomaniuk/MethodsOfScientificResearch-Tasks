using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Kolats_and_KrylovBogoliubov_methods
{
    class Program
    {
        private const double left = 0;
        private const double right = 1;
        private const int size = 100; // size > 5
        private const double h = (right - left) / (size - 1);
        private const double l2 = 540;

        private static Matrix CreateMatrix()
        {
            var matrix = new DenseMatrix(size, size);
            for (int i = 0; i < size; ++i)
            {
                matrix[i, i] = -2 / (h * h);
                if (i + 1 < size)
                {
                    matrix[i, i + 1] = 1 / (h * h);
                    matrix[i + 1, i] = 1 / (h * h);
                }
            }
            return matrix;
        }

        private static Vector GetInitialF()
        {
            var vector = new DenseVector(size);
            for (int i = 0; i < size; ++i)
            {
                vector[i] = 1;
                var val = left + (i * h);
                vector[i] = -((9 + val * val) * (9 + val * val));
            }
            return vector;
        }

        private static double Integrate(Vector u, Vector v)
        {
            double sum = 0;
            for (int i = 0; i < size; ++i)
            {
                sum += u[i] * v[i] * h;
            }
            return sum;
        }

        static void Main(string[] args)
        {
            Matrix T = CreateMatrix();
            List<Vector> f = new List<Vector>();
            f.Add(GetInitialF());

            //Console.WriteLine("T matrix:\t\t\t\t\t\t\t F vector:");
            //for (int i = 0; i < size; ++i)
            //{
            //    for (int j = 0; j < size; ++j)
            //    {
            //        Console.Write($"{T[i, j],-5}{(j == size - 1 ? "\t\t\t" + f[0][i].ToString() : string.Empty)}");
            //    }
            //    Console.WriteLine();
            //}
            //Console.WriteLine();

            //Iterations count
            int itersCount = 15;

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
