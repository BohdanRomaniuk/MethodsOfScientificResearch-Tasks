using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace Finding_the_initial_approximation
{
    using Helpers;
    using MathNet.Numerics.LinearAlgebra;

    class Program
    {
        public static double left;
        public static double right;
        public static int n;
        public static double h;
        public static double center;
        public static double rho;

        private static double F(double x)
        {
            return (-1) * ((9 + x * x) * (9 + x * x));
        }

        private static Vector<Complex> Vector()
        {
            Vector<Complex> vector = Vector<Complex>.Build.Dense(n);
            for (int i = 0; i < n; ++i)
            {
                vector[i] = F(left + i * h);
            }
            return vector;
        }

        private static Matrix<Complex> GetMatrixA()
        {
            return Matrix<Complex>.Build.DiagonalOfDiagonalVector(Vector());
        }

        private static Matrix<Complex> GetMatrixT()
        {
            var matrix = Matrix<Complex>.Build.Sparse(n, n);
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

        public static Matrix<Complex> LinearOperator(Complex lambda)
        {
            return GetMatrixT() - lambda * H2() * GetMatrixA();
        }

        public static Matrix<Complex> LinearOperatorDerivative()
        {
            return -H2() * GetMatrixA();
        }

        private static double H2()
        {
            return Math.Pow(h, 2);
        }

        public static void LU(Matrix<Complex> D, out Matrix<Complex> L, out Matrix<Complex> U)
        {
            L = Matrix<Complex>.Build.Sparse(n, n);
            U = Matrix<Complex>.Build.Sparse(n, n);

            for (int i = 0; i < n; ++i)
            {
                U[0, i] = D[0, i];

                for (int j = i; j < n; ++j)
                {
                    Complex summation = Complex.Zero;

                    for (int k = 0; k < i; ++k)
                    {
                        summation += L[i, k] * U[k, j];
                    }

                    U[i, j] = D[i, j] - summation;

                    if (i <= j)
                    {
                        summation = Complex.Zero;

                        for (int k = 0; k < i; ++k)
                        {
                            summation += L[j, k] * U[k, i];
                        }

                        L[j, i] = (D[j, i] - summation) / U[i, i];
                    }
                }
            }
        }

        public static void MULV(Matrix<Complex> B, out Matrix<Complex> M, Matrix<Complex> U, Matrix<Complex> L, out Matrix<Complex> V)
        {
            M = Matrix<Complex>.Build.Sparse(n, n);
            V = Matrix<Complex>.Build.Sparse(n, n);

            for (int i = 0; i < n; ++i)
            {
                V[0, i] = B[0, i];

                for (int j = i; j < n; ++j)
                {
                    Complex summation = Complex.Zero;

                    for (int k = 0; k < i; ++k)
                    {
                        summation += M[i, k] * U[k, j] + L[i, k] * V[k, j];
                    }

                    V[i, j] = B[i, j] - summation;

                    if (i <= j)
                    {
                        summation = Complex.Zero;

                        for (int k = 0; k < i; ++k)
                        {
                            summation += M[j, k] * U[k, i] + L[j, k] * V[k, i];
                        }

                        M[j, i] = (B[j, i] - summation - L[j, i] * V[i, i]) / U[i, i];
                    }
                }
            }
        }

        public static Complex Determinant(Matrix<Complex> U)
        {
            Complex result = Complex.One;
            for (int i = 0; i < U.RowCount; ++i)
            {
                result *= U[i, i];
            }
            return result;
        }

        public static Complex DeterminantDerivative(Matrix<Complex> V, Matrix<Complex> U)
        {
            var sum = Complex.Zero;
            for (int k = 0; k < V.RowCount; ++k)
            {
                Complex product = Complex.One;
                for (int i = 0; i < U.RowCount; ++i)
                {
                    if (i != k)
                    {
                        product *= U[i, i];
                    }
                }
                sum += V[k, k] * product;
            }

            return sum;
        }

        public static double RootsCount()
        {
            Complex summation = Complex.Zero;

            for (int i = 0; i < n; ++i)
            {
                var spectralRadius = SpectralRadius(i + 1, n);
                var lambda = center + spectralRadius;

                Matrix<Complex> D = LinearOperator(lambda);
                Matrix<Complex> B = LinearOperatorDerivative();
                Matrix<Complex> L, U, M, V;

                LU(D, out L, out U);
                MULV(B, out M, U, L, out V);

                var determinant = Determinant(U);
                var determinantDerivative = DeterminantDerivative(V, U);

                summation += spectralRadius * determinantDerivative / determinant;
            }

            return (summation / n).Magnitude;
        }

        public static double InitialApproximation(int j, int n)
        {
            return (center + SpectralRadius(j, n)).Magnitude;
        }

        private static Complex SpectralRadius(int j, int n)
        {
            return rho * Complex.Exp(Complex.ImaginaryOne * 2.0 * Math.PI * j / n);
        }

        static void Main(string[] args)
        {
            left = 0.0;
            right = 1.0;
            n = 50;
            h = (right - left) / (n + 1);
            center = 65;
            rho = 2;

            int parts = 4;
            for (int k = 0; k < parts; ++k)
            {
                left += k;
                right += k;
                double s0 = RootsCount();
                int rootsCount = (int)Math.Round(s0);
                Console.WriteLine($"[{left}, {right}] -> s0 = {s0}\nFound {rootsCount} root(s):");

                if (rootsCount > 0)
                {
                    for (int i = 0; i < rootsCount; ++i)
                    {
                        double initalApproximation = InitialApproximation(i + 1, rootsCount);
                        Console.WriteLine("{0}. lambda = {1}", i + 1, initalApproximation / 2 - 0.5);
                    }
                }
                else
                {
                    Console.WriteLine("-------------");
                }
                Console.WriteLine();
            }


            Console.ReadKey();
        }
    }
}
