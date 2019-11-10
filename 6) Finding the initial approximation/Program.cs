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
        public const double LowerBound = 0.0 ;
        public const double UpperBound = 1.0 ;
        public const int N = 50;
        public const double H = 1.0 / (N + 1);
        public const double Center = 65.0;
        public const double Rho = 2.0;

        private static double F(double x)
        {
            return (-1) * ((9 + x * x) * (9 + x * x));
        }

        private static Vector<Complex> Vector()
        {
            Vector<Complex> vector = Vector<Complex>.Build.Dense(N);
            for (int i = 0; i < N; ++i)
            {
                vector[i] = F(LowerBound + i * H);
            }
            return vector;
        }

        private static Matrix<Complex> GetMatrixA()
        {
            return Matrix<Complex>.Build.DiagonalOfDiagonalVector(Vector());
        }

        private static Matrix<Complex> GetMatrixT()
        {
            var matrix = Matrix<Complex>.Build.Sparse(N, N);
            for (int i = 0; i < N; ++i)
            {
                matrix[i, i] = -2;
                if (i + 1 < N)
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
            return Math.Pow(H, 2);
        }

        public static void LU(Matrix<Complex> D, out Matrix<Complex> L, out Matrix<Complex> U)
        {
            L = Matrix<Complex>.Build.Sparse(N, N);
            U = Matrix<Complex>.Build.Sparse(N, N);

            for (int i = 0; i < N; ++i)
            {
                U[0, i] = D[0, i];

                for (int j = i; j < N; ++j)
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
            M = Matrix<Complex>.Build.Sparse(N, N);
            V = Matrix<Complex>.Build.Sparse(N, N);

            for (int i = 0; i < N; ++i)
            {
                V[0, i] = B[0, i];

                for (int j = i; j < N; ++j)
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

            for (int i = 0; i < N; ++i)
            {
                var spectralRadius = SpectralRadius(i + 1, N);
                var lambda = Center + spectralRadius;

                Matrix<Complex> D = LinearOperator(lambda);
                Matrix<Complex> B = LinearOperatorDerivative();
                Matrix<Complex> L, U, M, V;

                LU(D, out L, out U);
                MULV(B, out M, U, L, out V);

                var determinant = Determinant(U);
                var determinantDerivative = DeterminantDerivative(V, U);

                summation += spectralRadius * determinantDerivative / determinant;
            }

            return (summation / N).Magnitude;
        }

        public static double InitialApproximation(int j, int n)
        {
            return (Center + SpectralRadius(j, n)).Magnitude;
        }

        private static Complex SpectralRadius(int j, int n)
        {
            return Rho * Complex.Exp(Complex.ImaginaryOne * 2.0 * Math.PI * j / n);
        }

        static void Main(string[] args)
        {
            double s0 = RootsCount();
            Console.WriteLine("s0 = {0}\n", s0);

            int rootsCount = (int)Math.Round(s0);
            if (rootsCount > 0)
            {
                Console.WriteLine("{0} root(s) found. Initial Approximation(s):\n", rootsCount);
                for (int i = 0; i < rootsCount; ++i)
                {
                    double initalApproximation = InitialApproximation(i + 1, rootsCount);

                    Console.WriteLine("{0}. lambda0 = {1}", i + 1, initalApproximation / 2);
                }
            }
            else
            {
                Console.WriteLine("No roots found.");
            }

            Console.ReadKey();
        }
    }
}
