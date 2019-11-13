using System;
using System.Numerics;

namespace Finding_the_initial_approximation
{
    using Helpers;
    class Program
    {
        public static double left;
        public static double right;
        public static int n;
        public static double h;
        public static double center;
        public static double radius;

        private static double F(double x)
        {
            return (-1) * ((9 + x * x) * (9 + x * x));
        }

        private static Complex[] GetDiscreteFunction()
        {
            var res = new Complex[n];
            for (int i = 0; i < n; ++i)
            {
                res[i] = F(left + i * h);
            }
            return res;
        }

        private static Complex[,] GetMatrixT()
        {
            var matrix = new Complex[n, n];
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

        private static Complex[,] GetMatrixA()
        {
            var f = GetDiscreteFunction();
            var diagonalMatrix = new Complex[n, n];
            for (int i = 0; i < n; ++i)
            {
                diagonalMatrix[i, i] = f[i];
            }
            return GetMatrixT().MultiplyComplex(diagonalMatrix);
        }

        private static Complex[,] GetMatrixD(Complex[,] A, Complex value)
        {
            var D = new Complex[n, n];
            D = D.FillDiagonalComplex(value);
            return A.SubComplex(D);
        }

        public static void LU(Complex[,] D, out Complex[,] L, out Complex[,] U)
        {
            L = new Complex[n, n];
            U = new Complex[n, n];

            for (int i = 0; i < n; ++i)
            {
                U[0, i] = D[0, i];
                for (int j = i; j < n; ++j)
                {
                    var sum = Complex.Zero;
                    for (int k = 0; k < i; ++k)
                    {
                        sum += L[i, k] * U[k, j];
                    }
                    U[i, j] = D[i, j] - sum;
                    if (i <= j)
                    {
                        sum = Complex.Zero;
                        for (int k = 0; k < i; ++k)
                        {
                            sum += L[j, k] * U[k, i];
                        }
                        L[j, i] = (D[j, i] - sum) / U[i, i];
                    }
                }
            }
        }

        public static void MV(Complex[,] B, out Complex[,] M, Complex[,] U, Complex[,] L, out Complex[,] V)
        {
            M = new Complex[n, n];
            V = new Complex[n, n];

            for (int i = 0; i < n; ++i)
            {
                V[0, i] = B[0, i];
                for (int j = i; j < n; ++j)
                {
                    var sum = Complex.Zero;
                    for (int k = 0; k < i; ++k)
                    {
                        sum += M[i, k] * U[k, j] + L[i, k] * V[k, j];
                    }
                    V[i, j] = B[i, j] - sum;
                    if (i <= j)
                    {
                        sum = Complex.Zero;
                        for (int k = 0; k < i; ++k)
                        {
                            sum += M[j, k] * U[k, i] + L[j, k] * V[k, i];
                        }
                        M[j, i] = (B[j, i] - sum - L[j, i] * V[i, i]) / U[i, i];
                    }
                }
            }
        }

        public static Complex Determinant(Complex[,] U)
        {
            var result = Complex.One;
            for (int i = 0; i < U.GetLength(0); ++i)
            {
                result *= U[i, i];
            }
            return result;
        }

        public static Complex DeterminantDerivative(Complex[,] V, Complex[,] U)
        {
            var sum = Complex.Zero;
            for (int k = 0; k < V.GetLength(0); ++k)
            {
                var product = Complex.One;
                for (int i = 0; i < U.GetLength(0); ++i)
                {
                    if (i == k)
                    {
                        continue;
                    }
                    product *= U[i, i];
                }
                sum += V[k, k] * product;
            }
            return sum;
        }

        //s0
        public static double RootsCount()
        {
            var sum = Complex.Zero;

            for (int i = 0; i < n; ++i)
            {
                var spectralRadius = Spec(i + 1, n);
                var lambda = center + spectralRadius;

                Complex[,] D = GetMatrixD(GetMatrixA(), lambda);
                Complex[,] B = new Complex[n, n];
                B.FillDiagonalComplex(Complex.One);

                Complex[,] L, U, M, V;

                LU(D, out L, out U);
                MV(B, out M, U, L, out V);

                ///////////////First approach	
                var determinant = Determinant(U);
                var determinantDerivative = DeterminantDerivative(V, U);
                sum += spectralRadius * determinantDerivative / determinant;

                ///////////////Second approach	
                //Complex prodSum = 0;	
                //for (int j = 0; j < n; ++j)	
                //{	
                //    prodSum += V[j, j] / U[j, j];	
                //}	
                //sum += spectralRadius * prodSum;
            }

            return Complex.Abs(sum / n);
        }

        public static double InitialApproximation(int j, int n)
        {
            return Complex.Abs(center + Spec(j, n));
        }

        private static Complex Spec(int j, int n)
        {
            return radius * Complex.Exp(2.0 * Math.PI * Complex.ImaginaryOne * j / n);
        }

        static void Main(string[] args)
        {
            left = 0.0;
            right = 1.0;
            n = 4;
            h = (right - left) / (n - 1);
            center = 32;
            radius = 1;

            double s0 = RootsCount();
            int rootsCount = (int)Math.Round(s0);
            Console.WriteLine($"[{center - radius}, {center + radius}] => s0 = {s0}\nFound {rootsCount} root(s):");

            if (rootsCount > 0)
            {
                for (int i = 0; i < rootsCount; ++i)
                {
                    double initalApproximation = InitialApproximation(i + 1, rootsCount);
                    Console.WriteLine("{0}. lambda = {1}", i + 1, initalApproximation);
                }
            }
            else
            {
                Console.WriteLine("-------------");
            }
            Console.ReadKey();
        }
    }
}
