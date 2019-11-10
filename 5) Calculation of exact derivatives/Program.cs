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

        private static double GetDeltaLambda(double[,] matrix, int n)
        {
            var b = new double[n, n];
            var c = new double[n, n];
            var d = matrix.Copy();
            b.FillDiagonal(1);

            var u = new double[n, n];
            var l = new double[n, n];

            var v = new double[n, n];
            var m = new double[n, n];
            var w = new double[n, n];
            var nn = new double[n, n];

            for (int r = 0; r < n; ++r)
            {
                int k = r;
                for (int i = r + 1; i < n; ++i, ++k)
                {
                    u[r, k] = d[r, k] - l.Row(r).Multiply(u.Col(k)).SumFirst(r);
                    l[i, r] = (d[i, r] - l.Row(i).Multiply(u.Col(r)).SumFirst(r)) / u[r, r];

                    v[r, k] = b[r, k] - m.Row(r).Multiply(u.Col(k)).Sum(l.Row(r).Multiply(v.Col(k))).SumFirst(r);
                    m[i, r] = (b[i, r] - m.Row(i).Multiply(u.Col(r)).Sum(l.Row(i).Multiply(v.Col(r))).SumFirst(r) - l[i, r] * v[r, r]) / u[r, r];

                    w[r, k] = c[r, k] - nn.Row(r).Multiply(u.Col(k)).Sum(m.Row(r).Multiply(v.Col(k), 2)).Sum(l.Row(r).Multiply(w.Col(k))).SumFirst(r);

                    nn[i, r] = (c[i, r] - nn.Row(i).Multiply(u.Col(r)).Sum(m.Row(i).Multiply(v.Col(r), 2)).Sum(l.Row(i).Multiply(w.Col(r)))
                        .SumFirst(r) - 2 * m[i, r] *v[r, r] - l[i, r] * w[r, r]) / u[r, r];
                }
                k = n - 1;

                u[r, k] = d[r, k] - l.Row(r).Multiply(u.Col(k)).SumFirst(r);
                v[r, k] = b[r, k] - m.Row(r).Multiply(u.Col(k)).Sum(l.Row(r).Multiply(v.Col(k))).SumFirst(r);
                w[r, k] = c[r, k] - nn.Row(r).Multiply(u.Col(k)).Sum(m.Row(r).Multiply(v.Col(k), 2)).Sum(l.Row(r).Multiply(w.Col(k))).SumFirst(r);
            }

            l.FillDiagonal(1);
            m.FillDiagonal(0);
            nn.FillDiagonal(0);

            double sum = 0;
            for(int i=0; i<n; ++i)
            {
                sum += v[i, i] / u[i, i];
            }

            //Delta lamdba
            return 1 / sum; 
        }

        private static double[,] SubDiagonal(double[,] matrix, double value, int n)
        {
            var sub = new double[n, n];
            sub = sub.FillDiagonal(value);
            return matrix.Sub(sub);
        }

        static void Main(string[] args)
        {
            left = 0;
            right = 1;
            n = 4;
            h = (right - left) / (n - 1);
            var deviation = 0.00001;
            var A = GetMatrixA();

            double lprev = 33;
            double lnext = 0;
            double deltaLambda = double.MaxValue;
            do
            {
                var matrix = SubDiagonal(A, lprev, n);
                deltaLambda = GetDeltaLambda(matrix, n);

                lnext = lprev + deltaLambda;
                lprev = lnext;
            } while (deltaLambda > deviation);

            Console.WriteLine($"First = {lnext}");
            var decA = new Accord.Math.Decompositions.EigenvalueDecomposition(A, false, true);
            var eigValsA = decA.RealEigenvalues.OrderBy(c => c);
            Console.WriteLine($"Real  = { string.Join("\t", eigValsA.Select(c => c))}");

            Console.ReadKey();
        }
    }
}
