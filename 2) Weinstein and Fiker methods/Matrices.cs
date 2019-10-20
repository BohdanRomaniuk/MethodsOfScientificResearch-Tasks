using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;

namespace Weinstein_and_Fiker_methods
{
    public static class Matrices
    {
        private static readonly double[] Coef = { 1, -2, 1 };
        private static Matrix derivMatrix;

        private static Matrix DerivMatrix
        {
            get
            {
                if (derivMatrix == null)
                {
                    derivMatrix = new DenseMatrix(Constants.N, Constants.N);

                    for (int i = 0; i < Constants.N; i++)
                    {
                        for (int j = 0; j < Coef.Count(); j++)
                        {
                            int k = i + Coef.Count() / 2 - j;
                            if (0 <= k && k < Constants.N)
                            {
                                derivMatrix[i, k] = Coef[j];
                            }
                        }
                    }
                    //

                    //
                }

                return derivMatrix;
            }
        }

        public static Matrix A
        {
            get
            {
                Matrix matrix = new DenseMatrix(Constants.N, Constants.N);
                double H = Constants.H;
                int N = Constants.N;
                double[] values = new double[3] { 1, -2, 1 };
                values[0] = (-1) * (1 / (H * H) + 1 / (2 * H));
                values[1] = 2 / (H * H);
                values[2] = (-1) * (1 / (H * H) - 1 / (2 * H));

                matrix[0, 0] = 1;

                matrix[N - 1, N - 1] = 1;
                for (int i = 1; i < N - 1; i++)
                {
                    int k = 0;
                    for (int j = i - 1; j <= i + 1; j++)
                    {
                        matrix[i, j] = values[k];
                        k++;
                    }
                }
                return matrix;
            }
        }

        public static Matrix B
        {
            get
            {
                Matrix matrix = new DenseMatrix(Constants.N, Constants.N);
                double H = Constants.H;
                int N = Constants.N;
                double[] values = new double[3] { 1, -2, 1 };
                values[0] = (-1) * (1 / (H * H) + 1 / (4 * H));
                values[1] = 2 / (H * H);
                values[2] = (-1) * (1 / (H * H) - 1 / (4 * H));

                matrix[0, 0] = 0.5;

                matrix[N - 1, N - 1] = 0.5;
                for (int i = 1; i < N - 1; i++)
                {
                    int k = 0;
                    for (int j = i - 1; j <= i + 1; j++)
                    {
                        matrix[i, j] = values[k];
                        k++;
                    }
                }
                return matrix;
            }
        }

        public static Matrix C
        {
            get { return (Matrix)(A - B); }
        }
    }
}
