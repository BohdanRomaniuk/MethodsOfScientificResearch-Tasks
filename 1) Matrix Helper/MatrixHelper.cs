using System;
using System.Numerics;

namespace Helpers
{
    public static class MatrixHelper
    {
        public static void Print(this double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < matrix.GetLength(1); ++j)
                {
                    Console.Write($"{matrix[i, j]:F10}\t");
                }
                Console.WriteLine();
            }
        }

        public static void Print(this double[] vector)
        {
            for (int i = 0; i < vector.Length; ++i)
            {
                Console.Write($"{vector[i]:F10}\t");
            }
            Console.WriteLine();
        }

        public static double[] Multiply(this double[] first, double[] second, double constant = 1.0)
        {
            var result = new double[first.Length];
            for (int i = 0; i < first.Length; ++i)
            {
                result[i] = first[i] * second[i] * constant;
            }
            return result;
        }

        public static double[,] Multiply(this double[,] first, double[,] second)
        {
            int rows = first.GetLength(0);
            int cols = first.GetLength(1);
            var result = new double[rows, cols];
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    for (int k = 0; k < cols; ++k)
                    {
                        result[i, j] += first[i, k] * second[k, j];
                    }
                }
            }
            return result;
        }

        public static Complex[,] MultiplyComplex(this Complex[,] first, Complex scalar)
        {
            int rows = first.GetLength(0);
            int cols = first.GetLength(1);
            var result = new Complex[rows, cols];
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    result[i, j] += first[i, j] * scalar;
                }
            }
            return result;
        }

        public static double[,] Sum(this double[,] first, double[,] second)
        {
            int rows = first.GetLength(0);
            int cols = first.GetLength(1);
            var result = new double[rows, cols];
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    result[i, j] += first[i, j] + second[i, j];
                }
            }
            return result;
        }

        public static double[,] Sub(this double[,] first, double[,] second)
        {
            int rows = first.GetLength(0);
            int cols = first.GetLength(1);
            var result = new double[rows, cols];
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    result[i, j] += first[i, j] - second[i, j];
                }
            }
            return result;
        }

        public static Complex[,] SubComplex(this Complex[,] first, Complex[,] second)
        {
            int rows = first.GetLength(0);
            int cols = first.GetLength(1);
            var result = new Complex[rows, cols];
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    result[i, j] += first[i, j] - second[i, j];
                }
            }
            return result;
        }

        public static double[] Sum(this double[] first, double[] second)
        {
            var result = new double[first.Length];
            for (int i = 0; i < first.Length; ++i)
            {
                result[i] = first[i] + second[i];
            }
            return result;
        }

        public static double[] Sub(this double[] first, double[] second)
        {
            var result = new double[first.Length];
            for (int i = 0; i < first.Length; ++i)
            {
                result[i] = first[i] - second[i];
            }
            return result;
        }

        public static double[,] FillDiagonal(this double[,] matrix, double value)
        {
            int n = matrix.GetLength(0);
            for (int i = 0; i < n; ++i)
            {
                matrix[i, i] = value;
            }
            return matrix;
        }

        public static double[,] CopyFrom(this double[,] matrix, double[,] second)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    matrix[i, j] = second[i, j];
                }
            }
            return matrix;
        }

        public static double[] Row(this double[,] matrix, int row)
        {
            int cols = matrix.GetLength(1);
            var res = new double[cols];
            for (int i = 0; i < cols; ++i)
            {
                res[i] = matrix[row, i];
            }
            return res;
        }

        public static double[] Col(this double[,] matrix, int column)
        {
            int rows = matrix.GetLength(0);
            var res = new double[rows];
            for (int i = 0; i < rows; ++i)
            {
                res[i] = matrix[i, column];
            }
            return res;
        }

        public static double SumFirst(this double[] vector, int size)
        {
            double res = 0;
            for (int i = 0; i < size; ++i)
            {
                res += vector[i];
            }
            return res;
        }
    }
}
